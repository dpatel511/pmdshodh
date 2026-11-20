//! HTML damage report from JSON run summaries (docs/report-proposal.md, agreed 2026-11-20).
//! Pure: parsed summaries in, one self-contained HTML string out (inline SVG, no JavaScript).
//! The same summaries always give the same bytes.

use std::fmt::Write as _;

use serde_json::Value;

/// A run summary read for the report. Counts are as in the JSON (summary.rs); rates are computed here.
#[derive(Debug, Clone)]
pub struct Sample {
    /// Set by the caller (see [`labels`]).
    pub label: String,
    /// The summary file as given.
    pub source: String,
    pub format_version: u64,
    pub tool_version: String,
    pub subcommand: String,
    pub command: String,
    pub input: String,
    /// `None`: not recorded (version 1).
    pub input_format: Option<String>,
    /// `None`: not recorded (version 1).
    pub reference: Option<Reference>,
    pub library: String,
    pub treatment: String,
    pub pmdtools_compat: bool,
    pub model: Vec<(&'static str, f64)>,
    pub threshold: Option<f64>,
    pub upper_threshold: Option<f64>,
    pub read_filters: Vec<(&'static str, String)>,
    pub reads: Reads,
    pub score_histogram: ScoreHistogram,
    /// (length, reads), by length.
    pub lengths_scored: Vec<(u64, u64)>,
    pub lengths_kept: Option<Vec<(u64, u64)>>,
    pub damage: Damage,
}

#[derive(Debug, Clone)]
pub struct Reference {
    pub path: Option<String>,
    pub bases_from: String,
    /// With `-r`: MD tags that agree, disagree, are absent.
    pub md_tags: Option<[u64; 3]>,
}

#[derive(Debug, Clone)]
pub struct Reads {
    pub input: u64,
    pub scored: u64,
    /// `None` unless the summary is from `filter`.
    pub kept: Option<u64>,
    pub removed_by_threshold: Option<u64>,
    pub no_informative_sites: u64,
    pub skipped: Vec<(String, u64)>,
    pub flags: Vec<(String, u64)>,
}

#[derive(Debug, Clone)]
pub struct ScoreHistogram {
    pub min: f64,
    pub bin_width: f64,
    pub counts: Vec<u64>,
    pub underflow: u64,
    pub overflow: u64,
}

#[derive(Debug, Clone)]
pub struct Damage {
    pub range: usize,
    /// `[5′, 3′][reference C, G][read base A, C, G, T][position]`.
    pub tables: [[[Vec<u64>; 4]; 2]; 2],
    /// First-base C→T rate, its standard error and the number of sites, as in the JSON.
    pub first: (Option<f64>, Option<f64>, u64),
}

const BASES: [&str; 4] = ["A", "C", "G", "T"];
const A: usize = 0;
const C: usize = 1;
const G: usize = 2;
const T: usize = 3;

/// Field access that names the file and field in every error.
struct Fields<'a> {
    file: &'a str,
    root: &'a Value,
}

impl<'a> Fields<'a> {
    fn err(&self, path: &str, what: &str) -> String {
        format!("{}: field {path} {what}", self.file)
    }
    fn get(&self, path: &str) -> Result<&'a Value, String> {
        let mut v = self.root;
        for key in path.split('.') {
            v = v.get(key).ok_or_else(|| self.err(path, "is missing"))?;
        }
        Ok(v)
    }
    fn u64(&self, path: &str) -> Result<u64, String> {
        self.get(path)?
            .as_u64()
            .ok_or_else(|| self.err(path, "is not a count"))
    }
    fn opt_u64(&self, path: &str) -> Result<Option<u64>, String> {
        match self.get(path)? {
            Value::Null => Ok(None),
            _ => self.u64(path).map(Some),
        }
    }
    fn f64(&self, path: &str) -> Result<f64, String> {
        self.get(path)?
            .as_f64()
            .ok_or_else(|| self.err(path, "is not a number"))
    }
    fn opt_f64(&self, path: &str) -> Result<Option<f64>, String> {
        match self.get(path)? {
            Value::Null => Ok(None),
            _ => self.f64(path).map(Some),
        }
    }
    fn str(&self, path: &str) -> Result<String, String> {
        self.get(path)?
            .as_str()
            .map(str::to_string)
            .ok_or_else(|| self.err(path, "is not a string"))
    }
    fn opt_str(&self, path: &str) -> Result<Option<String>, String> {
        match self.get(path)? {
            Value::Null => Ok(None),
            _ => self.str(path).map(Some),
        }
    }
    fn bool(&self, path: &str) -> Result<bool, String> {
        self.get(path)?
            .as_bool()
            .ok_or_else(|| self.err(path, "is not true or false"))
    }
    fn counts(&self, path: &str) -> Result<Vec<u64>, String> {
        let not = || self.err(path, "is not a list of counts");
        self.get(path)?
            .as_array()
            .ok_or_else(not)?
            .iter()
            .map(|v| v.as_u64().ok_or_else(not))
            .collect()
    }
    fn count_map(&self, path: &str) -> Result<Vec<(String, u64)>, String> {
        let not = || self.err(path, "is not a map of counts");
        self.get(path)?
            .as_object()
            .ok_or_else(not)?
            .iter()
            .map(|(k, v)| Ok((k.clone(), v.as_u64().ok_or_else(not)?)))
            .collect()
    }
    fn one_of(&self, path: &str, allowed: &[&str]) -> Result<String, String> {
        let v = self.str(path)?;
        if allowed.contains(&v.as_str()) {
            Ok(v)
        } else {
            Err(self.err(
                path,
                &format!("must be {}, not {v:?}", allowed.join(" or ")),
            ))
        }
    }
    /// Read-length histogram: keys are lengths, sorted by number.
    fn lengths(&self, path: &str) -> Result<Vec<(u64, u64)>, String> {
        let mut v = self
            .count_map(path)?
            .into_iter()
            .map(|(k, n)| {
                k.parse().map(|len| (len, n)).map_err(|_| {
                    self.err(path, &format!("has a length that is not a number: {k:?}"))
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        v.sort_unstable();
        Ok(v)
    }
}

/// Read one summary (`format_version` 1 or 2). `file` names it in errors, which also name the field.
pub fn parse(file: &str, json: &Value) -> Result<Sample, String> {
    let f = Fields { file, root: json };
    let version = f.u64("format_version")?;
    if !(1..=2).contains(&version) {
        return Err(format!(
            "{file}: format_version {version} is not supported (this pmdshodh reads versions 1 and 2)"
        ));
    }
    if f.str("tool.name")? != "pmdshodh" {
        return Err(f.err("tool.name", "is not pmdshodh"));
    }
    // Version 1 lacks these; version 2 requires them.
    let recorded = |path: &str| version >= 2 || f.get(path).is_ok();
    let input_format = if recorded("input_format") {
        Some(f.one_of("input_format", &["sam", "bam", "cram"])?)
    } else {
        None
    };
    let reference = if recorded("reference") {
        Some(Reference {
            path: f.opt_str("reference.path")?,
            bases_from: f.one_of("reference.bases_from", &["md", "reference"])?,
            md_tags: match f.get("reference.md_tags")? {
                Value::Null => None,
                _ => Some([
                    f.u64("reference.md_tags.agree")?,
                    f.u64("reference.md_tags.disagree")?,
                    f.u64("reference.md_tags.absent")?,
                ]),
            },
        })
    } else {
        None
    };

    let range = f.u64("damage.range")? as usize;
    if range == 0 {
        return Err(f.err("damage.range", "must be at least 1"));
    }
    let mut tables: [[[Vec<u64>; 4]; 2]; 2] = Default::default();
    for (e, end) in ["five_prime", "three_prime"].iter().enumerate() {
        for (r, refb) in ["C", "G"].iter().enumerate() {
            for (b, base) in BASES.iter().enumerate() {
                let path = format!("damage.{end}.{refb}.{base}");
                let col = f.counts(&path)?;
                if col.len() != range {
                    return Err(f.err(
                        &path,
                        &format!("has {} positions, not damage.range {range}", col.len()),
                    ));
                }
                tables[e][r][b] = col;
            }
        }
    }

    let bin_width = f.f64("score_histogram.bin_width")?;
    if !(bin_width.is_finite() && bin_width > 0.0) {
        return Err(f.err("score_histogram.bin_width", "must be above 0"));
    }
    let counts = f.counts("score_histogram.counts")?;
    if counts.is_empty() {
        return Err(f.err("score_histogram.counts", "is empty"));
    }
    let min = f.f64("score_histogram.min")?;

    Ok(Sample {
        label: String::new(),
        source: file.to_string(),
        format_version: version,
        tool_version: f.str("tool.version")?,
        subcommand: f.str("subcommand")?,
        command: f.str("command")?,
        input: f.str("input")?,
        input_format,
        reference,
        library: f.one_of("mode.library", &["ds", "ss"])?,
        treatment: f.one_of("mode.treatment", &["none", "half", "full"])?,
        pmdtools_compat: f.bool("mode.pmdtools_compat")?,
        model: vec![
            ("P (pmd_p)", f.f64("model.pmd_p")?),
            ("C (pmd_c)", f.f64("model.pmd_c")?),
            ("π ancient", f.f64("model.pi_ancient")?),
            ("π contamination", f.f64("model.pi_contamination")?),
            ("null background", f.f64("model.null_background")?),
        ],
        threshold: f.opt_f64("filters.threshold")?,
        upper_threshold: f.opt_f64("filters.upper_threshold")?,
        read_filters: vec![
            ("minimum MAPQ", f.u64("filters.min_mapq")?.to_string()),
            (
                "minimum base quality",
                f.u64("filters.min_baseq")?.to_string(),
            ),
            ("no soft clips", yes_no(f.bool("filters.no_clips")?)),
            ("no indels", yes_no(f.bool("filters.no_indels")?)),
        ],
        reads: Reads {
            input: f.u64("reads.input")?,
            scored: f.u64("reads.scored")?,
            kept: f.opt_u64("reads.kept")?,
            removed_by_threshold: f.opt_u64("reads.removed_by_threshold")?,
            no_informative_sites: f.u64("reads.no_informative_sites")?,
            skipped: f.count_map("reads.skipped")?,
            flags: f.count_map("reads.flags")?,
        },
        score_histogram: ScoreHistogram {
            min,
            bin_width,
            counts,
            underflow: f.u64("score_histogram.underflow")?,
            overflow: f.u64("score_histogram.overflow")?,
        },
        lengths_scored: f.lengths("length_histogram.scored")?,
        lengths_kept: match f.get("length_histogram.kept")? {
            Value::Null => None,
            _ => Some(f.lengths("length_histogram.kept")?),
        },
        damage: Damage {
            range,
            tables,
            first: (
                f.opt_f64("damage.first_position.c_to_t")?,
                f.opt_f64("damage.first_position.se")?,
                f.u64("damage.first_position.n")?,
            ),
        },
    })
}

fn yes_no(b: bool) -> String {
    if b { "yes" } else { "no" }.to_string()
}

/// `100 num / den`, or `None` when nothing was counted.
fn percent(num: u64, den: u64) -> Option<f64> {
    (den > 0).then(|| 100.0 * num as f64 / den as f64)
}

/// Damage rates by position from each end, in percent (index 0 = terminal base).
#[derive(Debug, Clone, PartialEq)]
pub struct Curves {
    pub five: Vec<Option<f64>>,
    pub three: Vec<Option<f64>>,
    /// The substitution shown at the 3′ end: "G→A" (ds) or "C→T" (ss).
    pub three_name: &'static str,
    /// All other substitutions at reference C and G, combined.
    pub five_other: Vec<Option<f64>>,
    pub three_other: Vec<Option<f64>>,
}

/// 5′ C→T; 3′ G→A for ds and C→T for ss, where damage shows as C→T at both ends (D6).
pub fn curves(s: &Sample) -> Curves {
    let t = &s.damage.tables;
    let sites = |e: usize, r: usize, i: usize| -> u64 { (0..4).map(|b| t[e][r][b][i]).sum() };
    // (reference row, read base) of the damage substitution at each end.
    let ss = s.library == "ss";
    let shown = [(0, T), if ss { (0, T) } else { (1, A) }];
    let damage = |e: usize| -> Vec<Option<f64>> {
        let (r, b) = shown[e];
        (0..s.damage.range)
            .map(|i| percent(t[e][r][b][i], sites(e, r, i)))
            .collect()
    };
    let other = |e: usize| -> Vec<Option<f64>> {
        (0..s.damage.range)
            .map(|i| {
                let mut n = 0;
                for (r, same) in [(0, C), (1, G)] {
                    for (b, col) in t[e][r].iter().enumerate() {
                        if b != same && (r, b) != shown[e] {
                            n += col[i];
                        }
                    }
                }
                percent(n, sites(e, 0, i) + sites(e, 1, i))
            })
            .collect()
    };
    Curves {
        five: damage(0),
        three: damage(1),
        three_name: if ss { "C→T" } else { "G→A" },
        five_other: other(0),
        three_other: other(1),
    }
}

/// Plot area: data ranges mapped to SVG pixels.
#[derive(Debug, Clone, Copy)]
pub struct Frame {
    pub x: (f64, f64),
    pub y: (f64, f64),
    /// Left, top, width and height of the plot area in pixels.
    pub area: (f64, f64, f64, f64),
}

impl Frame {
    /// A single x value is centred.
    fn px(&self, x: f64) -> f64 {
        let (left, _, w, _) = self.area;
        if self.x.1 > self.x.0 {
            left + (x - self.x.0) / (self.x.1 - self.x.0) * w
        } else {
            left + w / 2.0
        }
    }
    fn py(&self, y: f64) -> f64 {
        let (_, top, _, h) = self.area;
        if self.y.1 > self.y.0 {
            top + h - (y - self.y.0) / (self.y.1 - self.y.0) * h
        } else {
            top + h
        }
    }
}

/// `<polyline>`s for one series; a missing value breaks the line. A lone point is a circle.
pub fn series_svg(f: &Frame, points: &[(f64, Option<f64>)], color: &str) -> String {
    let mut out = String::new();
    let mut run: Vec<(f64, f64)> = Vec::new();
    let flush = |run: &mut Vec<(f64, f64)>, out: &mut String| {
        match run.len() {
            0 => {}
            1 => {
                let _ = write!(
                    out,
                    r#"<circle cx="{:.2}" cy="{:.2}" r="2.5" fill="{color}"/>"#,
                    run[0].0, run[0].1
                );
            }
            _ => {
                let pts: Vec<String> = run.iter().map(|(x, y)| format!("{x:.2},{y:.2}")).collect();
                let _ = write!(
                    out,
                    r#"<polyline fill="none" stroke="{color}" stroke-width="1.5" points="{}"/>"#,
                    pts.join(" ")
                );
            }
        }
        run.clear();
    };
    for &(x, y) in points {
        match y {
            Some(y) => run.push((f.px(x), f.py(y))),
            None => flush(&mut run, &mut out),
        }
    }
    flush(&mut run, &mut out);
    out
}

/// File name of a path written on any system.
fn file_name(p: &str) -> &str {
    p.rsplit(['/', '\\']).next().unwrap_or(p)
}

/// Labels for samples given as (summary path as given, input path recorded in it): the input's
/// file name (the summary's for stdin), then the summary's file name, then its path, until unique.
pub fn labels(samples: &[(&str, &str)]) -> Vec<String> {
    let base: Vec<String> = samples
        .iter()
        .map(|(summary, input)| file_name(if *input == "-" { summary } else { input }).to_string())
        .collect();
    let clashes = |v: &[String]| -> Vec<bool> {
        v.iter()
            .map(|l| v.iter().filter(|m| *m == l).count() > 1)
            .collect()
    };
    let mut out = base.clone();
    for (i, c) in clashes(&base).into_iter().enumerate() {
        if c {
            out[i] = format!("{} ({})", base[i], file_name(samples[i].0));
        }
    }
    for (i, c) in clashes(&out).into_iter().enumerate() {
        if c {
            out[i] = format!("{} ({})", base[i], samples[i].0);
        }
    }
    out
}

/// Escape text for HTML text and attribute values.
pub fn escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&#39;"),
            c => out.push(c),
        }
    }
    out
}

/// A number without trailing zeros (at most 6 decimals).
fn num(x: f64) -> String {
    let s = format!("{x:.6}");
    let s = s.trim_end_matches('0').trim_end_matches('.');
    if s == "-0" { "0".into() } else { s.into() }
}

/// Counts with thousands separators.
fn count(n: u64) -> String {
    let s = n.to_string();
    let mut out = String::new();
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

fn values_attr(values: &[Option<f64>]) -> String {
    values
        .iter()
        .map(|v| v.map_or("n/a".to_string(), |v| format!("{v:.3}")))
        .collect::<Vec<_>>()
        .join(" ")
}

/// Okabe–Ito colours, readable for colour-blind viewers and on light and dark backgrounds.
const PALETTE: [&str; 7] = [
    "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#999999",
];
const DAMAGE_COLOR: &str = "#D55E00";
const OTHER_COLOR: &str = "#999999";

/// A step of 1, 2 or 5 × 10^k giving about five intervals up to `max`.
fn nice_step(max: f64) -> f64 {
    let raw = max / 5.0;
    let mag = 10f64.powf(raw.log10().floor());
    let n = raw / mag;
    mag * if n <= 1.0 {
        1.0
    } else if n <= 2.0 {
        2.0
    } else if n <= 5.0 {
        5.0
    } else {
        10.0
    }
}

/// Multiples of `step` in `[lo, hi]`.
fn ticks(lo: f64, hi: f64, step: f64) -> Vec<f64> {
    let first = (lo / step).ceil() as i64;
    let last = (hi / step).floor() as i64;
    (first..=last).map(|i| i as f64 * step).collect()
}

struct Series {
    name: String,
    color: &'static str,
    faint: bool,
    points: Vec<(f64, Option<f64>)>,
}

/// Bars from `x0` to `x1` with height `y`.
struct Bars {
    name: String,
    color: &'static str,
    bars: Vec<(f64, f64, f64)>,
    /// Plotted values for scripts, e.g. "100:4".
    values: String,
}

struct Chart<'a> {
    title: &'a str,
    x_label: &'a str,
    y_label: &'a str,
    x: (f64, f64),
    x_ticks: Vec<f64>,
    /// Integer y axis (read counts).
    integer_y: bool,
    lines: Vec<Series>,
    bars: Vec<Bars>,
    /// Vertical lines: (x, label).
    marks: Vec<(f64, String)>,
}

const W: f64 = 440.0;
const H: f64 = 270.0;
const AREA: (f64, f64, f64, f64) = (60.0, 14.0, 364.0, 206.0);

fn chart(c: &Chart) -> String {
    let data_max = c
        .lines
        .iter()
        .flat_map(|s| s.points.iter().filter_map(|p| p.1))
        .chain(c.bars.iter().flat_map(|b| b.bars.iter().map(|b| b.2)))
        .fold(0.0f64, f64::max);
    let (y_max, step) = if data_max > 0.0 {
        let mut step = nice_step(data_max);
        if c.integer_y {
            step = step.max(1.0);
        }
        ((data_max / step).ceil() * step, step)
    } else if c.integer_y {
        (1.0, 1.0)
    } else {
        (1.0, 0.2)
    };
    let f = Frame {
        x: c.x,
        y: (0.0, y_max),
        area: AREA,
    };
    let (left, top, w, h) = AREA;
    let bottom = top + h;
    let mut s = String::new();
    let _ = writeln!(
        s,
        "<figure><svg viewBox=\"0 0 {W} {H}\" width=\"{W}\" height=\"{H}\" role=\"img\" aria-label=\"{}\">",
        escape(c.title)
    );
    for y in ticks(0.0, y_max, step) {
        let py = f.py(y);
        let _ = writeln!(
            s,
            "<line class=\"grid\" x1=\"{left}\" y1=\"{py:.2}\" x2=\"{:.2}\" y2=\"{py:.2}\"/><text x=\"{:.2}\" y=\"{:.2}\" text-anchor=\"end\">{}</text>",
            left + w,
            left - 6.0,
            py + 4.0,
            num(y)
        );
    }
    for x in &c.x_ticks {
        let px = f.px(*x);
        let _ = writeln!(
            s,
            "<line class=\"axis\" x1=\"{px:.2}\" y1=\"{bottom}\" x2=\"{px:.2}\" y2=\"{:.2}\"/><text x=\"{px:.2}\" y=\"{:.2}\" text-anchor=\"middle\">{}</text>",
            bottom + 4.0,
            bottom + 17.0,
            num(*x)
        );
    }
    let _ = writeln!(
        s,
        "<line class=\"axis\" x1=\"{left}\" y1=\"{bottom}\" x2=\"{:.2}\" y2=\"{bottom}\"/><line class=\"axis\" x1=\"{left}\" y1=\"{top}\" x2=\"{left}\" y2=\"{bottom}\"/>",
        left + w
    );
    let _ = writeln!(
        s,
        "<text x=\"{:.2}\" y=\"{:.2}\" text-anchor=\"middle\">{}</text><text transform=\"translate(14 {:.2}) rotate(-90)\" text-anchor=\"middle\">{}</text>",
        left + w / 2.0,
        H - 8.0,
        escape(c.x_label),
        top + h / 2.0,
        escape(c.y_label)
    );
    for b in &c.bars {
        let _ = write!(
            s,
            "<g data-series=\"{}\" data-values=\"{}\" fill=\"{}\" fill-opacity=\"0.75\">",
            escape(&b.name),
            b.values,
            b.color
        );
        for &(x0, x1, y) in &b.bars {
            let (a, z) = (f.px(x0), f.px(x1));
            let py = f.py(y);
            let _ = write!(
                s,
                "<rect x=\"{a:.2}\" y=\"{py:.2}\" width=\"{:.2}\" height=\"{:.2}\"/>",
                (z - a).max(1.0),
                bottom - py
            );
        }
        s.push_str("</g>\n");
    }
    for l in &c.lines {
        let values: Vec<Option<f64>> = l.points.iter().map(|p| p.1).collect();
        let _ = writeln!(
            s,
            "<g data-series=\"{}\" data-values=\"{}\"{}>{}</g>",
            escape(&l.name),
            values_attr(&values),
            if l.faint {
                " stroke-dasharray=\"4 3\" opacity=\"0.8\""
            } else {
                ""
            },
            series_svg(&f, &l.points, l.color)
        );
    }
    for (x, label) in &c.marks {
        let px = f.px(*x);
        let _ = writeln!(
            s,
            "<line class=\"mark\" data-threshold=\"{}\" x1=\"{px:.2}\" y1=\"{top}\" x2=\"{px:.2}\" y2=\"{bottom}\"/><text x=\"{:.2}\" y=\"{:.2}\">{}</text>",
            num(*x),
            px + 3.0,
            top + 10.0,
            escape(label)
        );
    }
    s.push_str("</svg>\n<figcaption>");
    let _ = write!(s, "<strong>{}</strong>", escape(c.title));
    let legend: Vec<(&str, &str)> = c
        .bars
        .iter()
        .map(|b| (b.name.as_str(), b.color))
        .chain(c.lines.iter().map(|l| (l.name.as_str(), l.color)))
        .collect();
    if legend.len() > 1 {
        s.push_str("<ul class=\"legend\">");
        for (name, color) in legend {
            let _ = write!(
                s,
                "<li><span class=\"key\" style=\"background:{color}\"></span>{}</li>",
                escape(name)
            );
        }
        s.push_str("</ul>");
    }
    s.push_str("</figcaption></figure>\n");
    s
}

fn position_ticks(range: usize) -> Vec<f64> {
    let step = nice_step(range as f64).max(1.0);
    let mut t = vec![1.0];
    t.extend(ticks(step, range as f64, step));
    t.dedup();
    t
}

fn positions(values: &[Option<f64>]) -> Vec<(f64, Option<f64>)> {
    values
        .iter()
        .enumerate()
        .map(|(i, v)| ((i + 1) as f64, *v))
        .collect()
}

fn damage_charts(s: &Sample) -> String {
    let c = curves(s);
    let range = s.damage.range;
    let panel = |end: &str, name: &str, values: &[Option<f64>], other: &[Option<f64>]| {
        chart(&Chart {
            title: &format!("{end} {name}"),
            x_label: &format!("Position from the {end} end"),
            y_label: "% of sites",
            x: (1.0, range as f64),
            x_ticks: position_ticks(range),
            integer_y: false,
            lines: vec![
                Series {
                    name: format!("{end} other substitutions"),
                    color: OTHER_COLOR,
                    faint: true,
                    points: positions(other),
                },
                Series {
                    name: format!("{end} {name}"),
                    color: DAMAGE_COLOR,
                    faint: false,
                    points: positions(values),
                },
            ],
            bars: vec![],
            marks: vec![],
        })
    };
    let mut out = String::from("<div class=\"charts\">\n");
    out += &panel("5′", "C→T", &c.five, &c.five_other);
    out += &panel("3′", c.three_name, &c.three, &c.three_other);
    out += "</div>\n";
    out
}

/// Percent of scored reads per score bin (`None` without scored reads).
fn score_percent(s: &Sample) -> Vec<Option<f64>> {
    s.score_histogram
        .counts
        .iter()
        .map(|&n| percent(n, s.reads.scored))
        .collect()
}

fn score_x(s: &Sample) -> (f64, f64) {
    let h = &s.score_histogram;
    (h.min, h.min + h.bin_width * h.counts.len() as f64)
}

fn x_ticks(x: (f64, f64)) -> Vec<f64> {
    ticks(x.0, x.1, nice_step(x.1 - x.0))
}

fn threshold_marks(s: &Sample, x: (f64, f64), notes: &mut Vec<String>) -> Vec<(f64, String)> {
    let mut marks = Vec::new();
    for (t, name) in [
        (s.threshold, "--threshold"),
        (s.upper_threshold, "--upper-threshold"),
    ] {
        let Some(t) = t else { continue };
        if t >= x.0 && t <= x.1 {
            marks.push((t, num(t)));
        } else {
            notes.push(format!("{name} {} is outside the plotted range.", num(t)));
        }
    }
    marks
}

fn score_chart(s: &Sample) -> String {
    let h = &s.score_histogram;
    let x = score_x(s);
    let mut notes = Vec::new();
    let marks = threshold_marks(s, x, &mut notes);
    let pct = score_percent(s);
    let mut out = String::from("<div class=\"figure-with-note\">\n");
    out += &chart(&Chart {
        title: "PMD score",
        x_label: "PMD score",
        y_label: "% of scored reads",
        x,
        x_ticks: x_ticks(x),
        integer_y: false,
        lines: vec![],
        bars: vec![Bars {
            name: "PMD score".into(),
            color: PALETTE[0],
            bars: pct
                .iter()
                .enumerate()
                .filter_map(|(i, p)| {
                    let x0 = h.min + i as f64 * h.bin_width;
                    p.map(|p| (x0, x0 + h.bin_width, p))
                })
                .collect(),
            values: values_attr(&pct),
        }],
        marks,
    });
    let _ = write!(
        out,
        "<p class=\"note\">Bins of {} from {} to {}. Below {}: {} reads; at or above {}: {} reads.",
        num(h.bin_width),
        num(x.0),
        num(x.1),
        num(x.0),
        count(h.underflow),
        num(x.1),
        count(h.overflow)
    );
    if s.reads.scored == 0 {
        out += " No scored reads.";
    }
    if let Some(kept) = s.reads.kept {
        let _ = write!(out, " {} kept.", kept_text(kept, s.reads.scored));
    }
    for n in notes {
        out += " ";
        out += &escape(&n);
    }
    out += "</p>\n</div>\n";
    out
}

fn kept_text(kept: u64, scored: u64) -> String {
    format!(
        "{} of {} scored reads ({})",
        count(kept),
        count(scored),
        percent(kept, scored).map_or("n/a".into(), |p| format!("{p:.1}%"))
    )
}

fn length_chart(s: &Sample) -> String {
    let all = s
        .lengths_scored
        .iter()
        .chain(s.lengths_kept.iter().flatten());
    let lo = all.clone().map(|l| l.0).min();
    let hi = all.map(|l| l.0).max();
    let (Some(lo), Some(hi)) = (lo, hi) else {
        return "<p class=\"note\">Read length: no scored reads.</p>\n".into();
    };
    let x = (lo as f64 - 0.5, hi as f64 + 0.5);
    let bars = |name: &str, color, lengths: &[(u64, u64)]| Bars {
        name: name.into(),
        color,
        bars: lengths
            .iter()
            .map(|&(l, n)| (l as f64 - 0.5, l as f64 + 0.5, n as f64))
            .collect(),
        values: lengths
            .iter()
            .map(|(l, n)| format!("{l}:{n}"))
            .collect::<Vec<_>>()
            .join(" "),
    };
    let mut b = vec![bars("scored reads", PALETTE[0], &s.lengths_scored)];
    if let Some(kept) = &s.lengths_kept {
        b.push(bars("kept reads", PALETTE[1], kept));
    }
    chart(&Chart {
        title: "Read length",
        x_label: "Read length (bases)",
        y_label: "Reads",
        x,
        x_ticks: ticks(x.0, x.1, nice_step(x.1 - x.0).max(1.0)),
        integer_y: true,
        lines: vec![],
        bars: b,
        marks: vec![],
    })
}

fn first_base_text(s: &Sample) -> String {
    match s.damage.first {
        (Some(r), Some(se), n) => format!(
            "{:.2}% ± {:.2}% (standard error), {} sites",
            100.0 * r,
            100.0 * se,
            count(n)
        ),
        _ => "n/a (no sites)".into(),
    }
}

fn row(out: &mut String, key: &str, value: &str) {
    let _ = writeln!(out, "<tr><th>{}</th><td>{value}</td></tr>", escape(key));
}

fn treatment_name(t: &str) -> &'static str {
    match t {
        "half" => "UDG-half",
        "full" => "full UDG",
        _ => "none",
    }
}

fn library_name(l: &str) -> &'static str {
    if l == "ss" {
        "single-stranded (ss)"
    } else {
        "double-stranded (ds)"
    }
}

fn sample_section(i: usize, s: &Sample) -> String {
    let mut out = String::new();
    let _ = writeln!(
        out,
        "<section id=\"sample-{}\" data-sample=\"{}\">\n<h2>{}</h2>",
        i + 1,
        escape(&s.label),
        escape(&s.label)
    );
    out += "<div class=\"tables\">\n<table>\n<caption>Run</caption>\n";
    let not_recorded = "<em>not recorded</em>".to_string();
    row(&mut out, "Summary file", &escape(&s.source));
    row(&mut out, "Input", &escape(&s.input));
    row(
        &mut out,
        "Input format",
        &s.input_format
            .as_deref()
            .map_or(not_recorded.clone(), escape),
    );
    row(
        &mut out,
        "Reference",
        &match &s.reference {
            None => not_recorded.clone(),
            Some(Reference { path: Some(p), .. }) => escape(p),
            Some(_) => "none (reference bases from MD tags)".into(),
        },
    );
    row(
        &mut out,
        "Command",
        &format!("<code>{}</code>", escape(&s.command)),
    );
    row(
        &mut out,
        "Written by",
        &format!(
            "pmdshodh {} ({}), summary format version {}",
            escape(&s.tool_version),
            escape(&s.subcommand),
            s.format_version
        ),
    );
    row(&mut out, "Library", library_name(&s.library));
    row(&mut out, "Treatment", treatment_name(&s.treatment));
    row(
        &mut out,
        "PMDtools compatibility",
        &yes_no(s.pmdtools_compat),
    );
    let model: Vec<String> = s
        .model
        .iter()
        .map(|(k, v)| format!("{k} {}", num(*v)))
        .collect();
    row(&mut out, "Model", &escape(&model.join(", ")));
    for (name, t) in [
        ("Threshold", s.threshold),
        ("Upper threshold", s.upper_threshold),
    ] {
        if let Some(t) = t {
            row(&mut out, name, &num(t));
        }
    }
    let filters: Vec<String> = s
        .read_filters
        .iter()
        .map(|(k, v)| format!("{k} {v}"))
        .collect();
    row(&mut out, "Read filters", &escape(&filters.join(", ")));
    out += "</table>\n<table>\n<caption>Reads</caption>\n";
    let r = &s.reads;
    row(&mut out, "Input", &count(r.input));
    row(&mut out, "Scored", &count(r.scored));
    if let Some(k) = r.kept {
        row(&mut out, "Kept", &kept_text(k, r.scored));
    }
    if let Some(n) = r.removed_by_threshold {
        row(&mut out, "Removed by threshold", &count(n));
    }
    row(
        &mut out,
        "No informative sites",
        &count(r.no_informative_sites),
    );
    for (reason, n) in &r.skipped {
        row(
            &mut out,
            &format!("Skipped: {}", reason.replace('_', " ")),
            &count(*n),
        );
    }
    for (flag, n) in &r.flags {
        row(
            &mut out,
            &format!("Flag: {}", flag.replace('_', " ")),
            &count(*n),
        );
    }
    if let Some(Reference {
        md_tags: Some([agree, disagree, absent]),
        ..
    }) = &s.reference
    {
        row(
            &mut out,
            "MD tags vs reference",
            &format!(
                "{} agree, {} disagree, {} absent",
                count(*agree),
                count(*disagree),
                count(*absent)
            ),
        );
    }
    out += "</table>\n</div>\n";
    let _ = writeln!(
        out,
        "<h3>Damage by position</h3>\n<p>First-base C→T: {}</p>",
        first_base_text(s)
    );
    out += &damage_charts(s);
    let mut note = String::from(
        "Percent of reference C (C→T) or G (G→A) sites in scored reads, in molecule orientation. \
         The dashed line is all other substitutions at reference C and G combined.",
    );
    if s.library == "ss" {
        note += " Single-stranded library: damage shows as C→T at both ends.";
    }
    match s.treatment.as_str() {
        "half" => {
            note += " UDG-half: only CpG sites, a C at the 5′ terminal base and a G at the 3′ \
                     terminal base are counted."
        }
        "full" => note += " Full UDG: only CpG sites are counted.",
        _ => {}
    }
    if s.pmdtools_compat && s.treatment != "none" {
        note += " With --pmdtools-compat, CpG context follows PMDtools' damage table.";
    }
    let _ = writeln!(out, "<p class=\"note\">{note}</p>");
    out += "<h3>PMD score and read length</h3>\n<div class=\"charts\">\n";
    out += &score_chart(s);
    out += &length_chart(s);
    out += "</div>\n</section>\n";
    out
}

fn overlays(samples: &[Sample]) -> String {
    let color = |i: usize| PALETTE[i % PALETTE.len()];
    let range = samples.iter().map(|s| s.damage.range).max().unwrap_or(1);
    let curves: Vec<Curves> = samples.iter().map(curves).collect();
    let damage = |title: &str, end: &str, pick: &dyn Fn(usize) -> (String, Vec<Option<f64>>)| {
        chart(&Chart {
            title,
            x_label: &format!("Position from the {end} end"),
            y_label: "% of sites",
            x: (1.0, range as f64),
            x_ticks: position_ticks(range),
            integer_y: false,
            lines: (0..samples.len())
                .map(|i| {
                    let (name, values) = pick(i);
                    Series {
                        name,
                        color: color(i),
                        faint: false,
                        points: positions(&values),
                    }
                })
                .collect(),
            bars: vec![],
            marks: vec![],
        })
    };
    let mut out =
        String::from("<section id=\"overlays\">\n<h2>All samples</h2>\n<div class=\"charts\">\n");
    out += &damage("5′ C→T", "5′", &|i| {
        (samples[i].label.clone(), curves[i].five.clone())
    });
    out += &damage("3′ damage", "3′", &|i| {
        (
            format!("{} ({})", samples[i].label, curves[i].three_name),
            curves[i].three.clone(),
        )
    });
    // Score histograms as lines at bin centres, over the union of the samples' ranges.
    let x = samples
        .iter()
        .map(score_x)
        .fold((f64::MAX, f64::MIN), |a, b| (a.0.min(b.0), a.1.max(b.1)));
    out += &chart(&Chart {
        title: "PMD score",
        x_label: "PMD score",
        y_label: "% of scored reads",
        x,
        x_ticks: x_ticks(x),
        integer_y: false,
        lines: samples
            .iter()
            .enumerate()
            .map(|(i, s)| {
                let h = &s.score_histogram;
                Series {
                    name: s.label.clone(),
                    color: color(i),
                    faint: false,
                    points: score_percent(s)
                        .into_iter()
                        .enumerate()
                        .map(|(j, p)| (h.min + (j as f64 + 0.5) * h.bin_width, p))
                        .collect(),
                }
            })
            .collect(),
        bars: vec![],
        marks: vec![],
    });
    out += "</div>\n<p class=\"note\">Each sample's 3′ curve shows G→A for double-stranded and C→T for single-stranded libraries.</p>\n</section>\n";
    out
}

const STYLE: &str = "
:root { --fg: #1a1a1a; --bg: #ffffff; --muted: #5f6368; --line: #d0d4d9; --head: #f3f4f6; --link: #0b57d0; }
@media (prefers-color-scheme: dark) {
  :root { --fg: #e8eaed; --bg: #16181c; --muted: #a0a4aa; --line: #3c4046; --head: #22252a; --link: #8ab4f8; }
}
a { color: var(--link); }
body { font: 14px/1.45 system-ui, sans-serif; color: var(--fg); background: var(--bg); margin: 0 auto; padding: 16px; max-width: 960px; }
h1 { font-size: 1.5em; } h2 { border-bottom: 1px solid var(--line); padding-bottom: 4px; margin-top: 2em; } h3 { margin-bottom: 0.3em; }
table { border-collapse: collapse; margin: 0 0 12px; }
caption { text-align: left; font-weight: 600; padding: 4px 0; }
th, td { text-align: left; vertical-align: top; padding: 3px 10px 3px 0; border-bottom: 1px solid var(--line); }
td.n, th.n { text-align: right; }
thead th { background: var(--head); padding-left: 4px; }
code { font-size: 0.9em; word-break: break-all; }
.tables { display: flex; flex-wrap: wrap; gap: 0 32px; }
.charts { display: flex; flex-wrap: wrap; gap: 8px 16px; align-items: flex-start; }
figure { margin: 0; max-width: 100%; }
.figure-with-note { max-width: 440px; }
svg { max-width: 100%; height: auto; font-size: 11px; fill: var(--fg); }
svg .axis { stroke: var(--muted); } svg .grid { stroke: var(--line); } svg .mark { stroke: var(--fg); stroke-dasharray: 2 2; }
figcaption { font-size: 0.9em; }
.legend { list-style: none; padding: 0; margin: 4px 0; display: flex; flex-wrap: wrap; gap: 2px 14px; }
.key { display: inline-block; width: 12px; height: 12px; margin-right: 5px; vertical-align: -1px; }
.note, footer { color: var(--muted); font-size: 0.9em; max-width: 880px; }
";

/// The report for one or more samples.
pub fn render(samples: &[Sample]) -> String {
    let mut out = String::from(
        "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n\
         <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n\
         <title>pmdshodh damage report</title>\n<style>",
    );
    out += STYLE;
    out += "</style>\n</head>\n<body>\n<h1>pmdshodh damage report</h1>\n";
    let filtered = samples.iter().any(|s| s.reads.kept.is_some());
    out += "<table>\n<thead><tr><th>Sample</th><th>Library</th><th>Treatment</th><th class=\"n\">Scored reads</th>";
    if filtered {
        out += "<th class=\"n\">Kept</th>";
    }
    out += "<th class=\"n\">First-base C→T</th></tr></thead>\n<tbody>\n";
    for (i, s) in samples.iter().enumerate() {
        let _ = write!(
            out,
            "<tr><td><a href=\"#sample-{}\">{}</a></td><td>{}</td><td>{}</td><td class=\"n\">{}</td>",
            i + 1,
            escape(&s.label),
            library_name(&s.library),
            treatment_name(&s.treatment),
            count(s.reads.scored)
        );
        if filtered {
            let kept = s
                .reads
                .kept
                .and_then(|k| percent(k, s.reads.scored))
                .map_or("n/a".into(), |p| format!("{p:.1}%"));
            let _ = write!(out, "<td class=\"n\">{kept}</td>");
        }
        let first = s
            .damage
            .first
            .0
            .map_or("n/a".into(), |r| format!("{:.2}%", 100.0 * r));
        let _ = writeln!(out, "<td class=\"n\">{first}</td></tr>");
    }
    out += "</tbody>\n</table>\n";
    if samples.len() > 1 {
        out += &overlays(samples);
    }
    for (i, s) in samples.iter().enumerate() {
        out += &sample_section(i, s);
    }
    out += "<footer><p>Written by pmdshodh report from the JSON run summaries, which remain the \
            source of truth. Rates are computed from their counts.</p></footer>\n</body>\n</html>\n";
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    /// A version 2 summary with `range` positions. 5′ ref C: 10 T of 100 at position 0;
    /// 3′ ref G: 5 A of 50, 3′ ref C: 8 T of 40 at position 0; other positions all unchanged.
    fn summary(library: &str, range: usize) -> Value {
        let col = |first: u64, rest: u64| {
            let mut v = vec![rest; range];
            v[0] = first;
            v
        };
        let zeros = vec![0u64; range];
        let five_c = json!({"A": zeros, "C": col(89, 100), "G": col(1, 0), "T": col(10, 0)});
        let five_g = json!({"A": col(2, 0), "C": zeros, "G": col(98, 100), "T": zeros});
        let three_c = json!({"A": zeros, "C": col(32, 40), "G": zeros, "T": col(8, 0)});
        let three_g = json!({"A": col(5, 0), "C": zeros, "G": col(45, 50), "T": zeros});
        let mut hist = vec![0u64; 120];
        hist[40] = 5; // [0, 0.5)
        hist[46] = 4; // [3, 3.5)
        json!({
            "format_version": 2,
            "tool": {"name": "pmdshodh", "version": "0.2.4"},
            "subcommand": "filter",
            "command": "pmdshodh filter -t 3 \"a<b>&c\".bam -o out.bam",
            "input": "C:\\data\\a<b>&c.bam",
            "input_format": "bam",
            "reference": {"path": null, "bases_from": "md", "md_tags": null},
            "mode": {"library": library, "treatment": "none", "pmdtools_compat": false},
            "model": {"pmd_p": 0.3, "pmd_c": 0.01, "pi_ancient": 0.001, "pi_contamination": 0.001, "null_background": 0.001},
            "filters": {"threshold": 3.0, "upper_threshold": null, "min_mapq": 0, "min_baseq": 0, "no_clips": false, "no_indels": false},
            "reads": {
                "input": 12, "scored": 10, "kept": 4, "removed_by_threshold": 6, "no_informative_sites": 1,
                "skipped": {"min_mapq": 0, "unmapped": 2},
                "flags": {"duplicate": 0, "paired": 0, "qc_fail": 0, "reverse": 5, "secondary": 0, "supplementary": 0}
            },
            "score_histogram": {"bin_width": 0.5, "min": -20.0, "max": 40.0, "underflow": 1, "overflow": 0, "counts": hist},
            "length_histogram": {"scored": {"9": 2, "10": 3, "100": 5}, "kept": {"100": 4}},
            "damage": {
                "range": range, "orientation": "molecule", "reads_used": "scored",
                "five_prime": {"C": five_c, "G": five_g},
                "three_prime": {"C": three_c, "G": three_g},
                "first_position": {"c_to_t": 0.1, "se": 0.03, "n": 100}
            }
        })
    }

    /// No NaN or infinity formatted into the page ("inf" followed by a delimiter, not "informative").
    fn all_finite(html: &str) -> bool {
        !html.contains("NaN")
            && !["inf\"", "inf,", "inf ", "inf<", "inf%"]
                .iter()
                .any(|p| html.contains(p))
    }

    fn sample(library: &str, range: usize) -> Sample {
        parse("run.json", &summary(library, range)).unwrap()
    }

    #[test]
    fn rates_from_counts_and_none_without_counts() {
        assert_eq!(percent(1, 4), Some(25.0));
        assert_eq!(percent(0, 7), Some(0.0));
        assert_eq!(percent(0, 0), None);
    }

    #[test]
    fn threshold_marks_inside_and_outside_the_plotted_range() {
        let mut j = summary("ds", 3);
        j["filters"]["threshold"] = json!(50.0);
        j["filters"]["upper_threshold"] = json!(-1.5);
        let html = render(&[parse("t.json", &j).unwrap()]);
        assert!(
            html.contains(r#"data-threshold="-1.5""#),
            "upper threshold drawn"
        );
        assert!(!html.contains(r#"data-threshold="50""#));
        assert!(html.contains("--threshold 50 is outside the plotted range."));
    }

    #[test]
    fn udg_notes_state_which_sites_are_counted() {
        let note = |treatment: &str, compat: bool| {
            let mut j = summary("ds", 3);
            j["mode"]["treatment"] = json!(treatment);
            j["mode"]["pmdtools_compat"] = json!(compat);
            render(&[parse("u.json", &j).unwrap()])
        };
        assert!(note("half", false).contains(
            "UDG-half: only CpG sites, a C at the 5′ terminal base and a G at the 3′ \
             terminal base are counted."
        ));
        assert!(note("full", false).contains("Full UDG: only CpG sites are counted."));
        assert!(note("full", true).contains("CpG context follows PMDtools"));
        let none = note("none", true);
        assert!(!none.contains("only CpG") && !none.contains("CpG context"));
    }

    #[test]
    fn ds_shows_3p_g_to_a_and_ss_shows_3p_c_to_t() {
        let ds = curves(&sample("ds", 3));
        assert_eq!(ds.five, vec![Some(10.0), Some(0.0), Some(0.0)]);
        assert_eq!(ds.three, vec![Some(10.0), Some(0.0), Some(0.0)]); // 5 of 50
        assert_eq!(ds.three_name, "G→A");
        let ss = curves(&sample("ss", 3));
        assert_eq!(ss.five, ds.five);
        assert_eq!(ss.three, vec![Some(20.0), Some(0.0), Some(0.0)]); // 8 of 40
        assert_eq!(ss.three_name, "C→T");
    }

    #[test]
    fn other_substitutions_exclude_the_plotted_one() {
        // 5′: C→G 1 and G→A 2 of 100 + 100 sites; 3′ ds: C→T 8 of 40 + 50 sites.
        let ds = curves(&sample("ds", 2));
        assert_eq!(ds.five_other, vec![Some(1.5), Some(0.0)]);
        assert_eq!(ds.three_other[0], Some(800.0 / 90.0));
        // 3′ ss: G→A 5 of 90.
        assert_eq!(curves(&sample("ss", 2)).three_other[0], Some(500.0 / 90.0));
    }

    #[test]
    fn positions_without_sites_are_gaps_not_zero() {
        let mut j = summary("ds", 2);
        for b in ["A", "C", "G", "T"] {
            j["damage"]["five_prime"]["C"][b][1] = json!(0);
        }
        let c = curves(&parse("run.json", &j).unwrap());
        assert_eq!(c.five, vec![Some(10.0), None]);
    }

    #[test]
    fn series_coordinates_and_gaps() {
        let f = Frame {
            x: (1.0, 3.0),
            y: (0.0, 10.0),
            area: (50.0, 10.0, 200.0, 100.0),
        };
        let svg = series_svg(
            &f,
            &[(1.0, Some(10.0)), (2.0, Some(5.0)), (3.0, None)],
            "#123",
        );
        assert!(
            svg.contains(r#"points="50.00,10.00 150.00,60.00""#),
            "{svg}"
        );
        assert_eq!(svg.matches("<polyline").count(), 1);
        let lone = series_svg(&f, &[(1.0, None), (2.0, Some(0.0)), (3.0, None)], "#123");
        assert!(
            lone.contains(r#"<circle cx="150.00" cy="110.00""#),
            "{lone}"
        );
        assert!(!lone.contains("<polyline"));
    }

    #[test]
    fn a_single_x_value_is_centred_not_nan() {
        let f = Frame {
            x: (1.0, 1.0),
            y: (0.0, 10.0),
            area: (50.0, 10.0, 200.0, 100.0),
        };
        let svg = series_svg(&f, &[(1.0, Some(5.0))], "#123");
        assert!(svg.contains(r#"cx="150.00" cy="60.00""#), "{svg}");
        assert!(all_finite(&render(&[sample("ds", 1)])), "range 1");
    }

    #[test]
    fn escaping() {
        assert_eq!(
            escape(r#"<a href="x">&'"#),
            "&lt;a href=&quot;x&quot;&gt;&amp;&#39;"
        );
        let html = render(&[sample("ds", 3)]);
        assert!(html.contains("a&lt;b&gt;&amp;c.bam"));
        assert!(html.contains("&quot;a&lt;b&gt;&amp;c&quot;.bam -o out.bam"));
        assert!(!html.contains("a<b>"));
    }

    #[test]
    fn labels_are_input_file_names_made_unique() {
        assert_eq!(
            labels(&[
                ("x/one.json", "/data/s1.bam"),
                ("two.json", "C:\\data\\s2.bam")
            ]),
            ["s1.bam", "s2.bam"]
        );
        // stdin input: the summary's file name.
        assert_eq!(labels(&[("x/run.json", "-")]), ["run.json"]);
        // A clash gets the summary's file name, then its path as given.
        assert_eq!(
            labels(&[("a/ds.json", "s.bam"), ("a/ss.json", "t/s.bam")]),
            ["s.bam (ds.json)", "s.bam (ss.json)"]
        );
        assert_eq!(
            labels(&[("a/run.json", "s.bam"), ("b/run.json", "s.bam")]),
            ["s.bam (a/run.json)", "s.bam (b/run.json)"]
        );
    }

    #[test]
    fn length_keys_sort_numerically() {
        assert_eq!(sample("ds", 3).lengths_scored, [(9, 2), (10, 3), (100, 5)]);
    }

    #[test]
    fn version_1_lacks_input_format_and_reference() {
        let mut j = summary("ds", 3);
        j["format_version"] = json!(1);
        let o = j.as_object_mut().unwrap();
        o.remove("input_format");
        o.remove("reference");
        let html = render(&[parse("old.json", &j).unwrap()]);
        assert!(html.contains("not recorded"));
        // Version 2 requires them.
        let mut j = summary("ds", 3);
        j.as_object_mut().unwrap().remove("reference");
        let e = parse("new.json", &j).unwrap_err();
        assert!(e.contains("new.json") && e.contains("reference"), "{e}");
    }

    #[test]
    fn bad_summaries_name_the_file_and_field() {
        let mut j = summary("ds", 3);
        j["format_version"] = json!(3);
        let e = parse("f.json", &j).unwrap_err();
        assert!(
            e.contains("f.json") && e.contains("format_version") && e.contains('3'),
            "{e}"
        );
        let mut j = summary("ds", 3);
        j["reads"]["scored"] = json!("ten");
        let e = parse("f.json", &j).unwrap_err();
        assert!(e.contains("reads.scored"), "{e}");
        let mut j = summary("ds", 3);
        j["damage"]["range"] = json!(4); // tables hold 3 positions
        let e = parse("f.json", &j).unwrap_err();
        assert!(e.contains("damage"), "{e}");
        let mut j = summary("ds", 3);
        j["mode"]["library"] = json!("xs");
        assert!(parse("f.json", &j).unwrap_err().contains("mode.library"));
        let e = parse("f.json", &json!({"a": 1})).unwrap_err();
        assert!(e.contains("format_version"), "{e}");
    }

    #[test]
    fn score_and_damage_summaries_have_no_kept_counts() {
        let mut j = summary("ds", 3);
        j["subcommand"] = json!("damage");
        j["reads"]["kept"] = Value::Null;
        j["reads"]["removed_by_threshold"] = Value::Null;
        j["length_histogram"]["kept"] = Value::Null;
        j["filters"]["threshold"] = Value::Null;
        let html = render(&[parse("d.json", &j).unwrap()]);
        assert!(!html.contains("Kept"), "no kept row for damage");
    }

    #[test]
    fn filter_report_states_the_kept_fraction_and_thresholds() {
        let html = render(&[sample("ds", 3)]);
        assert!(html.contains("4 of 10 scored reads (40.0%)"), "{html}");
        assert!(html.contains(r#"data-threshold="3""#));
    }

    #[test]
    fn empty_run_renders_without_nan() {
        let mut j = summary("ds", 3);
        j["reads"] = json!({"input": 0, "scored": 0, "kept": 0, "removed_by_threshold": 0,
            "no_informative_sites": 0, "skipped": {}, "flags": {}});
        j["score_histogram"]["counts"] = json!(vec![0; 120]);
        j["score_histogram"]["underflow"] = json!(0);
        j["length_histogram"] = json!({"scored": {}, "kept": {}});
        for end in ["five_prime", "three_prime"] {
            for r in ["C", "G"] {
                for b in ["A", "C", "G", "T"] {
                    j["damage"][end][r][b] = json!([0, 0, 0]);
                }
            }
        }
        j["damage"]["first_position"] = json!({"c_to_t": null, "se": null, "n": 0});
        let html = render(&[parse("e.json", &j).unwrap()]);
        assert!(all_finite(&html));
        assert!(html.contains("n/a"));
    }

    #[test]
    fn self_contained_and_deterministic() {
        let samples = [sample("ds", 30), sample("ss", 25)];
        let html = render(&samples);
        assert!(html.starts_with("<!DOCTYPE html>\n"));
        assert!(!html.contains("<script") && !html.contains("http") && !html.contains('\r'));
        assert_eq!(html, render(&samples));
    }

    #[test]
    fn overlays_only_with_several_samples() {
        let one = render(&[sample("ds", 30)]);
        assert!(!one.contains("id=\"overlays\""));
        let mut a = sample("ds", 30);
        a.label = "A".into();
        let mut b = sample("ss", 25);
        b.label = "B".into();
        let both = render(&[a, b]);
        assert!(both.contains("id=\"overlays\""));
        // Each 3′ curve is labelled with what it shows.
        assert!(
            both.contains("A (G→A)") && both.contains("B (C→T)"),
            "{both}"
        );
    }

    #[test]
    fn plotted_values_are_machine_readable() {
        let html = render(&[sample("ds", 3)]);
        assert!(
            html.contains(r#"data-series="5′ C→T" data-values="10.000 0.000 0.000""#),
            "{html}"
        );
    }
}
