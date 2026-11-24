//! What a run prints to stderr: the end-of-run summary block and, in a terminal, a live progress
//! line. Nothing here touches stdout, the JSON summary or the exit status.

use std::io::{IsTerminal, Read, Write};
use std::sync::OnceLock;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{Duration, Instant};

/// Bytes read so far from the input file, for the progress line's percentage.
static BYTES_READ: AtomicU64 = AtomicU64::new(0);
/// Size of the input file, if it is a regular file.
static INPUT_SIZE: OnceLock<u64> = OnceLock::new();

/// Counts the bytes read through it into `BYTES_READ`.
pub struct Counting<R>(pub R);

impl<R: Read> Read for Counting<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let n = self.0.read(buf)?;
        BYTES_READ.fetch_add(n as u64, Ordering::Relaxed);
        Ok(n)
    }
}

/// Record the input file's size, so progress can show how much of it has been read.
pub fn set_input_size(bytes: u64) {
    let _ = INPUT_SIZE.set(bytes);
}

/// Colour on stderr: only in a terminal that understands ANSI codes, and never with `NO_COLOR`.
#[derive(Clone, Copy)]
pub struct Style {
    color: bool,
}

impl Style {
    pub fn detect() -> Style {
        let no_color = std::env::var_os("NO_COLOR").is_some_and(|v| !v.is_empty());
        let dumb = std::env::var_os("TERM").is_some_and(|t| t == "dumb");
        // Windows consoles show raw escape codes unless the terminal enables them; Windows
        // Terminal does (WT_SESSION), as does any terminal that sets TERM.
        let ansi = cfg!(unix)
            || std::env::var_os("WT_SESSION").is_some()
            || std::env::var_os("TERM").is_some();
        Style {
            color: std::io::stderr().is_terminal() && ansi && !no_color && !dumb,
        }
    }

    #[cfg(test)]
    pub fn plain() -> Style {
        Style { color: false }
    }

    fn paint(self, code: &str, s: &str) -> String {
        if self.color {
            format!("\x1b[{code}m{s}\x1b[0m")
        } else {
            s.to_string()
        }
    }

    fn bold(self, s: &str) -> String {
        self.paint("1", s)
    }
    fn label(self, s: &str) -> String {
        self.paint("36", s)
    }
    fn good(self, s: &str) -> String {
        self.paint("32", s)
    }
    fn warn(self, s: &str) -> String {
        self.paint("33", s)
    }
}

/// An integer with thousands separators: 1234567 → "1,234,567".
pub fn grouped(n: u64) -> String {
    let digits = n.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

/// Elapsed time: "0.4 s", "12.3 s", "4:05", "1:02:03".
pub fn duration(d: Duration) -> String {
    let s = d.as_secs_f64();
    if s < 1.0 {
        return format!("{s:.2} s");
    }
    if s < 60.0 {
        return format!("{s:.1} s");
    }
    let t = d.as_secs();
    let (h, m, sec) = (t / 3600, t / 60 % 60, t % 60);
    if h > 0 {
        format!("{h}:{m:02}:{sec:02}")
    } else {
        format!("{m}:{sec:02}")
    }
}

/// Reads per second: "15,000 reads/s" or "1.92 M reads/s".
pub fn rate(reads: u64, d: Duration) -> String {
    let per_s = reads as f64 / d.as_secs_f64().max(1e-9);
    if per_s >= 1e6 {
        format!("{:.2} M reads/s", per_s / 1e6)
    } else {
        format!("{} reads/s", grouped(per_s.round() as u64))
    }
}

/// A PMD bound as given, without a trailing ".0": 3 → "3", 2.5 → "2.5".
fn bound(x: f64) -> String {
    let s = format!("{x}");
    s.strip_suffix(".0").map_or(s.clone(), str::to_string)
}

/// Everything the end-of-run summary shows.
pub struct RunReport {
    pub subcommand: String,
    pub input: String,
    pub format: &'static str,
    pub reference: Option<String>,
    pub library: &'static str,
    pub treatment: &'static str,
    pub compat: bool,
    pub threshold: Option<f64>,
    pub upper_threshold: Option<f64>,
    pub reads: u64,
    pub scored: u64,
    pub no_sites: u64,
    pub kept: Option<u64>,
    /// Skip reasons with their counts, in the summary's order; zero counts are left out here.
    pub skipped: Vec<(&'static str, u64)>,
    /// First-base 5′ C→T rate, its standard error, and the number of sites, if damage was counted.
    pub first_base: Option<(f64, f64, u64)>,
    pub md_disagreeing: u64,
    pub outputs: Vec<String>,
    pub elapsed: Duration,
}

impl RunReport {
    /// The summary block, one line per topic, labels aligned.
    pub fn render(&self, style: Style) -> String {
        let treatment = match self.treatment {
            "half" => "UDG-half",
            "full" => "full UDG",
            _ => "no UDG",
        };
        let mut head = format!(
            "pmdshodh {} {} · {} · {treatment}",
            self.subcommand,
            env!("CARGO_PKG_VERSION"),
            self.library
        );
        let keep = match (self.threshold, self.upper_threshold) {
            (Some(t), Some(u)) => Some(format!("keep {} ≤ PMD < {}", bound(t), bound(u))),
            (Some(t), None) => Some(format!("keep PMD ≥ {}", bound(t))),
            (None, Some(u)) => Some(format!("keep PMD < {}", bound(u))),
            (None, None) => None,
        };
        if let Some(k) = keep {
            head += &format!(" · {k}");
        }
        if self.compat {
            head += " · --pmdtools-compat";
        }

        let mut lines = vec![style.bold(&head)];
        let mut row = |label: &str, text: String| {
            lines.push(format!("  {} {text}", style.label(&format!("{label:<9}"))));
        };
        let bases = match &self.reference {
            Some(r) => format!("reference bases from {r}"),
            None => "reference bases from MD tags".to_string(),
        };
        row(
            "input",
            format!("{} ({}) · {bases}", self.input, self.format.to_uppercase()),
        );
        let skipped: u64 = self.skipped.iter().map(|(_, n)| n).sum();
        let mut reads = format!(
            "{} read · {} scored · {} skipped",
            grouped(self.reads),
            grouped(self.scored),
            grouped(skipped)
        );
        if self.no_sites > 0 {
            reads += &format!(
                " · {} scored without informative sites",
                grouped(self.no_sites)
            );
        }
        row("reads", reads);
        if skipped > 0 {
            let why: Vec<String> = self
                .skipped
                .iter()
                .map(|(r, n)| format!("{r} {}", grouped(*n)))
                .collect();
            row("skipped", style.warn(&why.join(" · ")));
        }
        if let Some(kept) = self.kept {
            let share = if self.scored > 0 {
                format!(" ({:.1}%)", 100.0 * kept as f64 / self.scored as f64)
            } else {
                String::new()
            };
            row(
                "kept",
                style.good(&format!(
                    "{} of {} scored{share}",
                    grouped(kept),
                    grouped(self.scored)
                )),
            );
        }
        if let Some((rate, se, n)) = self.first_base {
            row(
                "damage",
                format!(
                    "5′ C→T at the first base {:.1}% ± {:.*}% ({} site{})",
                    100.0 * rate,
                    // Large libraries have standard errors well under 0.05%.
                    if 100.0 * se < 0.05 { 2 } else { 1 },
                    100.0 * se,
                    grouped(n),
                    if n == 1 { "" } else { "s" }
                ),
            );
        }
        if self.md_disagreeing > 0 {
            row(
                "note",
                style.warn(&format!(
                    "{} scored reads have an MD tag that disagrees with the reference (written \
                     before realignment, or against N bases that the aligner replaced, as bwa \
                     does); they were scored from the reference",
                    grouped(self.md_disagreeing)
                )),
            );
        }
        if !self.outputs.is_empty() {
            row("written", self.outputs.join(" · "));
        }
        row(
            "time",
            format!(
                "{} · {}",
                duration(self.elapsed),
                rate(self.reads, self.elapsed)
            ),
        );
        lines.join("\n")
    }
}

/// A one-line progress display on stderr, in a terminal only: reads so far, rate, share of the
/// input file read, and elapsed time. It appears after a second, so short runs never show it,
/// redraws a few times a second, and is erased before anything else is printed.
pub struct Progress {
    label: &'static str,
    enabled: bool,
    start: Instant,
    last: Instant,
    width: usize,
}

impl Progress {
    const DELAY: Duration = Duration::from_secs(1);
    const EVERY: Duration = Duration::from_millis(250);

    pub fn new(label: &'static str, quiet: bool) -> Progress {
        let now = Instant::now();
        Progress {
            label,
            enabled: !quiet && std::io::stderr().is_terminal(),
            start: now,
            last: now,
            width: 0,
        }
    }

    pub fn elapsed(&self) -> Duration {
        self.start.elapsed()
    }

    /// Called after each batch of reads; cheap when there is nothing to draw.
    pub fn tick(&mut self, reads: u64) {
        if !self.enabled {
            return;
        }
        let now = Instant::now();
        if now.duration_since(self.start) < Self::DELAY
            || now.duration_since(self.last) < Self::EVERY
        {
            return;
        }
        self.last = now;
        let elapsed = now.duration_since(self.start);
        let mut text = format!(
            "  pmdshodh {} · {} reads · {}",
            self.label,
            grouped(reads),
            rate(reads, elapsed)
        );
        if let Some(&size) = INPUT_SIZE.get().filter(|&&s| s > 0) {
            let done = BYTES_READ.load(Ordering::Relaxed).min(size);
            text += &format!(" · {:.0}% of input", 100.0 * done as f64 / size as f64);
        }
        text += &format!(" · {}", duration(elapsed));
        let pad = self.width.saturating_sub(text.chars().count());
        self.width = text.chars().count();
        let mut err = std::io::stderr().lock();
        let _ = write!(err, "\r{text}{}", " ".repeat(pad));
        let _ = err.flush();
    }

    /// Erase the progress line, if one was drawn.
    pub fn clear(&mut self) {
        if self.width > 0 {
            let mut err = std::io::stderr().lock();
            let _ = write!(err, "\r{}\r", " ".repeat(self.width));
            let _ = err.flush();
            self.width = 0;
        }
    }
}

impl Drop for Progress {
    fn drop(&mut self) {
        self.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numbers_and_times_are_formatted_for_people() {
        assert_eq!(grouped(0), "0");
        assert_eq!(grouped(999), "999");
        assert_eq!(grouped(1_000), "1,000");
        assert_eq!(grouped(14_006_183), "14,006,183");
        assert_eq!(duration(Duration::from_millis(400)), "0.40 s");
        assert_eq!(duration(Duration::from_millis(12_340)), "12.3 s");
        assert_eq!(duration(Duration::from_secs(245)), "4:05");
        assert_eq!(duration(Duration::from_secs(3723)), "1:02:03");
        assert_eq!(rate(15_000, Duration::from_secs(1)), "15,000 reads/s");
        assert_eq!(rate(19_200_000, Duration::from_secs(10)), "1.92 M reads/s");
        assert_eq!(bound(3.0), "3");
        assert_eq!(bound(2.5), "2.5");
        assert_eq!(bound(-1000.0), "-1000");
    }

    fn report() -> RunReport {
        RunReport {
            subcommand: "filter".into(),
            input: "in.bam".into(),
            format: "bam",
            reference: None,
            library: "ds",
            treatment: "none",
            compat: false,
            threshold: Some(3.0),
            upper_threshold: None,
            reads: 1_500,
            scored: 1_500,
            no_sites: 0,
            kept: Some(524),
            skipped: vec![],
            first_base: Some((0.398, 0.025, 372)),
            md_disagreeing: 0,
            outputs: vec!["out.bam".into()],
            elapsed: Duration::from_millis(100),
        }
    }

    #[test]
    fn summary_block_for_a_filter_run() {
        let text = report().render(Style::plain());
        let expected = format!(
            "pmdshodh filter {} · ds · no UDG · keep PMD ≥ 3
  input     in.bam (BAM) · reference bases from MD tags
  reads     1,500 read · 1,500 scored · 0 skipped
  kept      524 of 1,500 scored (34.9%)
  damage    5′ C→T at the first base 39.8% ± 2.5% (372 sites)
  written   out.bam
  time      0.10 s · 15,000 reads/s",
            env!("CARGO_PKG_VERSION")
        );
        assert_eq!(text, expected);
    }

    #[test]
    fn summary_block_lists_skips_bounds_and_notes_only_when_present() {
        let mut r = report();
        r.subcommand = "score".into();
        r.kept = None;
        r.threshold = None;
        r.first_base = None;
        r.treatment = "half";
        r.library = "ss";
        r.compat = true;
        r.reference = Some("ref.fa".into());
        r.skipped = vec![("unmapped", 3), ("min_mapq", 1_200)];
        r.no_sites = 7;
        r.md_disagreeing = 2;
        r.first_base = Some((0.184, 0.00022, 3_097_363));
        assert!(
            r.render(Style::plain())
                .contains("5′ C→T at the first base 18.4% ± 0.02% (3,097,363 sites)")
        );
        r.first_base = Some((1.0, 0.0, 1));
        assert!(r.render(Style::plain()).contains("(1 site)"));
        r.first_base = None;
        let text = r.render(Style::plain());
        assert!(text.starts_with(&format!(
            "pmdshodh score {} · ss · UDG-half · --pmdtools-compat\n",
            env!("CARGO_PKG_VERSION")
        )));
        assert!(text.contains("in.bam (BAM) · reference bases from ref.fa"));
        assert!(text.contains(
            "reads     1,500 read · 1,500 scored · 1,203 skipped · 7 scored without informative sites"
        ));
        assert!(text.contains("skipped   unmapped 3 · min_mapq 1,200"));
        assert!(text.contains("note      2 scored reads have an MD tag that disagrees"));
        assert!(!text.contains("kept") && !text.contains("damage"));

        r.threshold = Some(3.0);
        r.upper_threshold = Some(10.0);
        assert!(r.render(Style::plain()).contains("keep 3 ≤ PMD < 10"));
        r.threshold = None;
        assert!(r.render(Style::plain()).contains("keep PMD < 10"));
    }

    #[test]
    fn colour_wraps_text_only_when_enabled() {
        let on = Style { color: true };
        assert_eq!(on.good("x"), "\x1b[32mx\x1b[0m");
        assert_eq!(Style::plain().good("x"), "x");
    }
}
