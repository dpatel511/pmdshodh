use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::{ArgGroup, Args, Parser, Subcommand, ValueEnum};
use noodles_util::alignment;

use noodles::sam::alignment::Record as _;
use pmdshodh::damage::DamageCounts;
use pmdshodh::model::Params;
use pmdshodh::reads::{Outcome, ReadFilters, SkipReason, process};
use pmdshodh::sites::ReadScore;
use pmdshodh::sites::{Library, Settings, Treatment};
use pmdshodh::summary::{RunInfo, Summary};

/// Per-read post-mortem damage (PMD) scores for ancient DNA alignments.
#[derive(Parser)]
#[command(name = "pmdshodh", version)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Keep reads by PMD score and write them as SAM or BAM.
    Filter(FilterArgs),
    /// Write a table of per-read PMD scores (TSV).
    Score(ScoreArgs),
    /// Write the JSON run summary only: read counts, score histogram and damage table.
    Damage(DamageArgs),
}

#[derive(Args)]
struct DamageArgs {
    /// Input SAM or BAM file ("-" for stdin).
    input: PathBuf,
    /// Output JSON file (default: stdout).
    #[arg(short, long)]
    output: Option<PathBuf>,
    /// Positions from each read end in the damage table (PMDtools --range).
    #[arg(long, default_value_t = 30)]
    range: usize,
    /// Stop after this many input reads (PMDtools -n). Exits normally.
    #[arg(short = 'n', long, visible_alias = "number")]
    max_reads: Option<u64>,
    #[command(flatten)]
    model: ModelArgs,
    /// Do not print the run summary to stderr.
    #[arg(long)]
    quiet: bool,
}

#[derive(Args)]
#[command(group(ArgGroup::new("bounds").args(["threshold", "upper_threshold"]).required(true).multiple(true)))]
struct FilterArgs {
    /// Input SAM or BAM file ("-" for stdin).
    input: PathBuf,
    /// Output file; format from the extension (.bam or .sam). "-" writes BAM to stdout.
    #[arg(short, long)]
    output: PathBuf,
    /// Output format for stdout ("-o -").
    #[arg(long, value_enum, default_value = "bam")]
    output_format: OutputFormat,
    /// Keep reads with PMD score >= this (PMDtools -t/--threshold).
    #[arg(short = 't', long, allow_negative_numbers = true)]
    threshold: Option<f64>,
    /// Keep reads with PMD score < this (PMDtools --upperthreshold).
    #[arg(long, visible_alias = "upperthreshold", allow_negative_numbers = true)]
    upper_threshold: Option<f64>,
    /// Add the PMD score to kept reads as a tag (overwrites an existing tag of that name).
    #[arg(long)]
    tag: bool,
    /// Tag name: "pm" (pm:f, 32-bit float) or "DS" (DS:Z, 3 decimals, as PMDtools --writesamfield).
    #[arg(long, value_enum, default_value = "pm", requires = "tag")]
    tag_name: TagName,
    /// Accepted for PMDtools compatibility and ignored: a header is always written.
    #[arg(long, hide = true)]
    header: bool,
    #[command(flatten)]
    summary: SummaryArgs,
    #[command(flatten)]
    model: ModelArgs,
    /// Do not print the run summary to stderr.
    #[arg(long)]
    quiet: bool,
}

#[derive(Clone, Copy, PartialEq, Eq, ValueEnum)]
enum OutputFormat {
    Bam,
    Sam,
}

#[derive(Clone, Copy, PartialEq, Eq, ValueEnum)]
enum TagName {
    #[value(name = "pm")]
    Pm,
    #[value(name = "DS")]
    Ds,
}

#[derive(Args)]
struct ScoreArgs {
    /// Input SAM or BAM file ("-" for stdin).
    input: PathBuf,
    /// Output TSV file (default: stdout).
    #[arg(short, long)]
    output: Option<PathBuf>,
    #[command(flatten)]
    summary: SummaryArgs,
    #[command(flatten)]
    model: ModelArgs,
    /// Do not print the run summary to stderr.
    #[arg(long)]
    quiet: bool,
}

#[derive(Clone, Copy, ValueEnum)]
enum LibraryArg {
    Ds,
    Ss,
}

#[derive(Clone, Copy, ValueEnum)]
enum TreatmentArg {
    None,
    Half,
    Full,
}

/// Scoring options shared by all subcommands. PMDtools option names are accepted as aliases.
#[derive(Args)]
#[command(group(ArgGroup::new("library_choice").args(["library", "ss"])))]
#[command(group(ArgGroup::new("treatment_choice").args(["treatment", "udg_minus", "udg_half", "udg_plus"])))]
struct ModelArgs {
    /// Library type: double-stranded (ds) or single-stranded (ss). Never guessed from the data.
    #[arg(long, value_enum, default_value = "ds")]
    library: LibraryArg,
    /// Same as --library ss (PMDtools option).
    #[arg(long)]
    ss: bool,
    /// UDG treatment: none (all C/G sites), half (CpG sites plus the terminal base), full (CpG only).
    #[arg(long, value_enum, default_value = "none")]
    treatment: TreatmentArg,
    /// Same as --treatment none (PMDtools option).
    #[arg(long = "UDGminus")]
    udg_minus: bool,
    /// Same as --treatment half (PMDtools option).
    #[arg(long = "UDGhalf")]
    udg_half: bool,
    /// Same as --treatment full (PMDtools --UDGplus / --CpG).
    #[arg(long = "UDGplus", visible_alias = "CpG")]
    udg_plus: bool,
    /// Geometric decay of damage away from the read end (PMDtools --PMDpparam).
    #[arg(long, visible_alias = "PMDpparam", default_value_t = 0.3)]
    pmd_p: f64,
    /// Damage added at every position (PMDtools --PMDconstant).
    #[arg(long, visible_alias = "PMDconstant", default_value_t = 0.01)]
    pmd_c: f64,
    /// Polymorphism rate in the ancient model (PMDtools --polymorphism_ancient).
    #[arg(long, visible_alias = "polymorphism_ancient", default_value_t = 0.001)]
    pi_ancient: f64,
    /// Polymorphism rate in the null model (PMDtools --polymorphism_contamination).
    #[arg(
        long,
        visible_alias = "polymorphism_contamination",
        default_value_t = 0.001
    )]
    pi_contamination: f64,
    /// Skip reads with mapping quality below this (PMDtools -m).
    #[arg(short = 'm', long, visible_alias = "requiremapq", default_value_t = 0)]
    min_mapq: u8,
    /// Ignore bases with base quality below this (PMDtools -q).
    #[arg(short = 'q', long, visible_alias = "requirebaseq", default_value_t = 0)]
    min_baseq: u8,
    /// Skip reads with soft or hard clips (PMDtools --noclips).
    #[arg(long, visible_alias = "noclips")]
    no_clips: bool,
    /// Skip reads with insertions or deletions (PMDtools --noindels).
    #[arg(long, visible_alias = "noindels")]
    no_indels: bool,
    /// Reproduce PMDtools 0.60 exactly, including its bugs: the reference rebuilt wrongly for reads
    /// combining an insertion or soft clip with a deletion or '=' CIGAR operations, CpG context
    /// wrapping around the read ends, and PMDtools' CpG rule for damage tables. For comparisons.
    #[arg(long)]
    pmdtools_compat: bool,
    /// Worker threads for scoring and BGZF (de)compression. Output is identical for any value.
    #[arg(short = '@', long, default_value_t = 1)]
    threads: usize,
}

impl ModelArgs {
    fn settings(&self) -> Result<(Settings, ReadFilters), String> {
        for (name, v) in [
            ("--pmd-p", self.pmd_p),
            ("--pmd-c", self.pmd_c),
            ("--pi-ancient", self.pi_ancient),
            ("--pi-contamination", self.pi_contamination),
        ] {
            if !(0.0..=1.0).contains(&v) {
                return Err(format!("{name} must be between 0 and 1, got {v}"));
            }
        }
        let library = if self.ss {
            Library::SingleStranded
        } else {
            match self.library {
                LibraryArg::Ds => Library::DoubleStranded,
                LibraryArg::Ss => Library::SingleStranded,
            }
        };
        let treatment = if self.udg_half {
            Treatment::Half
        } else if self.udg_plus {
            Treatment::Full
        } else if self.udg_minus {
            Treatment::None
        } else {
            match self.treatment {
                TreatmentArg::None => Treatment::None,
                TreatmentArg::Half => Treatment::Half,
                TreatmentArg::Full => Treatment::Full,
            }
        };
        let settings = Settings {
            library,
            treatment,
            params: Params {
                p: self.pmd_p,
                c: self.pmd_c,
                pi_ancient: self.pi_ancient,
                pi_contamination: self.pi_contamination,
            },
            min_baseq: self.min_baseq,
            pmdtools_compat: self.pmdtools_compat,
        };
        let filters = ReadFilters {
            min_mapq: self.min_mapq,
            no_clips: self.no_clips,
            no_indels: self.no_indels,
        };
        Ok((settings, filters))
    }
}

/// Exit status: 1 for runtime failures, 2 for invalid options or input (docs/design.md 4.5).
enum Failure {
    Input(String),
    Runtime(String),
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let result = match cli.command {
        Command::Filter(args) => run_filter(&args),
        Command::Score(args) => run_score(&args),
        Command::Damage(args) => run_damage(&args),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(Failure::Input(msg)) => {
            eprintln!("pmdshodh: error: {msg}");
            ExitCode::from(2)
        }
        Err(Failure::Runtime(msg)) => {
            eprintln!("pmdshodh: error: {msg}");
            ExitCode::from(1)
        }
    }
}

fn open_input(path: &PathBuf) -> Result<Box<dyn Read + Send>, Failure> {
    if path.as_os_str() == "-" {
        Ok(Box::new(io::stdin()))
    } else {
        File::open(path)
            .map(|f| Box::new(f) as Box<dyn Read + Send>)
            .map_err(|e| Failure::Input(format!("cannot open {}: {e}", path.display())))
    }
}

/// Options for the JSON run summary.
#[derive(Args)]
struct SummaryArgs {
    /// Also write a JSON run summary (counts, score histogram, damage table) to this file.
    #[arg(long)]
    summary: Option<PathBuf>,
    /// Positions from each read end in the damage table (PMDtools --range).
    #[arg(long, default_value_t = 30)]
    range: usize,
}

type InputReader = alignment::io::Reader<Box<dyn Read + Send>>;
type AlignmentRecord = alignment::Record;

fn open_reader(
    path: &Path,
    threads: usize,
) -> Result<(InputReader, noodles::sam::Header), Failure> {
    use std::io::BufRead;
    let cannot_read = |e: io::Error| Failure::Input(format!("cannot read {}: {e}", path.display()));
    let mut input = io::BufReader::new(open_input(&path.to_path_buf())?);
    let bgzf = input
        .fill_buf()
        .map_err(cannot_read)?
        .starts_with(&[0x1f, 0x8b]);
    let builder = alignment::io::reader::Builder::default();
    let source: Box<dyn Read + Send> = if bgzf {
        let workers = std::num::NonZero::new(threads).unwrap_or(std::num::NonZero::<usize>::MIN);
        Box::new(noodles::bgzf::io::MultithreadedReader::with_worker_count(
            workers, input,
        ))
    } else {
        Box::new(input)
    };
    let builder = if bgzf {
        builder.set_compression_method(None)
    } else {
        builder
    };
    let mut reader = builder.build_from_reader(source).map_err(cannot_read)?;
    let header = reader
        .read_header()
        .map_err(|e| Failure::Input(format!("invalid header: {e}")))?;
    Ok((reader, header))
}

/// Read every record (or the first `max_reads`), count it in `summary`, and pass scored reads to
/// `on_scored`, which returns whether the read was kept (`None` when there is no filtering).
/// Records read and scored per batch; bounds memory independently of the input size.
const BATCH: usize = 10_000;

/// Read every record (or the first `max_reads`), count it in `summary`, and pass scored reads to
/// `on_scored` in input order; it returns whether the read was kept (`None` without filtering).
/// Scoring runs on the rayon thread pool; results and summary are independent of the thread count.
fn scan(
    reader: &mut InputReader,
    header: &noodles::sam::Header,
    settings: &Settings,
    filters: &ReadFilters,
    summary: &mut Summary,
    max_reads: Option<u64>,
    mut on_scored: impl FnMut(&AlignmentRecord, &ReadScore) -> Result<Option<bool>, Failure>,
) -> Result<(), Failure> {
    use rayon::prelude::*;
    let range = summary.damage.range;
    let mut n_read = 0u64;
    loop {
        let mut batch: Vec<AlignmentRecord> = Vec::with_capacity(BATCH);
        while batch.len() < BATCH && max_reads.is_none_or(|m| n_read < m) {
            let mut record = AlignmentRecord::default();
            let n = reader
                .read_record(header, &mut record)
                .map_err(|e| Failure::Input(format!("record {}: {e}", n_read + 1)))?;
            if n == 0 {
                break;
            }
            n_read += 1;
            batch.push(record);
        }
        if batch.is_empty() {
            return Ok(());
        }
        let chunk = batch.len().div_ceil(rayon::current_num_threads()).max(1);
        let parts: Vec<(Vec<Result<Outcome, String>>, DamageCounts)> = batch
            .par_chunks(chunk)
            .map(|records| {
                let mut damage = DamageCounts::new(range);
                let outcomes = records
                    .iter()
                    .map(|r| process(r, settings, filters, Some(&mut damage)))
                    .collect();
                (outcomes, damage)
            })
            .collect();
        let first = n_read - batch.len() as u64;
        let outcomes = parts.iter().flat_map(|(o, _)| o.iter());
        for (k, (record, outcome)) in batch.iter().zip(outcomes).enumerate() {
            let flags = record.flags().map_err(|e| {
                Failure::Input(format!(
                    "record {}: invalid FLAG: {e}",
                    first + k as u64 + 1
                ))
            })?;
            summary.add_input(flags);
            match outcome {
                Ok(Outcome::Scored { score, length }) => {
                    let kept = on_scored(record, score)?;
                    summary.add_scored(score, *length, kept);
                }
                Ok(Outcome::Skipped(reason)) => summary.add_skipped(*reason),
                Err(msg) => {
                    let name = record.name().map_or_else(|| "*".into(), |n| n.to_string());
                    return Err(Failure::Input(format!("read {name}: {msg}")));
                }
            }
        }
        for (_, damage) in &parts {
            summary.damage.merge(damage);
        }
    }
}

/// Use `threads` worker threads for scoring (and BGZF, where used).
fn init_threads(threads: usize) -> Result<(), Failure> {
    if threads == 0 {
        return Err(Failure::Input("--threads must be at least 1".into()));
    }
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .map_err(|e| Failure::Runtime(format!("cannot start threads: {e}")))
}

fn command_line() -> String {
    std::env::args().collect::<Vec<_>>().join(" ")
}

/// Write JSON to a file (atomically) or to stdout (`None` or "-").
fn write_json(path: Option<&PathBuf>, value: &serde_json::Value) -> Result<(), Failure> {
    let text = serde_json::to_string_pretty(value)
        .map_err(|e| Failure::Runtime(format!("cannot format summary: {e}")))?
        + "\n";
    match path {
        Some(p) if p.as_os_str() != "-" => {
            let (atomic, mut f) = AtomicOutput::create(p)?;
            f.write_all(text.as_bytes())
                .map_err(|e| Failure::Runtime(format!("cannot write {}: {e}", p.display())))?;
            drop(f);
            atomic.commit()
        }
        _ => io::stdout()
            .write_all(text.as_bytes())
            .map_err(|e| Failure::Runtime(format!("cannot write output: {e}"))),
    }
}

/// Facts shared by the end of every subcommand.
struct Finish<'a> {
    subcommand: &'a str,
    input: &'a Path,
    settings: Settings,
    filters: ReadFilters,
    threshold: Option<f64>,
    upper_threshold: Option<f64>,
}

impl Finish<'_> {
    fn json(&self, summary: &Summary) -> serde_json::Value {
        summary.to_json(&RunInfo {
            subcommand: self.subcommand,
            command: command_line(),
            input: self.input.display().to_string(),
            settings: self.settings,
            filters: self.filters,
            threshold: self.threshold,
            upper_threshold: self.upper_threshold,
        })
    }

    /// One summary line on stderr, unless quiet.
    fn report(&self, json: &serde_json::Value, summary: &Summary, quiet: bool) {
        if quiet {
            return;
        }
        let r = &json["reads"];
        let kept = if r["kept"].is_null() {
            String::new()
        } else {
            format!(", {} kept", r["kept"])
        };
        eprintln!(
            "pmdshodh {}: {} reads read, {} scored{kept}, skipped: {}",
            self.subcommand,
            r["input"],
            r["scored"],
            skip_summary(summary.skipped())
        );
    }
}

fn run_score(args: &ScoreArgs) -> Result<(), Failure> {
    let (settings, filters) = args.model.settings().map_err(Failure::Input)?;
    init_threads(args.model.threads)?;
    let (mut reader, header) = open_reader(&args.input, args.model.threads)?;
    let out: Box<dyn Write> = match &args.output {
        Some(path) => Box::new(
            File::create(path)
                .map_err(|e| Failure::Runtime(format!("cannot create {}: {e}", path.display())))?,
        ),
        None => Box::new(io::stdout().lock()),
    };
    let mut out = BufWriter::new(out);
    let write_err = |e: io::Error| Failure::Runtime(format!("cannot write output: {e}"));

    writeln!(out, "read_name\tflag\tpmds\tn_sites").map_err(write_err)?;
    let mut summary = Summary::new(args.summary.range);
    scan(
        &mut reader,
        &header,
        &settings,
        &filters,
        &mut summary,
        None,
        |record, s| {
            let name = record
                .name()
                .map_or_else(|| "*".to_string(), |n| n.to_string());
            let flag = record.flags().map(u16::from).unwrap_or_default();
            writeln!(out, "{name}\t{flag}\t{}\t{}", s.pmds, s.n_sites).map_err(write_err)?;
            Ok(None)
        },
    )?;
    out.flush().map_err(write_err)?;
    let fin = Finish {
        subcommand: "score",
        input: &args.input,
        settings,
        filters,
        threshold: None,
        upper_threshold: None,
    };
    let json = fin.json(&summary);
    if let Some(p) = &args.summary.summary {
        write_json(Some(p), &json)?;
    }
    fin.report(&json, &summary, args.quiet);
    Ok(())
}

fn run_damage(args: &DamageArgs) -> Result<(), Failure> {
    let (settings, filters) = args.model.settings().map_err(Failure::Input)?;
    init_threads(args.model.threads)?;
    let (mut reader, header) = open_reader(&args.input, args.model.threads)?;
    let mut summary = Summary::new(args.range);
    scan(
        &mut reader,
        &header,
        &settings,
        &filters,
        &mut summary,
        args.max_reads,
        |_, _| Ok(None),
    )?;
    let fin = Finish {
        subcommand: "damage",
        input: &args.input,
        settings,
        filters,
        threshold: None,
        upper_threshold: None,
    };
    let json = fin.json(&summary);
    write_json(args.output.as_ref(), &json)?;
    fin.report(&json, &summary, args.quiet);
    Ok(())
}

/// Format a score as PMDtools --writesamfield does: Python's str(round(x, 3)).
fn pmdtools_ds_text(x: f64) -> String {
    let mut s = format!("{x:.3}");
    while s.ends_with('0') && !s.ends_with(".0") {
        s.pop();
    }
    s
}

fn skip_summary(skipped: &[u64]) -> String {
    let reasons: Vec<String> = SkipReason::ALL
        .iter()
        .zip(skipped)
        .filter(|(_, n)| **n > 0)
        .map(|(r, n)| format!("{}={n}", r.name()))
        .collect();
    if reasons.is_empty() {
        "none".to_string()
    } else {
        reasons.join(", ")
    }
}

/// Output file that only appears under its final name if the run succeeds.
struct AtomicOutput {
    tmp: PathBuf,
    dest: PathBuf,
}

impl AtomicOutput {
    fn create(dest: &Path) -> Result<(Self, File), Failure> {
        let mut name = dest.file_name().unwrap_or_default().to_os_string();
        name.push(".pmdshodh-tmp");
        let tmp = dest.with_file_name(name);
        let file = File::create(&tmp)
            .map_err(|e| Failure::Runtime(format!("cannot create {}: {e}", tmp.display())))?;
        Ok((
            AtomicOutput {
                tmp,
                dest: dest.to_path_buf(),
            },
            file,
        ))
    }

    fn commit(self) -> Result<(), Failure> {
        std::fs::rename(&self.tmp, &self.dest)
            .map_err(|e| Failure::Runtime(format!("cannot rename to {}: {e}", self.dest.display())))
    }
}

impl Drop for AtomicOutput {
    fn drop(&mut self) {
        // Removes the temporary file if commit() was not reached (no-op after a successful rename).
        let _ = std::fs::remove_file(&self.tmp);
    }
}

fn run_filter(args: &FilterArgs) -> Result<(), Failure> {
    use noodles::sam;
    use noodles::sam::alignment::RecordBuf;
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;
    use noodles::sam::header::record::value::{Map, map::program::tag as pg};
    use noodles_util::alignment::io::Format;

    let (settings, filters) = args.model.settings().map_err(Failure::Input)?;
    init_threads(args.model.threads)?;
    let to_stdout = args.output.as_os_str() == "-";
    let format = if to_stdout {
        match args.output_format {
            OutputFormat::Bam => Format::Bam,
            OutputFormat::Sam => Format::Sam,
        }
    } else {
        match args.output.extension().and_then(|e| e.to_str()) {
            Some("bam") => Format::Bam,
            Some("sam") => Format::Sam,
            _ => {
                return Err(Failure::Input(format!(
                    "cannot tell the output format of {}: use a .bam or .sam extension",
                    args.output.display()
                )));
            }
        }
    };

    let (mut reader, in_header) = open_reader(&args.input, args.model.threads)?;

    let mut out_header = in_header.clone();
    let program = Map::<sam::header::record::value::map::Program>::builder()
        .insert(pg::NAME, "pmdshodh")
        .insert(pg::VERSION, env!("CARGO_PKG_VERSION"))
        .insert(pg::COMMAND_LINE, command_line())
        .build()
        .map_err(|e| Failure::Runtime(format!("cannot build @PG line: {e}")))?;
    out_header
        .programs_mut()
        .add("pmdshodh", program)
        .map_err(|e| Failure::Runtime(format!("cannot add @PG line: {e}")))?;

    let (atomic, sink): (Option<AtomicOutput>, Box<dyn Write + Send>) = if to_stdout {
        (None, Box::new(io::stdout()))
    } else {
        let (a, f) = AtomicOutput::create(&args.output)?;
        (Some(a), Box::new(f))
    };
    let mut writer = OutWriter::new(format, sink, args.model.threads)?;
    let write_err = |e: io::Error| Failure::Runtime(format!("cannot write output: {e}"));
    writer.write_header(&out_header).map_err(write_err)?;

    let tag = match args.tag_name {
        TagName::Pm => Tag::new(b'p', b'm'),
        TagName::Ds => Tag::new(b'D', b'S'),
    };
    let mut summary = Summary::new(args.summary.range);
    scan(
        &mut reader,
        &in_header,
        &settings,
        &filters,
        &mut summary,
        None,
        |record, score| {
            let keep = args.threshold.is_none_or(|t| score.pmds >= t)
                && args.upper_threshold.is_none_or(|u| score.pmds < u);
            if !keep {
                return Ok(Some(false));
            }
            if args.tag {
                let mut buf = RecordBuf::try_from_alignment_record(&in_header, record)
                    .map_err(|e| Failure::Input(format!("cannot copy record: {e}")))?;
                let value = match args.tag_name {
                    TagName::Pm => Value::Float(score.pmds as f32),
                    TagName::Ds => Value::String(pmdtools_ds_text(score.pmds).into()),
                };
                buf.data_mut().insert(tag, value);
                writer.write_record(&out_header, &buf).map_err(write_err)?;
            } else {
                writer
                    .write_record(&out_header, record)
                    .map_err(write_err)?;
            }
            Ok(Some(true))
        },
    )?;
    writer.finish(&out_header).map_err(write_err)?;
    if let Some(a) = atomic {
        a.commit()?;
    }

    let fin = Finish {
        subcommand: "filter",
        input: &args.input,
        settings,
        filters,
        threshold: args.threshold,
        upper_threshold: args.upper_threshold,
    };
    let json = fin.json(&summary);
    if let Some(p) = &args.summary.summary {
        write_json(Some(p), &json)?;
    }
    fin.report(&json, &summary, args.quiet);
    Ok(())
}

/// Alignment output. BAM always goes through the multithreaded BGZF writer, so the bytes written
/// do not depend on the number of threads.
enum OutWriter {
    Sam(alignment::io::Writer<BufWriter<Box<dyn Write + Send>>>),
    Bam(noodles::bam::io::Writer<noodles::bgzf::io::MultithreadedWriter<Box<dyn Write + Send>>>),
}

impl OutWriter {
    fn new(
        format: noodles_util::alignment::io::Format,
        sink: Box<dyn Write + Send>,
        threads: usize,
    ) -> Result<Self, Failure> {
        use noodles_util::alignment::io::Format;
        Ok(match format {
            Format::Bam => {
                let workers =
                    std::num::NonZero::new(threads).unwrap_or(std::num::NonZero::<usize>::MIN);
                OutWriter::Bam(noodles::bam::io::Writer::from(
                    noodles::bgzf::io::MultithreadedWriter::with_worker_count(workers, sink),
                ))
            }
            _ => OutWriter::Sam(
                alignment::io::writer::Builder::default()
                    .set_format(format)
                    .build_from_writer(BufWriter::new(sink))
                    .map_err(|e| Failure::Runtime(format!("cannot open output: {e}")))?,
            ),
        })
    }

    fn write_header(&mut self, header: &noodles::sam::Header) -> io::Result<()> {
        match self {
            OutWriter::Sam(w) => w.write_header(header),
            OutWriter::Bam(w) => w.write_header(header),
        }
    }

    fn write_record<R: noodles::sam::alignment::Record>(
        &mut self,
        header: &noodles::sam::Header,
        record: &R,
    ) -> io::Result<()> {
        use noodles::sam::alignment::io::Write as _;
        match self {
            OutWriter::Sam(w) => w.write_record(header, record),
            OutWriter::Bam(w) => w.write_alignment_record(header, record),
        }
    }

    fn finish(&mut self, header: &noodles::sam::Header) -> io::Result<()> {
        match self {
            OutWriter::Sam(w) => w.finish(header),
            // Writes the BGZF end-of-file block and waits for the compression workers.
            OutWriter::Bam(w) => w.get_mut().finish().map(|_| ()),
        }
    }
}
