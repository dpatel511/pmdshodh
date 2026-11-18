use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::{ArgGroup, Args, Parser, Subcommand, ValueEnum};
use noodles_util::alignment;

use noodles::sam::alignment::Record as _;
use pmdshodh::damage::DamageCounts;
use pmdshodh::model::Params;
use pmdshodh::reads::{Outcome, ReadFilters, SkipReason, process, process_bam};
use pmdshodh::reference::Reference;
use pmdshodh::sites::ReadScore;
use pmdshodh::sites::{Library, Settings, Treatment};
use pmdshodh::summary::{RunInfo, Summary};

/// Per-read post-mortem damage (PMD) scores for ancient DNA alignments.
#[derive(Parser)]
#[command(name = "pmdshodh", bin_name = "pmdshodh", version)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Keep reads by PMD score and write them as SAM or BAM.
    #[command(after_help = FILTER_EXAMPLES)]
    Filter(FilterArgs),
    /// Write a table of per-read PMD scores (TSV).
    #[command(after_help = SCORE_EXAMPLES)]
    Score(ScoreArgs),
    /// Write the JSON run summary only: read counts, score histogram and damage table.
    #[command(after_help = DAMAGE_EXAMPLES)]
    Damage(DamageArgs),
}

const FILTER_EXAMPLES: &str = "\
Examples:
  # Keep reads with PMD score >= 3, as
  # samtools view -h in.bam | pmdtools --threshold 3 --header | samtools view -b -o out.bam -
  pmdshodh filter -t 3 -@ 8 in.bam -o out.bam

  # Single-stranded, UDG-half library; add the score as a pm tag and write a JSON summary
  pmdshodh filter -t 3 --library ss --treatment half --tag --summary run.json in.bam -o out.bam

  # In a pipe: SAM or BAM in on stdin, BAM out on stdout
  samtools view -h -q 30 in.bam | pmdshodh filter -t 3 - -o - | samtools sort -o out.bam";

const SCORE_EXAMPLES: &str = "\
Examples:
  # One line per read: read_name, flag, pmds, n_sites
  pmdshodh score -@ 8 in.bam -o scores.tsv

  # Scores as PMDtools 0.60 computes them, including its known bugs (for comparisons)
  pmdshodh score --pmdtools-compat in.bam -o scores.tsv";

const DAMAGE_EXAMPLES: &str = "\
Examples:
  # Substitution rates by position, score histogram and read counts as JSON
  pmdshodh damage -@ 8 in.bam -o damage.json

  # Only the first 100,000 reads (PMDtools -n)
  pmdshodh damage -n 100000 in.bam -o damage.json";

#[derive(Args)]
struct DamageArgs {
    /// Input SAM, BAM or CRAM file ("-" for stdin; CRAM needs -r).
    input: PathBuf,
    /// Output JSON file (default: stdout).
    #[arg(help_heading = "Output", short, long)]
    output: Option<PathBuf>,
    /// Positions from each read end in the damage table (PMDtools --range).
    #[arg(
        help_heading = "Summary",
        long,
        default_value_t = 30,
        allow_negative_numbers = true
    )]
    range: usize,
    /// Stop after this many input reads (PMDtools -n). Exits normally.
    #[arg(
        help_heading = "Input",
        short = 'n',
        long,
        visible_alias = "number",
        allow_negative_numbers = true
    )]
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
    /// Input SAM, BAM or CRAM file ("-" for stdin; CRAM needs -r).
    input: PathBuf,
    /// Output file; format from the extension (.bam or .sam). "-" writes BAM to stdout.
    #[arg(help_heading = "Output", short, long)]
    output: PathBuf,
    /// Output format for stdout ("-o -") and devices such as /dev/null [default: bam].
    #[arg(help_heading = "Output", long, value_enum)]
    output_format: Option<OutputFormat>,
    /// Keep reads with PMD score >= this (PMDtools -t/--threshold).
    #[arg(
        help_heading = "Filtering",
        short = 't',
        long,
        allow_negative_numbers = true
    )]
    threshold: Option<f64>,
    /// Keep reads with PMD score < this (PMDtools --upperthreshold). Reads with no informative
    /// sites carry no evidence and are not kept (PMDtools keeps them; so does --pmdtools-compat).
    #[arg(
        help_heading = "Filtering",
        long,
        visible_alias = "upperthreshold",
        allow_negative_numbers = true
    )]
    upper_threshold: Option<f64>,
    /// Add the PMD score to kept reads as a tag (overwrites an existing tag of that name).
    #[arg(help_heading = "Output", long)]
    tag: bool,
    /// Tag name: "pm" (pm:f, 32-bit float) or "DS" (DS:Z, 3 decimals, as PMDtools --writesamfield).
    #[arg(
        help_heading = "Output",
        long,
        value_enum,
        default_value = "pm",
        requires = "tag"
    )]
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
    /// Input SAM, BAM or CRAM file ("-" for stdin; CRAM needs -r).
    input: PathBuf,
    /// Output TSV file (default: stdout; "-" also writes to stdout).
    #[arg(help_heading = "Output", short, long)]
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
    #[arg(help_heading = "Scoring model", long, value_enum, default_value = "ds")]
    library: LibraryArg,
    /// Same as --library ss (PMDtools option).
    #[arg(help_heading = "Scoring model", long)]
    ss: bool,
    /// UDG treatment: none (all C/G sites), half (CpG sites plus the terminal base), full (CpG only).
    #[arg(
        help_heading = "Scoring model",
        long,
        value_enum,
        default_value = "none"
    )]
    treatment: TreatmentArg,
    /// Same as --treatment none (PMDtools option).
    #[arg(help_heading = "Scoring model", long = "UDGminus")]
    udg_minus: bool,
    /// Same as --treatment half (PMDtools option).
    #[arg(help_heading = "Scoring model", long = "UDGhalf")]
    udg_half: bool,
    /// Same as --treatment full (PMDtools --UDGplus / --CpG).
    #[arg(
        help_heading = "Scoring model",
        long = "UDGplus",
        visible_alias = "CpG"
    )]
    udg_plus: bool,
    /// Geometric decay of damage away from the read end (PMDtools --PMDpparam).
    #[arg(
        help_heading = "Scoring model",
        long,
        visible_alias = "PMDpparam",
        default_value_t = 0.3,
        allow_negative_numbers = true
    )]
    pmd_p: f64,
    /// Damage added at every position (PMDtools --PMDconstant).
    #[arg(
        help_heading = "Scoring model",
        long,
        visible_alias = "PMDconstant",
        default_value_t = 0.01,
        allow_negative_numbers = true
    )]
    pmd_c: f64,
    /// Polymorphism rate in the ancient model (PMDtools --polymorphism_ancient).
    #[arg(
        help_heading = "Scoring model",
        long,
        visible_alias = "polymorphism_ancient",
        default_value_t = 0.001,
        allow_negative_numbers = true
    )]
    pi_ancient: f64,
    /// Polymorphism rate in the null model (PMDtools --polymorphism_contamination).
    #[arg(
        help_heading = "Scoring model",
        long,
        visible_alias = "polymorphism_contamination",
        default_value_t = 0.001,
        allow_negative_numbers = true
    )]
    pi_contamination: f64,
    /// Skip reads with mapping quality below this (PMDtools -m).
    #[arg(
        help_heading = "Read filters",
        short = 'm',
        long,
        visible_alias = "requiremapq",
        default_value_t = 0,
        allow_negative_numbers = true
    )]
    min_mapq: u8,
    /// Ignore bases with base quality below this (PMDtools -q).
    #[arg(
        help_heading = "Read filters",
        short = 'q',
        long,
        visible_alias = "requirebaseq",
        default_value_t = 0,
        allow_negative_numbers = true
    )]
    min_baseq: u8,
    /// Skip reads with soft or hard clips (PMDtools --noclips).
    #[arg(help_heading = "Read filters", long, visible_alias = "noclips")]
    no_clips: bool,
    /// Skip reads with insertions or deletions (PMDtools --noindels).
    #[arg(help_heading = "Read filters", long, visible_alias = "noindels")]
    no_indels: bool,
    /// Reproduce PMDtools 0.60 exactly, including its bugs: the reference rebuilt wrongly for reads
    /// combining an insertion or soft clip with a deletion or '=' CIGAR operations, CpG context
    /// wrapping around the read ends, and PMDtools' CpG rule for damage tables. For comparisons.
    #[arg(help_heading = "Scoring model", long)]
    pmdtools_compat: bool,
    /// Threads to use: at most this many do work at once (reading, scoring, BGZF
    /// compression and writing). Output is identical for any value.
    #[arg(
        help_heading = "Performance",
        short = '@',
        long,
        default_value_t = 1,
        allow_negative_numbers = true
    )]
    threads: usize,
    /// Reference FASTA, indexed with samtools faidx (bgzip also needs its .gzi). Reference bases
    /// come from it instead of MD tags, so reads need no MD tag; MD tags are checked against it.
    /// Required for CRAM input.
    #[arg(help_heading = "Input", short = 'r', long)]
    reference: Option<PathBuf>,
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
        // The damage rate is largest at the read end, p + c; above 1 the scores are NaN.
        if pmdshodh::model::damage(1, self.pmd_p, self.pmd_c) > 1.0 {
            return Err(format!(
                "--pmd-p + --pmd-c is the damage rate at the read end and must be at most 1, got {} + {}",
                self.pmd_p, self.pmd_c
            ));
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
    /// The reader of stdout closed it (`| head`): exit quietly with status 141, as a process
    /// killed by SIGPIPE does (samtools behaves this way).
    BrokenPipe,
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
        Err(Failure::BrokenPipe) => ExitCode::from(141),
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
    #[arg(help_heading = "Summary", long)]
    summary: Option<PathBuf>,
    /// Positions from each read end in the damage table (PMDtools --range).
    #[arg(
        help_heading = "Summary",
        long,
        default_value_t = 30,
        allow_negative_numbers = true
    )]
    range: usize,
}

/// SAM and BAM go through noodles' reader; CRAM is decoded container by container on the pool.
enum InputReader {
    Alignment(alignment::io::Reader<Box<dyn Read + Send>>),
    Cram(pmdshodh::cram::Reader),
}
type AlignmentRecord = alignment::Record;

/// Input format, detected from the file contents.
#[derive(Clone, Copy)]
enum InputFormat {
    Sam,
    Bam,
    Cram,
}

impl InputFormat {
    fn name(self) -> &'static str {
        match self {
            InputFormat::Sam => "sam",
            InputFormat::Bam => "bam",
            InputFormat::Cram => "cram",
        }
    }
}

/// Open `-r` (if given) and the input, and check the input header against the reference before
/// any record is read. CRAM is decoded with `-r` only (docs/design.md section 11).
fn open_alignments(
    path: &Path,
    model: &ModelArgs,
) -> Result<
    (
        InputReader,
        noodles::sam::Header,
        Option<Reference>,
        InputFormat,
    ),
    Failure,
> {
    use std::io::BufRead;
    let reference_path = model.reference.as_deref();
    let mut reference = open_reference(reference_path)?;
    let cannot_read = |e: io::Error| Failure::Input(format!("cannot read {}: {e}", path.display()));
    let mut input = io::BufReader::new(open_input(&path.to_path_buf())?);
    let head = input.fill_buf().map_err(cannot_read)?;
    // No bytes at all is not an alignment file (samtools refuses it too); it usually means an
    // upstream step failed. A file with only a header is valid input.
    if head.is_empty() {
        return Err(Failure::Input(empty_input(path)));
    }
    let bgzf = head.starts_with(&[0x1f, 0x8b]);
    // BGZF sets FEXTRA (for its block size); plain gzip, as `gzip` writes it, does not.
    if bgzf && head.get(3).is_some_and(|flags| flags & 0x04 == 0) {
        return Err(Failure::Input(format!(
            "{} is gzip-compressed but not BGZF: recompress it with bgzip",
            path.display()
        )));
    }
    let mut source: Box<dyn BufRead + Send> = if bgzf {
        Box::new(io::BufReader::new(pmdshodh::bgzf::Reader::new(
            input,
            model.threads > 1,
        )))
    } else {
        Box::new(input)
    };
    let magic = source.fill_buf().map_err(cannot_read)?;
    if magic.is_empty() {
        return Err(Failure::Input(empty_input(path)));
    }
    let format = if magic.starts_with(b"CRAM") {
        InputFormat::Cram
    } else if magic.starts_with(b"BAM\x01") {
        InputFormat::Bam
    } else {
        if let Some(what) = not_sam(magic) {
            return Err(Failure::Input(format!(
                "{} is not SAM, BAM or CRAM: {what}",
                shown(path)
            )));
        }
        InputFormat::Sam
    };
    let source: Box<dyn Read + Send> = Box::new(source);
    let invalid_header = |e: io::Error| Failure::Input(format!("invalid header: {e}"));
    let (reader, header) = if let InputFormat::Cram = format {
        let (Some(r), Some(p)) = (reference.as_mut(), reference_path) else {
            return Err(Failure::Input(format!(
                "{} is CRAM, which needs its reference: pass it with -r FILE",
                path.display()
            )));
        };
        let mut reader = pmdshodh::cram::Reader::new(source);
        let header = reader.read_header().map_err(invalid_header)?;
        r.enable_cram(p, &header)
            .map_err(|e| reference_error(p, e))?;
        (InputReader::Cram(reader), header)
    } else {
        use noodles_util::alignment::io::Format;
        let mut reader = alignment::io::reader::Builder::default()
            .set_compression_method(None)
            .set_format(match format {
                InputFormat::Bam => Format::Bam,
                _ => Format::Sam,
            })
            .build_from_reader(source)
            .map_err(cannot_read)?;
        let header = reader.read_header().map_err(invalid_header)?;
        (InputReader::Alignment(reader), header)
    };
    if let (Some(r), Some(p)) = (reference.as_ref(), reference_path) {
        check_reference(r, p, &header)?;
    }
    Ok((reader, header, reference, format))
}

/// The input as the user named it ("-" is standard input).
fn shown(path: &Path) -> String {
    if path.as_os_str() == "-" {
        "standard input".into()
    } else {
        path.display().to_string()
    }
}

fn empty_input(path: &Path) -> String {
    format!(
        "{} is empty (no data): expected SAM, BAM or CRAM",
        shown(path)
    )
}

/// Why the start of an input that is neither BAM nor CRAM cannot be SAM, if it cannot.
fn not_sam(start: &[u8]) -> Option<&'static str> {
    let start = &start[..start.len().min(512)];
    if start
        .iter()
        .any(|&b| (b < 0x20 && !matches!(b, b'\t' | b'\n' | b'\r')) || b == 0x7f)
    {
        return Some("it is binary data");
    }
    match start.first() {
        Some(b'>') => Some("it looks like FASTA (sequences, not alignments)"),
        // SAM header lines are @HD, @SQ, @RG, @PG or @CO followed by a tab; FASTQ starts with @name.
        Some(b'@')
            if start.len() >= 4
                && !(matches!(&start[1..3], b"HD" | b"SQ" | b"RG" | b"PG" | b"CO")
                    && start[3] == b'\t') =>
        {
            Some("it looks like FASTQ (reads that are not aligned yet)")
        }
        _ => None,
    }
}

/// A record names its reference and its mate's by @SQ line; BAM stores the line's index, so a
/// name the header lacks cannot be written.
fn check_reference_ids(
    record: &AlignmentRecord,
    header: &noodles::sam::Header,
) -> Result<(), String> {
    let n = header.reference_sequences().len();
    for (what, id) in [
        ("reference", record.reference_sequence_id(header)),
        (
            "mate's reference",
            record.mate_reference_sequence_id(header),
        ),
    ] {
        match id {
            Some(Ok(i)) if i < n => {}
            None => {}
            _ => {
                return Err(format!(
                    "its {what} is not in the header (no @SQ line for it)"
                ));
            }
        }
    }
    Ok(())
}

/// Read every record (or the first `max_reads`), count it in `summary`, and pass scored reads to
/// `on_scored`, which returns whether the read was kept (`None` when there is no filtering).
/// Records read and scored per batch; bounds memory independently of the input size.
const BATCH: usize = 10_000;

/// Read every record (or the first `max_reads`), count it in `summary`, and pass scored reads to
/// `on_scored` in input order; it returns whether the read was kept (`None` without filtering).
/// Scoring runs on the rayon thread pool; results and summary are independent of the thread count.
///
/// With more than one thread, the stages overlap: the pool scores up to [`SCORED_AHEAD`] batches
/// while one pool thread reads the next and consumes the oldest. That thread does all serial work
/// (reading, consuming, writing) and runs other pool jobs (scoring, BGZF) while it waits. Batches are consumed in input
/// order, so output and summary are the same as the sequential loop used with one thread.
#[allow(clippy::too_many_arguments)]
fn scan(
    reader: &mut InputReader,
    header: &noodles::sam::Header,
    reference: &mut Option<Reference>,
    settings: &Settings,
    filters: &ReadFilters,
    summary: &mut Summary,
    count_damage: bool,
    max_reads: Option<u64>,
    threads: usize,
    mut on_scored: impl FnMut(&AlignmentRecord, &ReadScore) -> Result<Option<bool>, Failure> + Send,
) -> Result<(), Failure> {
    // Damage tables are only written to the JSON summary; counting them is a quarter of the work.
    let range = count_damage.then_some(summary.damage.range);
    if threads == 1 {
        let mut n_read = 0u64;
        loop {
            let mut batch = read_batch(reader, header, reference, max_reads, &mut n_read)?;
            if batch.is_empty() {
                return Ok(());
            }
            batch.decode(header)?;
            let parts = score_batch(&batch, settings, filters, range);
            let first = n_read - batch.len() as u64;
            consume_batch(&batch.records, &parts, first, summary, &mut on_scored)?;
        }
    }

    // The serial loop runs on a pool thread, so the pool's `threads` threads are all there is;
    // the calling thread just waits.
    rayon::scope(|s| {
        let spawn = |mut batch: Batch| {
            let (tx, rx) = std::sync::mpsc::sync_channel(1);
            s.spawn(move |_| {
                let parts = batch
                    .decode(header)
                    .map(|()| score_batch(&batch, settings, filters, range));
                let _ = tx.send((batch, parts));
            });
            rx
        };
        // A task that panicked drops its sender; the scope then re-raises that panic.
        let join = |rx: std::sync::mpsc::Receiver<(Batch, Result<Scored, Failure>)>| {
            pmdshodh::bgzf::wait_for(&rx).expect("scoring task panicked")
        };
        let mut n_read = 0u64;
        let mut first = 0u64;
        let mut in_flight = std::collections::VecDeque::new();
        let mut consume_oldest = |in_flight: &mut std::collections::VecDeque<_>| {
            let (batch, parts) = join(in_flight.pop_front().unwrap());
            consume_batch(&batch.records, &parts?, first, summary, &mut on_scored)?;
            first += batch.len() as u64;
            Ok::<(), Failure>(())
        };
        loop {
            match read_batch(reader, header, reference, max_reads, &mut n_read) {
                Ok(batch) if !batch.is_empty() => in_flight.push_back(spawn(batch)),
                Ok(_) => {
                    while !in_flight.is_empty() {
                        consume_oldest(&mut in_flight)?;
                    }
                    return Ok(());
                }
                Err(e) => {
                    // Earlier reads come first: report their errors before the read error.
                    while !in_flight.is_empty() {
                        consume_oldest(&mut in_flight)?;
                    }
                    return Err(e);
                }
            }
            if in_flight.len() > SCORED_AHEAD {
                consume_oldest(&mut in_flight)?;
            }
        }
    })
}

/// Batches scored on the pool while the serial thread consumes an older one (several threads).
const SCORED_AHEAD: usize = 2;

/// Outcomes and damage counts (if counted) of one batch, in input order, one entry per scoring
/// chunk.
type Scored = Vec<(Vec<Result<Outcome, String>>, Option<DamageCounts>)>;

/// Records read together, and with `-r` each mapped record's reference window.
struct Batch {
    records: Vec<AlignmentRecord>,
    /// All windows back to back; `windows[i]` is record `i`'s range in it.
    bases: Vec<u8>,
    windows: Vec<Option<std::ops::Range<usize>>>,
    /// CRAM containers still to decode into `records`, `bases` and `windows`.
    cram: Option<CramBatch>,
}

/// CRAM containers read on the serial thread and decoded with the batch's scoring, on the pool,
/// so the serial thread can load the next batch's reference sequences meanwhile.
struct CramBatch {
    containers: Vec<noodles::cram::io::reader::Container>,
    repository: noodles::fasta::Repository,
    /// `M5` checks of the reference sequences this batch loaded; all must pass before its records
    /// are used.
    checks: std::sync::Mutex<Vec<pmdshodh::reference::Md5Check>>,
    /// Records kept (fewer than decoded when `max_reads` ends in this batch).
    keep: usize,
    /// Records before this batch.
    first: u64,
}

impl Batch {
    fn len(&self) -> usize {
        self.records.len() + self.cram.as_ref().map_or(0, |c| c.keep)
    }

    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Decode CRAM containers, if any (on the pool).
    fn decode(&mut self, header: &noodles::sam::Header) -> Result<(), Failure> {
        let Some(c) = self.cram.take() else {
            return Ok(());
        };
        let error = |e: io::Error| {
            Failure::Input(format!(
                "CRAM records {} to {}: {e}",
                c.first + 1,
                c.first + c.keep as u64
            ))
        };
        for check in c.checks.into_inner().unwrap() {
            pmdshodh::bgzf::wait_for(&check)
                .unwrap_or_else(|| Err(io::Error::other("M5 check stopped")))
                .map_err(error)?;
        }
        let decoded =
            pmdshodh::cram::decode(&c.containers, &c.repository, header).map_err(error)?;
        self.bases = decoded.bases;
        self.windows = decoded.windows;
        self.windows.truncate(c.keep);
        self.records = decoded
            .records
            .into_iter()
            .take(c.keep)
            .map(AlignmentRecord::Cram)
            .collect();
        Ok(())
    }

    fn window(&self, i: usize) -> Option<&[u8]> {
        self.windows.get(i)?.clone().map(|r| &self.bases[r])
    }
}

/// Read up to [`BATCH`] records (CRAM: one container per pool thread); an empty batch means the
/// end of the input (or of `max_reads`). With a reference, fetch each mapped record's window
/// here, on the serial thread (CRAM: on the pool, while decoding).
fn read_batch(
    reader: &mut InputReader,
    header: &noodles::sam::Header,
    reference: &mut Option<Reference>,
    max_reads: Option<u64>,
    n_read: &mut u64,
) -> Result<Batch, Failure> {
    let reader = match reader {
        InputReader::Alignment(reader) => reader,
        InputReader::Cram(reader) => {
            return read_cram_batch(reader, header, reference, max_reads, n_read);
        }
    };
    let mut batch = Batch {
        records: Vec::with_capacity(BATCH),
        bases: Vec::new(),
        windows: Vec::new(),
        cram: None,
    };
    while batch.len() < BATCH && max_reads.is_none_or(|m| *n_read < m) {
        let mut record = AlignmentRecord::default();
        let n = reader
            .read_record(header, &mut record)
            .map_err(|e| Failure::Input(format!("record {}: {e}", *n_read + 1)))?;
        if n == 0 {
            break;
        }
        *n_read += 1;
        if let Some(reference) = reference.as_mut() {
            let window = fetch_window(reference, header, &record, &mut batch.bases)
                .map_err(|e| Failure::Input(format!("record {}: {e}", *n_read)))?;
            batch.windows.push(window);
        }
        batch.records.push(record);
    }
    Ok(batch)
}

/// Read one CRAM container per pool thread and load their reference sequences; they are decoded
/// on the pool (`Batch::decode`). Records are counted from the container headers.
fn read_cram_batch(
    reader: &mut pmdshodh::cram::Reader,
    header: &noodles::sam::Header,
    reference: &mut Option<Reference>,
    max_reads: Option<u64>,
    n_read: &mut u64,
) -> Result<Batch, Failure> {
    let error = |e: io::Error| Failure::Input(format!("CRAM records after record {n_read}: {e}"));
    let mut batch = Batch {
        records: Vec::new(),
        bases: Vec::new(),
        windows: Vec::new(),
        cram: None,
    };
    if max_reads.is_some_and(|m| *n_read >= m) {
        return Ok(batch);
    }
    let containers = reader
        .read_containers(rayon::current_num_threads())
        .map_err(error)?;
    let ids: Vec<usize> = containers
        .iter()
        .filter_map(pmdshodh::cram::reference_sequence_id)
        .collect();
    let (repository, checks) = reference
        .as_mut()
        .expect("CRAM input is opened with -r")
        .cram_batch(header, &ids)
        .map_err(error)?;
    let records: u64 = containers
        .iter()
        .map(|c| c.header().record_count() as u64)
        .sum();
    let keep = max_reads.map_or(records, |m| records.min(m - *n_read));
    if keep > 0 {
        batch.cram = Some(CramBatch {
            containers,
            repository,
            checks: checks.into(),
            keep: keep as usize,
            first: *n_read,
        });
        *n_read += keep;
    }
    Ok(batch)
}

/// Append a mapped record's reference bases (alignment start, CIGAR reference span) to `bases`
/// and return their range; `None` for an unmapped record.
fn fetch_window(
    reference: &mut Reference,
    header: &noodles::sam::Header,
    record: &AlignmentRecord,
    bases: &mut Vec<u8>,
) -> io::Result<Option<std::ops::Range<usize>>> {
    if record.flags()?.is_unmapped() {
        return Ok(None);
    }
    let invalid = |msg: &str| io::Error::new(io::ErrorKind::InvalidData, msg.to_string());
    let id = record
        .reference_sequence_id(header)
        .transpose()?
        .ok_or_else(|| invalid("mapped read has no reference sequence"))?;
    let start = record
        .alignment_start()
        .transpose()?
        .ok_or_else(|| invalid("mapped read has no position"))?;
    let span = record.alignment_span().transpose()?.unwrap_or(0);
    let (name, _) = header
        .reference_sequences()
        .get_index(id)
        .ok_or_else(|| invalid("reference sequence ID is not in the header"))?;
    let name = std::str::from_utf8(name).map_err(|_| invalid("reference name is not UTF-8"))?;
    let from = bases.len();
    reference.fetch(name, usize::from(start) - 1, span, bases)?;
    Ok(Some(from..bases.len()))
}

/// Open `-r`; a missing file or index is an input error.
fn open_reference(path: Option<&Path>) -> Result<Option<Reference>, Failure> {
    let Some(path) = path else {
        return Ok(None);
    };
    if !path.exists() {
        return Err(Failure::Input(format!(
            "cannot open reference {}: no such file",
            path.display()
        )));
    }
    Reference::open(path)
        .map(Some)
        .map_err(|e| reference_error(path, e))
}

fn reference_error(path: &Path, e: io::Error) -> Failure {
    Failure::Input(format!(
        "cannot read reference {} or its index: {e} (index it with samtools faidx {})",
        path.display(),
        path.display()
    ))
}

/// Check that every header sequence is in the reference with the same length (docs/design.md
/// section 11).
fn check_reference(
    reference: &Reference,
    path: &Path,
    header: &noodles::sam::Header,
) -> Result<(), Failure> {
    let mut differ = Vec::new();
    for (name, seq) in header.reference_sequences() {
        let name = name.to_string();
        let len = usize::from(seq.length());
        match reference.length(&name) {
            Some(l) if l == len => {}
            Some(l) => differ.push(format!(
                "{name} has {len} bases in the header, {l} in the reference"
            )),
            None => differ.push(format!("{name} is not in the reference")),
        }
    }
    if let Some(first) = differ.first() {
        return Err(Failure::Input(format!(
            "reference {} does not match the input header ({} of {} sequences differ; {first})",
            path.display(),
            differ.len(),
            header.reference_sequences().len()
        )));
    }
    Ok(())
}

fn score_batch(
    batch: &Batch,
    settings: &Settings,
    filters: &ReadFilters,
    range: Option<usize>,
) -> Scored {
    use rayon::prelude::*;
    let chunk = batch.len().div_ceil(rayon::current_num_threads()).max(1);
    batch
        .records
        .par_chunks(chunk)
        .enumerate()
        .map(|(c, records)| {
            let mut damage = range.map(DamageCounts::new);
            let outcomes = records
                .iter()
                .enumerate()
                .map(|(k, r)| {
                    let window = batch.window(c * chunk + k);
                    match r {
                        AlignmentRecord::Bam(b) => {
                            process_bam(b, settings, filters, window, damage.as_mut())
                        }
                        _ => process(r, settings, filters, window, damage.as_mut()),
                    }
                })
                .collect();
            (outcomes, damage)
        })
        .collect()
}

/// Count a scored batch in `summary` and pass its scored reads to `on_scored`, in order.
/// `first` is the number of records before this batch (for error messages).
fn consume_batch(
    batch: &[AlignmentRecord],
    parts: &Scored,
    first: u64,
    summary: &mut Summary,
    on_scored: &mut impl FnMut(&AlignmentRecord, &ReadScore) -> Result<Option<bool>, Failure>,
) -> Result<(), Failure> {
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
            Ok(Outcome::Scored { score, length, md }) => {
                let kept = on_scored(record, score)?;
                summary.add_scored(score, *length, kept);
                summary.add_md_check(*md);
            }
            Ok(Outcome::Skipped(reason)) => summary.add_skipped(*reason),
            Err(msg) => {
                let name = record.name().map_or_else(|| "*".into(), |n| n.to_string());
                return Err(Failure::Input(format!("read {name}: {msg}")));
            }
        }
    }
    for damage in parts.iter().filter_map(|(_, d)| d.as_ref()) {
        summary.damage.merge(damage);
    }
    Ok(())
}

/// Largest --range: far beyond ancient DNA read lengths; each position costs memory per table.
const MAX_RANGE: usize = 1000;
/// Largest -@: far beyond any node, while still failing fast instead of starting threads for minutes.
const MAX_THREADS: usize = 1024;

/// The damage table needs at least the terminal position.
fn check_range(range: usize) -> Result<(), Failure> {
    if range == 0 {
        return Err(Failure::Input("--range must be at least 1".into()));
    }
    if range > MAX_RANGE {
        return Err(Failure::Input(format!(
            "--range must be at most {MAX_RANGE}, got {range}"
        )));
    }
    Ok(())
}

/// At most `threads` threads do work at once. With one, this thread does everything (the pool
/// thread only runs scoring while this thread waits). With more, the pool's `threads` threads do
/// all the work, one of them the serial part (scan, bgzf.rs), while this thread waits.
fn init_threads(threads: usize) -> Result<(), Failure> {
    if threads == 0 {
        return Err(Failure::Input("--threads must be at least 1".into()));
    }
    if threads > MAX_THREADS {
        return Err(Failure::Input(format!(
            "--threads must be at most {MAX_THREADS}, got {threads}"
        )));
    }
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .map_err(|e| Failure::Runtime(format!("cannot start threads: {e}")))
}

fn command_line() -> String {
    std::env::args().collect::<Vec<_>>().join(" ")
}

/// Write JSON to an output opened before the run, and commit it.
fn write_json(
    (output, mut w): (Output, Box<dyn Write + Send>),
    value: &serde_json::Value,
) -> Result<(), Failure> {
    let text = serde_json::to_string_pretty(value)
        .map_err(|e| Failure::Runtime(format!("cannot format summary: {e}")))?
        + "\n";
    let write_err = output.write_err();
    w.write_all(text.as_bytes()).map_err(&write_err)?;
    w.flush().map_err(&write_err)?;
    drop(w);
    output.commit()
}

/// Facts shared by the end of every subcommand.
struct Finish<'a> {
    subcommand: &'a str,
    input: &'a Path,
    input_format: InputFormat,
    settings: Settings,
    filters: ReadFilters,
    threshold: Option<f64>,
    upper_threshold: Option<f64>,
    reference: Option<&'a Path>,
}

impl Finish<'_> {
    fn json(&self, summary: &Summary) -> serde_json::Value {
        summary.to_json(&RunInfo {
            subcommand: self.subcommand,
            command: command_line(),
            input: self.input.display().to_string(),
            input_format: self.input_format.name(),
            settings: self.settings,
            filters: self.filters,
            threshold: self.threshold,
            upper_threshold: self.upper_threshold,
            reference: self.reference.map(|p| p.display().to_string()),
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
        if summary.md_disagreeing() > 0 {
            eprintln!(
                "pmdshodh {}: {} scored reads have an MD tag that disagrees with the reference \
                 (usually stale after realignment); they were scored from the reference",
                self.subcommand,
                summary.md_disagreeing()
            );
        }
    }
}

fn run_score(args: &ScoreArgs) -> Result<(), Failure> {
    let (settings, filters) = args.model.settings().map_err(Failure::Input)?;
    check_range(args.summary.range)?;
    init_threads(args.model.threads)?;
    check_outputs(
        &args.input,
        args.model.reference.as_deref(),
        &[args.output.as_deref(), args.summary.summary.as_deref()],
    )?;
    let (mut reader, header, mut reference, input_format) =
        open_alignments(&args.input, &args.model)?;
    // Both outputs are created before any record is read, so a bad path fails at once.
    let (output, out) = Output::open(args.output.as_deref())?;
    let summary_out = args
        .summary
        .summary
        .as_deref()
        .map(|p| Output::open(Some(p)))
        .transpose()?;
    let mut out = BufWriter::new(out);
    let write_err = output.write_err();

    writeln!(out, "read_name\tflag\tpmds\tn_sites").map_err(&write_err)?;
    let mut summary = Summary::new(args.summary.range);
    scan(
        &mut reader,
        &header,
        &mut reference,
        &settings,
        &filters,
        &mut summary,
        args.summary.summary.is_some(),
        None,
        args.model.threads,
        |record, s| {
            let name = record
                .name()
                .map_or_else(|| "*".to_string(), |n| n.to_string());
            let flag = record.flags().map(u16::from).unwrap_or_default();
            writeln!(out, "{name}\t{flag}\t{}\t{}", s.pmds, s.n_sites).map_err(&write_err)?;
            Ok(None)
        },
    )?;
    out.flush().map_err(&write_err)?;
    drop(out);
    output.commit()?;
    let fin = Finish {
        subcommand: "score",
        input: &args.input,
        input_format,
        settings,
        filters,
        threshold: None,
        upper_threshold: None,
        reference: args.model.reference.as_deref(),
    };
    let json = fin.json(&summary);
    if let Some(s) = summary_out {
        write_json(s, &json)?;
    }
    fin.report(&json, &summary, args.quiet);
    Ok(())
}

fn run_damage(args: &DamageArgs) -> Result<(), Failure> {
    let (settings, filters) = args.model.settings().map_err(Failure::Input)?;
    check_range(args.range)?;
    init_threads(args.model.threads)?;
    check_outputs(
        &args.input,
        args.model.reference.as_deref(),
        &[args.output.as_deref()],
    )?;
    let (mut reader, header, mut reference, input_format) =
        open_alignments(&args.input, &args.model)?;
    let json_out = Output::open(args.output.as_deref())?;
    let mut summary = Summary::new(args.range);
    scan(
        &mut reader,
        &header,
        &mut reference,
        &settings,
        &filters,
        &mut summary,
        true,
        args.max_reads,
        args.model.threads,
        |_, _| Ok(None),
    )?;
    let fin = Finish {
        subcommand: "damage",
        input: &args.input,
        input_format,
        settings,
        filters,
        threshold: None,
        upper_threshold: None,
        reference: args.model.reference.as_deref(),
    };
    let json = fin.json(&summary);
    write_json(json_out, &json)?;
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

/// Where a command writes. A regular file, or a path that does not exist yet, is written under a
/// temporary name and renamed when the run succeeds, so a failed run leaves nothing under the
/// final name; a symlink is followed and the file it points to is replaced. Anything else that
/// exists (a FIFO, a device such as /dev/null, /dev/stdout, a process substitution) is written in
/// place: renaming over it would replace it.
struct Output {
    /// The output as the user named it, for messages.
    shown: String,
    /// Temporary file and final name, when written by rename.
    rename: Option<(PathBuf, PathBuf)>,
    to_stdout: bool,
}

impl Output {
    /// `None` or "-" is stdout.
    fn open(path: Option<&Path>) -> Result<(Self, Box<dyn Write + Send>), Failure> {
        let Some(path) = path.filter(|p| p.as_os_str() != "-") else {
            let out = Output {
                shown: "standard output".into(),
                rename: None,
                to_stdout: true,
            };
            return Ok((out, Box::new(io::stdout())));
        };
        let shown = path.display().to_string();
        let cannot = |e: io::Error| Failure::Runtime(format!("cannot create {shown}: {e}"));
        let dest = match std::fs::metadata(path) {
            _ if is_windows_device(path) => None,
            Ok(m) if m.is_dir() => {
                return Err(Failure::Input(format!("output {shown} is a directory")));
            }
            Ok(m) if m.is_file() => Some(std::fs::canonicalize(path).map_err(cannot)?),
            Ok(_) => None,
            // A dangling symlink: creating the file creates its target.
            Err(_) if path.symlink_metadata().is_ok() => None,
            Err(_) => Some(path.to_path_buf()),
        };
        let (file, rename) = match dest {
            Some(dest) => {
                let mut name = dest.file_name().unwrap_or_default().to_os_string();
                name.push(".pmdshodh-tmp");
                let tmp = dest.with_file_name(name);
                (File::create(&tmp).map_err(cannot)?, Some((tmp, dest)))
            }
            None => (File::create(path).map_err(cannot)?, None),
        };
        let out = Output {
            shown,
            rename,
            to_stdout: false,
        };
        Ok((out, Box::new(file)))
    }

    /// Error while writing; a reader that closed stdout ends the run quietly.
    fn write_err(&self) -> impl Fn(io::Error) -> Failure + Send + Sync + use<> {
        let (shown, to_stdout) = (self.shown.clone(), self.to_stdout);
        move |e| {
            if to_stdout && e.kind() == io::ErrorKind::BrokenPipe {
                Failure::BrokenPipe
            } else {
                Failure::Runtime(format!("cannot write {shown}: {e}"))
            }
        }
    }

    /// Give a renamed output its final name (call after the writer is closed).
    fn commit(mut self) -> Result<(), Failure> {
        if let Some((tmp, dest)) = self.rename.take()
            && let Err(e) = std::fs::rename(&tmp, &dest)
        {
            let _ = std::fs::remove_file(&tmp);
            return Err(Failure::Runtime(format!(
                "cannot write {}: {e}",
                self.shown
            )));
        }
        Ok(())
    }
}

impl Drop for Output {
    fn drop(&mut self) {
        // Removes the temporary file if commit() was not reached.
        if let Some((tmp, _)) = &self.rename {
            let _ = std::fs::remove_file(tmp);
        }
    }
}

/// An existing output that is not a regular file or directory: a FIFO or a device.
fn is_special(path: &Path) -> bool {
    is_windows_device(path) || std::fs::metadata(path).is_ok_and(|m| !m.is_file() && !m.is_dir())
}

/// Windows device names (NUL, CON, COM1, ...), which have no file metadata; a temporary file
/// next to one (NUL.pmdshodh-tmp) is the device again, so it cannot be renamed into place.
fn is_windows_device(path: &Path) -> bool {
    if !cfg!(windows) {
        return false;
    }
    let Some(name) = path.file_name().and_then(|n| n.to_str()) else {
        return false;
    };
    let name = name.to_ascii_uppercase();
    matches!(name.as_str(), "CON" | "PRN" | "AUX" | "NUL")
        || (name.len() == 4
            && (name.starts_with("COM") || name.starts_with("LPT"))
            && matches!(name.as_bytes()[3], b'1'..=b'9'))
}

/// Whether two paths name the same existing regular file (also through links or `./`).
fn same_file(a: &Path, b: &Path) -> bool {
    let (Ok(ma), Ok(mb)) = (std::fs::metadata(a), std::fs::metadata(b)) else {
        return false;
    };
    if !ma.is_file() || !mb.is_file() {
        return false;
    }
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        ma.dev() == mb.dev() && ma.ino() == mb.ino()
    }
    #[cfg(not(unix))]
    {
        std::fs::canonicalize(a).ok() == std::fs::canonicalize(b).ok()
    }
}

/// Refuse to write over the input or the reference, or two outputs to one file: renaming the
/// output into place would silently replace the file being read.
fn check_outputs(
    input: &Path,
    reference: Option<&Path>,
    outputs: &[Option<&Path>],
) -> Result<(), Failure> {
    let outputs: Vec<&Path> = outputs
        .iter()
        .flatten()
        .copied()
        .filter(|p| p.as_os_str() != "-")
        .collect();
    for (i, out) in outputs.iter().enumerate() {
        for (what, other) in [("input", Some(input)), ("reference", reference)] {
            if other.is_some_and(|o| same_file(out, o)) {
                return Err(Failure::Input(format!(
                    "output {} is the {what} file: write to a different file",
                    out.display()
                )));
            }
        }
        if outputs[i + 1..]
            .iter()
            .any(|later| same_file(out, later) || (out == later && !is_special(out)))
        {
            return Err(Failure::Input(format!(
                "{} is given for two outputs",
                out.display()
            )));
        }
    }
    Ok(())
}

fn run_filter(args: &FilterArgs) -> Result<(), Failure> {
    use noodles::sam;
    use noodles::sam::alignment::RecordBuf;
    use noodles::sam::alignment::record::data::field::Tag;
    use noodles::sam::alignment::record_buf::data::field::Value;
    use noodles::sam::header::record::value::{Map, map::program::tag as pg};

    let (settings, filters) = args.model.settings().map_err(Failure::Input)?;
    for (name, bound) in [
        ("--threshold", args.threshold),
        ("--upper-threshold", args.upper_threshold),
    ] {
        if bound.is_some_and(|b| !b.is_finite()) {
            return Err(Failure::Input(format!(
                "{name} must be a finite number, got {}",
                bound.unwrap()
            )));
        }
    }
    check_range(args.summary.range)?;
    init_threads(args.model.threads)?;
    let format = output_format(&args.output, args.output_format)?;
    check_outputs(
        &args.input,
        args.model.reference.as_deref(),
        &[Some(&args.output), args.summary.summary.as_deref()],
    )?;

    let (mut reader, in_header, mut reference, input_format) =
        open_alignments(&args.input, &args.model)?;

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

    // Both outputs are created before any record is read, so a bad path fails at once.
    let (output, sink) = Output::open(Some(&args.output))?;
    let summary_out = args
        .summary
        .summary
        .as_deref()
        .map(|p| Output::open(Some(p)))
        .transpose()?;
    let mut writer = OutWriter::new(format, sink, args.model.threads)?;
    let write_err = output.write_err();
    writer.write_header(&out_header).map_err(&write_err)?;

    let tag = match args.tag_name {
        TagName::Pm => Tag::new(b'p', b'm'),
        TagName::Ds => Tag::new(b'D', b'S'),
    };
    let mut summary = Summary::new(args.summary.range);
    scan(
        &mut reader,
        &in_header,
        &mut reference,
        &settings,
        &filters,
        &mut summary,
        args.summary.summary.is_some(),
        None,
        args.model.threads,
        |record, score| {
            let keep = args.threshold.is_none_or(|t| score.pmds >= t)
                && args.upper_threshold.is_none_or(|u| {
                    score.pmds < u && (score.n_sites > 0 || settings.pmdtools_compat)
                });
            if !keep {
                return Ok(Some(false));
            }
            // score and damage read SAM without a header (`samtools view | ...`); a kept record
            // can only be written if its references are @SQ lines.
            check_reference_ids(record, &in_header).map_err(|msg| {
                let name = record.name().map_or_else(|| "*".into(), |n| n.to_string());
                Failure::Input(format!("read {name}: {msg}"))
            })?;
            if args.tag {
                let mut buf = RecordBuf::try_from_alignment_record(&in_header, record)
                    .map_err(|e| Failure::Input(format!("cannot copy record: {e}")))?;
                let value = match args.tag_name {
                    TagName::Pm => Value::Float(score.pmds as f32),
                    TagName::Ds => Value::String(pmdtools_ds_text(score.pmds).into()),
                };
                buf.data_mut().insert(tag, value);
                writer.write_record(&out_header, &buf).map_err(&write_err)?;
            } else {
                writer
                    .write_record(&out_header, record)
                    .map_err(&write_err)?;
            }
            Ok(Some(true))
        },
    )?;
    writer.finish(&out_header).map_err(&write_err)?;
    drop(writer);
    output.commit()?;

    let fin = Finish {
        subcommand: "filter",
        input: &args.input,
        input_format,
        settings,
        filters,
        threshold: args.threshold,
        upper_threshold: args.upper_threshold,
        reference: args.model.reference.as_deref(),
    };
    let json = fin.json(&summary);
    if let Some(s) = summary_out {
        write_json(s, &json)?;
    }
    fin.report(&json, &summary, args.quiet);
    Ok(())
}

/// The format of `filter`'s output: from a .bam or .sam extension (any case). Stdout and devices
/// such as /dev/null have none, so --output-format decides there, BAM by default.
fn output_format(
    path: &Path,
    flag: Option<OutputFormat>,
) -> Result<noodles_util::alignment::io::Format, Failure> {
    use noodles_util::alignment::io::Format;
    let name = |f: OutputFormat| match f {
        OutputFormat::Bam => "bam",
        OutputFormat::Sam => "sam",
    };
    let to_stdout = path.as_os_str() == "-";
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .map(str::to_ascii_lowercase);
    let from_ext = match ext.as_deref() {
        _ if to_stdout => None,
        Some("bam") => Some(OutputFormat::Bam),
        Some("sam") => Some(OutputFormat::Sam),
        Some("cram") => {
            return Err(Failure::Input(format!(
                "cannot write {}: CRAM output is not supported; use .bam or .sam",
                path.display()
            )));
        }
        _ => None,
    };
    let format = match (from_ext, flag) {
        (Some(e), Some(f)) if e != f => {
            return Err(Failure::Input(format!(
                "--output-format {} contradicts the .{} extension of {}",
                name(f),
                name(e),
                path.display()
            )));
        }
        (Some(e), _) => e,
        (None, f) if to_stdout || is_special(path) => f.unwrap_or(OutputFormat::Bam),
        (None, _) => {
            return Err(Failure::Input(format!(
                "cannot tell the output format of {}: use a .bam or .sam extension",
                path.display()
            )));
        }
    };
    Ok(match format {
        OutputFormat::Bam => Format::Bam,
        OutputFormat::Sam => Format::Sam,
    })
}

/// Alignment output. BAM goes through pmdshodh's BGZF writer, whose bytes do not depend on the
/// number of threads.
enum OutWriter {
    Sam(alignment::io::Writer<BufWriter<Box<dyn Write + Send>>>),
    Bam(noodles::bam::io::Writer<pmdshodh::bgzf::Writer<Box<dyn Write + Send>>>),
}

impl OutWriter {
    fn new(
        format: noodles_util::alignment::io::Format,
        sink: Box<dyn Write + Send>,
        threads: usize,
    ) -> Result<Self, Failure> {
        use noodles_util::alignment::io::Format;
        Ok(match format {
            Format::Bam => OutWriter::Bam(noodles::bam::io::Writer::from(
                pmdshodh::bgzf::Writer::new(sink, threads > 1),
            )),
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
            OutWriter::Bam(w) => w.get_mut().finish(),
        }
    }
}
