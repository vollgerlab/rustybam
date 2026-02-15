use clap::IntoApp;
use clap::{AppSettings, Parser, Subcommand};

#[derive(Parser, Debug)]
#[clap(
    author,
    version,
    about,
    propagate_version = true,
    subcommand_required = true,
    infer_subcommands = true,
    arg_required_else_help = true,
    help_expected = true
)]
#[clap(global_setting(AppSettings::DeriveDisplayOrder))]
pub struct Cli {
    /// Threads for decompression.
    #[clap(short, long, default_value_t = 8)]
    pub threads: usize,

    /// Logging level [-v: Info, -vv: Debug, -vvv: Trace].
    #[clap(short, long, parse(from_occurrences), help_heading = "DEBUG")]
    pub verbose: usize,

    #[clap(subcommand)]
    pub command: Option<Commands>,
}

///
/// This structure contains all the subcommands for rustybam and their help descriptions.
///
/// Because of naming conventions for rust enums the commands names have
/// different capitalization than on the command line.
/// For example, the `Liftover` enum is invoked using `rustybam liftover`
/// and the `TrimPaf` command with `rustybam trim-paf`.
///
#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Get percent identity stats from a sam/bam/cram or PAF.
    /// Requires =/X operations in the CIGAR string!
    ///
    /// ## output column descriptions:
    /// ### perID_by_matches is calculated as:
    ///  `matches / (matches + mismatches)`
    /// ### perID_by_events is calculated as:
    ///  `matches / (matches + mismatches + insertion events + deletion events)`
    /// ### perID_by_all is calculated as:
    ///  `matches / (matches + mismatches + insertion bases + deletion bases)`
    Stats {
        /// Input sam/bam/cram/file.
        #[clap(default_value = "-")]
        bam: String,
        /// Print query coordinates first
        #[clap(short, long)]
        qbed: bool,
        /// Specify that the input is paf format,
        /// (must have cg tag with extended cigar).
        #[clap(short, long)]
        paf: bool,
    },
    /// Count the number of bases in a bed file.
    #[clap(visible_aliases = &["bedlen", "bl", "bedlength"])]
    BedLength {
        /// Input bed file.
        #[clap(default_value = "-")]
        bed: Vec<String>,
        /// Make the output human readable (Mbp).
        #[clap(short, long)]
        readable: bool,
        /// Count bases for each category in this column <COLUMN>.
        #[clap(short, long)]
        column: Option<u8>,
    },
    /// Filter PAF records in various ways.
    Filter {
        /// PAF file from minimap2 or unimap. Must have the cg tag, and n matches will be zero unless the cigar uses =X.
        #[clap(default_value = "-")]
        paf: String,
        /// Minimum number of aligned bases across all alignments between a target and query in order to keep those records.
        #[clap(short, long, default_value_t = 0)]
        paired_len: u64,
        /// Minimum alignment length.
        #[clap(short, long, default_value_t = 0)]
        aln: u64,
        /// Minimum query length.
        #[clap(short, long, default_value_t = 0)]
        query: u64,
    },
    /// Invert the target and query sequences in a PAF along with the CIGAR string.
    Invert {
        /// PAF file from minimap2 or unimap. Must have the cg tag, and n matches will be zero unless the cigar uses =X.
        #[clap(default_value = "-")]
        paf: String,
    },
    /// Liftover target sequence coordinates onto query sequence using a PAF.
    ///
    /// This is a function for lifting over coordinates from a reference (<BED>) to a query using a PAF file from minimap2 or unimap (note, you can use `paftools.js sam2paf` to convert SAM data to PAF format).
    /// The returned file is a PAF file that is trimmed to the regions in the bed file. Even the cigar in the returned PAF file is trimmed so it can be used downstream! Additionally, a tag with the format `id:Z:<>` is added to the PAF where `<>` is either the 4th column of the input bed file or if not present `chr_start_end`.
    #[clap(visible_aliases=&["lo"], aliases = &["william-t-harvey", "wth"])]
    Liftover {
        /// PAF file from minimap2 or unimap run with -c and --eqx [i.e. the PAF file must have the cg tag and use extended CIGAR opts (=/X)].
        #[clap(default_value = "-")]
        paf: String,
        /// BED file of reference regions to liftover to the query.
        #[clap(short, long)]
        bed: String,
        /// Specifies that the BED file contains query coordinates to be lifted onto the reference (reverses direction of liftover).
        ///
        /// Note, that this will make the query in the input `PAF` the target in the output `PAF`.
        #[clap(short, long)]
        qbed: bool,
        /// If multiple records overlap the same region in the <bed> return only the largest liftover. The default is to return all liftovers.
        #[clap(short, long)]
        largest: bool,
    },
    /// Trim PAF records that overlap in query sequence to find an optimal splitting point using dynamic programming.
    ///
    /// Note, this can be combined with `rb invert` to also trim the target sequence.
    ///
    /// This idea is to mimic some of the trimming that happens in PAV to improve breakpoint detection. Starts with the largest overlap and iterates until no query overlaps remain.
    #[clap(visible_aliases = &["trim", "tp"])]
    TrimPaf {
        /// PAF file from minimap2 or unimap. Must have the cg tag, and n matches will be zero unless the cigar uses =X.
        #[clap(default_value = "-")]
        paf: String,
        /// Value added for a matching base.
        #[clap(short, long, default_value_t = 1)]
        match_score: i32,
        /// Value subtracted for a mismatching base.
        #[clap(short, long, default_value_t = 1)]
        diff_score: i32,
        /// Value subtracted for an insertion or deletion.
        #[clap(short, long, default_value_t = 1)]
        indel_score: i32,
        /// Remove contained alignments as well as overlaps.
        #[clap(short, long)]
        remove_contained: bool,
    },
    /// Orient paf records so that most of the bases are in the forward direction.
    ///
    /// Optionally scaffold the queriers so that there is one query per target.
    Orient {
        /// PAF file from minimap2 or unimap. Must have the cg tag, and n matches will be zero unless the cigar uses =X.
        #[clap(default_value = "-")]
        paf: String,
        /// Generate ~1 query per target that scaffolds together all the records that map to that target sequence.
        ///
        /// The order of the scaffold will be determined by the average target position across all the queries, and the name of the new scaffold will be.
        #[clap(short, long)]
        scaffold: bool,
        /// Space to add between records.
        #[clap(short, long, default_value_t = 1_000_000)]
        insert: u64,
    },
    /// Break PAF records with large indels into multiple records (useful for SafFire).
    #[clap(visible_aliases = &["breakpaf", "bp"])]
    BreakPaf {
        /// PAF file from minimap2 or unimap. Must have the cg tag, and n matches will be zero unless the cigar uses =X.
        #[clap(default_value = "-")]
        paf: String,
        /// Maximum indel size to keep in the paf.
        #[clap(short, long, default_value_t = 100)]
        max_size: u32,
    },
    /// Convert a PAF file into a SAM file. Warning, all alignments will be marked as primary!
    #[clap(visible_aliases = &["paftosam", "p2s", "paf2sam"])]
    PafToSam {
        /// PAF file from minimap2 or unimap. Must have a CIGAR tag.
        #[clap(default_value = "-")]
        paf: String,
        /// Optional query fasta file (with index) to populate the query seq field.
        #[clap(short, long)]
        fasta: Option<String>,
    },
    /// Splits fastx from stdin into multiple files.
    ///
    /// Specifically it reads fastx format (fastq, fasta, or mixed) from stdin and divides the records across multiple output files. Output files can be compressed by adding `.gz`, and the input can also be compressed or uncompressed.
    #[clap(visible_aliases = &["fxs", "fasta-split", "fastq-split" ,"fa-split", "fq-split"])]
    FastxSplit {
        /// List of fastx files to write to.
        fastx: Vec<String>,
    },
    /// Mimic bedtools getfasta but allow for bgzip in both bed and fasta inputs.
    #[clap(visible_aliases = &["getfasta", "gf"])]
    GetFasta {
        /// Fasta file to extract sequences from.
        #[clap(short, long, default_value = "-")]
        fasta: String,
        /// BED file of regions to extract.
        #[clap(short, long)]
        bed: String,
        /// Reverse complement the sequence if the strand is "-".
        #[clap(short, long)]
        strand: bool,
        /// Add the name (4th column) to the header of the fasta output.
        #[clap(short, long)]
        name: bool,
    },
    /// Get the frequencies of each bp at each position.
    Nucfreq {
        /// Input sam/bam/cram/file.
        #[clap(default_value = "-")]
        bam: String,
        /// Print nucfreq info from the input region e.g "chr1:1-1000".
        #[clap(short, long)]
        region: Option<String>,
        /// Print nucfreq info from regions in the bed file
        /// output is optionally tagged using the 4th column.
        #[clap(short, long)]
        bed: Option<String>,
        /// Smaller output format.
        #[clap(short, long)]
        small: bool,
    },
    /// Report the longest exact repeat length at every position in a fasta.
    Repeat {
        /// Input fasta file.
        #[clap(default_value = "-")]
        fasta: String,
        /// The smallest repeat length to report.
        #[clap(short, long, default_value_t = 21)]
        min: usize,
    },
    /// Extract the intervals in a genome (fasta) that are made up of SUNs.
    Suns {
        /// Input fasta file with the genome.
        #[clap(short, long, default_value = "-")]
        fasta: String,
        /// The size of the required unique kmer.
        #[clap(short, long, default_value_t = 21)]
        kmer_size: usize,
        /// The maximum size SUN interval to report.
        #[clap(short, long, default_value_t = std::usize::MAX)]
        max_size: usize,
        /// Confirm all the SUNs (very slow) only for debugging.
        #[clap(short, long)]
        validate: bool,
    },
    /// Add RG lines from a source BAM file to the BAM from stdin to the BAM going to stdout
    AddRg {
        /// Source BAM file to read RG lines from
        #[clap(required = true)]
        source: String,

        /// Number of threads to use
        #[clap(short, long, default_value_t = 8)]
        threads: usize,

        /// Write uncompressed output
        #[clap(short = 'u', long)]
        uncompressed: bool,

        /// Add this string as the sample name (SM) to each read group in the output BAM
        #[clap(short = 's', long)]
        sample: Option<String>,
    },
    /// Calculate summary statistics from fasta/q, sam, bam, or bed files. e.g. N50, mean, quantiles.
    SeqStats {
        /// Input files (fast{a,q}(.gz), sam, bam, bed)
        #[clap(required = true)]
        infiles: Vec<String>,

        /// Number of threads to use
        #[clap(short, long, default_value_t = 4)]
        threads: usize,

        /// Print human-readable output
        #[clap(short = 'r', long)]
        human: bool,

        /// Quantiles to calculate
        #[clap(short, long, default_value = "0.5")]
        quantiles: Vec<f64>,

        /// Genome size for NG50 calculation
        #[clap(short = 'g', long)]
        genome_size: Option<usize>,
    },
}

pub fn make_cli_parse() -> Cli {
    Cli::parse()
}

pub fn make_cli_app() -> clap::Command<'static> {
    Cli::command()
}
