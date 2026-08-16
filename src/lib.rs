//! # Command line interface for rustybam
//! [rustybam command line interface, subcommands, and options.](cli::Commands)
//! # README for rustybam
#![doc = include_str!("../README.md")]
/// Annotation of bed files.
pub mod annotate;
/// Calculate stats from sam/bam/cram and paf files.
pub mod bamstats;
/// Bed file utilities.
pub mod bed;
/// Command line interface for rustybam.
pub mod cli;
/// Functions for fastx files.
pub mod fastx;
/// Mimic of `bedtools getfasta` that works with bgzipped files.
pub mod getfasta;
/// Liftover between assemblies using a PAF file.
pub mod liftover;
/// Module for automatically reading a writing compressed or uncompressed files.
pub mod myio;
/// Make statistics for NucFreq plotting.
pub mod nucfreq;
/// PAF file utilities.
pub mod paf;
/// Find SUNs within a fasta.
pub mod suns;
/// Functions for trimming query overlaps from a PAF file.
pub mod trim_overlap;

/// bed summary functions
pub mod bed_stats;

//pub mod test;
pub mod add_rg;
pub mod seq_content;
pub mod seq_stats;
