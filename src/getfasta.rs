use super::*;
use bio::alphabets::dna::revcomp;
use bio_types::strand::Strand::{Forward, Reverse};
use rust_htslib::faidx;
use std::str;

pub fn fetch_fasta(reader: &faidx::Reader, chrom: &str, start: usize, end: usize) -> Vec<u8> {
    // rust-htslib 1.0 fetch_seq returns an owned Vec and frees the
    // htslib buffer itself. No leak, no to_owned copy.
    reader.fetch_seq(chrom, start, end).unwrap()
}
/// # Example
/// ```
/// use rustybam::getfasta::get_fasta;
/// get_fasta(".test/test.fa", ".test/getfasta.bed", true, true);
/// get_fasta(".test/test.fa", ".test/getfasta.bed", false, true);
/// get_fasta(".test/test.fa", ".test/getfasta.bed", true, false);
/// get_fasta(".test/test.fa", ".test/getfasta.bed", false, false);
/// get_fasta(".test/test.fa.gz", ".test/getfasta.bed.gz", true, false);
/// ```
pub fn get_fasta(path: &str, bed: &str, add_name: bool, use_strand: bool) {
    let bed_recs = bed::parse_bed(bed);
    let reader = faidx::Reader::from_path(path).unwrap();

    // iterate over be records
    for rec in bed_recs {
        let mut name = format!("{}:{}-{}", rec.name, rec.st, rec.en);
        let mut seq = fetch_fasta(&reader, &rec.name, rec.st as usize, rec.en as usize);
        // add name information if present
        if add_name {
            if let Some(n) = rec.record.name() {
                name = format!("{}::{}", n, name)
            }
        }

        // add the stand information if present
        if use_strand {
            match rec.record.strand() {
                Some(strand) => {
                    match strand {
                        Reverse => {
                            seq = revcomp(seq);
                            name.push_str("(-)");
                        }
                        Forward => {
                            name.push_str("(+)");
                        }
                        _ => {
                            name.push_str("(.)");
                        }
                    };
                }
                None => name.push_str("(.)"),
            };
        }
        println!(">{}\n{}", name, str::from_utf8(&seq).unwrap());
    }
}
