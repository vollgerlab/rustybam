use super::bed;
use needletail::parse_fastx_file;
use num_format::{Locale, ToFormattedString};
use rayon::prelude::*;
use rust_htslib::bam::{self, Read};
use std::fs;
use std::io::{self, BufRead};
use std::path::Path;

fn read_bam(file: &str, threads: usize) -> Option<Vec<usize>> {
    let mut lengths = Vec::new();

    let mut bam = bam::Reader::from_path(file).ok()?;
    bam.set_threads(threads).ok()?;
    for record in bam.records() {
        let rec = record.ok()?;
        if rec.is_unmapped() || (!rec.is_secondary() && !rec.is_supplementary()) {
            lengths.push(rec.seq().len());
        }
    }
    eprintln!("SAM/BAM read: {}", file);
    Some(lengths)
}

fn read_bed(file: &str) -> Option<Vec<usize>> {
    let regions = bed::parse_bed(file);
    if regions.is_empty() {
        return None;
    }
    Some(regions.iter().map(|r| (r.en - r.st) as usize).collect())
}

/// Read sequence lengths straight out of a `.fai` index, without touching the sequence data.
fn read_fai(fai_path: &str) -> Option<Vec<usize>> {
    let file = fs::File::open(fai_path).ok()?;
    let reader = io::BufReader::new(file);
    let mut lengths = Vec::new();
    for line in reader.lines() {
        let line = line.ok()?;
        let len: usize = line.split('\t').nth(1)?.parse().ok()?;
        lengths.push(len);
    }
    if lengths.is_empty() {
        None
    } else {
        Some(lengths)
    }
}

fn read_fasta_or_fastx(file: &str) -> Option<Vec<usize>> {
    let fai_path = format!("{}.fai", file);
    if Path::new(&fai_path).exists() {
        if let Some(lengths) = read_fai(&fai_path) {
            eprintln!("Index read: {}", file);
            return Some(lengths);
        }
    }

    let mut reader = parse_fastx_file(file).ok()?;
    let mut lengths = Vec::new();
    while let Some(record) = reader.next() {
        match record {
            Ok(rec) => lengths.push(rec.num_bases()),
            Err(e) => log::warn!("Error reading record in {}: {}", file, e),
        }
    }
    eprintln!("Fastx read: {}", file);
    Some(lengths)
}

fn get_lengths(file: &str, threads: usize) -> Option<Vec<usize>> {
    if file.ends_with(".bam") || file.ends_with(".sam") || file.ends_with(".cram") {
        log::info!("Reading BAM/SAM/CRAM file: {}", file);
        read_bam(file, threads)
    } else if file.ends_with(".bed") || file.ends_with(".bed.gz") {
        log::info!("Reading BED file: {}", file);
        read_bed(file)
    } else {
        log::info!("Reading fasta/fastq file: {}", file);
        read_fasta_or_fastx(file)
    }
}

/// Linear-interpolation quantile, matching numpy's default `np.quantile` method.
fn quantile(sorted_asc: &[usize], q: f64) -> f64 {
    let n = sorted_asc.len();
    if n == 0 {
        return 0.0;
    }
    let h = q * (n - 1) as f64;
    let lo = h.floor() as usize;
    let hi = (h.ceil() as usize).min(n - 1);
    let lo_v = sorted_asc[lo] as f64;
    lo_v + (h - lo as f64) * (sorted_asc[hi] as f64 - lo_v)
}

fn calc_stats(
    lengths: &[usize],
    quantiles: &[f64],
    genome_size: Option<usize>,
) -> (usize, usize, f64, Vec<f64>, usize, usize, usize, f64) {
    let n = lengths.len();
    let total: usize = genome_size.unwrap_or_else(|| lengths.iter().sum());
    let mut asc = lengths.to_vec();
    asc.sort_unstable();

    let min = *asc.first().unwrap_or(&0);
    let max = *asc.last().unwrap_or(&0);
    let mean = total as f64 / n as f64;

    let au_n: f64 = asc.iter().map(|&x| (x * x) as f64).sum::<f64>() / total as f64;

    let quantile_values: Vec<f64> = quantiles.iter().map(|&q| quantile(&asc, q)).collect();

    let mut cumulative = 0;
    let mut n50 = 0;
    for &len in asc.iter().rev() {
        cumulative += len;
        if cumulative >= total / 2 {
            n50 = len;
            break;
        }
    }

    (total, n, mean, quantile_values, min, max, n50, au_n)
}

pub fn h_fmt<T>(num: T) -> String
where
    T: Into<f64> + Copy,
{
    let mut num: f64 = num.into();
    for unit in ["", "Kbp", "Mbp"] {
        if num < 1000.0 {
            return format!("{:.2}{}", num, unit);
        }
        num /= 1000.0;
    }
    format!("{:.2}{}", num, "Gbp")
}

fn quantile_header_label(q: f64) -> String {
    // Round to one decimal so 0.333 prints as 33.3%, not 33.300000000000004%.
    let pct = (q * 1000.0).round() / 10.0;
    format!("{}%", pct)
}

fn process_file(
    file: &str,
    threads: usize,
    human_readable: bool,
    quantiles: &[f64],
    genome_size: Option<usize>,
) -> Option<String> {
    let lengths = get_lengths(file, threads);

    let Some(lengths) = lengths else {
        eprintln!("Skipping file: {}", file);
        return None;
    };

    let (total, n, mean, quantile_values, min, max, n50, au_n) =
        calc_stats(&lengths, quantiles, genome_size);

    let quantile_str = quantile_values
        .iter()
        .map(|q| {
            if human_readable {
                h_fmt(*q)
            } else {
                format!("{:.2}", q)
            }
        })
        .collect::<Vec<_>>()
        .join("\t");

    let line = if human_readable {
        format!(
            "{}\t{}\t{}\t{:}\t{}\t{}\t{}\t{}\t{:}\n",
            file,
            h_fmt(total as f64),
            n.to_formatted_string(&Locale::en),
            h_fmt(mean),
            quantile_str,
            h_fmt(min as f64),
            h_fmt(max as f64),
            h_fmt(n50 as f64),
            h_fmt(au_n)
        )
    } else {
        format!(
            "{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{:.2}\n",
            file, total, n, mean, quantile_str, min, max, n50, au_n
        )
    };
    Some(line)
}

pub fn seq_stats(
    infiles: &[String],
    threads: usize,
    human_readable: bool,
    quantiles: &[f64],
    genome_size: Option<usize>,
) {
    let infiles: Vec<&String> = infiles
        .iter()
        .filter(|f| {
            // Keep FIFOs and process substitutions. Skip only files
            // that are missing, or regular and empty.
            let exists_nonempty = fs::metadata(f.as_str())
                .map(|m| !m.is_file() || m.len() > 0)
                .unwrap_or(false);
            if !exists_nonempty {
                eprintln!("Skipping, because missing or empty: {}", f);
            }
            exists_nonempty
        })
        .collect();

    let quantile_headers = quantiles
        .iter()
        .map(|&q| quantile_header_label(q))
        .collect::<Vec<_>>()
        .join("\t");
    let mut output = format!(
        "file\ttotalBp\tnSeqs\tmean\t{}\tmin\tmax\tN50\tauN\n",
        quantile_headers
    );

    let rows: Vec<String> = infiles
        .par_iter()
        .filter_map(|file| process_file(file, threads, human_readable, quantiles, genome_size))
        .collect();

    for row in rows {
        output.push_str(&row);
    }

    print!("{}", output);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quantile_matches_numpy_linear() {
        assert_eq!(quantile(&[1, 2, 3, 4], 0.5), 2.5);
        assert_eq!(quantile(&[1, 2, 3, 4], 0.0), 1.0);
        assert_eq!(quantile(&[1, 2, 3, 4], 1.0), 4.0);
    }

    #[test]
    fn reads_fasta_via_fai_fast_path() {
        let lengths = read_fasta_or_fastx(".test/test.fa").expect("should read test.fa");
        assert_eq!(lengths.len(), 2);
    }

    #[test]
    fn reads_fasta_without_fai_via_needletail() {
        let dir = std::env::temp_dir();
        let path = dir.join("rustybam_seq_stats_test.fa");
        fs::write(&path, b">a\nACGT\n>b\nACGTACGT\n").unwrap();
        let lengths =
            read_fasta_or_fastx(path.to_str().unwrap()).expect("should read fasta via needletail");
        assert_eq!(lengths, vec![4, 8]);
        fs::remove_file(&path).ok();
    }
}
