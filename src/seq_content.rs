use super::bed;
use rust_htslib::faidx;
use std::collections::HashMap;

/// Represents k-mer counts for a specific region category
#[derive(Debug, Clone)]
pub struct KmerCounts {
    pub name: String,
    pub kmer_counts: HashMap<String, u64>,
}

/// Generate reverse complement of a DNA sequence
fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => c,
        })
        .collect()
}

/// Get the canonical form of a k-mer (lexicographically smaller of k-mer and its reverse complement)
fn canonical_kmer(kmer: &str) -> String {
    let rc = reverse_complement(kmer);
    if kmer <= rc.as_str() {
        kmer.to_string()
    } else {
        rc
    }
}

/// Generate all possible canonical k-mers of given length using DNA alphabet
fn generate_all_kmers(k: usize) -> Vec<String> {
    let bases = ['A', 'C', 'G', 'T'];
    let mut kmers = vec!["".to_string()];

    for _ in 0..k {
        let mut new_kmers = Vec::new();
        for kmer in &kmers {
            for &base in &bases {
                new_kmers.push(format!("{}{}", kmer, base));
            }
        }
        kmers = new_kmers;
    }

    // Convert to canonical form and deduplicate
    let mut canonical_kmers: Vec<String> = kmers
        .into_iter()
        .map(|kmer| canonical_kmer(&kmer))
        .collect();

    canonical_kmers.sort();
    canonical_kmers.dedup();
    canonical_kmers
}

/// Extract k-mers from a DNA sequence using canonical form
fn extract_kmers(sequence: &[u8], k: usize) -> HashMap<String, u64> {
    let mut kmer_counts = HashMap::new();

    if sequence.len() < k {
        return kmer_counts;
    }

    for i in 0..=(sequence.len() - k) {
        let kmer_bytes = &sequence[i..i + k];

        // Convert to uppercase and check if all bases are valid DNA
        let kmer_str: Result<String, _> = kmer_bytes
            .iter()
            .map(|&b| match b.to_ascii_uppercase() {
                b'A' | b'C' | b'G' | b'T' => Ok(b.to_ascii_uppercase() as char),
                _ => Err(()),
            })
            .collect();

        if let Ok(kmer) = kmer_str {
            // Convert to canonical form (merge with reverse complement)
            let canonical = canonical_kmer(&kmer);
            *kmer_counts.entry(canonical).or_insert(0) += 1;
        }
        // Skip k-mers with invalid characters (N, etc.)
    }

    kmer_counts
}

/// Count k-mers in genomic regions specified by BED file
pub fn count_kmers_in_regions(
    fasta_path: &str,
    bed_path: &str,
    k: usize,
) -> Result<Vec<KmerCounts>, Box<dyn std::error::Error>> {
    // Parse BED file
    let bed_regions = bed::parse_bed(bed_path);

    // Open FASTA file
    let fasta_reader = faidx::Reader::from_path(fasta_path)?;

    // Group regions by chromosome, then by region name (4th column)
    let mut regions_by_chrom: HashMap<String, Vec<&bed::Region>> = HashMap::new();
    for region in &bed_regions {
        regions_by_chrom
            .entry(region.name.clone())
            .or_default()
            .push(region);
    }

    // Initialize results map for each unique region name
    let mut results_map: HashMap<String, HashMap<String, u64>> = HashMap::new();
    for region in &bed_regions {
        results_map.entry(region.id.clone()).or_default();
    }

    // Process each chromosome
    for (chrom, regions) in regions_by_chrom {
        log::info!("Processing chromosome: {}", chrom);

        // Fetch each region on its own. This avoids one whole-chromosome
        // buffer per contig.
        for region in regions {
            let start = region.st as usize;
            let end = region.en as usize;

            if start >= end {
                log::warn!("Invalid region coordinates: {}:{}-{}", chrom, start, end);
                continue;
            }

            // fetch_seq takes an inclusive end position. rust-htslib 1.0
            // returns an owned Vec and frees the htslib buffer itself.
            let region_sequence = fasta_reader.fetch_seq(&chrom, start, end - 1)?;

            // Count k-mers in this region
            let region_kmers = extract_kmers(&region_sequence, k);

            // Add to combined counts for this region name
            let combined_counts = results_map.get_mut(&region.id).unwrap();
            for (kmer, count) in region_kmers {
                *combined_counts.entry(kmer).or_insert(0) += count;
            }
        }
    }

    // Convert results map to vector
    let results: Vec<KmerCounts> = results_map
        .into_iter()
        .map(|(name, kmer_counts)| KmerCounts { name, kmer_counts })
        .collect();

    Ok(results)
}

/// Print k-mer counts in a tabular format
pub fn print_kmer_results(results: &[KmerCounts], k: usize) {
    if results.is_empty() {
        return;
    }

    // Generate all possible k-mers for consistent column ordering
    let all_kmers = generate_all_kmers(k);

    // Sort results by region name alphabetically
    let mut sorted_results = results.to_vec();
    sorted_results.sort_by(|a, b| a.name.cmp(&b.name));

    // Print header
    print!("name\ttotal_kmers");
    for kmer in &all_kmers {
        print!("\t{}", kmer);
    }
    println!();

    // Print results for each region name in alphabetical order
    for result in &sorted_results {
        // Calculate total k-mer count
        let total_count: u64 = result.kmer_counts.values().sum();

        print!("{}\t{}", result.name, total_count);
        for kmer in &all_kmers {
            let count = result.kmer_counts.get(kmer).unwrap_or(&0);
            print!("\t{}", count);
        }
        println!();
    }
}

/// Main function to run seq-content analysis
pub fn run_seq_content(
    fasta_path: &str,
    bed_path: &str,
    k: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    log::info!("Analyzing {}-mer content in regions from {}", k, bed_path);
    log::info!("Using FASTA file: {}", fasta_path);

    let results = count_kmers_in_regions(fasta_path, bed_path, k)?;

    log::info!("Found {} unique region names", results.len());

    print_kmer_results(&results, k);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("ATCG"), "CGAT");
        assert_eq!(reverse_complement("AAAA"), "TTTT");
        assert_eq!(reverse_complement("GCGC"), "GCGC");
    }

    #[test]
    fn test_canonical_kmer() {
        assert_eq!(canonical_kmer("AAT"), "AAT"); // AAT < ATT
        assert_eq!(canonical_kmer("ATT"), "AAT"); // ATT > AAT, so return AAT
        assert_eq!(canonical_kmer("GGG"), "CCC"); // GGG > CCC, so return CCC
        assert_eq!(canonical_kmer("ACG"), "ACG"); // ACG < CGT, so return ACG
        assert_eq!(canonical_kmer("CGT"), "ACG"); // CGT > ACG, so return ACG
    }

    #[test]
    fn test_generate_all_kmers() {
        let kmers_2 = generate_all_kmers(2);
        assert_eq!(kmers_2.len(), 10); // 16 - 6 pairs merged (AA, AC, AG, AT, CC, CG remain)
        assert!(kmers_2.contains(&"AA".to_string()));
        assert!(!kmers_2.contains(&"TT".to_string())); // TT should be merged with AA
        assert!(kmers_2.contains(&"AC".to_string()));

        let kmers_1 = generate_all_kmers(1);
        assert_eq!(kmers_1.len(), 2); // A and C (T->A, G->C)
        assert!(kmers_1.contains(&"A".to_string()));
        assert!(kmers_1.contains(&"C".to_string()));
        assert!(!kmers_1.contains(&"G".to_string())); // G should be merged with C
        assert!(!kmers_1.contains(&"T".to_string())); // T should be merged with A
    }

    #[test]
    fn test_extract_kmers() {
        let sequence = b"ATCGATCG";
        let kmers = extract_kmers(sequence, 3);

        // With canonical k-mers:
        // ATC (canonical: ATC, since ATC < GAT), TCG (canonical: CGA, since TCG > CGA),
        // CGA (canonical: CGA), GAT (canonical: ATC, since GAT > ATC),
        // ATC (canonical: ATC), TCG (canonical: CGA)
        assert_eq!(kmers.get("ATC").unwrap_or(&0), &3); // ATC appears 2 times + GAT->ATC 1 time
        assert_eq!(kmers.get("CGA").unwrap_or(&0), &3); // CGA appears 1 time + TCG->CGA 2 times
    }

    #[test]
    fn test_extract_kmers_with_invalid_chars() {
        let sequence = b"ATCNATCG";
        let kmers = extract_kmers(sequence, 3);

        // Should skip k-mers containing N
        assert_eq!(kmers.get("TCN").unwrap_or(&0), &0);
        assert_eq!(kmers.get("CNA").unwrap_or(&0), &0);
        assert_eq!(kmers.get("NAT").unwrap_or(&0), &0);

        // Should count valid k-mers in canonical form
        // ATC appears 2 times (positions 0 and 5)
        // TCG appears 1 time (position 6) -> canonical CGA
        assert_eq!(kmers.get("ATC").unwrap_or(&0), &2);
        assert_eq!(kmers.get("CGA").unwrap_or(&0), &1);
    }
}
