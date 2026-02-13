use super::paf::*;
use log;
use rust_htslib::bam::record::Cigar::*;
use std::cmp::{max, min};

/// Walk compressed CIGAR to compute per-position scores for query range [q_start, q_end).
/// O(n_cigar_ops + overlap_length) — no expansion to 1-bp resolution.
fn compute_overlap_scores(
    paf: &PafRecord,
    q_start: u64,
    q_end: u64,
    match_score: i32,
    diff_score: i32,
    indel_score: i32,
) -> Vec<i32> {
    let overlap_len = (q_end - q_start) as usize;
    let mut scores = vec![0i32; overlap_len];

    let mut q_pos = if paf.strand == '+' { paf.q_st } else { paf.q_en };

    for op in paf.cigar.into_iter() {
        let op_len = op.len() as u64;
        let moves_q = consumes_query(op);

        if !moves_q {
            continue;
        }

        // Compute absolute query range [q_lo, q_hi) for this op
        let (q_lo, q_hi) = if paf.strand == '+' {
            let lo = q_pos;
            q_pos += op_len;
            (lo, q_pos)
        } else {
            q_pos -= op_len;
            (q_pos, q_pos + op_len)
        };

        // Compute overlap with [q_start, q_end)
        let ovl_st = max(q_lo, q_start);
        let ovl_en = min(q_hi, q_end);
        if ovl_st >= ovl_en {
            continue;
        }

        // Determine score for this op type
        let score = match op {
            Equal(_) => match_score,
            Ins(_) | Del(_) => -indel_score,
            _ => -diff_score, // Match, Diff
        };

        // Fill scores for the overlapping portion
        let idx_start = (ovl_st - q_start) as usize;
        let idx_end = (ovl_en - q_start) as usize;
        for s in scores[idx_start..idx_end].iter_mut() {
            *s = score;
        }
    }

    scores
}

/// Example
/// ```
/// use rustybam::trim_overlap::*;
/// use rustybam::paf::*;
///
/// let mut left = PafRecord::new("Q 10 0 10 + T 20 0 10 3 9 60 cg:Z:7=1X2=").unwrap();
/// let mut right = PafRecord::new("Q 10 5 10 - T 20 10 15 3 9 60 cg:Z:3=1X1=").unwrap();
///
/// trim_overlapping_pafs(&mut left, &mut right, 1 ,1 ,1);
///
/// assert_eq!(left.cigar.to_string(), "7=");
/// assert_eq!(right.cigar.to_string(), "3=");
/// ```
pub fn trim_overlapping_pafs(
    left: &mut PafRecord,
    right: &mut PafRecord,
    match_score: i32,
    diff_score: i32,
    indel_score: i32,
) {
    let st_ovl = max(left.q_st, right.q_st);
    let en_ovl = min(left.q_en, right.q_en);
    log::info!("Number of overlapping bases {}", en_ovl - st_ovl);

    let l_scores = compute_overlap_scores(left, st_ovl, en_ovl, match_score, diff_score, indel_score);
    let r_scores = compute_overlap_scores(right, st_ovl, en_ovl, match_score, diff_score, indel_score);

    // Build cumulative scores: l_score has a leading 0, r_score has a trailing 0
    let mut l_score = vec![0i32];
    l_score.extend_from_slice(&l_scores);
    let mut r_score = r_scores;
    r_score.push(0);

    // make cumulative from left
    l_score.iter_mut().fold(0, |acc, x| {
        *x += acc;
        *x
    });
    // make cumulative from right
    r_score.iter_mut().rev().fold(0, |acc, x| {
        *x += acc;
        *x
    });

    // determine the split point
    let mut max_idx: u64 = 0;
    let mut max = 0;
    for (idx, (l, r)) in l_score.iter().zip(r_score.iter()).enumerate() {
        if l + r > max {
            max = l + r;
            max_idx = idx as u64;
        }
    }
    left.truncate_record_by_query(left.q_st, st_ovl + max_idx);
    right.truncate_record_by_query(st_ovl + max_idx, right.q_en);

    log::info!(
        "Split position was found to be {} bases after the overlap start ({}) with a score of {}.",
        max_idx,
        st_ovl,
        max
    );
}

/// Example
/// ```
/// use rustybam::trim_overlap::*;
/// use rustybam::paf::*;
///
/// ```

#[cfg(test)]
mod tests {
    use super::*;

    fn pretty_print_paf(paf: &Paf) {
        println!("\nPAF RECORDS");
        for rec in &paf.records {
            println!(
                "  {}:{}-{} ({}) -> {}:{}-{} cigar={}",
                rec.q_name, rec.q_st, rec.q_en, rec.strand,
                rec.t_name, rec.t_st, rec.t_en, rec.cigar
            );
        }
    }

    #[test]
    fn test_inversion_trimming() {
        let mut left = PafRecord::new("Q 20 0 10 + T 20 0 10 3 9 60 cg:Z:7=1X2=").unwrap();
        left.check_integrity().unwrap();

        let mut center = PafRecord::new("Q 20 4 15 - T 20 5 16 3 9 60 cg:Z:3=1X3=1M1X2=").unwrap();
        center.check_integrity().unwrap();

        let mut right =
            PafRecord::new("Q 20 10 20 + T 20 10 20 3 9 60 cz:Z:10= cg:Z:2=2X2=2X2=").unwrap();
        right.check_integrity().unwrap();

        let mut paf = Paf::new();
        paf.records = vec![left, center, right];

        pretty_print_paf(&paf);
        paf.overlapping_paf_recs(1, 1, 1, false);
        pretty_print_paf(&paf);

        let expected_cigars = vec!["7=", "2=1X3=1M", "2=2X2="];
        for (i, rec) in paf.records.iter().enumerate() {
            assert_eq!(rec.cigar.to_string(), expected_cigars[i]);
        }
    }
}
