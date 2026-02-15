use super::bed;
use super::paf::*;
use colored::Colorize;
use itertools::Itertools;
use rayon::iter::ParallelBridge;
use rayon::prelude::*;
use rust_htslib::bam::record::Cigar::*;
use rust_htslib::bam::record::CigarString;
use std::cmp;
pub enum Error {
    PafParseCigar { msg: String },
    PafParseCS { msg: String },
    ParseIntError { msg: String },
    ParsePafColumn {},
}
//type LiftoverResult<T> = Result<T, crate::liftover::Error>;

/// Trim a PAF record to a target region by walking compressed CIGAR ops directly.
/// Complexity: O(n_cigar_ops) per call.
pub fn trim_paf_rec_to_rgn_fast(rgn: &bed::Region, paf: &PafRecord) -> Option<PafRecord> {
    let mut trimmed_paf = paf.small_copy();
    trimmed_paf.id = rgn.id.clone();

    // check if we can return right away
    if paf.t_st > rgn.st && paf.t_en < rgn.en {
        let mut tmp_rec = paf.clone();
        tmp_rec.id = rgn.id.clone();
        return Some(tmp_rec);
    }

    let trim_t_st = cmp::max(rgn.st, paf.t_st);
    let trim_t_en = cmp::min(rgn.en, paf.t_en);

    if trim_t_st >= trim_t_en {
        return None;
    }

    let cs_ref = &paf.cs_ops;
    let mut new_cs_ops: Option<CsOps> = cs_ref.as_ref().map(|cd| CsOps {
        ops: Vec::new(),
        seq_data: cd.seq_data.clone(),
    });

    let mut t_pos = paf.t_st;
    let mut q_before: u64 = 0;
    let mut q_in: u64 = 0;
    let mut new_cigar_ops: Vec<rust_htslib::bam::record::Cigar> = Vec::new();
    let mut started = false;
    let mut finished = false;
    let mut new_t_st = trim_t_st;
    let mut new_t_en = trim_t_en;

    for (ci, op) in paf.cigar.iter().enumerate() {
        if finished {
            break;
        }
        let op_len = op.len() as u64;
        let moves_t = consumes_reference(op);
        let moves_q = consumes_query(op);

        if !moves_t {
            if started {
                new_cigar_ops.push(*op);
                if let (Some(ref mut ncs), Some(ref cs)) = (&mut new_cs_ops, cs_ref) {
                    ncs.ops.push(cs.ops[ci]);
                }
                if moves_q {
                    q_in += op_len;
                }
            } else if moves_q {
                q_before += op_len;
            }
            continue;
        }

        let op_t_end = t_pos + op_len;
        if op_t_end <= trim_t_st {
            if moves_q {
                q_before += op_len;
            }
            t_pos = op_t_end;
            continue;
        }
        if t_pos >= trim_t_en {
            break;
        }

        let overlap_start = cmp::max(t_pos, trim_t_st);
        let overlap_end = cmp::min(op_t_end, trim_t_en);
        let skip_before = overlap_start - t_pos;
        let keep = overlap_end - overlap_start;

        if !started {
            started = true;
            new_t_st = overlap_start;
            if moves_q {
                q_before += skip_before;
            }
        }

        new_t_en = overlap_end;
        new_cigar_ops.push(update_cigar_opt_len(op, keep as u32));
        if let (Some(ref mut ncs), Some(ref cs)) = (&mut new_cs_ops, cs_ref) {
            ncs.ops.push(cs.ops[ci].trim(skip_before as u32, keep as u32));
        }
        if moves_q {
            q_in += keep;
        }

        if overlap_end >= trim_t_en {
            finished = true;
        }
        t_pos = op_t_end;
    }

    if new_cigar_ops.is_empty() {
        return None;
    }

    // Remove leading indels and adjust coordinates
    let mut lead_t_adjust: u64 = 0;
    let mut lead_q_adjust: u64 = 0;
    while !new_cigar_ops.is_empty() {
        let first = new_cigar_ops[0];
        if matches!(first, Match(_) | Equal(_) | Diff(_)) {
            break;
        }
        if consumes_reference(&first) {
            lead_t_adjust += first.len() as u64;
        }
        if consumes_query(&first) {
            lead_q_adjust += first.len() as u64;
        }
        new_cigar_ops.remove(0);
        if let Some(ref mut ncs) = new_cs_ops {
            ncs.ops.remove(0);
        }
    }
    new_t_st += lead_t_adjust;
    q_before += lead_q_adjust;
    q_in -= lead_q_adjust;

    // Remove trailing indels and adjust coordinates
    let mut trail_t_adjust: u64 = 0;
    let mut trail_q_adjust: u64 = 0;
    while !new_cigar_ops.is_empty() {
        let last = *new_cigar_ops.last().unwrap();
        if matches!(last, Match(_) | Equal(_) | Diff(_)) {
            break;
        }
        if consumes_reference(&last) {
            trail_t_adjust += last.len() as u64;
        }
        if consumes_query(&last) {
            trail_q_adjust += last.len() as u64;
        }
        new_cigar_ops.pop();
        if let Some(ref mut ncs) = new_cs_ops {
            ncs.ops.pop();
        }
    }
    new_t_en -= trail_t_adjust;
    q_in -= trail_q_adjust;

    if new_cigar_ops.is_empty() {
        return None;
    }

    // Check for match ops
    let has_match = new_cigar_ops
        .iter()
        .any(|op| matches!(op, Match(_) | Equal(_) | Diff(_)));
    if !has_match {
        return None;
    }

    // Set coordinates
    trimmed_paf.t_st = new_t_st;
    trimmed_paf.t_en = new_t_en;
    if paf.strand == '+' {
        trimmed_paf.q_st = paf.q_st + q_before;
        trimmed_paf.q_en = paf.q_st + q_before + q_in;
    } else {
        trimmed_paf.q_en = paf.q_en - q_before;
        trimmed_paf.q_st = paf.q_en - q_before - q_in;
    }

    trimmed_paf.cigar = CigarString(new_cigar_ops);
    trimmed_paf.cs_ops = new_cs_ops;

    if trimmed_paf.q_st > trimmed_paf.q_en || trimmed_paf.t_st > trimmed_paf.t_en {
        eprintln!(
            "Warning: liftover of {} failed. {} > {} or {} > {}",
            rgn, trimmed_paf.q_st, trimmed_paf.q_en, trimmed_paf.t_st, trimmed_paf.t_en
        );
        return None;
    }

    // check that this is still a valid paf record
    if let Err(e) = trimmed_paf.check_integrity() {
        eprintln!("WARNING: {:?}", e);
        return None;
    };

    Some(trimmed_paf)
}

pub fn trim_helper(name: &str, recs: &[PafRecord], rgns: &[bed::Region]) -> Vec<PafRecord> {
    let cur_recs: Vec<&PafRecord> = recs
        .par_iter()
        .filter(|rec| rec.t_name == name)
        .collect();
    let cur_rgns: Vec<&bed::Region> = rgns
        .par_iter()
        .filter(|rgn| rgn.name == name)
        .collect();

    let cur_trimmed_paf: Vec<PafRecord> = cur_recs
        .iter()
        .cartesian_product(cur_rgns) // make all pairwise combs
        .par_bridge()
        .filter(|(paf, rgn)| paf.paf_overlaps_rgn(rgn)) //filter to overlapping pairs
        .filter_map(|(paf, rgn)| trim_paf_rec_to_rgn_fast(rgn, paf))
        .collect();

    cur_trimmed_paf
}

pub fn trim_paf_by_rgns(
    rgns: &[bed::Region],
    paf_recs: &[PafRecord],
    invert_query: bool,
) -> Vec<PafRecord> {
    // swap query and ref if needed.
    let mut newvec = Vec::new();
    let recs: &[PafRecord] = if invert_query {
        for rec in paf_recs.iter() {
            newvec.push(paf_swap_query_and_target(rec));
        }
        &newvec
    } else {
        paf_recs
    };

    // define the unique contig names to looks at
    let names: Vec<&String> = recs.iter().map(|rec| &rec.t_name).unique().collect();
    let mut trimmed_paf = Vec::new();

    // iterate over one contigs name at a time
    for (idx, name) in names.iter().enumerate() {
        eprint!(
            "\rProcessing contig {}   {}/{}  ",
            name.bright_green().bold(),
            idx + 1,
            names.len()
        );
        let mut tmp = trim_helper(name, recs, rgns);
        trimmed_paf.append(&mut tmp);
    }
    eprintln!();
    trimmed_paf
}

/// # Example
/// ```
/// use rustybam::paf;
/// use rustybam::liftover;
/// let mut rec = paf::PafRecord::new("Q 10 0 10 + T 15 0 15 9 15 60 cg:Z:5=5D5=").unwrap();
/// let mut rec = paf::PafRecord::new("Q 15 0 15 + T 10 0 10 9 15 60 cg:Z:5=5I5=").unwrap();
/// let rec = paf::PafRecord::new("Q 15 0 15 - T 10 0 10 9 15 60 cg:Z:5=5I5=").unwrap();
/// for paf in liftover::break_paf_on_indels(&rec, 0){
///     eprintln!("{}", paf);
///     assert!(paf.t_en - paf.t_st == 5, "Incorrect size.");
/// }
/// ```
pub fn break_paf_on_indels(paf: &PafRecord, break_length: u32) -> Vec<PafRecord> {
    let mut rtn = Vec::new();
    let mut cur_tpos = paf.t_st;
    let mut pre_tpos = paf.t_st;
    for opt in paf.cigar.into_iter() {
        let opt_len = opt.len();
        if opt_len > break_length && matches!(opt, Del(_i) | Ins(_i)) {
            if cur_tpos > pre_tpos {
                let rgn = bed::Region {
                    name: paf.t_name.clone(),
                    st: pre_tpos,
                    en: cur_tpos,
                    id: paf.id.clone(),
                    ..Default::default()
                };
                if let Some(mut x) = trim_paf_rec_to_rgn_fast(&rgn, paf) {
                    x.check_integrity().unwrap();
                    rtn.push(x)
                }
            }
            pre_tpos = cur_tpos;
            if consumes_reference(opt) {
                pre_tpos += opt_len as u64;
            }
        }
        if consumes_reference(opt) {
            cur_tpos += opt_len as u64;
        }
    }
    if cur_tpos > pre_tpos {
        let rgn = bed::Region {
            name: paf.t_name.clone(),
            st: pre_tpos,
            en: cur_tpos,
            id: paf.id.clone(),
            ..Default::default()
        };
        if let Some(x) = trim_paf_rec_to_rgn_fast(&rgn, paf) {
            rtn.push(x)
        }
    }
    rtn
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bed::Region;

    #[test]
    fn test_liftover() {
        println!(
            "
            /// Example alignment
            /// 14-18         XXXXX
            /// 0123456789012345567890....
            /// ACTGACTGAAACTGAC-TAGA
            /// ------------||||I|D||
            ///             TGACGT-AC
            ///           01234567789 (forward)
            ///               XXXXX
            ///             98765433210 (reverse)
            "
        );
        let f_paf = PafRecord::new("Q 10 2 10 + T 40 12 20 3 9 60 cg:Z:4M1I1=1D2=").unwrap();
        let r_paf = PafRecord::new("Q 10 2 10 - T 40 12 20 3 9 60 cg:Z:4M1I1=1D2=").unwrap();

        let rgns = vec![
            Region {
                name: "T".to_string(),
                st: 14,
                en: 15,
                id: "None".to_string(),
                ..Default::default()
            },
            Region {
                name: "T".to_string(),
                st: 14,
                en: 18,
                id: "".to_string(),
                ..Default::default()
            },
            Region {
                name: "T".to_string(),
                st: 12,
                en: 20,
                id: "".to_string(),
                ..Default::default()
            },
            Region {
                name: "T".to_string(),
                st: 12,
                en: 30,
                id: "".to_string(),
                ..Default::default()
            },
            Region {
                name: "T".to_string(),
                st: 5,
                en: 20,
                id: "".to_string(),
                ..Default::default()
            },
            Region {
                name: "T".to_string(),
                st: 5,
                en: 30,
                id: "".to_string(),
                ..Default::default()
            },
        ];

        let sts = vec![4, 7, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2];
        let ens = vec![5, 8, 8, 8, 10, 10, 10, 10, 10, 10, 10, 10];
        let mut idx = 0;
        for r in &rgns {
            let trim = trim_paf_rec_to_rgn_fast(r, &f_paf).unwrap();
            eprintln!("fwd: {}", trim);
            assert_eq!(
                trim.q_st, sts[idx],
                "fwd q_st mismatch for rgn {}:{}-{} (idx {})",
                r.name, r.st, r.en, idx
            );
            assert_eq!(
                trim.q_en, ens[idx],
                "fwd q_en mismatch for rgn {}:{}-{} (idx {})",
                r.name, r.st, r.en, idx
            );
            idx += 1;

            let trim = trim_paf_rec_to_rgn_fast(r, &r_paf).unwrap();
            eprintln!("rev: {}", trim);
            assert_eq!(
                trim.q_st, sts[idx],
                "rev q_st mismatch for rgn {}:{}-{} (idx {})",
                r.name, r.st, r.en, idx
            );
            assert_eq!(
                trim.q_en, ens[idx],
                "rev q_en mismatch for rgn {}:{}-{} (idx {})",
                r.name, r.st, r.en, idx
            );
            idx += 1;
        }
    }

    #[test]
    fn test_break_paf() {
        let rec =
            PafRecord::new("Q 15 0 15 + T 10 0 10 9 15 60 cg:Z:5=5I5=").unwrap();
        let broken = break_paf_on_indels(&rec, 0);
        for paf in &broken {
            eprintln!("{}", paf);
            assert_eq!(paf.t_en - paf.t_st, 5, "Incorrect target size.");
        }
        assert_eq!(broken.len(), 2, "Should break into 2 records");

        // Reverse strand
        let rec =
            PafRecord::new("Q 15 0 15 - T 10 0 10 9 15 60 cg:Z:5=5I5=").unwrap();
        let broken = break_paf_on_indels(&rec, 0);
        for paf in &broken {
            eprintln!("{}", paf);
            assert_eq!(paf.t_en - paf.t_st, 5, "Incorrect target size.");
        }
        assert_eq!(broken.len(), 2, "Should break into 2 records");

        // Test with deletion break
        let rec =
            PafRecord::new("Q 10 0 10 + T 15 0 15 9 15 60 cg:Z:5=5D5=").unwrap();
        let broken = break_paf_on_indels(&rec, 0);
        for paf in &broken {
            eprintln!("{}", paf);
            assert_eq!(paf.t_en - paf.t_st, 5, "Incorrect target size.");
        }
        assert_eq!(broken.len(), 2, "Should break into 2 records");
    }
}
