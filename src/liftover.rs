use super::bed;
use super::paf::*;
use colored::Colorize;
use itertools::Itertools;
use rayon::iter::ParallelBridge;
use rayon::prelude::*;
use rust_htslib::bam::record::Cigar::*;
use std::cmp;
pub enum Error {
    PafParseCigar { msg: String },
    PafParseCS { msg: String },
    ParseIntError { msg: String },
    ParsePafColumn {},
}
//type LiftoverResult<T> = Result<T, crate::liftover::Error>;

pub fn trim_paf_rec_to_rgn(rgn: &bed::Region, paf: &PafRecord) -> Option<PafRecord> {
    // initialize a trimmed paf record
    let mut trimmed_paf = paf.small_copy();
    trimmed_paf.id = rgn.id.clone();

    // check if we can return right away
    if paf.t_st > rgn.st && paf.t_en < rgn.en {
        let mut tmp_rec = paf.clone();
        tmp_rec.id = rgn.id.clone();
        return Some(tmp_rec);
    }

    // index at the start of trimmed alignment
    trimmed_paf.t_st = cmp::max(rgn.st, paf.t_st);
    let start_idx = match paf.tpos_to_idx_match(trimmed_paf.t_st, true) {
        Ok(idx) => idx,
        Err(_) => panic!(
            "\nProblem getting index in cigar:\n{}\n{}\n{}\n",
            trimmed_paf.t_st, rgn, paf
        ),
    };

    // index at the end of trimmed alignment
    trimmed_paf.t_en = cmp::min(rgn.en, paf.t_en);
    // the end index is not inclusive so minus 1
    let end_idx = match paf.tpos_to_idx_match(trimmed_paf.t_en - 1, false) {
        Ok(idx) => idx,
        Err(error_idx) => panic!(
            "\nProblem getting index in cigar:\n{}\n{:?}\n{}\n{}\n",
            trimmed_paf.t_en - 1,
            &paf.tpos_aln[error_idx - 1],
            rgn,
            paf
        ),
    };

    // if this is the case then the whole internal regions is an indel
    if start_idx > end_idx {
        return None;
    }

    // update the start and end positions
    trimmed_paf.t_st = paf.tpos_aln[start_idx];
    trimmed_paf.q_st = paf.qpos_aln[start_idx];
    trimmed_paf.t_en = paf.tpos_aln[end_idx];
    trimmed_paf.q_en = paf.qpos_aln[end_idx];

    // get the cigar opts
    trimmed_paf.cigar = PafRecord::collapse_long_cigar(&paf.subset_cigar(start_idx, end_idx));

    // make sure there are opts in the cigar that are not indels
    let mut no_match_opts = true;
    for opt in &trimmed_paf.cigar {
        if matches!(opt, Match(_) | Equal(_) | Diff(_)) {
            no_match_opts = false;
            break;
        }
    }
    if no_match_opts {
        return None;
    }

    if paf.strand == '-' {
        std::mem::swap(&mut trimmed_paf.q_en, &mut trimmed_paf.q_st);
    }
    // make end index not inclusive
    trimmed_paf.t_en += 1;
    trimmed_paf.q_en += 1;

    // trim indels from start and end
    trimmed_paf.remove_trailing_indels();

    if trimmed_paf.cigar.len() == 0 {
        return None;
    }
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
    let mut cur_recs: Vec<PafRecord> = recs
        .into_par_iter()
        .filter(|rec| rec.t_name == name)
        .map(|paf| (*paf).clone()) // clone the records so we can change them
        .collect();
    let cur_rgns: Vec<&bed::Region> = rgns
        .into_par_iter()
        .filter(|rgn| rgn.name == name)
        .collect();

    // calculated base by base liftover chain for each paf record
    cur_recs
        .par_iter_mut()
        .for_each(|paf| (paf).aligned_pairs());

    let cur_trimmed_paf: Vec<PafRecord> = cur_recs
        .iter()
        .cartesian_product(cur_rgns) // make all pairwise combs
        .par_bridge()
        .filter(|(paf, rgn)| paf.paf_overlaps_rgn(rgn)) //filter to overlapping pairs
        .filter_map(|(paf, rgn)| trim_paf_rec_to_rgn(rgn, paf))
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
/// let mut rec = paf::PafRecord::new("Q 15 0 15 - T 10 0 10 9 15 60 cg:Z:5=5I5=").unwrap();
/// rec.aligned_pairs();
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
                //rtn.push(trim_paf_rec_to_rgn(&rgn, paf));
                if let Some(mut x) = trim_paf_rec_to_rgn(&rgn, paf) {
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
        //eprintln!("{} {} {:?}", pre_tpos, cur_tpos, opt);
    }
    if cur_tpos > pre_tpos {
        let rgn = bed::Region {
            name: paf.t_name.clone(),
            st: pre_tpos,
            en: cur_tpos,
            id: paf.id.clone(),
            ..Default::default()
        };
        if let Some(x) = trim_paf_rec_to_rgn(&rgn, paf) {
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
    fn test_aln_pair_liftover() {
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
        let mut f_paf = PafRecord::new("Q 10 2 10 + T 40 12 20 3 9 60 cg:Z:4M1I1=1D2=").unwrap();
        f_paf.aligned_pairs();
        let mut r_paf = PafRecord::new("Q 10 2 10 - T 40 12 20 3 9 60 cg:Z:4M1I1=1D2=").unwrap();
        r_paf.aligned_pairs();

        let rgn = Region {
            name: "T".to_string(),
            st: 14,
            en: 15,
            id: "None".to_string(),
            ..Default::default()
        };
        let rgn2 = Region {
            name: "T".to_string(),
            st: 14,
            en: 18,
            id: "".to_string(),
            ..Default::default()
        };
        let rgn3 = Region {
            name: "T".to_string(),
            st: 12,
            en: 20,
            id: "".to_string(),
            ..Default::default()
        };
        // test right extend
        let rgn4 = Region {
            name: "T".to_string(),
            st: 12,
            en: 30,
            id: "".to_string(),
            ..Default::default()
        };
        // test left extend
        let rgn5 = Region {
            name: "T".to_string(),
            st: 5,
            en: 20,
            id: "".to_string(),
            ..Default::default()
        };
        // test both extend
        let rgn6 = Region {
            name: "T".to_string(),
            st: 5,
            en: 30,
            id: "".to_string(),
            ..Default::default()
        };

        let sts = vec![4, 7, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2];
        let ens = vec![5, 8, 8, 8, 10, 10, 10, 10, 10, 10, 10, 10];
        let mut idx = 0;
        for r in [rgn, rgn2, rgn3, rgn4, rgn5, rgn6] {
            let trim = trim_paf_rec_to_rgn(&r, &f_paf).unwrap();
            eprintln!("{}", trim);
            eprintln!("{:?}", f_paf.tpos_aln);
            eprintln!("{:?}", f_paf.qpos_aln);
            eprintln!("{:?}", f_paf.long_cigar.to_string());
            assert_eq!(trim.q_st, sts[idx]);
            assert_eq!(trim.q_en, ens[idx]);
            idx += 1;

            eprintln!();
            let trim = trim_paf_rec_to_rgn(&r, &r_paf).unwrap();
            eprintln!("{}", trim);
            eprintln!("{}", trim_paf_rec_to_rgn(&r, &r_paf).unwrap());
            eprintln!("{:?}", r_paf.tpos_aln);
            eprintln!("{:?}", r_paf.qpos_aln);
            eprintln!("{:?}", f_paf.long_cigar.to_string());
            assert_eq!(trim.q_st, sts[idx]);
            assert_eq!(trim.q_en, ens[idx]);
            idx += 1;

            eprintln!("\n");
        }
    }

    #[test]
    /// this function tests that if we subset the ref and then invert that paf
    /// we can get back the original ref subset.
    fn check_invertible() {
        /*
        let mut paf = make_fake_paf_rec();
        paf.aligned_pairs();
        let rgn = Region {
            name: "T".to_string(),
            st: 14,
            en: 17,
            id: "None".to_string(),
        };
        let trim_paf = trim_paf_rec_to_rgn(&rgn, &paf);
        eprintln!("{}", trim_paf);
        let q_rgn = Region {
            name: trim_paf.q_name.clone(),
            st: trim_paf.q_st,
            en: trim_paf.q_en,
            id: "None".to_string(),
        };

        let mut trim_paf_swap = paf_swap_query_and_target(&trim_paf);
        trim_paf_swap.aligned_pairs();
        let ref_paf = paf_swap_query_and_target(&trim_paf_rec_to_rgn(&q_rgn, &trim_paf_swap));
        // make sure we recreated what we wanted
        assert_eq!(ref_paf.t_st, trim_paf.t_st);
        assert_eq!(ref_paf.t_en, trim_paf.t_en);
        assert_eq!(ref_paf.q_st, trim_paf.q_st);
        assert_eq!(ref_paf.q_en, trim_paf.q_en);
        assert_eq!(ref_paf.cigar, trim_paf.cigar);
        eprintln!("{}", ref_paf);
        */
    }
}
