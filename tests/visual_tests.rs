//! Visual test cases for rustybam's liftover, trim, and break-paf operations.
//!
//! Each test contains an ASCII diagram showing the alignment, the region being
//! operated on, and the expected result. This makes it easy to verify correctness
//! by visual inspection.

use rustybam::bed::Region;
use rustybam::liftover::{break_paf_on_indels, trim_paf_rec_to_rgn_fast};
use rustybam::paf::PafRecord;

// =============================================================================
// TEST ALIGNMENT 1: Forward strand with insertion and deletion
//
// PAF: Q 12 0 12 + T 12 0 12 10 12 60 cg:Z:3=2I3=2D4=
//
// Query:   0  1  2  3  4  5  6  7           8  9  10 11
//          A  C  T  g  a  G  A  C           T  A  G  A
// Target:  0  1  2        3  4  5  6  7     8  9  10 11
//          A  C  T        G  A  C  c  g     T  A  G  A
// CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
//          |--3=-|  |2I|  |--3=-|  |2D|    |----4=-----|
//
// Block 1: query [0,3)   -> target [0,3)    3=
// Ins:     query [3,5)   -> target ---       2I
// Block 2: query [5,8)   -> target [3,6)    3=
// Del:     query ---      -> target [6,8)    2D
// Block 3: query [8,12)  -> target [8,12)   4=
// =============================================================================

fn fwd_paf() -> PafRecord {
    PafRecord::new("Q 12 0 12 + T 12 0 12 10 12 60 cg:Z:3=2I3=2D4=").unwrap()
}

fn rgn(st: u64, en: u64) -> Region {
    Region {
        name: "T".to_string(),
        st,
        en,
        id: "".to_string(),
        ..Default::default()
    }
}

// =============================================================================
// FORWARD STRAND: Trim to exact block boundaries
// =============================================================================

#[test]
fn fwd_trim_block1_exact() {
    // Query:   0  1  2  3  4  5  6  7           8  9  10 11
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //          [======)
    //          T:[0,3) -> keep 3=
    //          Q:[0,3)
    let paf = fwd_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(0, 3), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "3=");
    assert_eq!((trim.t_st, trim.t_en), (0, 3));
    assert_eq!((trim.q_st, trim.q_en), (0, 3));
}

#[test]
fn fwd_trim_block2_exact() {
    // Query:   0  1  2  3  4  5  6  7           8  9  10 11
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //                         [======)
    //                         T:[3,6) -> keep 3=
    //                         Q:[5,8)
    let paf = fwd_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(3, 6), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "3=");
    assert_eq!((trim.t_st, trim.t_en), (3, 6));
    assert_eq!((trim.q_st, trim.q_en), (5, 8));
}

#[test]
fn fwd_trim_block3_exact() {
    // Query:   0  1  2  3  4  5  6  7           8  9  10 11
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //                                           [=========)
    //                                           T:[8,12) -> keep 4=
    //                                           Q:[8,12)
    let paf = fwd_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(8, 12), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "4=");
    assert_eq!((trim.t_st, trim.t_en), (8, 12));
    assert_eq!((trim.q_st, trim.q_en), (8, 12));
}

// =============================================================================
// FORWARD STRAND: Trim spanning gaps
// =============================================================================

#[test]
fn fwd_trim_span_insertion() {
    // Query:   0  1  2  3  4  5  6  7           8  9  10 11
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //          [=================)
    //          T:[0,6) -> keep 3=2I3=
    //          Q:[0,8)
    //
    // The insertion (2I) is "inside" the target range because
    // it sits between target positions 2 and 3.
    let paf = fwd_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(0, 6), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "3=2I3=");
    assert_eq!((trim.t_st, trim.t_en), (0, 6));
    assert_eq!((trim.q_st, trim.q_en), (0, 8));
}

#[test]
fn fwd_trim_span_deletion() {
    // Query:   0  1  2  3  4  5  6  7           8  9  10 11
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //                         [======================)
    //                         T:[3,10) -> keep 3=2D2=
    //                         Q:[5,10)
    //
    // Spans the deletion gap: target [6,8) has no query bases
    // but the deletion is kept in the CIGAR.
    let paf = fwd_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(3, 10), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "3=2D2=");
    assert_eq!((trim.t_st, trim.t_en), (3, 10));
    assert_eq!((trim.q_st, trim.q_en), (5, 10));
}

#[test]
fn fwd_trim_span_all_blocks() {
    // Query:   0  1  2  3  4  5  6  7           8  9  10 11
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //          [==============================================)
    //          T:[0,12) -> keep entire CIGAR
    //          Q:[0,12)
    let paf = fwd_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(0, 12), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "3=2I3=2D4=");
    assert_eq!((trim.t_st, trim.t_en), (0, 12));
    assert_eq!((trim.q_st, trim.q_en), (0, 12));
}

// =============================================================================
// FORWARD STRAND: Trim cutting into blocks (partial ops)
// =============================================================================

#[test]
fn fwd_trim_cut_into_block1() {
    // Query:   0  1  2  3  4  5  6  7           8  9  10 11
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //             [==)
    //             T:[1,2) -> keep 1=
    //             Q:[1,2)
    let paf = fwd_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(1, 2), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "1=");
    assert_eq!((trim.t_st, trim.t_en), (1, 2));
    assert_eq!((trim.q_st, trim.q_en), (1, 2));
}

#[test]
fn fwd_trim_cut_block1_into_block2() {
    // Query:   0  1  2  3  4  5  6  7           8  9  10 11
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //             [=============)
    //             T:[1,5) -> keep 2=2I2=
    //             Q:[1,7)
    //
    // Cuts 1bp off the start of block1 and 1bp off the end of block2.
    let paf = fwd_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(1, 5), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "2=2I2=");
    assert_eq!((trim.t_st, trim.t_en), (1, 5));
    assert_eq!((trim.q_st, trim.q_en), (1, 7));
}

#[test]
fn fwd_trim_middle_of_alignment() {
    // Query:   0  1  2  3  4  5  6  7           8  9  10 11
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //                      [==================)
    //                      T:[4,9) -> keep 2=2D1=
    //                      Q:[6,9)
    //
    // Cuts into block2 (skip 1bp), includes the 2D, cuts into block3 (keep 1bp).
    let paf = fwd_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(4, 9), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "2=2D1=");
    assert_eq!((trim.t_st, trim.t_en), (4, 9));
    assert_eq!((trim.q_st, trim.q_en), (6, 9));
}

#[test]
fn fwd_trim_both_sides() {
    // Query:   0  1  2  3  4  5  6  7           8  9  10 11
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //             [=============================)
    //             T:[1,10) -> keep 2=2I3=2D2=
    //             Q:[1,10)
    //
    // Trims 1bp off block1 start and 2bp off block3 end.
    let paf = fwd_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(1, 10), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "2=2I3=2D2=");
    assert_eq!((trim.t_st, trim.t_en), (1, 10));
    assert_eq!((trim.q_st, trim.q_en), (1, 10));
}

// =============================================================================
// FORWARD STRAND: Trim to deletion-only region returns None
// =============================================================================

#[test]
fn fwd_trim_deletion_only_returns_none() {
    // Query:   0  1  2  3  4  5  6  7           8  9  10 11
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //                                  [==)
    //                                  T:[6,8) -> only deletion, no aligned bases
    //
    // A region covering only the deletion gap produces no matches.
    let paf = fwd_paf();
    assert!(trim_paf_rec_to_rgn_fast(&rgn(6, 8), &paf).is_none());
}

// =============================================================================
// FORWARD STRAND: Trim ending at deletion boundary strips trailing del
// =============================================================================

#[test]
fn fwd_trim_ending_at_deletion() {
    // Query:   0  1  2  3  4  5  6  7           8  9  10 11
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //          [===========================)
    //          T:[0,8) -> CIGAR would be 3=2I3=2D, but trailing 2D is stripped
    //          -> Result: 3=2I3=, T:[0,6), Q:[0,8)
    //
    // Trailing indels are always removed and coordinates adjusted.
    let paf = fwd_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(0, 8), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "3=2I3=");
    assert_eq!((trim.t_st, trim.t_en), (0, 6));
    assert_eq!((trim.q_st, trim.q_en), (0, 8));
}

// =============================================================================
// FORWARD STRAND: Trim starting at deletion boundary strips leading del
// =============================================================================

#[test]
fn fwd_trim_starting_at_deletion() {
    // Query:   0  1  2  3  4  5  6  7           8  9  10 11
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //                                  [===============)
    //                                  T:[6,12) -> CIGAR would be 2D4=
    //                                  -> Leading 2D stripped: 4=, T:[8,12), Q:[8,12)
    let paf = fwd_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(6, 12), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "4=");
    assert_eq!((trim.t_st, trim.t_en), (8, 12));
    assert_eq!((trim.q_st, trim.q_en), (8, 12));
}

// =============================================================================
// FORWARD STRAND: Superset region returns original
// =============================================================================

#[test]
fn fwd_trim_superset_region() {
    // Query:   0  1  2  3  4  5  6  7           8  9  10 11
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //      [===================================================)
    //      T:[0,20) -> region larger than alignment, returns full record
    //      Q:[0,12)
    let paf = fwd_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(0, 20), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "3=2I3=2D4=");
    assert_eq!((trim.t_st, trim.t_en), (0, 12));
    assert_eq!((trim.q_st, trim.q_en), (0, 12));
}

// =============================================================================
// FORWARD STRAND: Non-overlapping region returns None
// =============================================================================

#[test]
fn fwd_trim_no_overlap_returns_none() {
    // Query:   0  1  2  3  4  5  6  7           8  9  10 11
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //                                                         [======)
    //                                                         T:[15,20)
    let paf = fwd_paf();
    assert!(trim_paf_rec_to_rgn_fast(&rgn(15, 20), &paf).is_none());
}

// =============================================================================
// TEST ALIGNMENT 2: Reverse strand with same CIGAR
//
// PAF: Q 12 0 12 - T 12 0 12 10 12 60 cg:Z:3=2I3=2D4=
//
// CIGAR processes target left-to-right, query right-to-left:
//
// Query:  11 10  9  8  7  6  5  4           3  2  1  0
//          A  C  T  g  a  G  A  C           T  A  G  A
// Target:  0  1  2        3  4  5  6  7     8  9  10 11
//          A  C  T        G  A  C  c  g     T  A  G  A
// CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
//          |--3=-|  |2I|  |--3=-|  |2D|    |----4=-----|
//
// Block 1: query [9,12)  -> target [0,3)   3=  (q goes 12→9)
// Ins:     query [7,9)   -> target ---     2I  (q goes 9→7)
// Block 2: query [4,7)   -> target [3,6)   3=  (q goes 7→4)
// Del:     query ---      -> target [6,8)   2D
// Block 3: query [0,4)   -> target [8,12)  4=  (q goes 4→0)
// =============================================================================

fn rev_paf() -> PafRecord {
    PafRecord::new("Q 12 0 12 - T 12 0 12 10 12 60 cg:Z:3=2I3=2D4=").unwrap()
}

// =============================================================================
// REVERSE STRAND: Trim to exact block boundaries
// =============================================================================

#[test]
fn rev_trim_block1_exact() {
    // Query:  11 10  9  8  7  6  5  4           3  2  1  0
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //          [======)
    //          T:[0,3) -> keep 3=
    //          Q:[9,12) (reverse: q_en=12-0=12, q_st=12-3=9)
    let paf = rev_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(0, 3), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "3=");
    assert_eq!((trim.t_st, trim.t_en), (0, 3));
    assert_eq!((trim.q_st, trim.q_en), (9, 12));
}

#[test]
fn rev_trim_block2_exact() {
    // Query:  11 10  9  8  7  6  5  4           3  2  1  0
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //                         [======)
    //                         T:[3,6) -> keep 3=
    //                         Q:[4,7) (reverse: q_en=12-5=7, q_st=7-3=4)
    let paf = rev_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(3, 6), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "3=");
    assert_eq!((trim.t_st, trim.t_en), (3, 6));
    assert_eq!((trim.q_st, trim.q_en), (4, 7));
}

#[test]
fn rev_trim_block3_exact() {
    // Query:  11 10  9  8  7  6  5  4           3  2  1  0
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //                                           [=========)
    //                                           T:[8,12) -> keep 4=
    //                                           Q:[0,4) (reverse: q_en=12-8=4, q_st=4-4=0)
    let paf = rev_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(8, 12), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "4=");
    assert_eq!((trim.t_st, trim.t_en), (8, 12));
    assert_eq!((trim.q_st, trim.q_en), (0, 4));
}

// =============================================================================
// REVERSE STRAND: Trim spanning gaps
// =============================================================================

#[test]
fn rev_trim_span_insertion() {
    // Query:  11 10  9  8  7  6  5  4           3  2  1  0
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //          [=================)
    //          T:[0,6) -> keep 3=2I3=
    //          Q:[4,12) (reverse: q_en=12-0=12, q_st=12-8=4)
    let paf = rev_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(0, 6), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "3=2I3=");
    assert_eq!((trim.t_st, trim.t_en), (0, 6));
    assert_eq!((trim.q_st, trim.q_en), (4, 12));
}

#[test]
fn rev_trim_span_deletion() {
    // Query:  11 10  9  8  7  6  5  4           3  2  1  0
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //                         [======================)
    //                         T:[3,10) -> keep 3=2D2=
    //                         Q:[2,7) (reverse: q_en=12-5=7, q_st=7-5=2)
    let paf = rev_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(3, 10), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "3=2D2=");
    assert_eq!((trim.t_st, trim.t_en), (3, 10));
    assert_eq!((trim.q_st, trim.q_en), (2, 7));
}

#[test]
fn rev_trim_span_all_blocks() {
    // Query:  11 10  9  8  7  6  5  4           3  2  1  0
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //          [==============================================)
    //          T:[0,12) -> keep entire CIGAR
    //          Q:[0,12)
    let paf = rev_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(0, 12), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "3=2I3=2D4=");
    assert_eq!((trim.t_st, trim.t_en), (0, 12));
    assert_eq!((trim.q_st, trim.q_en), (0, 12));
}

// =============================================================================
// REVERSE STRAND: Partial trims
// =============================================================================

#[test]
fn rev_trim_both_sides() {
    // Query:  11 10  9  8  7  6  5  4           3  2  1  0
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //             [=============================)
    //             T:[1,10) -> keep 2=2I3=2D2=
    //             Q:[2,11) (reverse: q_en=12-1=11, q_st=11-9=2)
    let paf = rev_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(1, 10), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "2=2I3=2D2=");
    assert_eq!((trim.t_st, trim.t_en), (1, 10));
    assert_eq!((trim.q_st, trim.q_en), (2, 11));
}

// =============================================================================
// TEST ALIGNMENT 3: Forward strand with mismatches
//
// PAF: Q 10 0 10 + T 10 0 10 8 10 60 cg:Z:3=1X2=1X3=
//
// Query:   0  1  2  3  4  5  6  7  8  9
//          A  C  T  a  G  A  g  T  A  C
// Target:  0  1  2  3  4  5  6  7  8  9
//          A  C  T  G  G  A  A  T  A  C
// CIGAR:   =  =  =  X  =  =  X  =  =  =
//          |-3=--|  X  |-2=-|  X  |--3=-|
// =============================================================================

fn mismatch_paf() -> PafRecord {
    PafRecord::new("Q 10 0 10 + T 10 0 10 8 10 60 cg:Z:3=1X2=1X3=").unwrap()
}

#[test]
fn mismatch_trim_around_first_snp() {
    // Query:   0  1  2  3  4  5  6  7  8  9
    // Target:  0  1  2  3  4  5  6  7  8  9
    // CIGAR:   =  =  =  X  =  =  X  =  =  =
    //                [=====)
    //                T:[2,5) -> keep 1=1X1=
    //                Q:[2,5)
    let paf = mismatch_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(2, 5), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "1=1X1=");
    assert_eq!((trim.t_st, trim.t_en), (2, 5));
    assert_eq!((trim.q_st, trim.q_en), (2, 5));
}

#[test]
fn mismatch_trim_between_snps() {
    // Query:   0  1  2  3  4  5  6  7  8  9
    // Target:  0  1  2  3  4  5  6  7  8  9
    // CIGAR:   =  =  =  X  =  =  X  =  =  =
    //                   [========)
    //                   T:[3,7) -> keep 1X2=1X
    //
    // But leading/trailing indel stripping doesn't apply to X/=.
    // X is a Diff which is a "match" type — it stays.
    //                   Q:[3,7)
    let paf = mismatch_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(3, 7), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "1X2=1X");
    assert_eq!((trim.t_st, trim.t_en), (3, 7));
    assert_eq!((trim.q_st, trim.q_en), (3, 7));
}

// =============================================================================
// TEST ALIGNMENT 4: CS tag preservation through trimming
//
// PAF: Q 12 0 12 + T 12 0 12 10 12 60 cs:Z::3+ga:3-cg:4
//
// This has no cg tag, so cs is parsed to generate CIGAR:
//   :3 → 3=, +ga → 2I, :3 → 3=, -cg → 2D, :4 → 4=
//
// Same alignment as Test 1, but with cs tag input.
// Query:   0  1  2  3  4  5  6  7           8  9  10 11
// Target:  0  1  2        3  4  5  6  7     8  9  10 11
// CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
// CS:     |--:3--|  +ga  |--:3--|  -cg    |----:4------|
// =============================================================================

fn cs_paf() -> PafRecord {
    PafRecord::new("Q 12 0 12 + T 12 0 12 10 12 60 cs:Z::3+ga:3-cg:4").unwrap()
}

#[test]
fn cs_parsed_cigar_correct() {
    // Verify the cs tag was correctly parsed into CIGAR ops
    let paf = cs_paf();
    assert_eq!(paf.cigar.to_string(), "3=2I3=2D4=");
    assert!(paf.cs_ops.is_some());
}

#[test]
fn cs_trim_block1_preserves_cs() {
    // Query:   0  1  2  3  4  5  6  7           8  9  10 11
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CS:     |--:3--|  +ga  |--:3--|  -cg    |----:4------|
    //          [======)
    //          T:[0,3) -> CIGAR: 3=, CS: :3
    let paf = cs_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(0, 3), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "3=");
    assert!(trim.cs_ops.is_some());
    let cs_str: String = trim.cs_ops.unwrap().to_cs_string();
    assert_eq!(cs_str, ":3");
}

#[test]
fn cs_trim_span_insertion_preserves_cs() {
    // Query:   0  1  2  3  4  5  6  7           8  9  10 11
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CS:     |--:3--|  +ga  |--:3--|  -cg    |----:4------|
    //          [=================)
    //          T:[0,6) -> CIGAR: 3=2I3=, CS: :3+ga:3
    let paf = cs_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(0, 6), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "3=2I3=");
    assert!(trim.cs_ops.is_some());
    let cs_str: String = trim.cs_ops.unwrap().to_cs_string();
    assert_eq!(cs_str, ":3+ga:3");
}

#[test]
fn cs_trim_partial_preserves_cs() {
    // Query:   0  1  2  3  4  5  6  7           8  9  10 11
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CS:     |--:3--|  +ga  |--:3--|  -cg    |----:4------|
    //             [=============================)
    //             T:[1,10) -> CIGAR: 2=2I3=2D2=
    //                         CS:    :2+ga:3-cg:2
    //
    // The first :3 is trimmed to :2 (skip 1, keep 2)
    // The last :4 is trimmed to :2 (skip 0, keep 2)
    let paf = cs_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(1, 10), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "2=2I3=2D2=");
    assert!(trim.cs_ops.is_some());
    let cs_str: String = trim.cs_ops.unwrap().to_cs_string();
    assert_eq!(cs_str, ":2+ga:3-cg:2");
}

// =============================================================================
// TEST ALIGNMENT 5: CS tag with extended form (=ACGT and *xy)
//
// PAF: Q 8 0 8 + T 8 0 8 7 8 60 cs:Z:=ACT*ag=ACGT
//
// CIGAR: 3=1X4=
//
// Query:   0  1  2  3  4  5  6  7
//          A  C  T  g  A  C  G  T
// Target:  0  1  2  3  4  5  6  7
//          A  C  T  a  A  C  G  T
// CIGAR:   =  =  =  X  =  =  =  =
// CS:     |=ACT-| *ag |--=ACGT---|
// =============================================================================

fn cs_extended_paf() -> PafRecord {
    PafRecord::new("Q 8 0 8 + T 8 0 8 7 8 60 cs:Z:=ACT*ag=ACGT").unwrap()
}

#[test]
fn cs_extended_parsed_correctly() {
    let paf = cs_extended_paf();
    assert_eq!(paf.cigar.to_string(), "3=1X4=");
    assert!(paf.cs_ops.is_some());
}

#[test]
fn cs_extended_trim_cuts_matchseq() {
    // Query:   0  1  2  3  4  5  6  7
    // Target:  0  1  2  3  4  5  6  7
    // CS:     |=ACT-| *ag |--=ACGT---|
    //                [=====)
    //                T:[2,5) -> CIGAR: 1=1X1=
    //                CS: =T*ag=A
    //
    // =ACT trimmed: skip 2, keep 1 → =T
    // *ag: kept entirely (length 1)
    // =ACGT trimmed: skip 0, keep 1 → =A
    let paf = cs_extended_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(2, 5), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "1=1X1=");
    assert!(trim.cs_ops.is_some());
    let cs_str: String = trim.cs_ops.unwrap().to_cs_string();
    assert_eq!(cs_str, "=T*ag=A");
}

#[test]
fn cs_extended_trim_full_preserves() {
    // Query:   0  1  2  3  4  5  6  7
    // Target:  0  1  2  3  4  5  6  7
    // CS:     |=ACT-| *ag |--=ACGT---|
    //          [===========================)
    //          T:[0,8) -> full alignment preserved
    //          CS: =ACT*ag=ACGT
    let paf = cs_extended_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(0, 8), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "3=1X4=");
    let cs_str: String = trim.cs_ops.unwrap().to_cs_string();
    assert_eq!(cs_str, "=ACT*ag=ACGT");
}

// =============================================================================
// BREAK-PAF TESTS
//
// Using Alignment 1:
// PAF: Q 12 0 12 + T 12 0 12 10 12 60 cg:Z:3=2I3=2D4=
//
// Query:   0  1  2  3  4  5  6  7           8  9  10 11
// Target:  0  1  2        3  4  5  6  7     8  9  10 11
// CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
//          |--3=-|  |2I|  |--3=-|  |2D|    |----4=-----|
//
// break_length=0: break at ALL indels (len > 0)
//   -> Break at 2I (between t:2 and t:3): region1 T:[0,3), region2 starts at T:[3,...)
//   -> Break at 2D (between t:5 and t:8): region2 T:[3,6), region3 T:[8,12)
//
//   Record 1: T:[0,3),  Q:[0,3),  CIGAR: 3=
//   Record 2: T:[3,6),  Q:[5,8),  CIGAR: 3=
//   Record 3: T:[8,12), Q:[8,12), CIGAR: 4=
// =============================================================================

#[test]
fn break_fwd_at_all_indels() {
    // Query:   0  1  2  3  4  5  6  7           8  9  10 11
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //          [======)       [======)          [=========)
    //          T:[0,3)        T:[3,6)           T:[8,12)
    //          3=             3=                4=
    let paf = fwd_paf();
    let broken = break_paf_on_indels(&paf, 0);
    assert_eq!(broken.len(), 3, "Should break into 3 records");

    // Record 1
    assert_eq!(broken[0].cigar.to_string(), "3=");
    assert_eq!((broken[0].t_st, broken[0].t_en), (0, 3));
    assert_eq!((broken[0].q_st, broken[0].q_en), (0, 3));

    // Record 2
    assert_eq!(broken[1].cigar.to_string(), "3=");
    assert_eq!((broken[1].t_st, broken[1].t_en), (3, 6));
    assert_eq!((broken[1].q_st, broken[1].q_en), (5, 8));

    // Record 3
    assert_eq!(broken[2].cigar.to_string(), "4=");
    assert_eq!((broken[2].t_st, broken[2].t_en), (8, 12));
    assert_eq!((broken[2].q_st, broken[2].q_en), (8, 12));
}

#[test]
fn break_rev_at_all_indels() {
    // Query:  11 10  9  8  7  6  5  4           3  2  1  0
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //          [======)       [======)          [=========)
    //          T:[0,3)        T:[3,6)           T:[8,12)
    //          Q:[9,12)       Q:[4,7)           Q:[0,4)
    let paf = rev_paf();
    let broken = break_paf_on_indels(&paf, 0);
    assert_eq!(broken.len(), 3, "Should break into 3 records");

    // Record 1
    assert_eq!(broken[0].cigar.to_string(), "3=");
    assert_eq!((broken[0].t_st, broken[0].t_en), (0, 3));
    assert_eq!((broken[0].q_st, broken[0].q_en), (9, 12));

    // Record 2
    assert_eq!(broken[1].cigar.to_string(), "3=");
    assert_eq!((broken[1].t_st, broken[1].t_en), (3, 6));
    assert_eq!((broken[1].q_st, broken[1].q_en), (4, 7));

    // Record 3
    assert_eq!(broken[2].cigar.to_string(), "4=");
    assert_eq!((broken[2].t_st, broken[2].t_en), (8, 12));
    assert_eq!((broken[2].q_st, broken[2].q_en), (0, 4));
}

#[test]
fn break_no_break_when_threshold_too_high() {
    // Query:   0  1  2  3  4  5  6  7           8  9  10 11
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //
    // break_length=2: only break at indels > 2bp
    // Both indels are 2bp, so 2 > 2 is false. No breaks.
    let paf = fwd_paf();
    let broken = break_paf_on_indels(&paf, 2);
    assert_eq!(broken.len(), 1, "No breaks when threshold too high");
    assert_eq!(broken[0].cigar.to_string(), "3=2I3=2D4=");
}

#[test]
fn break_selective_threshold() {
    // PAF: Q 14 0 14 + T 14 0 14 10 14 60 cg:Z:3=1I3=5D4=3I
    //
    // Query:   0  1  2  3  4  5  6  7                8  9  10 11 12 13
    //          A  C  T  g  G  A  C  -  -  -  -  -    T  A  G  A  a  c  a
    // Target:  0  1  2     3  4  5  6  7  8  9  10  11 12 13
    //          A  C  T     G  A  C  g  a  t  c  g    T  A  G  A
    // CIGAR:   =  =  =  I  =  =  =  D  D  D  D  D   =  =  =  =  I  I  I
    //          |-3=--|  1I |-3=--|  |-----5D------|  |----4=------|  |-3I-|
    //
    // break_length=2: break at indels > 2bp
    //   1I (len=1): 1 > 2? No.
    //   5D (len=5): 5 > 2? Yes! Break here.
    //   3I (len=3): 3 > 2? Yes! Break here.
    //
    //   But 3I is at the end, so after the last break, the final region
    //   T:[11,14) produces 4= (trailing 3I stripped).
    //
    //   Record 1: T:[0,6), Q:[0,7), CIGAR: 3=1I3=    (trailing indels stripped: none needed)
    //   Record 2: T:[11,14), Q:[8,14)... wait.
    //
    // Actually let me retrace. The 5D consumes reference [6,11).
    // The code sets pre_tpos = cur_tpos = 6, then adds 5 (Del) → pre_tpos = 11.
    //
    //   Region 1: T:[0,6)
    //   After 5D break: pre_tpos=11
    //
    //   3I at end: len=3 > 2. cur_tpos=14 (from 4=). 14 > 11.
    //   Region 2: T:[11,14+1)... wait, let me re-check cur_tpos. By the time
    //   we hit 3I, cur_tpos = 0+3+3+5+4 = 15. That's wrong.
    //   Actually: 3= moves ref 3, 1I doesn't, 3= moves ref 3, 5D moves ref 5, 4= moves ref 4, 3I doesn't.
    //   target consumed = 3+3+5+4 = 15.
    //   But wait, q consumed = 3+1+3+4+3 = 14. t consumed = 15.
    //   So PAF would be: Q 14 0 14 + T 15 0 15 ...
    //
    // Let me fix this example to be consistent.
    //
    // With CIGAR 3=1I3=5D4=:
    //   q consumed = 3+1+3+4 = 11
    //   t consumed = 3+3+5+4 = 15
    //   PAF: Q 11 0 11 + T 15 0 15 ...
    //
    //   break_length=2: 1I (1>2? No), 5D (5>2? Yes)
    //   Record 1: T:[0,6), CIGAR: 3=1I3=
    //   Record 2: T:[11,15), CIGAR: 4=
    let paf = PafRecord::new("Q 11 0 11 + T 15 0 15 10 15 60 cg:Z:3=1I3=5D4=").unwrap();
    let broken = break_paf_on_indels(&paf, 2);
    assert_eq!(broken.len(), 2, "Should break at 5D only");

    assert_eq!(broken[0].cigar.to_string(), "3=1I3=");
    assert_eq!((broken[0].t_st, broken[0].t_en), (0, 6));

    assert_eq!(broken[1].cigar.to_string(), "4=");
    assert_eq!((broken[1].t_st, broken[1].t_en), (11, 15));
}

// =============================================================================
// TRUNCATE BY QUERY: Test via trim_overlapping_pafs
//
// Alignment for truncation testing:
// PAF: Q 12 0 12 + T 12 0 12 10 12 60 cg:Z:3=2I3=2D4=
//
// Query:   0  1  2  3  4  5  6  7           8  9  10 11
// Target:  0  1  2        3  4  5  6  7     8  9  10 11
// CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
// =============================================================================

#[test]
fn truncate_by_query_first_half() {
    // Query:   0  1  2  3  4  5  6  7           8  9  10 11
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //          [=================)
    //          Q:[0,8) -> keep 3=2I3=
    //          T:[0,6)
    let mut paf = fwd_paf();
    paf.truncate_record_by_query(0, 8);
    assert_eq!(paf.cigar.to_string(), "3=2I3=");
    assert_eq!((paf.t_st, paf.t_en), (0, 6));
    assert_eq!((paf.q_st, paf.q_en), (0, 8));
}

#[test]
fn truncate_by_query_second_half() {
    // Query:   0  1  2  3  4  5  6  7           8  9  10 11
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //                            [========================)
    //                            Q:[5,12) -> keep 3=2D4=
    //                            T:[3,12)
    let mut paf = fwd_paf();
    paf.truncate_record_by_query(5, 12);
    assert_eq!(paf.cigar.to_string(), "3=2D4=");
    assert_eq!((paf.t_st, paf.t_en), (3, 12));
    assert_eq!((paf.q_st, paf.q_en), (5, 12));
}

#[test]
fn truncate_by_query_middle() {
    // Query:   0  1  2  3  4  5  6  7           8  9  10 11
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //                      [=============)
    //                      Q:[4,9) -> skip_before cuts 1bp off 2I, starts in I
    //
    // Walking CIGAR with query coords:
    //   3= (q:0→3): skip. t_before=3.
    //   2I (q:3→5): overlap [4,5), skip 1, keep 1 → 1I. Started.
    //   3= (q:5→8): keep 3 → 3=. t_in=3.
    //   2D: keep → 2D. t_in += 2.
    //   4= (q:8→12): overlap [8,9), keep 1 → 1=. t_in += 1.
    //
    // Leading indel (1I) stripped → lead_q_adjust=1, q_st bumps to 5.
    // Result: 3=2D1=, t consumed = 3+2+1 = 6, T:[3,9), Q:[5,9)
    let mut paf = fwd_paf();
    paf.truncate_record_by_query(4, 9);
    assert_eq!(paf.cigar.to_string(), "3=2D1=");
    assert_eq!((paf.t_st, paf.t_en), (3, 9));
    assert_eq!((paf.q_st, paf.q_en), (5, 9));
}

// =============================================================================
// REVERSE STRAND: Truncate by query
// =============================================================================

#[test]
fn truncate_rev_by_query_high_end() {
    // Query:  11 10  9  8  7  6  5  4           3  2  1  0
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //
    // Truncate to Q:[4,12) (keep high query coordinates)
    // In forward query order: Q:[4,12) means keep right side of query
    // For reverse strand, this corresponds to the left side of target:
    //   Block 1: Q:[9,12) -> T:[0,3)  -- inside [4,12)
    //   Ins:     Q:[7,9)  -> ---      -- inside [4,12)
    //   Block 2: Q:[4,7)  -> T:[3,6)  -- inside [4,12)
    //   Block 3: Q:[0,4)  -> T:[8,12) -- outside [4,12)
    //
    // Result: 3=2I3=, T:[0,6), Q:[4,12)
    let mut paf = rev_paf();
    paf.truncate_record_by_query(4, 12);
    assert_eq!(paf.cigar.to_string(), "3=2I3=");
    assert_eq!((paf.t_st, paf.t_en), (0, 6));
    assert_eq!((paf.q_st, paf.q_en), (4, 12));
}

#[test]
fn truncate_rev_by_query_low_end() {
    // Query:  11 10  9  8  7  6  5  4           3  2  1  0
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    //
    // Truncate to Q:[0,7) (keep low query coordinates)
    // For reverse strand, this corresponds to the right side of target:
    //   Block 2: Q:[4,7) -> T:[3,6)  -- inside [0,7)
    //   Block 3: Q:[0,4) -> T:[8,12) -- inside [0,7)
    //
    // Result: 3=2D4=, T:[3,12), Q:[0,7)
    let mut paf = rev_paf();
    paf.truncate_record_by_query(0, 7);
    assert_eq!(paf.cigar.to_string(), "3=2D4=");
    assert_eq!((paf.t_st, paf.t_en), (3, 12));
    assert_eq!((paf.q_st, paf.q_en), (0, 7));
}

// =============================================================================
// TEST ALIGNMENT 6: Existing test case from liftover.rs
//
// PAF: Q 10 2 10 + T 40 12 20 3 9 60 cg:Z:4M1I1=1D2=
//
// This uses the M operator (match/mismatch), unlike = (sequence match).
//
// Query (fwd):  2  3  4  5  6  7  -  8  9
// Target:      12 13 14 15 16 17 18 19
//
// CIGAR:   M  M  M  M  I  =  D  =  =
//          |----4M---| 1I 1= 1D |-2=-|
//
// Offset: q_st=2, t_st=12
// Block 1: query [2,6) -> target [12,16) 4M
// Ins:     query [6,7) -> target ---     1I
// Block 2: query [7,8) -> target [16,17) 1=
// Del:     query ---   -> target [17,18) 1D
// Block 3: query [8,10) -> target [18,20) 2=
//
// ## visual alignment at absolute coords
//
//  14-18         XXXXX
//  0123456789012345567890....
//  ACTGACTGAAACTGAC-TAGA
//  ------------||||I|D||
//              TGACGT-AC
//            01234567789 (forward)
//                XXXXX
//              98765433210 (reverse)
// =============================================================================

#[test]
fn existing_test_fwd_trim_14_15() {
    // Target:      12 13 14 15 16 17 18 19
    // CIGAR:        M  M  M  M  I  =  D  =  =
    //                     [==)
    //                     T:[14,15) -> M trimmed to 1M (skip 2, keep 1)
    //                     Q:[4,5) (fwd: q_st=2+2=4, q_en=4+1=5)
    let paf = PafRecord::new("Q 10 2 10 + T 40 12 20 3 9 60 cg:Z:4M1I1=1D2=").unwrap();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(14, 15), &paf).unwrap();
    assert_eq!((trim.q_st, trim.q_en), (4, 5));
    assert_eq!((trim.t_st, trim.t_en), (14, 15));
}

#[test]
fn existing_test_fwd_trim_14_18() {
    // Target:      12 13 14 15 16 17 18 19
    // CIGAR:        M  M  M  M  I  =  D  =  =
    //                     [===========)
    //                     T:[14,18) -> 2M1I1=1D stripped → 2M1I1=
    //                     Q:[4,8) ... wait, 2M=2 query + 1I=1 query + 1==1 query = 4 query
    //                     Actually trailing 1D is stripped → T becomes [14,17)
    //
    // Let me retrace:
    //   4M (t:12→16): overlap [14,16), skip=2, keep=2. Push 2M. q_in=2.
    //   1I: Push 1I. q_in += 1 = 3.
    //   1= (t:16→17): overlap [16,17), keep=1. Push 1=. q_in += 1 = 4.
    //   1D (t:17→18): overlap [17,18), keep=1. Push 1D. 18>=18 → done.
    //   Trailing: 1D stripped. trail_t_adjust=1. new_t_en=18-1=17.
    //   Result: 2M1I1=, T:[14,17), Q:[4,8)
    let paf = PafRecord::new("Q 10 2 10 + T 40 12 20 3 9 60 cg:Z:4M1I1=1D2=").unwrap();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(14, 18), &paf).unwrap();
    assert_eq!((trim.q_st, trim.q_en), (4, 8));
    assert_eq!((trim.t_st, trim.t_en), (14, 17));
    assert_eq!(trim.cigar.to_string(), "2M1I1=");
}

#[test]
fn existing_test_rev_trim_14_15() {
    // PAF: Q 10 2 10 - T 40 12 20 ...
    // For reverse strand:
    //   CIGAR walks target left-to-right but query right-to-left.
    //   q_pos starts at q_en=10
    //   4M (t:12→16): q goes 10→6, so query [6,10)
    //   1I (t:---): q goes 6→5, so query [5,6)
    //   1= (t:16→17): q goes 5→4, so query [4,5)
    //   1D (t:17→18): no query
    //   2= (t:18→20): q goes 4→2, so query [2,4)
    //
    // Trim T:[14,15):
    //   4M (t:12→16): overlap [14,15), skip=2, keep=1. q_before=2. Push 1M. q_in=1.
    //   Done.
    //   Result: 1M, T:[14,15)
    //   q_en = 10 - 2 = 8. q_st = 8 - 1 = 7. Q:[7,8)
    let paf = PafRecord::new("Q 10 2 10 - T 40 12 20 3 9 60 cg:Z:4M1I1=1D2=").unwrap();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(14, 15), &paf).unwrap();
    assert_eq!((trim.q_st, trim.q_en), (7, 8));
    assert_eq!((trim.t_st, trim.t_en), (14, 15));
}

// =============================================================================
// EDGE CASE: Single base alignment
// =============================================================================

#[test]
fn single_base_trim() {
    // PAF: Q 1 0 1 + T 1 0 1 1 1 60 cg:Z:1=
    //
    // Query:   0
    // Target:  0
    // CIGAR:   =
    //          [)
    //          T:[0,1) -> 1=
    let paf = PafRecord::new("Q 1 0 1 + T 1 0 1 1 1 60 cg:Z:1=").unwrap();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(0, 1), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "1=");
    assert_eq!((trim.t_st, trim.t_en), (0, 1));
    assert_eq!((trim.q_st, trim.q_en), (0, 1));
}

// =============================================================================
// EDGE CASE: Alignment is entirely an insertion (should return None)
// =============================================================================

#[test]
fn all_insertion_returns_none() {
    // PAF: Q 5 0 5 + T 0 0 0 0 0 60 cg:Z:5I
    // This is degenerate but the function should handle it gracefully.
    // Actually, t_st=0, t_en=0, so any region will have trim_t_st >= trim_t_en.
    let paf = PafRecord::new("Q 5 0 5 + T 5 0 0 0 5 60 cg:Z:5I");
    // This may fail to parse or return None on trim
    if let Ok(paf) = paf {
        assert!(trim_paf_rec_to_rgn_fast(&rgn(0, 5), &paf).is_none());
    }
}

// =============================================================================
// COMPLEX: Multiple insertions and deletions
//
// PAF: Q 15 0 15 + T 17 0 17 11 17 60 cg:Z:2=3I2=1D2=2D1=1I3=
//
// Walk:
//   2=:  q[0,2)   t[0,2)
//   3I:  q[2,5)   t---
//   2=:  q[5,7)   t[2,4)
//   1D:  q---     t[4,5)
//   2=:  q[7,9)   t[5,7)
//   2D:  q---     t[7,9)
//   1=:  q[9,10)  t[9,10)
//   1I:  q[10,11) t---
//   3=:  q[11,14) t[10,13)... wait that's wrong.
//
// Let me recount: t consumed = 2+2+1+2+2+1+3 = 13. But I need 17.
// q consumed = 2+3+2+2+1+1+3 = 14. But I need 15.
// Let me adjust: cg:Z:2=3I2=1D2=2D1=1I4=
//   t consumed = 2+2+1+2+2+1+4 = 14, q consumed = 2+3+2+2+1+1+4 = 15. Still not 17 for t.
// Let me use: Q 15 0 15 + T 16 0 16 11 16 60 cg:Z:2=3I2=1D2=2D2=1I2=
//   t consumed = 2+2+1+2+2+2+2 = 13. Nope.
// OK let me just be careful:
// cg:Z:2=3I3=1D2=2D1=1I3=
//   q = 2+3+3+2+1+1+3 = 15 ✓
//   t = 2+3+1+2+2+1+3 = 14
// PAF: Q 15 0 15 + T 14 0 14 ...
//
// Query:   0  1  2  3  4  5  6  7     8  9     10    11 12 13 14
// Target:  0  1        2  3  4  5  6  7  8  9  10    11 12 13
// CIGAR:   =  =  I  I  I  =  =  =  D  =  =  D  D    =  I  =  =  =
//
// Actually this is getting complex. Let me simplify.
// =============================================================================

#[test]
fn complex_multiple_indels() {
    // PAF: Q 15 0 15 + T 14 0 14 10 14 60 cg:Z:2=3I3=1D2=2D1=1I3=
    //
    // Step-by-step CIGAR walk:
    //   2=:  q[0,2)    t[0,2)     Match
    //   3I:  q[2,5)    t---       Insertion
    //   3=:  q[5,8)    t[2,5)     Match
    //   1D:  q---      t[5,6)     Deletion
    //   2=:  q[8,10)   t[6,8)     Match
    //   2D:  q---      t[8,10)    Deletion
    //   1=:  q[10,11)  t[10,11)   Match
    //   1I:  q[11,12)  t---       Insertion
    //   3=:  q[12,15)  t[11,14)   Match
    //
    // Query:   0  1  2  3  4  5  6  7     8  9        10    11 12 13 14
    // Target:  0  1        2  3  4  5     6  7  8  9  10    11 12 13
    // CIGAR:   =  =  I  I  I  =  =  =  D  =  =  D  D  =  I  =  =  =
    //
    // Trim T:[2,11):
    //   2= (t:0→2): skip. q_before=2.
    //   3I: not started. q_before += 3.
    //   3= (t:2→5): overlap [2,5), keep=3. Push 3=. q_in=3. t_pos=5.
    //   1D (t:5→6): overlap [5,6), keep=1. Push 1D.
    //   2= (t:6→8): overlap [6,8), keep=2. Push 2=. q_in += 2 = 5.
    //   2D (t:8→10): overlap [8,10), keep=2. Push 2D.
    //   1= (t:10→11): overlap [10,11), keep=1. Push 1=. q_in += 1 = 6. 11>=11.
    //   Result: 3=1D2=2D1=, T:[2,11), Q:[5,11)
    let paf =
        PafRecord::new("Q 15 0 15 + T 14 0 14 10 14 60 cg:Z:2=3I3=1D2=2D1=1I3=").unwrap();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(2, 11), &paf).unwrap();
    assert_eq!(trim.cigar.to_string(), "3=1D2=2D1=");
    assert_eq!((trim.t_st, trim.t_en), (2, 11));
    assert_eq!((trim.q_st, trim.q_en), (5, 11));
}

// =============================================================================
// CS TAG + BREAK-PAF: Verify cs ops are preserved through break
// =============================================================================

#[test]
fn cs_break_preserves_cs_ops() {
    // PAF with cs tag: Q 12 0 12 + T 12 0 12 10 12 60 cs:Z::3+ga:3-cg:4
    //
    // Query:   0  1  2  3  4  5  6  7           8  9  10 11
    // Target:  0  1  2        3  4  5  6  7     8  9  10 11
    // CIGAR:   =  =  =  I  I  =  =  =  D  D    =  =  =  =
    // CS:     |--:3--|  +ga  |--:3--|  -cg    |----:4------|
    //
    // break_length=0:
    //   Record 1: T:[0,3), CS: :3
    //   Record 2: T:[3,6), CS: :3
    //   Record 3: T:[8,12), CS: :4
    let paf = cs_paf();
    let broken = break_paf_on_indels(&paf, 0);
    assert_eq!(broken.len(), 3);

    // All records should have cs_ops
    for rec in &broken {
        assert!(rec.cs_ops.is_some(), "cs_ops should be preserved through break-paf");
    }

    let cs1: String = broken[0].cs_ops.as_ref().unwrap().to_cs_string();
    let cs2: String = broken[1].cs_ops.as_ref().unwrap().to_cs_string();
    let cs3: String = broken[2].cs_ops.as_ref().unwrap().to_cs_string();
    assert_eq!(cs1, ":3");
    assert_eq!(cs2, ":3");
    assert_eq!(cs3, ":4");
}

// =============================================================================
// DELETION BOUNDARY TESTS: Behavior when trim boundary falls inside a deletion
//
// This tests the code path that causes a 1bp difference vs paftools.js liftover
// in --qbed mode. When --qbed is used, the PAF is swapped (I↔D) and trimmed
// by what were originally query coordinates (now target coordinates after swap).
// An insertion in the original PAF becomes a deletion after swap. Trimming
// into this deletion, then stripping trailing indels, produces the correct
// result. paftools.js reports 1bp more in these cases.
//
// We test directly with a deletion-containing PAF to exercise the same
// code path without needing to swap:
//
// PAF: Q 100000 0 100000 + T 130000 0 130000 100000 130000 60 cg:Z:50000M30000D50000M
//
// Query:  0 ............... 49999 |                              | 50000 ............. 99999
//         |------ 50000M --------|                                |------ 50000M ----------|
// Target: 0 ............... 49999 | 50000 ............... 79999  | 80000 ............. 129999
//                                   |------- 30000D (no query) --|
//
// Block 1: query [0, 50000)     -> target [0, 50000)      50000M
// Del:     query ---             -> target [50000, 80000)  30000D (no query movement)
// Block 2: query [50000, 100000) -> target [80000, 130000) 50000M
// =============================================================================

fn big_del_paf() -> PafRecord {
    PafRecord::new(
        "Q 100000 0 100000 + T 130000 0 130000 100000 130000 60 cg:Z:50000M30000D50000M",
    )
    .unwrap()
}

#[test]
fn trim_boundary_before_deletion() {
    // Target region ends before the deletion — no ambiguity.
    //
    // Target: [0 .......... 39999]
    //         |--- 50000M --|--- 30000D ---|--- 50000M ---|
    // Query:  [0 .......... 39999]
    //
    // Expected: T:[0, 40000), Q:[0, 40000), CIGAR: 40000M
    let paf = big_del_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(0, 40000), &paf).unwrap();
    assert_eq!((trim.t_st, trim.t_en), (0, 40000));
    assert_eq!((trim.q_st, trim.q_en), (0, 40000));
    assert_eq!(trim.cigar.to_string(), "40000M");
}

#[test]
fn trim_boundary_at_deletion_start() {
    // Target region ends exactly where the deletion starts.
    //
    // Target: [0 ................. 49999]
    //         |------- 50000M ---------|--- 30000D ---|--- 50000M ---|
    // Query:  [0 ................. 49999]
    //
    // Expected: T:[0, 50000), Q:[0, 50000), CIGAR: 50000M
    let paf = big_del_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(0, 50000), &paf).unwrap();
    assert_eq!((trim.t_st, trim.t_en), (0, 50000));
    assert_eq!((trim.q_st, trim.q_en), (0, 50000));
    assert_eq!(trim.cigar.to_string(), "50000M");
}

#[test]
fn trim_boundary_inside_deletion() {
    // Target region ends inside the deletion (t=70000 is within [50000, 80000)).
    // The deletion doesn't consume query bases, so the last aligned query
    // position is still 49999.
    //
    // Target: [0 ................. 49999 | 50000 ..... 69999]
    //         |------- 50000M ---------|  |-- partial 30000D --|
    // Query:  [0 ................. 49999]  (no query movement)
    //
    // Expected: T:[0, 50000), Q:[0, 50000), CIGAR: 50000M
    // The partial deletion is stripped as a trailing indel.
    //
    // NOTE: In --qbed mode, this corresponds to the case where paftools.js
    // reports a +1bp overcount. rb correctly strips the trailing deletion.
    let paf = big_del_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(0, 70000), &paf).unwrap();
    assert_eq!((trim.t_st, trim.t_en), (0, 50000));
    assert_eq!((trim.q_st, trim.q_en), (0, 50000));
    assert_eq!(trim.cigar.to_string(), "50000M");
}

#[test]
fn trim_boundary_at_deletion_end() {
    // Target region ends at the end of the deletion (t=80000).
    // Same as "inside" — the deletion is stripped as trailing.
    //
    // Target: [0 ................. 49999 | 50000 ............... 79999]
    //         |------- 50000M ---------|  |-------- 30000D ----------|
    // Query:  [0 ................. 49999]  (no query movement)
    //
    // Expected: T:[0, 50000), Q:[0, 50000), CIGAR: 50000M
    let paf = big_del_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(0, 80000), &paf).unwrap();
    assert_eq!((trim.t_st, trim.t_en), (0, 50000));
    assert_eq!((trim.q_st, trim.q_en), (0, 50000));
    assert_eq!(trim.cigar.to_string(), "50000M");
}

#[test]
fn trim_boundary_after_deletion() {
    // Target region spans past the deletion into the second match block.
    //
    // Target: [0 ............. 49999 | 50000 ....... 79999 | 80000 .... 89999]
    //         |------ 50000M -------|  |---- 30000D ------|  |- 10000M -|
    // Query:  [0 ............. 49999]                        [50000 . 59999]
    //
    // Expected: T:[0, 90000), Q:[0, 60000), CIGAR: 50000M30000D10000M
    let paf = big_del_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(0, 90000), &paf).unwrap();
    assert_eq!((trim.t_st, trim.t_en), (0, 90000));
    assert_eq!((trim.q_st, trim.q_en), (0, 60000));
    assert_eq!(trim.cigar.to_string(), "50000M30000D10000M");
}

#[test]
fn trim_start_inside_deletion() {
    // Target region starts inside the deletion.
    // The leading partial deletion is stripped.
    //
    // Target:                                 [60000 . 79999 | 80000 ..... 89999]
    //         |------ 50000M -------|  |---- 30000D --------|  |- 10000M -|
    // Query:                                                   [50000 . 59999]
    //
    // Expected: T:[80000, 90000), Q:[50000, 60000), CIGAR: 10000M
    let paf = big_del_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(60000, 90000), &paf).unwrap();
    assert_eq!((trim.t_st, trim.t_en), (80000, 90000));
    assert_eq!((trim.q_st, trim.q_en), (50000, 60000));
    assert_eq!(trim.cigar.to_string(), "10000M");
}

#[test]
fn trim_entirely_inside_deletion() {
    // Target region falls entirely within the deletion.
    // No query bases are covered — should return None.
    //
    // Target:                          [55000 ......... 75000]
    //         |------ 50000M -------|  |---- 30000D --------|  |--- 50000M ---|
    // Query:                           (no query movement)
    //
    // Expected: None (no aligned bases)
    let paf = big_del_paf();
    let trim = trim_paf_rec_to_rgn_fast(&rgn(55000, 75000), &paf);
    assert!(
        trim.is_none(),
        "Region entirely within deletion should produce no output"
    );
}
