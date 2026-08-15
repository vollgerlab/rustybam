use rustybam::paf::PafRecord;

/// A line with fewer than 12 columns must return an error.
#[test]
fn short_line_returns_err() {
    assert!(PafRecord::new("").is_err());
    assert!(PafRecord::new("A 1 2 3 + B 1 2 3 10 11").is_err());
}

/// A tag token that is too short must return an error.
#[test]
fn short_tag_returns_err() {
    assert!(PafRecord::new("A 1 2 3 + B 1 2 3 10 11 60 x").is_err());
}

/// A malformed tag without the two colons must return an error.
#[test]
fn malformed_tag_returns_err() {
    assert!(PafRecord::new("A 1 2 3 + B 1 2 3 10 11 60 badtag").is_err());
}

/// A truncated cs or cg tag without a value must return an error.
#[test]
fn truncated_cs_and_cg_tags_return_err() {
    assert!(PafRecord::new("A 1 2 3 + B 1 2 3 10 11 60 cs:Z").is_err());
    assert!(PafRecord::new("A 1 2 3 + B 1 2 3 10 11 60 cg:Z").is_err());
}

/// A well formed line with tags must parse correctly.
#[test]
fn well_formed_line_parses() {
    let rec = PafRecord::new("A 1 2 3 + B 1 2 3 10 11 60 tp:A:P cg:Z:1=").unwrap();
    assert_eq!(rec.q_name, "A");
    assert_eq!(rec.cigar.to_string(), "1=");
}

/// check_integrity must not count mismatches (X) as matches in column 10.
#[test]
fn nmatch_excludes_mismatches() {
    let mut rec = PafRecord::new("Q 11 0 11 + T 20 0 11 10 11 60 cg:Z:5=1X5=").unwrap();
    rec.check_integrity().unwrap();
    assert_eq!(rec.nmatch, 10);
    assert_eq!(rec.aln_len, 11);
}
