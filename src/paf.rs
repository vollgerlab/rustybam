use super::bed;
use super::getfasta;
use super::myio;
use super::trim_overlap::trim_overlapping_pafs;
use bio::alphabets::dna::revcomp;
use core::fmt;
use itertools::Itertools;
use lazy_static::lazy_static;
use natord;
use regex::Regex;
use rust_htslib::bam::record::Cigar::*;
use rust_htslib::bam::record::CigarString;
use rust_htslib::bam::record::*;
use rust_htslib::faidx;
use std::collections::{HashMap, HashSet};
use std::convert::TryFrom;
use std::io::BufRead;
use std::str::FromStr;

lazy_static! {
    static ref PAF_TAG: Regex = Regex::new("(..):(.):(.*)").unwrap();
}

#[derive(Debug)]
pub enum Error {
    PafParseCigar { msg: String },
    PafParseCS { msg: String },
    ParseIntError { msg: String },
    ParsePafColumn {},
}
type PafResult<T> = Result<T, crate::paf::Error>;

#[derive(Debug)]
pub struct Paf {
    pub records: Vec<PafRecord>,
    //pub records_by_contig: HashMap<String, Vec<&'a PafRecord>>,
}

impl Default for Paf {
    fn default() -> Self {
        Self::new()
    }
}

impl Paf {
    pub fn new() -> Paf {
        Paf {
            records: Vec::new(),
            //records_by_contig: HashMap::new(),
        }
    }
    /// read in the paf from a file pass "-" for stdin
    /// # Example
    /// ```
    /// use rustybam::paf;
    /// use std::fs::File;
    /// use std::io::*;
    /// let mut paf = paf::Paf::from_file(".test/asm_small.paf");
    /// assert_eq!(paf.records.len(), 249);
    ///
    /// ```
    pub fn from_file(file_name: &str) -> Paf {
        let paf_file = myio::reader(file_name);
        let mut paf = Paf::new();
        // read the paf recs into a vector
        for (index, line) in paf_file.lines().enumerate() {
            log::trace!("{:?}", line);
            match PafRecord::new(&line.unwrap()) {
                Ok(mut rec) => {
                    rec.check_integrity().unwrap();
                    paf.records.push(rec);
                }
                Err(_) => eprintln!("\nUnable to parse PAF record. Skipping line {}", index + 1),
            }
            log::debug!("Read PAF record number: {}", index + 1);
        }
        paf
    }

    pub fn query_name_map(&mut self) -> HashMap<String, Vec<&PafRecord>> {
        let mut records_by_contig = HashMap::new();
        for rec in &self.records {
            (*records_by_contig
                .entry(rec.q_name.clone())
                .or_insert_with(Vec::new))
            .push(rec);
        }
        records_by_contig
    }

    pub fn filter_aln_pairs(&mut self, paired_len: u64) {
        let mut dict = HashMap::new();
        for rec in self.records.iter_mut() {
            let aln_bp = dict
                .entry((rec.t_name.clone(), rec.q_name.clone()))
                .or_insert(0_u64);
            *aln_bp += rec.t_en - rec.t_st;
        }
        self.records.retain(|rec| {
            paired_len < *dict.get(&(rec.t_name.clone(), rec.q_name.clone())).unwrap()
        });
    }

    pub fn filter_query_len(&mut self, min_query_len: u64) {
        self.records.retain(|rec| rec.q_len > min_query_len);
    }

    /// Filter on alignment length
    pub fn filter_aln_len(&mut self, min_aln_len: u64) {
        self.records.retain(|rec| rec.t_en - rec.t_st > min_aln_len);
    }

    /// orient queries relative to their target (inverts if most bases are aligned rc).
    pub fn orient(&mut self) {
        let mut orient_order_dict = HashMap::new();
        let mut t_names = HashSet::new();
        // calculate whether a contig is mostly forward or reverse strand
        // and determine the middle alignment position with respect to the target
        for rec in &self.records {
            let (orient, total_bp, order) = orient_order_dict
                .entry((rec.t_name.clone(), rec.q_name.clone()))
                .or_insert((0_i64, 0_u64, 0_u64));
            // set the orientation of the query relative to the target
            if rec.strand == '-' {
                *orient -= (rec.q_en - rec.q_st) as i64;
            } else {
                *orient += (rec.q_en - rec.q_st) as i64;
            }
            // set a number that will determine the order of the contig
            let weight = rec.t_en - rec.t_st;
            *total_bp += weight;
            *order += weight * (rec.t_st + rec.t_en) / 2;
            // make a list of targets
            t_names.insert(rec.t_name.clone());
        }

        // set the order and orientation of records
        for rec in &mut self.records {
            // set the order of the records
            let (orient, total_bp, order) = orient_order_dict
                .get(&(rec.t_name.clone(), rec.q_name.clone()))
                .unwrap();
            rec.order = *order / *total_bp;

            // reverse record if it is mostly on the rc
            if *orient < 0 {
                rec.q_name = format!("{}-", rec.q_name);
                let new_st = rec.q_len - rec.q_en;
                let new_en = rec.q_len - rec.q_st;
                rec.q_st = new_st;
                rec.q_en = new_en;
                rec.strand = if rec.strand == '+' { '-' } else { '+' };
            } else {
                rec.q_name = format!("{}+", rec.q_name);
            }
        }
    }

    /// scaffold oriented contigs into one fake super contig
    pub fn scaffold(&mut self, spacer_size: u64) {
        // sort the records by their target name and order
        self.records.sort_by(|a, b| {
            a.t_name
                .cmp(&b.t_name) // group by target
                .then(a.order.cmp(&b.order)) // order query by position in target
                .then(a.q_st.cmp(&b.q_st)) // order by position in query
        });

        // group by t_name
        for (_t_name, t_recs) in &self.records.iter_mut().group_by(|rec| rec.t_name.clone()) {
            let mut t_recs: Vec<&mut PafRecord> = t_recs.collect();
            // sort recs by order
            t_recs.sort_by(|a, b| {
                a.order
                    .cmp(&b.order) // order query by position in target
                    .then(a.q_st.cmp(&b.q_st)) // order by position in query
            });

            // new scaffold name
            let scaffold_name = t_recs
                .iter()
                .map(|rec| rec.q_name.clone())
                .unique()
                .collect::<Vec<String>>()
                .join("::");

            let mut scaffold_len = 0_u64;
            for (_q_name, q_recs) in &t_recs.iter_mut().group_by(|rec| rec.q_name.clone()) {
                let q_recs: Vec<&mut &mut PafRecord> = q_recs.collect();
                let q_min = q_recs.iter().map(|rec| rec.q_st).min().unwrap_or(0);
                let q_max = q_recs.iter().map(|rec| rec.q_en).max().unwrap_or(0);
                let added_q_bases = q_max - q_min;
                for rec in q_recs {
                    rec.q_st = rec.q_st - q_min + scaffold_len;
                    rec.q_en = rec.q_en - q_min + scaffold_len;
                }
                scaffold_len += added_q_bases + spacer_size;
            }
            // remove padding insert on the end of rec
            scaffold_len -= spacer_size;

            for rec in t_recs {
                rec.q_name = scaffold_name.clone();
                rec.q_len = scaffold_len;
            }
        }
    }

    /// Identify overlapping pairs in Paf set
    pub fn overlapping_paf_recs(
        &mut self,
        match_score: i32,
        diff_score: i32,
        indel_score: i32,
        remove_contained: bool,
    ) {
        // remove trailing indels in all cases
        for rec in &mut self.records {
            rec.remove_trailing_indels();
        }

        let mut overlap_pairs = Vec::new();
        self.records.sort_by_key(|rec| rec.q_name.clone());
        let mut contained_indexes = vec![false; self.records.len()];

        // check if there are enough records to even try this operation
        if self.records.len() < 2 {
            return;
        }

        for i in 0..(self.records.len() - 1) {
            let rec1 = &self.records[i];
            let rgn1 = rec1.get_query_as_region();
            let mut j = i + 1;
            while j < self.records.len() && rec1.q_name == self.records[j].q_name {
                let rec2 = &self.records[j];
                let rgn2 = rec2.get_query_as_region();
                // count overlap
                let overlap = bed::get_overlap(&rgn1, &rgn2);
                // check if rec2 is contained
                if overlap < 1 {
                    j += 1;
                    continue;
                } else if overlap == (rec2.q_en - rec2.q_st) {
                    contained_indexes[j] = true;
                    log::debug!("{}\n^is contained in another alignment", rec2);
                } else if overlap == (rec1.q_en - rec1.q_st) {
                    contained_indexes[i] = true;
                    log::debug!("{}\n^is contained in another alignment", rec1);
                } else {
                    // put recs in left, right order
                    if rec1.q_st <= rec2.q_st {
                        overlap_pairs.push((overlap, i, j));
                    } else {
                        overlap_pairs.push((overlap, j, i));
                    }
                }
                // go to next
                j += 1;
            }
        }
        overlap_pairs.sort_by_key(|rec| u64::MAX - rec.0);
        log::debug!("{} overlapping pairs found", overlap_pairs.len());
        let mut q_seen: HashSet<String> = HashSet::new();
        let mut unseen = 0;
        for (_overlap, i, j) in overlap_pairs {
            let mut left = self.records[i].clone();
            let mut right = self.records[j].clone();
            let q_name = left.q_name.clone();
            // if we have not seen the q_name before it cannot be
            // in conflict with previous trimming steps
            if !q_seen.contains(&q_name) {
                trim_overlapping_pafs(&mut left, &mut right, match_score, diff_score, indel_score);
                log::trace!("{}", left);
                log::trace!("{}", right);
                self.records[i] = left;
                self.records[j] = right;
                q_seen.insert(q_name);
            } else {
                unseen += 1;
            }
        }

        if unseen > 0 {
            // recursively call for next overlap
            self.overlapping_paf_recs(match_score, diff_score, indel_score, remove_contained);
        } else if remove_contained {
            let n_to_remove = contained_indexes.iter().filter(|&x| *x).count();
            log::info!("Removing {} contained alignments.", n_to_remove);
            log::info!("{} total alignments.", self.records.len());
            let mut new_records = Vec::new();
            assert!(self.records.len() == contained_indexes.len());
            for (i, rec) in self.records.iter().enumerate() {
                if !contained_indexes[i] {
                    new_records.push(rec.clone());
                }
            }
            self.records = new_records;
            log::info!("{} total alignments.", self.records.len());
            // remove contained records
            //self.records.retain(|rec| !rec.contained);
        }
    }

    /// Make a SAM header from a Paf
    /// # Example
    /// ```
    /// use rustybam::paf;
    /// use std::fs::File;
    /// use std::io::*;
    /// let mut paf = paf::Paf::from_file(".test/asm_small.paf");
    /// let header = paf.sam_header();
    /// assert_eq!(header[0..3], "@HD".to_string());
    /// assert_eq!(header.split("\n").count(), 5);
    /// ```
    pub fn sam_header(&self) -> String {
        /*
        @HD	VN:1.6	SO:coordinate
        @SQ	SN:chr1	LN:248387497
        ...
        @SQ	SN:chrM	LN:16569
        @SQ	SN:chrY	LN:57227415
        @PG	ID:unimap	PN:unimap	VN:0.1-r41	CL:
        */
        let mut header = "@HD\tVN:1.6\n".to_string();

        // sort names naturally
        let mut names: Vec<(String, u64)> = self
            .records
            .iter()
            .map(|rec| (rec.t_name.clone(), rec.t_len))
            .unique()
            .collect();

        names.sort_by(|a, b| natord::compare(&a.0, &b.0));
        for (name, length) in names {
            header.push_str(&format!("@SQ\tSN:{}\tLN:{}\n", name, length));
        }
        header.push_str("@PG\tID:rustybam\tPN:rustybam");
        header
    }
}

/// A single cs-tag operation, stored in 1:1 correspondence with CIGAR ops.
#[derive(Debug, Clone)]
pub enum CsOp {
    /// `:N` — N matching bases (compact form, no sequence)
    Matches(u32),
    /// `=ACGT` — matching bases with explicit sequence
    MatchSeq(Vec<u8>),
    /// `*xy` — single-base mismatch (ref base, query base)
    Mismatch(u8, u8),
    /// `+acgt` — insertion of query bases
    Insertion(Vec<u8>),
    /// `-acgt` — deletion of reference bases
    Deletion(Vec<u8>),
}

impl CsOp {
    /// Trim this op to keep only bases at [skip..skip+keep).
    pub fn trim(&self, skip: u32, keep: u32) -> CsOp {
        match self {
            CsOp::Matches(_) => CsOp::Matches(keep),
            CsOp::MatchSeq(s) => {
                CsOp::MatchSeq(s[skip as usize..(skip + keep) as usize].to_vec())
            }
            CsOp::Mismatch(r, q) => {
                debug_assert!(skip == 0 && keep == 1);
                CsOp::Mismatch(*r, *q)
            }
            CsOp::Insertion(s) => {
                CsOp::Insertion(s[skip as usize..(skip + keep) as usize].to_vec())
            }
            CsOp::Deletion(s) => {
                CsOp::Deletion(s[skip as usize..(skip + keep) as usize].to_vec())
            }
        }
    }
}

impl fmt::Display for CsOp {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            CsOp::Matches(n) => write!(f, ":{}", n),
            CsOp::MatchSeq(s) => write!(f, "={}", std::str::from_utf8(s).unwrap()),
            CsOp::Mismatch(r, q) => write!(f, "*{}{}", *r as char, *q as char),
            CsOp::Insertion(s) => write!(f, "+{}", std::str::from_utf8(s).unwrap()),
            CsOp::Deletion(s) => write!(f, "-{}", std::str::from_utf8(s).unwrap()),
        }
    }
}

#[derive(Debug, Clone)]
pub struct PafRecord {
    pub q_name: String,
    pub q_len: u64,
    pub q_st: u64,
    pub q_en: u64,
    pub strand: char,
    pub t_name: String,
    pub t_len: u64,
    pub t_st: u64,
    pub t_en: u64,
    pub nmatch: u64,
    pub aln_len: u64,
    pub mapq: u64,
    pub cigar: CigarString,
    /// Parsed cs-tag ops, in 1:1 correspondence with `cigar`.
    /// `None` when input was CIGAR (cg tag); `Some` when input was cs tag.
    pub cs_ops: Option<Vec<CsOp>>,
    pub tags: String,
    pub id: String,
    pub order: u64,
    pub contained: bool,
}

impl PafRecord {
    /// # Example
    /// ```
    /// use rustybam::paf;
    /// let _paf = paf::PafRecord::new("A 1 2 3 + B 1 2 3 10 11 60").unwrap();
    /// let rec = paf::make_fake_paf_rec();
    /// assert_eq!("4M1I1D3=", rec.cigar.to_string());
    ///
    /// ```
    pub fn new(line: &str) -> PafResult<PafRecord> {
        let t: Vec<&str> = line.split_ascii_whitespace().collect();
        assert!(t.len() >= 12); // must have all required columns
                                // collect all the tags if any
        let mut tags = "".to_string();
        // find the cigar if it is there
        let mut cigar = CigarString(vec![]);
        let mut cs_ops: Option<Vec<CsOp>> = None;
        for token in t.iter().skip(12) {
            assert!(PAF_TAG.is_match(token));
            let caps = PAF_TAG.captures(token).unwrap();
            let tag = &caps[1];
            let value = &caps[3];
            if tag == "cg" {
                // cg tag takes precedence — parse CIGAR, discard any cs
                log::trace!("parsing cigar of length: {}", value.len());
                cigar =
                    CigarString::try_from(value.as_bytes()).expect("Unable to parse cigar string.");
                cs_ops = None;
            } else if tag == "cs" && cigar.is_empty() {
                // Only use cs tag when no cg tag is present
                let parsed = parse_cs_string(value)?;
                cigar = parsed.0;
                cs_ops = Some(parsed.1);
            } else {
                tags.push('\t');
                tags.push_str(token);
            }
        }

        // make the record
        let rec = PafRecord {
            q_name: t[0].to_string(),
            q_len: t[1].parse::<u64>().map_err(|_| Error::ParsePafColumn {})?,
            q_st: t[2].parse::<u64>().map_err(|_| Error::ParsePafColumn {})?,
            q_en: t[3].parse::<u64>().map_err(|_| Error::ParsePafColumn {})?,
            strand: t[4].parse::<char>().map_err(|_| Error::ParsePafColumn {})?,
            t_name: t[5].to_string(),
            t_len: t[6].parse::<u64>().map_err(|_| Error::ParsePafColumn {})?,
            t_st: t[7].parse::<u64>().map_err(|_| Error::ParsePafColumn {})?,
            t_en: t[8].parse::<u64>().map_err(|_| Error::ParsePafColumn {})?,
            nmatch: t[9].parse::<u64>().map_err(|_| Error::ParsePafColumn {})?,
            aln_len: t[10].parse::<u64>().map_err(|_| Error::ParsePafColumn {})?,
            mapq: t[11].parse::<u64>().map_err(|_| Error::ParsePafColumn {})?,
            cigar,
            cs_ops,
            tags,
            id: "".to_string(),
            order: 0,
            contained: false,
        };
        Ok(rec)
    }

    #[must_use]
    pub fn small_copy(&self) -> PafRecord {
        PafRecord {
            q_name: self.q_name.clone(),
            q_len: self.q_len,
            q_st: self.q_st,
            q_en: self.q_en,
            strand: self.strand,
            t_name: self.t_name.clone(),
            t_len: self.t_len,
            t_st: self.t_st,
            t_en: self.t_en,
            nmatch: self.nmatch,
            aln_len: self.aln_len,
            mapq: self.mapq,
            cigar: CigarString(Vec::new()),
            cs_ops: None,
            tags: self.tags.clone(),
            id: self.id.clone(),
            order: self.order,
            contained: self.contained,
        }
    }

    /// This function returns a region type from the paf query
    pub fn get_query_as_region(&self) -> bed::Region {
        bed::Region {
            name: self.q_name.clone(),
            st: self.q_st,
            en: self.q_en,
            ..Default::default()
        }
    }

    /// This function returns a region type from the paf target
    /// Example:
    /// ```
    /// use rustybam::bed::*;
    /// use rustybam::paf::*;
    /// let paf = PafRecord::new("Q 10 0 10 + T 20 12 20 3 9 60 cg:Z:7=1X2=").unwrap().get_target_as_region();
    /// let rgn = Region {name: "T".to_string(), st: 12, en: 20, ..Default::default()};
    /// assert_eq!((rgn.name, rgn.st, rgn.en),
    ///             (paf.name, paf.st, paf.en)
    /// );
    /// ```
    pub fn get_target_as_region(&self) -> bed::Region {
        bed::Region {
            name: self.t_name.clone(),
            st: self.t_st,
            en: self.t_en,
            ..Default::default()
        }
    }

    pub fn collapse_long_cigar(cigar: &CigarString) -> CigarString {
        let mut rtn = Vec::new();
        let mut pre_opt = cigar.0[0];
        let mut pre_len = 1;
        let mut idx = 1;
        while idx < cigar.len() {
            let cur_opt = cigar.0[idx];
            if std::mem::discriminant(&cur_opt) == std::mem::discriminant(&pre_opt) {
                pre_len += 1;
            } else {
                rtn.push(update_cigar_opt_len(&pre_opt, pre_len));
                pre_opt = cur_opt;
                pre_len = 1;
            }
            idx += 1;
        }
        rtn.push(update_cigar_opt_len(&pre_opt, pre_len));
        CigarString(rtn)
    }

    pub fn paf_overlaps_rgn(&self, rgn: &bed::Region) -> bool {
        if self.t_name != rgn.name {
            return false;
        }
        self.t_en > rgn.st && self.t_st < rgn.en
    }

    /// Return a tuple with the n of bases in the query and
    /// target inferred from the cigar string
    pub fn infer_n_bases(&mut self) -> (u64, u64, u64, u64) {
        let mut t_bases = 0;
        let mut q_bases = 0;
        let mut n_matches = 0;
        let mut aln_len = 0;
        for opt in self.cigar.into_iter() {
            if consumes_reference(opt) {
                t_bases += opt.len()
            }
            if consumes_query(opt) {
                q_bases += opt.len()
            }
            if is_match(opt) {
                n_matches += opt.len()
            }
            aln_len += opt.len();
        }
        (
            t_bases as u64,
            q_bases as u64,
            n_matches as u64,
            aln_len as u64,
        )
    }

    pub fn remove_trailing_indels(&mut self) {
        // if we check integrity, here it may fail due to indels on the ends
        // self.check_integrity().unwrap();

        let cigar_len = self.cigar.len();

        // find start to trim
        let mut st_opt = *self.cigar.first().unwrap();
        let mut remove_st_t = 0;
        let mut remove_st_q = 0;
        let mut remove_st_opts = 0;
        let mut removed_st_opts = Vec::new();
        while matches!(st_opt, Ins(_) | Del(_)) {
            if matches!(st_opt, Del(_)) {
                // consumes reference
                remove_st_t += st_opt.len();
                // TODO learn why I need this
                remove_st_q += 1;
            } else {
                remove_st_q += st_opt.len();
                // TODO learn why I need this
                // TODO handle the case when it is a del and then and ins
                // remove_st_t += 1;
            }
            remove_st_opts += 1;
            removed_st_opts.push(st_opt);
            if remove_st_opts < cigar_len {
                st_opt = self.cigar[remove_st_opts];
            } else {
                break;
            }
        }

        // remove extra counts put in my the case of Del followed by Ins
        if removed_st_opts.len() > 1 {
            for i in 0..(removed_st_opts.len() - 1) {
                let pre_opt = removed_st_opts[i];
                let cur_opt = removed_st_opts[i + 1];
                if (matches!(pre_opt, Del(_)) && matches!(cur_opt, Ins(_)))
                    || (matches!(pre_opt, Ins(_)) && matches!(cur_opt, Del(_)))
                {
                    remove_st_t += 1;
                    remove_st_q -= 1;
                }
            }
        }

        // find ends to trim
        let mut en_opt = *self.cigar.last().unwrap();
        let mut remove_en_t = 0;
        let mut remove_en_q = 0;
        let mut remove_en_opts = 0;
        let mut removed_en_opts = Vec::new();
        while matches!(en_opt, Ins(_) | Del(_)) {
            if matches!(en_opt, Del(_)) {
                // consumes reference
                remove_en_t += en_opt.len();
            } else {
                remove_en_q += en_opt.len();
            }
            remove_en_opts += 1;
            removed_en_opts.push(en_opt);
            if cigar_len - remove_en_opts > 0 {
                en_opt = self.cigar[cigar_len - 1 - remove_en_opts];
            } else {
                break;
            }
        }

        // log that we did something
        if remove_en_opts > 0 || remove_st_opts > 0 {
            self.id += &format!(
                "_TO.{}.{}",
                CigarString(removed_st_opts.clone()),
                CigarString(removed_en_opts.clone())
            );
        }

        // some logging.
        if remove_en_opts > 0 || remove_st_opts > 0 {
            log::debug!(
                "\nRemoved {} leading and {} trailing indels:\n{}\n{}\ntarget changes:{},{}\nquery changes: {},{}\n{}:{}-{}\n{}:{}-{}",
                remove_st_opts,
                remove_en_opts,
                CigarString(removed_st_opts),
                CigarString(removed_en_opts),
                remove_st_t,
                remove_en_t,
                remove_st_q,
                remove_en_q,
                self.q_name,
                self.q_st,
                self.q_en,
                self.t_name,
                self.t_st,
                self.t_en,
            );
        }

        // update the cigar string (and cs_ops if present)
        self.cigar = CigarString(self.cigar.0[remove_st_opts..].to_vec());
        self.cigar.0.truncate(self.cigar.len() - remove_en_opts);
        if let Some(ref mut cs) = self.cs_ops {
            *cs = cs[remove_st_opts..].to_vec();
            cs.truncate(cs.len() - remove_en_opts);
        }

        // update the target coordinates
        self.t_st += remove_st_t as u64;
        self.t_en -= remove_en_t as u64;

        // update the query coordinates if rc
        if self.strand == '-' {
            std::mem::swap(&mut remove_st_q, &mut remove_en_q);
        }
        // fix the query positions that need to be
        self.q_st += remove_st_q as u64;
        self.q_en -= remove_en_q as u64;

        // check we removed the indels
        if !self.cigar.is_empty() {
            let st_opt = *self.cigar.first().unwrap();
            let en_opt = *self.cigar.last().unwrap();
            if matches!(st_opt, Ins(_) | Del(_)) || matches!(en_opt, Ins(_) | Del(_)) {
                eprintln!("Why are there still indels?\n{}", self);
                //self.remove_trailing_indels();
            }
        }

        // make sure we did not break the cigar
        self.check_integrity().unwrap();
    }

    /// Truncate this record to keep only the portion within [new_q_st, new_q_en)
    /// in query coordinates. Walks compressed CIGAR ops directly — O(n_cigar_ops).
    pub fn truncate_record_by_query(&mut self, new_q_st: u64, new_q_en: u64) {
        // checks
        assert!(new_q_st >= self.q_st, "New start is less than old start.");
        assert!(new_q_en <= self.q_en, "New end is greater than old end.");

        if new_q_st == self.q_st && new_q_en == self.q_en {
            return;
        }

        let cs_vec = self.cs_ops.take(); // take ownership for parallel processing
        let mut new_cs_ops: Option<Vec<CsOp>> = cs_vec.as_ref().map(|_| Vec::new());

        let mut q_pos = if self.strand == '+' { self.q_st } else { self.q_en };
        let mut new_cigar_ops: Vec<Cigar> = Vec::new();
        let mut t_before: u64 = 0;
        let mut q_consumed_before: u64 = 0;
        let mut q_consumed_in: u64 = 0;
        let mut started = false;
        let mut finished = false;

        for (ci, op) in self.cigar.into_iter().enumerate() {
            if finished {
                break;
            }
            let op_len = op.len() as u64;
            let moves_t = consumes_reference(op);
            let moves_q = consumes_query(op);

            if !moves_q {
                // Deletion or ref skip — doesn't advance query
                if started {
                    new_cigar_ops.push(*op);
                    if let (Some(ref mut ncs), Some(ref cs)) = (&mut new_cs_ops, &cs_vec) {
                        ncs.push(cs[ci].clone());
                    }
                } else if moves_t {
                    t_before += op_len;
                }
                continue;
            }

            // Query-consuming op: compute absolute query range [q_lo, q_hi)
            let (q_lo, q_hi) = if self.strand == '+' {
                let lo = q_pos;
                q_pos += op_len;
                (lo, q_pos)
            } else {
                q_pos -= op_len;
                (q_pos, q_pos + op_len)
            };

            // Check overlap with [new_q_st, new_q_en)
            let overlap_start = std::cmp::max(q_lo, new_q_st);
            let overlap_end = std::cmp::min(q_hi, new_q_en);
            if overlap_start >= overlap_end {
                // No overlap
                if !started {
                    if moves_t {
                        t_before += op_len;
                    }
                    q_consumed_before += op_len;
                }
                continue;
            }

            // Bases to skip at the CIGAR-start side of this op
            let skip_before_cigar = if self.strand == '+' {
                overlap_start - q_lo
            } else {
                q_hi - overlap_end
            };
            let keep = overlap_end - overlap_start;

            if !started {
                started = true;
                if moves_t {
                    t_before += skip_before_cigar;
                }
                q_consumed_before += skip_before_cigar;
            }

            q_consumed_in += keep;
            new_cigar_ops.push(update_cigar_opt_len(op, keep as u32));
            if let (Some(ref mut ncs), Some(ref cs)) = (&mut new_cs_ops, &cs_vec) {
                ncs.push(cs[ci].trim(skip_before_cigar as u32, keep as u32));
            }

            // Check if we've reached the far boundary
            if (self.strand == '+' && overlap_end >= new_q_en)
                || (self.strand == '-' && overlap_start <= new_q_st)
            {
                finished = true;
            }
        }

        if new_cigar_ops.is_empty() {
            self.cs_ops = new_cs_ops;
            return;
        }

        // Remove leading indels and adjust coordinates
        let mut lead_t: u64 = 0;
        let mut lead_q: u64 = 0;
        while !new_cigar_ops.is_empty() {
            let first = new_cigar_ops[0];
            if matches!(first, Match(_) | Equal(_) | Diff(_)) {
                break;
            }
            if consumes_reference(&first) {
                lead_t += first.len() as u64;
            }
            if consumes_query(&first) {
                lead_q += first.len() as u64;
            }
            new_cigar_ops.remove(0);
            if let Some(ref mut ncs) = new_cs_ops {
                ncs.remove(0);
            }
        }
        t_before += lead_t;
        q_consumed_before += lead_q;
        q_consumed_in -= lead_q;

        // Remove trailing indels
        let mut trail_q: u64 = 0;
        while !new_cigar_ops.is_empty() {
            let last = *new_cigar_ops.last().unwrap();
            if matches!(last, Match(_) | Equal(_) | Diff(_)) {
                break;
            }
            if consumes_query(&last) {
                trail_q += last.len() as u64;
            }
            new_cigar_ops.pop();
            if let Some(ref mut ncs) = new_cs_ops {
                ncs.pop();
            }
        }
        q_consumed_in -= trail_q;

        if new_cigar_ops.is_empty() {
            self.cs_ops = new_cs_ops;
            return;
        }

        // Compute new target coordinates
        self.t_st += t_before;
        let t_bases: u64 = new_cigar_ops
            .iter()
            .filter(|op| consumes_reference(op))
            .map(|op| op.len() as u64)
            .sum();
        self.t_en = self.t_st + t_bases;

        // Set query coordinates
        if self.strand == '+' {
            self.q_st += q_consumed_before;
            self.q_en = self.q_st + q_consumed_in;
        } else {
            self.q_en -= q_consumed_before;
            self.q_st = self.q_en - q_consumed_in;
        }

        self.cigar = CigarString(new_cigar_ops);
        self.cs_ops = new_cs_ops;

        // check integrity and update aln_len and nmatch
        self.check_integrity().unwrap();
    }

    pub fn check_integrity(&mut self) -> PafResult<()> {
        let (t_bases, q_bases, nmatch, aln_len) = self.infer_n_bases();
        if self.t_en - self.t_st != t_bases {
            return Err(Error::PafParseCigar {
                msg: format!(
                    "target bases {} from cigar does not equal {}-{}={}\n{}\n",
                    t_bases,
                    self.t_en,
                    self.t_st,
                    self.t_en - self.t_st,
                    self
                ),
            });
        }
        if self.q_en - self.q_st != q_bases {
            return Err(Error::PafParseCigar {
                msg: format!(
                    "query bases {} from cigar does not equal {}-{}={}\n{}\n",
                    q_bases,
                    self.q_en,
                    self.q_st,
                    self.q_en - self.q_st,
                    self
                ),
            });
        }

        // update other fields
        self.nmatch = nmatch;
        self.aln_len = aln_len;

        Ok(())
    }

    /// Print the paf record as a SAM record
    /// Example:
    /// ```
    /// use rustybam::bed::*;
    /// use rustybam::paf::*;
    /// let paf = PafRecord::new("Q 10 0 10 + T 20 12 20 3 9 60 cg:Z:7=1X2=").unwrap();
    /// let sam = paf.to_sam_string(None);
    /// ```
    pub fn to_sam_string(&self, reader: Option<&faidx::Reader>) -> String {
        /*
        m64062_190807_194840/133628256/ccs	0	chr1	1	60	396=	*	0	0   *   *
        */
        let mut clip_char = 'H';
        let seq = match reader {
            Some(reader) => {
                let seq = getfasta::fetch_fasta(
                    reader,
                    &self.q_name,
                    0, //self.q_st as usize,
                    self.q_len as usize,
                );
                clip_char = 'S';
                let seq = if self.strand == '-' {
                    revcomp(seq)
                } else {
                    seq
                };
                std::str::from_utf8(&seq).unwrap().to_string()
            }
            None => "*".to_string(),
        };
        let qual = "*".to_string();
        let flag = if self.strand == '-' { 16 } else { 0 };
        let mut leading_clip = if self.q_st > 0 {
            format!("{}{}", self.q_st, clip_char)
        } else {
            "".to_string()
        };
        let mut trailing_clip = if self.q_len - self.q_en > 0 {
            format!("{}{}", self.q_len - self.q_en, clip_char)
        } else {
            "".to_string()
        };
        if self.strand == '-' {
            std::mem::swap(&mut leading_clip, &mut trailing_clip);
        }
        let o_cigar = format!("{}{}{}", leading_clip, self.cigar, trailing_clip);
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.q_name,
            flag,
            self.t_name,
            self.t_st + 1,
            self.mapq,
            o_cigar,
            "*",
            0,
            0,
            seq,
            qual
        )
    }
}

impl fmt::Display for PafRecord {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tid:Z:{}\tcg:Z:{}",
            self.q_name,
            self.q_len,
            self.q_st,
            self.q_en,
            self.strand,
            self.t_name,
            self.t_len,
            self.t_st,
            self.t_en,
            self.nmatch,
            self.aln_len,
            self.mapq,
            self.id,
            self.cigar,
        )?;
        // Also output cs tag if present
        if let Some(ref cs_ops) = self.cs_ops {
            write!(f, "\tcs:Z:")?;
            for op in cs_ops {
                write!(f, "{}", op)?;
            }
        }
        Ok(())
    }
}

pub fn consumes_reference(cigar_opt: &Cigar) -> bool {
    matches!(
        cigar_opt,
        Match(_i) | Del(_i) | RefSkip(_i) | Diff(_i) | Equal(_i)
    )
}
/// # Example
/// ```
/// use rustybam::paf;
/// use rust_htslib::bam::record::Cigar::*;
/// assert!(paf::consumes_query(&Diff(5)));
/// ```
pub fn consumes_query(cigar_opt: &Cigar) -> bool {
    matches!(
        cigar_opt,
        Match(_i) | Ins(_i) | SoftClip(_i) | Diff(_i) | Equal(_i)
    )
}

/// # Example
/// ```
/// use rustybam::paf;
/// use rust_htslib::bam::record::Cigar::*;
/// assert!(paf::is_match(&Match(5)));
/// assert!(paf::is_match(&Diff(5)));
/// assert!(paf::is_match(&Equal(5)));
/// ```
pub fn is_match(cigar_opt: &Cigar) -> bool {
    matches!(cigar_opt, Match(_i) | Diff(_i) | Equal(_i))
}

/// # Example
/// ```
/// use rustybam::paf;
/// use rust_htslib::bam::record::Cigar::*;
/// assert_eq!(Diff(5), paf::update_cigar_opt_len(&Diff(10), 5));
/// assert_eq!(Diff(10), paf::update_cigar_opt_len(&Diff(1), 10));
/// ```
pub fn update_cigar_opt_len(opt: &Cigar, new_opt_len: u32) -> Cigar {
    match opt {
        Match(_) => Match(new_opt_len),
        Ins(_) => Ins(new_opt_len),
        Del(_) => Del(new_opt_len),
        RefSkip(_) => RefSkip(new_opt_len),
        HardClip(_) => HardClip(new_opt_len),
        SoftClip(_) => SoftClip(new_opt_len),
        Pad(_) => Pad(new_opt_len),
        Equal(_) => Equal(new_opt_len),
        Diff(_) => Diff(new_opt_len),
    }
}

/// Create a CigarString from given str.
/// # Example
/// ```
/// use rustybam::paf;
/// use rust_htslib::bam::record::*;
/// use rust_htslib::bam::record::CigarString;
/// use rust_htslib::bam::record::Cigar::*;
/// use std::convert::TryFrom;
/// use std::str::FromStr;
/// let cigars = vec!["10M4D100I1102=", "100000M20=5P10X4M"];
/// for cigar_str in cigars{
///     let my_parse = paf::cigar_from_str(cigar_str).expect("Unable to parse cigar");
///     let hts_parse = CigarString::try_from(cigar_str).expect("Unable to parse cigar");
///     assert_eq!(my_parse, hts_parse);
/// }
/// ```
pub fn cigar_from_str(text: &str) -> PafResult<CigarString> {
    let bytes = text.as_bytes();
    let mut inner = Vec::new();
    let mut i = 0;
    let text_len = text.len();
    while i < text_len {
        let mut j = i;
        while j < text_len && bytes[j].is_ascii_digit() {
            j += 1;
        }
        let n = u32::from_str(&text[i..j]).map_err(|_| Error::PafParseCigar {
            msg: "expected integer".to_owned(),
        })?;
        let op = &text[j..j + 1];
        inner.push(match op {
            "M" => Cigar::Match(n),
            "I" => Cigar::Ins(n),
            "D" => Cigar::Del(n),
            "N" => Cigar::RefSkip(n),
            "H" => Cigar::HardClip(n),
            "S" => Cigar::SoftClip(n),
            "P" => Cigar::Pad(n),
            "=" => Cigar::Equal(n),
            "X" => Cigar::Diff(n),
            op => {
                return Err(Error::PafParseCigar {
                    msg: format!("Cannot parse opt: {}", op),
                })
            }
        });
        i = j + 1;
    }
    Ok(CigarString(inner))
}

/// Basically swaps the query and the reference in a cigar
pub fn cigar_swap_target_query(cigar: &CigarString, strand: char) -> CigarString {
    // flip cigar
    let mut new_cigar = Vec::new();
    for opt in cigar.into_iter() {
        let new_opt = match opt {
            Ins(l) => Del(*l),
            Del(l) => Ins(*l),
            _ => *opt,
        };
        new_cigar.push(new_opt);
    }
    if strand == '-' {
        new_cigar.reverse();
    }
    CigarString(new_cigar)
}

/// Swaps the query and reference and inverts the cigar sting
pub fn paf_swap_query_and_target(paf: &PafRecord) -> PafRecord {
    let mut flipped = paf.clone();
    // flip target
    flipped.t_name = paf.q_name.clone();
    flipped.t_len = paf.q_len;
    flipped.t_st = paf.q_st;
    flipped.t_en = paf.q_en;
    // flip query
    flipped.q_name = paf.t_name.clone();
    flipped.q_len = paf.t_len;
    flipped.q_st = paf.t_st;
    flipped.q_en = paf.t_en;

    // flip the cigar
    flipped.cigar = cigar_swap_target_query(&paf.cigar, paf.strand);

    // flip cs_ops: swap Insertion ↔ Deletion
    flipped.cs_ops = paf.cs_ops.as_ref().map(|ops| {
        let mut new_ops: Vec<CsOp> = ops
            .iter()
            .map(|op| match op {
                CsOp::Insertion(s) => CsOp::Deletion(s.clone()),
                CsOp::Deletion(s) => CsOp::Insertion(s.clone()),
                other => other.clone(),
            })
            .collect();
        if paf.strand == '-' {
            new_ops.reverse();
        }
        new_ops
    });

    flipped
}

pub fn make_fake_paf_rec() -> PafRecord {
    PafRecord::new("Q 10 2 10 - T 20 12 20 3 9 60 cg:Z:4M1I1D3=").unwrap()
}

/// Parse a cs-tag string into both a CigarString and a Vec<CsOp>.
/// Supports both long form (`=ACGT`, `*at`) and short form (`:10`).
///
/// # Example
/// ```
/// use rust_htslib::bam::record::Cigar::*;
/// use rustybam::paf;
/// let (cigar, cs_ops) = paf::parse_cs_string(":10=ACGTN+acgtn-acgtn*at=A").unwrap();
/// assert_eq!(cigar[0], Equal(10));
/// assert_eq!(cigar[1], Equal(5));
/// assert_eq!(cigar[2], Ins(5));
/// assert_eq!(cigar[3], Del(5));
/// assert_eq!(cigar[4], Diff(1));
/// assert_eq!(cigar[5], Equal(1));
/// assert_eq!(cs_ops.len(), cigar.len());
/// ```
pub fn parse_cs_string(cs: &str) -> PafResult<(CigarString, Vec<CsOp>)> {
    let bytes = cs.as_bytes();
    let length = bytes.len();
    let mut i = 0;
    let mut cigar = vec![];
    let mut cs_ops = vec![];
    while i < length {
        let cs_opt = bytes[i];
        i += 1; // past the operator
        match cs_opt {
            b'=' => {
                let start = i;
                while i < length && matches!(bytes[i], b'A' | b'C' | b'G' | b'T' | b'N') {
                    i += 1;
                }
                let seq = bytes[start..i].to_vec();
                let l = seq.len() as u32;
                cigar.push(Cigar::Equal(l));
                cs_ops.push(CsOp::MatchSeq(seq));
            }
            b':' => {
                let start = i;
                while i < length && bytes[i].is_ascii_digit() {
                    i += 1;
                }
                let l = u32::from_str(&cs[start..i]).map_err(|_| Error::ParseIntError {
                    msg: "Expected integer in cs :N operator".to_string(),
                })?;
                cigar.push(Cigar::Equal(l));
                cs_ops.push(CsOp::Matches(l));
            }
            b'*' => {
                let ref_base = bytes[i];
                let query_base = bytes[i + 1];
                i += 2;
                cigar.push(Cigar::Diff(1));
                cs_ops.push(CsOp::Mismatch(ref_base, query_base));
            }
            b'+' => {
                let start = i;
                while i < length && bytes[i].is_ascii_lowercase() {
                    i += 1;
                }
                let seq = bytes[start..i].to_vec();
                let l = seq.len() as u32;
                cigar.push(Cigar::Ins(l));
                cs_ops.push(CsOp::Insertion(seq));
            }
            b'-' => {
                let start = i;
                while i < length && bytes[i].is_ascii_lowercase() {
                    i += 1;
                }
                let seq = bytes[start..i].to_vec();
                let l = seq.len() as u32;
                cigar.push(Cigar::Del(l));
                cs_ops.push(CsOp::Deletion(seq));
            }
            b'~' => {
                return Err(Error::PafParseCS {
                    msg: "Splice operations not yet supported.".to_string(),
                });
            }
            _ => {
                return Err(Error::PafParseCS {
                    msg: format!("Unexpected operator in the cs string: {}", cs_opt as char),
                });
            }
        }
    }
    Ok((CigarString(cigar), cs_ops))
}

/// Parse a cs-tag string into a CigarString (convenience wrapper).
///
/// # Example
/// ```
/// use rust_htslib::bam::record::Cigar::*;
/// use rustybam::paf;
/// let cigar = paf::cs_to_cigar(":10=ACGTN+acgtn-acgtn*at=A").unwrap();
/// assert_eq!(cigar[0], Equal(10));
/// assert_eq!(cigar[1], Equal(5));
/// assert_eq!(cigar[2], Ins(5));
/// assert_eq!(cigar[3], Del(5));
/// assert_eq!(cigar[4], Diff(1));
/// assert_eq!(cigar[5], Equal(1));
/// ```
pub fn cs_to_cigar(cs: &str) -> PafResult<CigarString> {
    let (cigar, _) = parse_cs_string(cs)?;
    Ok(cigar)
}
