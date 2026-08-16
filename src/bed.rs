use crate::myio;
use bio::io::bed;
use lazy_static::lazy_static;
use regex::Regex;
use std::cmp::{max, min};
use std::fmt;
use std::str;

lazy_static! {
    static ref BED_RE: Regex = Regex::new(r"([^\s]+)\t([0-9]+)\t([0-9]+)\t?([^\s]+)?.*").unwrap();
    static ref RGN_RE: Regex = Regex::new(r"(.+):([0-9]+)-([0-9]+)").unwrap();
}

#[derive(Default)]
pub struct Region {
    pub name: String,
    pub st: u64,
    pub en: u64,
    pub id: String,
    pub record: bed::Record,
}

impl Region {
    pub fn get_column(&self, column: u8) -> String {
        match column {
            1 => self.name.clone(),
            2 => self.st.to_string(),
            3 => self.en.to_string(),
            4 => self.record.name().unwrap_or("no-value").to_string(),
            5 => self.record.score().unwrap_or("no-value").to_string(),
            6 => self
                .record
                .strand()
                .unwrap_or(bio_types::strand::Strand::Unknown)
                .to_string(),
            _ => "no-value".to_string(),
        }
    }
}

impl fmt::Display for Region {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}:{}-{}", self.name, self.st + 1, self.en)
    }
}

/// Checks if two regions overlap
/// # Example
/// ```
/// use rustybam::bed::*;
/// let rgn1 = &parse_bed_rec("chr1\t10\t15");
/// let rgn2 = &parse_bed_rec("chr1\t15\t20");
/// let rgn3 = &parse_bed_rec("chr1\t5\t10");
/// let big = &parse_bed_rec("chr1\t0\t20");
/// let small = &parse_bed_rec("chr1\t11\t12");
/// let left = &parse_bed_rec("chr1\t8\t12");
/// let right = &parse_bed_rec("chr1\t14\t16");
///
/// assert_eq!(has_overlap(rgn1, rgn2), false);
/// assert_eq!(has_overlap(rgn1, rgn3), false);
/// assert_eq!(has_overlap(rgn1, big), true);
/// assert_eq!(has_overlap(rgn1, small), true);
/// assert_eq!(has_overlap(rgn1, left), true);
/// assert_eq!(has_overlap(rgn1, left), true);
/// ```
pub fn has_overlap(rgn1: &Region, rgn2: &Region) -> bool {
    if rgn1.name != rgn2.name {
        return false;
    }
    rgn1.en > rgn2.st && rgn1.st < rgn2.en
}

/// return the # of overlapping bases
pub fn get_overlap(rgn1: &Region, rgn2: &Region) -> u64 {
    //max(0, min(rgn1.en, rgn2.en) - max(rgn1.st, rgn2.st))
    if rgn1.name != rgn2.name {
        return 0;
    }
    let my_min = min(rgn1.en, rgn2.en);
    let my_max = max(rgn1.st, rgn2.st);
    if my_min < my_max {
        return 0;
    }
    my_min - my_max
}

/// parse region strings
/// # Example
/// ```
/// let rgn = rustybam::bed::parse_region("chr1:1-1000");
/// assert_eq!("chr1", rgn.name);
/// assert_eq!(0, rgn.st);
/// assert_eq!(1000, rgn.en);
///
/// let rgn2 = rustybam::bed::parse_region("chr1:2-2000:1-1000");
/// assert_eq!("chr1:2-2000", rgn2.name);
/// ```
pub fn parse_region(region: &str) -> Region {
    let caps = RGN_RE
        .captures(region)
        .expect("Failed to parse region string.");

    let name = caps.get(1).unwrap().as_str().to_string();
    let st = caps.get(2).unwrap().as_str().parse::<u64>().unwrap() - 1;
    let en = caps.get(3).unwrap().as_str().parse().unwrap_or(4294967295); //this is 2^32-1
    let id = caps
        .get(4)
        .map_or(format!("{}:{}-{}", name, st + 1, en), |m| {
            m.as_str().to_string()
        });

    assert!(
        (st <= en),
        "Region start must be less than end.\n{}",
        region
    );

    Region {
        name,
        st,
        en,
        id,
        ..Default::default()
    }
}

/// parse bed strings
/// # Example
/// ```
/// let rgn = rustybam::bed::parse_bed_rec("chr1\t0\t1000\tid");
/// assert_eq!("chr1", rgn.name);
/// assert_eq!(0, rgn.st);
/// assert_eq!(1000, rgn.en);
/// assert_eq!("id", rgn.id);
///
/// let rgn2 = rustybam::bed::parse_bed_rec("chr1\t2\t2000");
/// assert_eq!("chr1", rgn2.name);
/// assert_eq!("chr1:3-2000", rgn2.id);
/// ```
pub fn parse_bed_rec(rec: &str) -> Region {
    let mut reader = bed::Reader::new(rec.as_bytes());
    let record = reader.records().next().unwrap().unwrap();
    parse_bed_record(record)
}

pub fn parse_bed_record(record: bed::Record) -> Region {
    let name = record.chrom().to_owned();
    let st = record.start();
    let en = record.end();
    let id = match record.name() {
        Some(x) => x.to_owned(),
        _ => format!("{}:{}-{}", name, st + 1, en),
    };
    Region {
        name,
        st,
        en,
        id,
        record,
    }
}

/// parse bed file
/// # Example
/// ```
/// use rustybam::bed::*;
/// let vec = parse_bed(".test/asm_small.bed");
/// assert_eq!(vec.len(), 10);
/// let vec = parse_bed(".test/asm_small.bed.gz");
/// assert_eq!(vec.len(), 10);
/// ```
pub fn parse_bed(filename: &str) -> Vec<Region> {
    let mut vec = Vec::new();
    let reader = myio::reader(filename);
    let mut records = bed::Reader::new(reader);
    let mut rec_num = 0;
    for (idx, rec) in records.records().enumerate() {
        match rec {
            Ok(r) => {
                let rgn = parse_bed_record(r);
                vec.push(rgn);
            }
            Err(e) => log::warn!(
                "Unable to parse bed at line {}, skipping. Reason: {}",
                idx + 1,
                e
            ),
        }
        rec_num += 1;
        log::debug!("Read bed record number {}", rec_num);
    }
    vec
}

/// # Example
/// ```
/// use rust_htslib::{bam, bam::Read};
/// let mut bam = bam::IndexedReader::from_path(".test/test_nucfreq.bam").unwrap();
/// let rgn = rustybam::bed::Region {
///     name : "CHROMOSOME_I".to_string(),
///     st :  0,
///     en : 95,
///     id : "None".to_string(),
///     ..Default::default()
/// };
/// let small_rgns = rustybam::bed::split_region(&rgn, 10);
/// assert_eq!(small_rgns[0].st, 0);
/// assert_eq!(small_rgns[0].en, 10);
/// assert_eq!(small_rgns[9].st, 90);
/// assert_eq!(small_rgns[9].en, 95);
/// let small_rgns_2 = rustybam::bed::split_region(&rgn, 100);
/// assert_eq!(small_rgns_2[0].st, 0);
/// assert_eq!(small_rgns_2[0].en, 95);
/// ```
pub fn split_region(rgn: &Region, window: u64) -> Vec<Region> {
    // make many smaller regions
    let mut start = rgn.st;
    let mut small_rgns = Vec::new();
    while start < rgn.en {
        let mut end = start + window;
        if end > rgn.en {
            end = rgn.en;
        }
        let tmprgn = Region {
            name: rgn.name.clone(),
            st: start,
            en: end,
            id: rgn.id.clone(),
            ..Default::default()
        };
        small_rgns.push(tmprgn);
        start = end;
    }
    small_rgns
}
