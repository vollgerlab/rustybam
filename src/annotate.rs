use bio::data_structures::interval_tree::{Entry, IntervalTree};
use bio::io::*;
use itertools::Itertools;

pub trait IntervalTreeExt<N: Ord + Clone, D> {
    fn find_bed_overlaps(&self, rec: &bed::Record) -> Vec<Entry<'_, u64, bed::Record>>;
}

impl IntervalTreeExt<u64, &bed::Record> for IntervalTree<u64, bed::Record> {
    fn find_bed_overlaps(&self, rec: &bed::Record) -> Vec<Entry<'_, u64, bed::Record>> {
        self.find(rec.start()..rec.end())
            .filter(|entry| entry.data().chrom() == rec.chrom())
            .collect_vec()
    }
}

pub fn interval_tree_from_bed_file(path: &str) -> IntervalTree<u64, bed::Record> {
    let mut tree = IntervalTree::new();
    let mut reader = bed::Reader::from_file(path).unwrap();
    for rec in reader.records() {
        let rec = rec.unwrap();
        tree.insert(rec.start()..rec.end(), rec.clone());
    }
    tree
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio::utils::Interval;

    #[test]
    fn test_chrom_specific_overlaps() {
        let mut tree = IntervalTree::new();
        let mut rec = bed::Record::new();
        rec.set_chrom("Range_1");
        rec.set_start(11);
        rec.set_end(20);

        let mut rec2 = bed::Record::new();
        rec2.set_chrom("Range_2");
        rec2.set_start(10);
        rec2.set_end(30);

        tree.insert(11..20, rec.clone());
        tree.insert(10..30, rec2.clone());

        for r in tree.find_bed_overlaps(&rec) {
            assert_eq!(r.interval(), &(Interval::from(rec.start()..rec.end())));
            assert_eq!(r.data().chrom(), "Range_1");
            eprintln!("found filtered range");
        }
        for r in tree.find(15..25) {
            if r.data().chrom() == "Range_1" {
                assert_eq!(r.interval(), &(Interval::from(rec.start()..rec.end())));
                assert_eq!(r.data().chrom(), "Range_1");
            } else if r.data().chrom() == "Range_2" {
                assert_eq!(r.interval(), &(Interval::from(rec2.start()..rec2.end())));
                assert_eq!(r.data().chrom(), "Range_2");
            }
        }
    }
}
