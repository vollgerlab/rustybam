#!/usr/bin/env bash
set -euo pipefail

# Compare rb liftover vs paftools.js liftover coordinates
#
# KEY INSIGHT: paftools.js liftover takes query.bed and lifts to target space.
# rb liftover (default) takes target BED and lifts to query space.
# So we use rb --qbed to match paftools behavior: both take query coords, output target coords.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
RUSTYBAM_DIR="$(dirname "$SCRIPT_DIR")"
PROJECT_DIR="$(dirname "$RUSTYBAM_DIR")"

PAF="$PROJECT_DIR/T2Tv2.0-to-GRCh38.no.eqx.gt50k.paf"
BED="$PROJECT_DIR/paper/benchmarks/results/liftover_regions.bed"
TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

echo "PAF: $PAF"
echo "BED: $BED ($(wc -l < "$BED") regions)"
echo "TMPDIR: $TMPDIR"
echo ""

# ---- Create BED with standardized 4th column for region IDs ----
awk -F'\t' 'NF >= 3 { print $1 "\t" $2 "\t" $3 "\t" $1 "_" $2 "_" $3 }' "$BED" > "$TMPDIR/regions.bed"
echo "Prepared $(wc -l < "$TMPDIR/regions.bed") regions with standardized IDs"

# ---- Run paftools.js liftover ----
# paftools interprets BED as QUERY coordinates, outputs TARGET coordinates
echo ""
echo "Running paftools.js liftover (query BED -> target coords)..."
paftools.js liftover -l 0 -q 0 "$PAF" "$TMPDIR/regions.bed" > "$TMPDIR/paftools_raw.bed" 2>/dev/null
echo "  paftools output: $(wc -l < "$TMPDIR/paftools_raw.bed") records"

# paftools output: target_name  t_start  t_end  region_id[_tN]  score  strand
awk -F'\t' '{
    id = $4
    sub(/_t[0-9]+$/, "", id)
    print id "\t" $1 "\t" $2 "\t" $3 "\t" $6
}' "$TMPDIR/paftools_raw.bed" | sort -k1,1 -k2,2 -k3,3n > "$TMPDIR/paftools.tsv"

# ---- Run rb liftover with --qbed (same direction as paftools) ----
# --qbed: interprets BED as QUERY coordinates, lifts to target space
# Output PAF swaps query/target: output query = original target (lifted coords)
echo ""
echo "Running rb liftover --qbed (query BED -> target coords)..."
cargo run --release --bin rb -- liftover --qbed --bed "$TMPDIR/regions.bed" "$PAF" > "$TMPDIR/rb_qbed_raw.paf" 2>/dev/null
echo "  rb --qbed output: $(wc -l < "$TMPDIR/rb_qbed_raw.paf") records"

# With --qbed, the output PAF has swapped query/target:
#   col 1-4: original target (= lifted coordinates)
#   col 6-9: original query (= input BED space)
# So columns 1,3,4 give us the lifted target coords (same as paftools output)
awk -F'\t' '{
    lifted_name = $1; lifted_st = $3; lifted_en = $4; strand = $5
    id = ""
    for (i = 13; i <= NF; i++) {
        if (substr($i, 1, 5) == "id:Z:") {
            id = substr($i, 6)
            break
        }
    }
    if (id == "") next
    print id "\t" lifted_name "\t" lifted_st "\t" lifted_en "\t" strand
}' "$TMPDIR/rb_qbed_raw.paf" | sort -k1,1 -k2,2 -k3,3n > "$TMPDIR/rb_qbed.tsv"

# ---- Also run rb liftover --qbed --largest ----
echo ""
echo "Running rb liftover --qbed --largest..."
cargo run --release --bin rb -- liftover --qbed --largest --bed "$TMPDIR/regions.bed" "$PAF" > "$TMPDIR/rb_qbed_largest.paf" 2>/dev/null
echo "  rb --qbed --largest output: $(wc -l < "$TMPDIR/rb_qbed_largest.paf") records"

awk -F'\t' '{
    lifted_name = $1; lifted_st = $3; lifted_en = $4; strand = $5
    id = ""
    for (i = 13; i <= NF; i++) {
        if (substr($i, 1, 5) == "id:Z:") {
            id = substr($i, 6)
            break
        }
    }
    if (id == "") next
    print id "\t" lifted_name "\t" lifted_st "\t" lifted_en "\t" strand
}' "$TMPDIR/rb_qbed_largest.paf" | sort -k1,1 -k2,2 -k3,3n > "$TMPDIR/rb_qbed_largest.tsv"

echo ""
echo "=== Record counts ==="
echo "  paftools:             $(wc -l < "$TMPDIR/paftools.tsv")"
echo "  rb --qbed all:        $(wc -l < "$TMPDIR/rb_qbed.tsv")"
echo "  rb --qbed --largest:  $(wc -l < "$TMPDIR/rb_qbed_largest.tsv")"

# ---- Show sample output for sanity check ----
echo ""
echo "=== Sanity check: first region from each tool ==="
echo "paftools (first 3):"
head -3 "$TMPDIR/paftools.tsv" | column -t
echo "rb --qbed (first 3):"
head -3 "$TMPDIR/rb_qbed.tsv" | column -t

# ---- Compare ----
echo ""
echo "=== Coordinate comparison ==="

python3 - "$TMPDIR/paftools.tsv" "$TMPDIR/rb_qbed.tsv" "$TMPDIR/rb_qbed_largest.tsv" <<'PYEOF'
import sys
from collections import defaultdict

def load_records(path):
    """Load records as dict: (region, chrom, strand) -> list of (start, end)"""
    recs = defaultdict(list)
    with open(path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            key = (parts[0], parts[1], parts[4])
            recs[key].append((int(parts[2]), int(parts[3])))
    for k in recs:
        recs[k].sort()
    return recs

def compare(name_a, recs_a, name_b, recs_b):
    keys_a = set(recs_a.keys())
    keys_b = set(recs_b.keys())
    shared = keys_a & keys_b

    print(f"  {name_a} keys: {len(keys_a)},  {name_b} keys: {len(keys_b)},  Shared: {len(shared)}")
    print(f"  Only {name_a}: {len(keys_a - keys_b)},  Only {name_b}: {len(keys_b - keys_a)}")
    print()

    bins = {0: 0, 1: 0, 10: 0, 100: 0, 1000: 0, 'big': 0}
    total = 0
    unmatched = 0
    diffs = []
    bad = []

    for key in sorted(shared):
        a_list = recs_a[key]
        b_list = recs_b[key]

        b_used = set()
        for a_st, a_en in a_list:
            best_dist = float('inf')
            best_idx = -1
            for i, (b_st, b_en) in enumerate(b_list):
                if i in b_used:
                    continue
                dist = abs(a_st - b_st) + abs(a_en - b_en)
                if dist < best_dist:
                    best_dist = dist
                    best_idx = i

            if best_idx < 0:
                unmatched += 1
                continue

            b_used.add(best_idx)
            b_st, b_en = b_list[best_idx]
            d_st = abs(a_st - b_st)
            d_en = abs(a_en - b_en)
            max_d = max(d_st, d_en)
            diffs.append(max_d)
            total += 1

            if max_d == 0:
                bins[0] += 1
            elif max_d <= 1:
                bins[1] += 1
            elif max_d <= 10:
                bins[10] += 1
            elif max_d <= 100:
                bins[100] += 1
            elif max_d <= 1000:
                bins[1000] += 1
            else:
                bins['big'] += 1
                if len(bad) < 10:
                    bad.append((key, (a_st, a_en), (b_st, b_en), d_st, d_en))

    print(f"  Pairs compared: {total}  ({unmatched} unmatched)")
    if total == 0:
        return
    print(f"  Exact (0bp):     {bins[0]:6d}  ({100*bins[0]/total:.1f}%)")
    print(f"  <=1bp:           {bins[1]:6d}  ({100*bins[1]/total:.1f}%)")
    print(f"  <=10bp:          {bins[10]:6d}  ({100*bins[10]/total:.1f}%)")
    print(f"  <=100bp:         {bins[100]:6d}  ({100*bins[100]/total:.1f}%)")
    print(f"  <=1000bp:        {bins[1000]:6d}  ({100*bins[1000]/total:.1f}%)")
    print(f"  >1000bp:         {bins['big']:6d}  ({100*bins['big']/total:.1f}%)")

    if diffs:
        sd = sorted(diffs)
        print(f"\n  Max coord diff:  median={sd[len(sd)//2]}, mean={sum(diffs)/len(diffs):.1f}, max={max(diffs)}")
        for p in [90, 95, 99]:
            idx = int(len(sd) * p / 100)
            print(f"  {p}th percentile: {sd[min(idx, len(sd)-1)]}")

    if bad:
        print(f"\n  Examples of >1000bp differences:")
        for key, a, b, dst, den in bad:
            rgn, chrom, strand = key
            print(f"    {rgn} {chrom}({strand})  "
                  f"{name_a}:[{a[0]:,}-{a[1]:,}]  "
                  f"{name_b}:[{b[0]:,}-{b[1]:,}]  "
                  f"d_st={dst:,} d_en={den:,}")

paftools = load_records(sys.argv[1])
rb_all = load_records(sys.argv[2])
rb_largest = load_records(sys.argv[3])

print("--- paftools vs rb --qbed (all records) ---")
compare("paftools", paftools, "rb_qbed", rb_all)

print()
print("--- paftools vs rb --qbed --largest ---")
compare("paftools", paftools, "rb_largest", rb_largest)
PYEOF

echo ""
echo "Done."
