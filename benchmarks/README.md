# rustybam Benchmarks

Reproducible benchmarks for the rustybam & SafFire application note.

## Prerequisites

- `rb` binary: `cargo build --release` in `../../rustybam/`
- `hyperfine`: `brew install hyperfine`
- `paftools.js`: `conda install -c bioconda minimap2`
- `samply` (optional, for flamegraph profiling): `cargo install samply`
- Input PAF in `input-pafs/` (compressed with bgzip)

## Running

```bash
# Run everything (benchmarks + profiles + summary)
snakemake -j4 -s Snakefile --forceall

# Run subsets
snakemake -j4 -s Snakefile bench      # just hyperfine benchmarks
snakemake -j4 -s Snakefile profile    # just macOS sample + samply profiles
snakemake -j4 -s Snakefile samply     # just samply flamegraph profiles
snakemake -j4 -s Snakefile summary    # just collect results into summary table
```

## Outputs

- `results/summary.tsv`: Aggregated benchmark results
- `results/summary.md`: Markdown summary with two tables (standalone + paftools comparison)

All intermediate files (derived PAFs, bench JSONs, profile outputs) are marked `temp()` and cleaned up automatically.

## What is benchmarked

| Benchmark | Description |
|-----------|-------------|
| `trim_paf` | Resolve overlapping alignments (1x PAF, 1,460 records) |
| `orient` | Orient alignments to forward strand (10x PAF) |
| `filter` | Filter by paired alignment length (10x PAF) |
| `invert` | Invert alignment coordinates (10x PAF) |
| `stats_paf` | Compute alignment statistics (10x PAF) |
| `break_paf_m{1000,5000}` | Split at indels of various sizes (10x PAF) |
| `stats_comparison` | rb stats vs paftools.js stat (1x and 10x) |
| `liftover_comparison` | rb liftover vs paftools.js liftover |
