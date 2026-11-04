# C++ host mutation caller

`mutation_caller.cpp` is the CPU-side prototype for the SCOPE v2 pipeline. It
does not run on the FPGA. It consumes a normalized basecalled read plus the
completed FPGA traceback result, extracts per-read evidence for a versioned
catalogue, and applies the locked coverage/allele-fraction thresholds.

Build with CMake:

```bash
cmake -S host -B build/host
cmake --build build/host
```

Or compile directly:

```bash
c++ -std=c++17 -O2 -Wall -Wextra -pedantic host/mutation_caller.cpp -o tb_mutation_caller
```

Run:

```bash
./tb_mutation_caller references.tsv catalogue.tsv alignments.tsv calls.tsv
```

All input files are tab-separated, with comment lines beginning with `#`.
There are no third-party dependencies.

## `references.tsv`

```text
# amplicon_id  name  gene_oriented_reference  genome_start_1based  genome_step  reference_version
0\trpoB_RRDR\tACGT...\t761155\t1\tNC_000962.3
```

`genome_step` is `1` or `-1`. Use `-1` for a gene-oriented window that runs
opposite to the H37Rv genome coordinate direction.

## `catalogue.tsv`

```text
# catalogue_id  amplicon_id  reference_offset  reference  alternate  gene  annotation  catalogue_version
WHO-EXAMPLE\t0\t24\tC\tT\trpoB\tp.Example\t<pinned-WHO-release>
```

Offsets are zero-based within the gene-oriented amplicon reference. Use `-`
for the inserted/deleted side of an indel. The exact WHO release, amplicon
coordinates, and indel-normalization rules still need to be supplied by the
Bio workstream before this is used on real data.

## `alignments.tsv`

```text
# read_id  amplicon_id  aligned  best_score  query_begin  query_end  reference_begin  reference_end  cigar  normalized_read_sequence
read-001\t0\t1\t4102\t0\t1500\t0\t1500\t1500M\tACGT...
```

The board/DMA adapter must produce a **complete** `M`/`I`/`D` CIGAR and
zero-based half-open query/reference endpoints. A score alone cannot identify
a mutation. Reads must already be reverse-complemented to the selected
gene-oriented reference; SCOPE v2 assigns primer matching and strand handling
to the FPGA path.

## Output

One output row is emitted per catalogue variant. It includes reference/alternate
counts, eligible coverage (`REF + ALT`), allele fraction, H37Rv position, and
one of:

- `INSUFFICIENT_COVERAGE`: fewer than 50 eligible reads;
- `RESISTANT_SUPPORTED`: alternate fraction >= 0.70;
- `MIXED`: alternate fraction from 0.10 through <0.70;
- `NOT_SUPPORTED`: alternate fraction <0.10.

This is a benchmark research component. It is not clinical diagnostic software.
