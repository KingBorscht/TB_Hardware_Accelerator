# SCOPE v2

Status: DRAFT v2.0 (supersedes v0.1)
Last updated: 2026-10-10

## Goal

Implement a banded, affine-gap Smith-Waterman alignment kernel with traceback for long, noisy nanopore amplicon reads (1-3 kb) covering *M. tuberculosis* drug-resistance regions, on a PYNQ-Z2 (Zynq-7020). Characterize latency, throughput, energy per read and memory use against CPU, GPU and same-board ARM baselines that run the same algorithm, scoring and precision. Validate that the hardware matches the software specification model bit-for-bit, and that the pipeline recovers known WHO-catalogue resistance mutations in samples whose true content is independently known.

This is a kernel benchmark and memory-architecture study, not a full genotyping instrument. Basecalling happens upstream of the accelerator.

## Pipeline

Sample DNA -> PCR amplicons -> nanopore sequencing and Dorado basecalling -> basecalled reads (FASTQ) into the FPGA -> primer match (amplicon + strand) -> banded affine-gap alignment with traceback -> per-read differences at catalogued positions -> consensus mutation call on the host ARM core.

## Decisions

| # | Decision | Choice | Focus | Status |
|---|----------|--------|-------|--------|
| 1 | Panel | 6-amplicon wet-lab core: rpoB RRDR, katG codon 315 region, inhA/fabG1 promoter, pncA, embB, gyrA. Hardware parameter: up to 32 amplicons | Bio | PROVISIONAL: pncA, embB, gyrA subject to change |
| 2 | Organism / sample source | No live *M. tuberculosis*. Non-infectious material only (synthetic controls, surrogate organism, or purified DNA)\* | Bio | Possible surrogate organism: *M. smegmatis* |
| 3 | Reference | H37Rv, NC_000962.3 (second genome added if a surrogate organism is used) | Bio | Locked |
| 4 | Reference windows | One window per amplicon (up to 2.5 kb), gene orientation, so positions match WHO codon numbering | Bio / Algorithm | PROVISIONAL: coordinates TBD with primers |
| 5 | Sequencing input | Nanopore, Dorado HAC basecalling. FPGA input is basecalled FASTQ, not raw signal. Assumed flow cell R10.x | Bio | Locked |
| 6 | Read length | Target 1-3 kb. L_max = 3,072 bp. 5 kb is a stretch needing checkpoint-and-recompute | Algorithm / RTL | Locked |
| 7 | Read identification | Match first 25 bp to primers (up to 3 mismatches) to get amplicon and strand; flip reverse reads | Algorithm / RTL | Locked |
| 8 | Algorithm | Banded local alignment (Smith-Waterman), single-piece affine gap, with traceback | Algorithm | Locked |
| 9 | Scoring | match +2, mismatch -4, gap open 4, extend 2 (minimap2 map-ont preset) | Algorithm | Locked |
| 10 | Score width | 13 bits unsigned (max score 2 x 3,072 = 6,144 < 8,192) | Algorithm | Locked for L_max 3,072 |
| 11 | Band width | 128 (+/-64) | Algorithm | PROVISIONAL |
| 12 | Traceback storage | Baseline: store all pointers in BRAM, 4 bits/cell (~1.57 Mbit, ~30% of BRAM). Alternative: checkpoint-and-recompute\* | RTL | Switch triggers at first synthesis |
| 13 | Array and reads in flight | One 128-PE array. Target: double-buffered (one read in forward pass, one in traceback) | RTL / Memory-DMA | PROVISIONAL |
| 14 | Clock and latency | 100 MHz target (150 MHz stretch). Hardware-only latency target under 100 us per 3 kb read | Bring-up / RTL | PROVISIONAL: unsynthesized estimate |
| 15 | Mutation calling | Host ARM core, Python majority vote. Coverage at least 50 reads per amplicon (100 target). Resistant at allele fraction 70% or more; mixed 10-70%; ignore below 10% | Bio / Algorithm | Locked |
| 16 | RTL language | Amaranth HDL generating Verilog. | RTL | Locked |
| 17 | Target board | PYNQ-Z2 (XC7Z020-1CLG400C) with PYNQ overlay | Bring-up | Locked |
| 18 | Baselines | Tiered: Python (correctness only), optimized CPU, same-board ARM, GPU, minimap2 as labeled context\* | Benchmarking | Locked structure |
| 19 | Benchmark reporting | Latency and throughput separate; hardware-only and end-to-end separate; GPU launch overhead counted; energy per read | Benchmarking | Locked |
| 20 | Power measurement | INA219 on supply rail plus board power monitor; document method including idle subtraction | Benchmarking | Locked |
| 21 | Ground truth | Independently known sample content | Bio / Verification | PROVISIONAL |

## Locked vs. parameterized

**Locked:** banded affine-gap local alignment with traceback; H37Rv as the TB reference; the target board; bit-exact verification against the specification model; fair baselines; no learned or predictive components; no live *M. tuberculosis* handled by the club; basecalled reads as FPGA input.

**Parameters:** panel size (up to 32), L_max (3,072), band width (128), scoring values, score width (13), clock, PE count (128), buffered reads, primer-match window, coverage and calling thresholds, traceback strategy.

## Resource estimates (arithmetic, not measured)

- Traceback: 3,072 x 128 x 4 bits = 1.57 Mbit of ~5 Mbit BRAM (630 KB): ~30% single-buffered, ~60% double-buffered.
- Reference storage: 32 x 2.5 kb x 2 bits = 160 kbit (~3%).
- Checkpoint-and-recompute: roughly 5% of BRAM per read in flight (rough estimate).
- Array logic: a 128-PE affine array is guessed at 20-32k of 53,200 LUTs. **This is the number most likely to force a change**; confirm at first synthesis.
- Throughput: the array should far outpace a single MinION. The host (DMA, Python) will likely limit end-to-end rate, so report hardware-only and end-to-end separately.

## Out of scope

- A full genotyping instrument or pipeline (read filtering, variant reporting).
- Basecalling on the FPGA; raw-signal methods (subsequence DTW, read-until).
- Whole-genome sequencing, seed-and-chain mapping, whole-genome alignment.
- Hotspot-local traceback and other approximate traceback schemes.
- Two-piece affine gaps in hardware.
- Handling live *M. tuberculosis*.
- Clinical claims beyond sensitivity and specificity against the WHO catalogue.

## Tasks

| Focus | Details |
|-------|---------|
| Algorithm / quantization | Affine banded design, precision study on long noisy reads, optimization evaluation, specification model and full SW reference |
| RTL / kernel | 128-PE array, affine cell, traceback walker, primer-match logic, testbenches |
| Memory / DMA | AXI DMA, PYNQ block design, cache coherency, read double-buffering |
| Verification | Simulation, bit-exact validation, coverage including indels near amplicon and primer boundaries |
| Benchmarking | Tiered baselines, fair-comparison method, power measurement, hardware-only vs end-to-end reporting |
| Bring-up / integration | Timing closure, constraints, bitstream, debugging, first-synthesis resource check |
| Bio | Primers and PCR, controls, sequencing and Dorado, ground truth, calling thresholds, WHO catalogue |

## Open questions

1. Organism and sample source with the supervising adult and school. *M. smegmatis* mutants selected on antibiotic plates would give real mutations, but they are *M. smegmatis* mutations, not new TB mutations.
2. Whether novel-mutation discovery is still a goal, and whether surrogate-organism mutants satisfy it.
3. Which deterministic optimizations are adopted, with accuracy thresholds defined before implementation.
4. Whether ksw2, Parasail or SSW support the exact scoring, banding and traceback of the spec.
5. Whether two-piece affine gaps are worth a stretch goal.
