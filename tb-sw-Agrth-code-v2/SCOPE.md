# SCOPE

Status: **DRAFT v0.2 — proposal for club ratification**
Last updated: 2026-10-05

Rows marked **PROPOSED** are new in v0.2 and come from the Algorithm focus. They are not agreed yet.

## Goal

Implement a banded Smith-Waterman alignment kernel with traceback for the *M. tuberculosis* drug-resistance hotspots below, on a PYNQ-Z2 (Zynq-7020). Characterize latency, throughput, and power against CPU and GPU baselines that use the same scoring scheme and precision. Validate that the hardware's calls match a software reference bit-for-bit, and that the alignments recover known resistance mutations from the WHO catalogue.

The project's claim is **cost-effectiveness, not accuracy**. Accuracy is a fixed requirement (results equal to the reference). The numbers we compete on are:

- reads per second per dollar of hardware
- energy per read
- latency per read

This is a kernel benchmark, not a full genotyping instrument.

## Architecture (PROPOSED, working name "Silicon Brain")

```
read --> 1. diagonal vote --> 2. gate --+-- pass --> output (SNPs only, no SW)
                                        |
                                        +-- fail --> 3. sliding-window SW + traceback --> output
```

1. **Diagonal vote.** Compare the read with the reference window base by base (2-bit XNOR), count matches on each diagonal, keep the best one.
2. **Gate.** If the ungapped score on that diagonal is high enough, report the mismatch positions and skip SW.
3. **Sliding-window SW.** A window of W cells centred on the voted diagonal. Each step shifts one read base in and reuses the same W cells. This is banded SW with band width W.

The design borrows four ideas from neural hardware: sparse activity, event-driven updates, low precision, and local connections. It does not use trained networks, asynchronous logic, or learning. In writing, call it "neural-inspired".

## Decisions

Anything marked PROVISIONAL is an assumption to be confirmed by the owner. Do not build on it until it is confirmed.

| # | Decision | Choice | Rationale | Focus | Status / revisit |
|---|----------|--------|-----------|-------|------------------|
| 1 | Loci | rpoB RRDR (Rv0667, codons 426-452, M. tuberculosis numbering); katG codon 315 region (Rv1908c); promoter region upstream of the fabG1/inhA operon (contains the commonly called inhA -15C>T) | Covers the main rifampicin and isoniazid resistance hotspots, consistent with WHO-endorsed diagnostics (Xpert MTB/RIF, MTBDRplus) | Bio | PROVISIONAL: confirm exact hotspot list against the WHO catalogue, week 2 |
| 2 | Reference | H37Rv, GenBank NC_000962.3 | Standard reference genome | Bio | Locked |
| 3 | Reference window | Each hotspot plus ~30-50 bp of flank on both sides; roughly 100-180 bp per window | Sets DP matrix size. Extract exact coordinates from the NC_000962.3 annotation | Bio | PROVISIONAL: coordinates TBD, week 2 |
| 4 | Strand handling | katG is on the reverse strand (Rv1908c). Decide whether to store windows in genome orientation or gene orientation, and how reads are oriented | Affects reference storage and how reads map. The diagonal vote would have to run for both orientations | Algorithm | OPEN: decide week 3 |
| 5 | Read length | 100-150 bp assumed | Typical short-read Illumina data. Confirm read structure of the CRyPTIC data actually used | Bio | PROVISIONAL: confirm week 2 |
| 6 | Algorithm | Local alignment (Smith-Waterman), banded, linear gap penalty, with traceback. The band is a sliding window of W cells (decision 16) placed by a diagonal vote (decision 15) | Linear gap keeps one state value per cell, so the PE stays small. Affine gaps are more realistic for indels and are a documented limitation and stretch goal | Algorithm | Locked family; window form PROPOSED |
| 7 | Scoring values | Match +2, mismatch -1, gap -2 | Starting point only. Not yet studied | Algorithm | PROVISIONAL: revisit week 4 |
| 8 | Score width | 9 bits unsigned for absolute scores (max 300). If cells store the difference to a neighbouring cell instead, 3 bits are enough (measured range -2 to +4). First RTL uses absolute scores | Difference storage shrinks each cell, so more parallel lanes fit. Range measured on simulated reads only | Algorithm | PROVISIONAL: confirm on real reads; RTL to decide when to switch |
| 9 | Traceback storage | 2-bit pointer per cell (diag / up / left / stop), stored only for cells inside the window: read length x W x 2 bits. For W = 17 that is 150 x 17 x 2 = 5,100 bits | v0.1 sized this for the full matrix (54 kbit). Banding cuts it about tenfold | RTL | PROVISIONAL: revisit week 7 |
| 10 | Software reference | Plain Python (`reference/`), later NumPy. Hardware output is diffed bit-exactly against it. The reference must model the same architecture as the hardware (vote, gate, window), not full SW | Correctness baseline for everything else | Verification | Locked; architecture-matched model PROPOSED |
| 11 | RTL language | Amaranth HDL (Python), generating Verilog for Vivado | Same language as the reference model, so the test flow is easy to share. Risk: an extra translation layer when debugging timing | RTL | OPEN: confirm by end of week 3. Do not mix Amaranth and hand-written Verilog |
| 12 | Target board | PYNQ-Z2 (XC7Z020-1CLG400C), PYNQ overlay for host-side streaming | Python/Jupyter host flow reduces harness engineering | Bring-up | Locked |
| 13 | Baselines | CPU and GPU implementing the same algorithm, scoring, and precision. Add one low-cost CPU board (for example Raspberry Pi 5) | Fair comparison is a central claim of the paper. A cost-effectiveness claim needs a cheap baseline, not only a laptop | Benchmarking | Locked; low-cost board PROPOSED |
| 14 | Power measurement | INA219 on the supply rail plus the board's own power monitoring where available. Document exactly how it is measured | Reviewers check power-measurement methods. Board-level power includes the ARM cores and static power, which may hide savings in the fabric | Benchmarking | PROVISIONAL: validate weeks 17-19 |
| 15 | Band placement | Diagonal vote: the diagonal with the most matching bases becomes the window centre | With the vote, a fixed band matched full SW whenever the half-width was at least the indel length (about 0.3% exceptions) | Algorithm | PROPOSED |
| 16 | Window length W | Parameter. Measured: W = 5 handles indels up to 2 bp, W = 9 up to 3 bp, W = 17 up to 8 bp. A 4-cell window handles only 1 bp | W sets the number of cells and the longest indel supported | Algorithm, Bio | PROPOSED: Bio to give the longest indel that must be supported |
| 17 | Gate | Reads whose ungapped score passes a threshold skip SW | 99.3% of SNP-only reads already had the full-SW score without SW. Risk: 3.6% of indel reads looked gap-free because the indel sat within 14 bp of the end of the overlap | Algorithm | PROPOSED: threshold and end-of-read rule OPEN |
| 18 | Event gating | A cell updates only when it deviates from the default (one more match, or still zero) | In simulation 63-78% of steps had no such cell. Whether this saves measurable power is unknown | RTL, Benchmarking | PROPOSED stretch goal |
| 19 | Correctness criteria | (a) Hardware equals the architecture-matched Python reference bit-for-bit. (b) Report how often that reference agrees with full SW, for score and for mutation call | The gate and the band are heuristics, so (a) alone does not show the calls are right | Verification, Algorithm | PROPOSED |
| 20 | Second, low-cost target | Port the same RTL to a small FPGA board as a stretch goal. Candidates: Tang Nano 9K (GW1NR-9, 8,640 LUT4, about $16-39), or an iCE40UP5K board (5,280 LUT4, 75 uW typical static power, chip about $11) | A sizing sketch needs about 3,500 LUT4 for W = 17 plus the vote (pessimistic mapper, no host interface, no place-and-route). PYNQ-Z2 is $129 through the AMD University Program. Prices as of October 2026 | Bring-up, RTL | PROPOSED stretch goal; host link (UART / SPI) undecided |

## Evidence behind the v0.2 proposals

All of it comes from simulated 150 bp reads drawn from the three H37Rv regions (0-2 SNPs, 0.5% substitution errors, one insertion or deletion of 1-12 bp, forward strand only). Nothing here has been run on CRyPTIC data or on hardware.

| Finding | Result | Script |
|---------|--------|--------|
| Fixed band after the diagonal vote | Equal to full SW when half-width >= indel length; 5-18% of the matrix computed | `reference/band_study.py` |
| Tier 1 alone | Equal to full SW for 99.3% of SNP-only reads | `reference/band_study.py` |
| Band that steers toward the row maximum | Worse than a fixed band of the same width | `reference/band_study.py` |
| Difference storage | 3 bits (neighbouring cell) or 2 bits (same cell, next step) | `reference/ideas_study.py`, `reference/window_study.py` |
| Per-segment vote ("attention" band) | One third of the cells of a half-width-16 band, but misses 7-13% of 5-12 bp indels near the end of the overlap | `reference/ideas_study.py` |
| Sliding window | Results by W as in decision 16; 9-12% of cells deviate from the default update | `reference/window_study.py` |
| Area sketch | Window core W = 17: 2,731 LUT4, 234 FF. Diagonal vote: 781 LUT4, 776 FF. Sketch matched the Python model on 30 simulated reads | `estimate/area_estimate.py`, `estimate/check_window_core.py` |

## Locked vs. parameterized

**Locked:** the three loci, H37Rv reference, local alignment with traceback and linear gaps, the target board, the bit-exact verification approach, and the fair-baseline requirement.

**Parameters (config values, not frozen numbers):** read length, window size, window length W, scoring values, score bit width, gate threshold, clock frequency, number of parallel lanes.

## Out of scope

- A full genotyping instrument or complete pipeline (read filtering, variant calling, reporting).
- Whole-genome alignment.
- Affine gap penalties in hardware (documented as a limitation and stretch goal).
- Clinical claims beyond sensitivity/specificity against the WHO catalogue.
- Trained neural networks (CNN, SNN, or attention models) inside the kernel. Results would stop being bit-exact.
- Asynchronous or analog circuits. An FPGA cannot do them.
- A band that moves during alignment. Tested and worse than a fixed band.
- Replacing SW with an event-driven (wavefront) algorithm. A rough estimate suggests about ten times fewer steps, but it changes the locked algorithm; record it as future work.
- A custom PCB or a chip tape-out. Possible follow-ups after the main project.

## Tasks

| Focus | Details |
|-------|---------|
| Algorithm / quantization | Architecture-matched Python reference (vote, gate, window, traceback) with test-vector generation; gate threshold and end-of-read rule; scoring study; strand handling; rerun all studies on real reads |
| RTL / kernel | Diagonal-vote block; window of W cells with traceback pointers; best-score tracker; later: difference storage, event gating, parallel lanes |
| Memory / DMA | AXI DMA architecture, PYNQ block design, cache-coherency handling. Reads are short, so batch many per transfer; define the input and output record formats with Algorithm |
| Verification | Simulation infra; bit-exact validation against the architecture-matched reference; coverage of indels near read ends, both strands, reads that fail the gate |
| Benchmarking | GPU/CPU baselines with the same two-stage algorithm; low-cost CPU board baseline; power-measurement setup; report reads/s per dollar and energy per read |
| Bring-up / integration | Timing closure, constraints, bitstream generation, hardware debugging. Stretch: second low-cost board |
| Bio | References, catalogue, test reads, validation. Needed soon: longest indel to support (sets W), window coordinates, CRyPTIC read structure |

## Open questions

- Exact hotspot list and coordinates (see decisions 1 and 3).
- Which CRyPTIC subset to use, and its read structure.
- Longest indel that must be supported at these loci (sets W, decision 16).
- Gate threshold, and how to catch indels near the end of a read (decision 17).
- Whether affine gaps matter for the mutations in scope (many resistance-relevant variants are SNPs, but indels near hotspot boundaries are a required test case).
- Does event gating reduce power that can be measured at the board supply (decisions 14 and 18)?
- How reads reach a small FPGA board without an on-board processor (decision 20).
- Possibility of leaning towards nanopore sequencing; requires restructuring of large pillars of the project.
- Exact loci used in the project (may expand from 3 to entire WHO drug-resistance gene panel).
- Pursue the raw-signal angle (adaptive sampling / sDTW) instead of basecalled-sequence alignment? This becomes a different algorithm entirely from Smith-Waterman.

## Change log

| Date | Change | Reason |
|------|--------|--------|
| <date> | v0.1 draft | Initial proposal |
| 2026-10-05 | v0.2 draft: goal restated around cost-effectiveness; architecture section; decisions 6, 8, 9, 10, 13 revised; decisions 15-20 proposed; tasks, out of scope, and open questions updated | Algorithm focus studies on simulated reads |
