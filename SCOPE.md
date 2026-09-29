# SCOPE

Status: **DRAFT v0.1 — proposal for club ratification**
Last updated: <date>

## Goal

Implement a banded Smith-Waterman alignment kernel with traceback for the *M. tuberculosis* drug-resistance hotspots below, on a PYNQ-Z2 (Zynq-7020). Characterize latency, throughput, and power against a CPU and GPU baseline that uses the same scoring scheme and precision. Validate that the hardware's calls match a software reference bit-for-bit, and that the alignments recover known resistance mutations from the WHO catalogue.

This is a kernel benchmark, not a full genotyping instrument.

## Decisions

Anything marked PROVISIONAL is an assumption to be confirmed by the owner. Do not build on it until it is confirmed.

| # | Decision | Choice | Rationale | Owner | Status / revisit |
|---|----------|--------|-----------|-------|------------------|
| 1 | Loci | rpoB RRDR (Rv0667, codons 426-452, M. tuberculosis numbering); katG codon 315 region (Rv1908c); promoter region upstream of the fabG1/inhA operon (contains the commonly called inhA -15C>T) | Covers the main rifampicin and isoniazid resistance hotspots, consistent with WHO-endorsed diagnostics (Xpert MTB/RIF, MTBDRplus) | Bio lead | PROVISIONAL: confirm exact hotspot list against the WHO catalogue, week 2 |
| 2 | Reference | H37Rv, GenBank NC_000962.3 | Standard reference genome | Bio lead | Locked |
| 3 | Reference window | Each hotspot plus ~30-50 bp of flank on both sides; roughly 100-180 bp per window | Sets DP matrix size. Extract exact coordinates from the NC_000962.3 annotation | Bio lead | PROVISIONAL: coordinates TBD, week 2 |
| 4 | Strand handling | katG is on the reverse strand (Rv1908c). Decide whether to store windows in genome orientation or gene orientation, and how reads are oriented | Affects reference storage and how reads map | Algorithm lead | OPEN: decide week 3 |
| 5 | Read length | 100-150 bp assumed | Typical short-read Illumina data. Confirm read structure of the CRyPTIC data actually used | Bio lead | PROVISIONAL: confirm week 2 |
| 6 | Algorithm | Local alignment (Smith-Waterman), banded, linear gap penalty, with traceback | Linear gap keeps one state value per cell, so the PE stays small. Affine gaps are more realistic for indels and are a documented limitation and stretch goal | Algorithm lead | Locked family; parameters below |
| 7 | Scoring values | Match +2, mismatch -1, gap -2 | Starting point only. Quantization study (weeks 3-4) decides final values and score bit width | Algorithm lead | PROVISIONAL: revisit week 4 |
| 8 | Score width | Max score is 2 x read length (about 300 for 150 bp), so at least 9 bits unsigned | Precision study will test whether narrower widths still give correct calls | Algorithm lead | PROVISIONAL: revisit week 4 |
| 9 | Traceback storage | 2-bit pointer per cell (diag / up / left / stop) | 150 x 180 cells x 2 bits is about 54 kbit, which fits in Zynq-7020 block RAM. Compare against a recompute approach in the write-up | RTL lead | PROVISIONAL: revisit week 7 |
| 10 | Software reference | Plain Python (`reference/`), later NumPy. Hardware output is diffed bit-exactly against it | Correctness baseline for everything else | Verification lead | Locked |
| 11 | RTL language | Amaranth HDL (Python), generating Verilog for Vivado | Same language as the reference model, so the test flow is easy to share. Risk: an extra translation layer when debugging timing | RTL lead | OPEN: confirm by end of week 3. Do not mix Amaranth and hand-written Verilog |
| 12 | Target board | PYNQ-Z2 (XC7Z020-1CLG400C), PYNQ overlay for host-side streaming | Python/Jupyter host flow reduces harness engineering | Bring-up lead | Locked |
| 13 | Baselines | CPU and GPU implementing the same algorithm, scoring, and precision | Fair comparison is a central claim of the paper | Benchmarking lead | Locked |
| 14 | Power measurement | INA219 on the supply rail plus the board's own power monitoring where available. Document exactly how it is measured | Reviewers check power-measurement methods | Benchmarking lead | PROVISIONAL: validate weeks 17-19 |

## Locked vs. parameterized

**Locked:** the three loci, H37Rv reference, local alignment with traceback and linear gaps, the target board, the bit-exact verification approach, and the fair-baseline requirement.

**Parameters (config values, not frozen numbers):** read length, window size, band width, scoring values, score bit width, clock frequency, array size.

## Out of scope

- A full genotyping instrument or complete pipeline (read filtering, variant calling, reporting).
- Whole-genome alignment.
- Affine gap penalties in hardware (documented as a limitation and stretch goal).
- Clinical claims beyond sensitivity/specificity against the WHO catalogue.

## Tasks

| Role |
|------|
| Algorithm / quantization lead |
| RTL / kernel lead |
| Memory / DMA lead |
| Verification lead |
| Benchmarking lead |
| Bring-up / integration lead |
| Bio lead (references, catalogue, test reads, validation) |


## Open questions

- Exact hotspot list and coordinates (see decisions 1 and 3).
- Which CRyPTIC subset to use, and its read structure.
- Whether affine gaps matter for the mutations in scope (many resistance-relevant variants are SNPs, but indels near hotspot boundaries are a required test case).
- Possibility of leaning towards nanopore sequencing; requires restructuring of large pillars of the project
- Exact loci used in the project (may expand from 3 to entire WHO drug-resistance gene panel)
- Pursue the raw-signal angle (adaptive sampling / sDTW) instead of basecalled-sequence alignment? this becomes a different algorithm entirely from Smith Waterman

## Change log

| Date | Change | Reason |
|------|--------|--------|
| <date> | v0.1 draft | Initial proposal |
