# TB Hardware Accelerator

A banded Smith-Waterman alignment kernel for *Mycobacterium tuberculosis* drug-resistance genotyping, implemented on a PYNQ-Z2 (Zynq-7020) FPGA and benchmarked against CPU/GPU baselines for latency, throughput, and power.

This is a **kernel benchmark**, not a full diagnostic instrument. See [`SCOPE.md`](./SCOPE.md) for the exact technical scope, locked decisions, and open questions.

## What this project does

We align sequencing reads covering the *rpoB*, *katG*, and *inhA* resistance loci against the H37Rv reference genome, in custom FPGA hardware, and compare it against equivalent CPU/GPU implementations. The goal is to show, with real measured data, how much power and latency you can save by moving this workload onto custom silicon instead of a CPU or GPU, and to map out exactly where those gains hold up.

Background precedent: [Gliddon et al. 2021](https://doi.org/10.1128/Spectrum.00610-21), a targeted isothermal amplification + nanopore sequencing workflow for the same three loci.

## Repo layout

```
reference/   Software (Python) Smith-Waterman reference implementation
rtl/         Amaranth HDL / hardware kernel source
bench/       CPU/GPU baselines and benchmarking scripts
bio/         Reference sequences, WHO mutation catalogue extracts, test data notes
data/        Small test vectors only (large/raw datasets are gitignored — see below)
docs/        Write-up drafts, figures, notes
SCOPE.md     Living document of technical decisions, rationale, and owners
```

## Getting started

**Requirements**
- [Git](https://git-scm.com/)
- [Python 3.10+](https://www.python.org/)
- [Amaranth HDL](https://amaranth-lang.org/) (`pip install amaranth`)
- [AMD Vivado](https://www.xilinx.com/support/download.html) (Standard edition is free) — only needed for synthesis/bitstream generation
- [VS Code](https://code.visualstudio.com/) (recommended, not required)

**Clone the repo**
```bash
git clone https://github.com/<your-username>/TB_Hardware_Accelerator.git
cd TB_Hardware_Accelerator
```

**Workflow**
1. Read [`SCOPE.md`](./SCOPE.md) before starting any work — it defines what's locked vs. still open.
2. Create a branch for your work: `git checkout -b your-feature-name`
3. Commit small, working changes with clear messages.
4. Open a pull request into `main` for review before merging.

## Core focus areas

- Algorithm / quantization
- RTL / kernel design
- Memory / DMA
- Verification
- Benchmarking
- Bring-up / integration
- Biology / reference data

## Status

Active development — see [`SCOPE.md`](./SCOPE.md) change log for the latest decisions.

## License

Not yet decided. All rights reserved until specified otherwise.
