"""Sliding-window model: one read base in, one out, a short memory in between.

The window holds W cells, one per diagonal around the diagonal picked by the
vote. Each step shifts in the next read base and reuses the same W cells, so
the hardware needs W cells instead of the whole band laid out in space.

Three questions, on the same simulated reads as band_study.py:
1. For each memory length W, how often does the result equal full SW?
2. How many cells do anything other than the default update in a step?
   Default = "one more match on my diagonal" (score + 2) or "still zero".
   In the spike picture a default cell fires exactly on schedule, so only the
   other cells carry new information.
3. How many bits does a cell need if it stores differences?

Run:  python reference/window_study.py          (from the project root)
"""

import argparse
import math
import random

import numpy as np

from band_study import (GAP, INDEL_LENGTHS, LOCI, MATCH, MISMATCH, READ_LEN, ROOT,
                        diagonal_vote, encode, read_fasta, simulate_read, sw_score)

WINDOWS = [4, 5, 9, 17]     # memory length in cells (bases); 4 cells = 8 bits of sequence


def window_offsets(W):
    """Diagonal offsets covered by a W-cell window, relative to the vote."""
    lo = -((W - 1) // 2)
    return lo, lo + W - 1


def sliding_window(read, ref, d0, W):
    """Banded SW computed as a W-cell window sliding one read base per step.

    Returns (H, mask): score matrix with zeros outside the window, and the
    boolean mask of cells that were inside it.
    """
    m, n = len(read), len(ref)
    lo_off, hi_off = window_offsets(W)
    ramp = GAP * np.arange(n + 1)
    H = np.zeros((m + 1, n + 1), dtype=np.int64)
    mask = np.zeros((m + 1, n + 1), dtype=bool)
    for i in range(m):
        a, b = max(i + d0 + lo_off, 0), min(i + d0 + hi_off, n - 1)
        if a > b:
            continue
        sub = np.where(ref == read[i], MATCH, MISMATCH)
        e = np.zeros(n + 1, dtype=np.int64)
        e[1:] = np.maximum(0, np.maximum(H[i, :-1] + sub, H[i, 1:] + GAP))
        e[: a + 1] = 0
        e[b + 2:] = 0
        cur = np.maximum.accumulate(e - ramp) + ramp
        cur[: a + 1] = 0
        cur[b + 2:] = 0
        H[i + 1] = cur
        mask[i + 1, a + 1:b + 2] = True
    return H, mask


def activity(H, mask):
    """Cells that deviate from the default update, and steps with none at all.

    Returns (active cells, cells in window, quiet steps, steps with a window).
    """
    prev = H[:-1, :-1]                       # same diagonal, previous step
    cur = H[1:, 1:]
    inside = mask[1:, 1:]
    default = (cur == prev + MATCH) | ((cur == 0) & (prev == 0))
    active = inside & ~default
    steps = inside.any(axis=1)
    quiet = steps & ~active.any(axis=1)
    return int(active.sum()), int(inside.sum()), int(quiet.sum()), int(steps.sum())


def diff_range(H, mask):
    """Range of differences between neighbouring cells that are both in the window."""
    side = (H[:, 1:] - H[:, :-1])[mask[:, 1:] & mask[:, :-1]]      # adjacent cells, same step
    time = (H[1:, 1:] - H[:-1, :-1])[mask[1:, 1:] & mask[:-1, :-1]]  # same cell, next step
    return side, time


def bits_for(lo, hi):
    return max(1, math.ceil(math.log2(hi - lo + 1)))


def run(n_per_cell, seed):
    rng = random.Random(seed)
    regions = read_fasta(ROOT / "data" / "regions.fasta")
    conditions = [("none", 0)] + [(k, L) for k in ("del", "ins") for L in INDEL_LENGTHS]
    rows = []
    ranges = {W: {"side": [0, 0], "time": [0, 0]} for W in WINDOWS}

    for locus, (centre, half) in LOCI.items():
        region = regions[locus]
        ref = encode(region[centre - half:centre + half])
        for kind, length in conditions:
            for _ in range(n_per_cell):
                read_s, _ = simulate_read(rng, region, centre, half, kind, length)
                read = encode(read_s)
                full, _ = sw_score(read, ref)
                d0, _ = diagonal_vote(read, ref)
                row = {"indel_len": length}
                for W in WINDOWS:
                    H, mask = sliding_window(read, ref, d0, W)
                    act, cells, quiet, steps = activity(H, mask)
                    row[f"exact{W}"] = int(H.max()) == full
                    row[f"active{W}"] = act
                    row[f"cells{W}"] = cells
                    row[f"quiet{W}"] = quiet
                    row[f"steps{W}"] = steps
                    side, time = diff_range(H, mask)
                    for name, d in (("side", side), ("time", time)):
                        if d.size:
                            ranges[W][name][0] = min(ranges[W][name][0], int(d.min()))
                            ranges[W][name][1] = max(ranges[W][name][1], int(d.max()))
                rows.append(row)
    return rows, ranges


def report(rows, ranges):
    lengths = [0] + INDEL_LENGTHS
    print("1. Reads whose score equals full SW, by memory length W (cells)")
    print("   " + "  ".join(f"{h:>7s}" for h in ["indel"] + [f"W={W}" for W in WINDOWS]))
    for L in lengths:
        sub = [r for r in rows if r["indel_len"] == L]
        print("   " + "  ".join([f"{L:>7d}"] + [f"{np.mean([r[f'exact{W}'] for r in sub]):>7.3f}" for W in WINDOWS]))

    print("\n2. Activity (all reads)")
    print("   " + "  ".join(f"{h:>22s}" for h in ["", "active cells / cells", "fully quiet steps", "active cells per read"]))
    for W in WINDOWS:
        act = sum(r[f"active{W}"] for r in rows)
        cells = sum(r[f"cells{W}"] for r in rows)
        quiet = sum(r[f"quiet{W}"] for r in rows)
        steps = sum(r[f"steps{W}"] for r in rows)
        print("   " + "  ".join([f"{'W=' + str(W):>22s}", f"{act / cells:>22.3f}",
                                 f"{quiet / steps:>22.3f}", f"{act / len(rows):>22.1f}"]))
    print("   by indel length, W=17: active cells / cells")
    for L in lengths:
        sub = [r for r in rows if r["indel_len"] == L]
        print(f"   {L:>7d}  {sum(r['active17'] for r in sub) / sum(r['cells17'] for r in sub):.3f}")

    print("\n3. Bits per cell if differences are stored")
    for W in WINDOWS:
        s, t = ranges[W]["side"], ranges[W]["time"]
        print(f"   W={W:<2d}  neighbour cell: {s[0]:+d}..{s[1]:+d} -> {bits_for(*s)} bits   "
              f"same cell, next step: {t[0]:+d}..{t[1]:+d} -> {bits_for(*t)} bits")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=60, help="reads per locus per condition")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    report(*run(args.n, args.seed))
