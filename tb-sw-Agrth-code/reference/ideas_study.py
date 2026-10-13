"""Three follow-up experiments on the same simulated reads as band_study.py.

1. Difference encoding (quantization): how many bits does a cell need if it
   stores the difference to its neighbour instead of the absolute score?
2. Attention-guided band: vote the diagonal separately for each segment of the
   read, and place a narrow band that follows the segment diagonals.
3. Event-driven work estimate: how many steps does a wavefront (O(ND)) edit
   distance alignment take, compared with the cells of a banded SW?

Run:  python reference/ideas_study.py           (from the project root)
"""

import argparse
import math
import random

import numpy as np

from band_study import (GAP, INDEL_LENGTHS, LOCI, MATCH, MISMATCH, READ_LEN, ROOT,
                        diagonal_vote, encode, read_fasta, simulate_read, sw_score)
from sw_reference import sw_align

SEGMENT = 25          # read bases per attention segment
RADIUS = 16           # segment votes look this far either side of the global diagonal
MIN_VOTE_FRAC = 0.7   # a segment vote counts only if this fraction of bases match


# ------------------------------------------------- 1. difference encoding

def sw_matrix(read, ref):
    """Full score matrix, (m+1) x (n+1), same recurrence as sw_score."""
    m, n = len(read), len(ref)
    ramp = GAP * np.arange(n + 1)
    H = np.zeros((m + 1, n + 1), dtype=np.int64)
    for i in range(m):
        sub = np.where(ref == read[i], MATCH, MISMATCH)
        e = np.zeros(n + 1, dtype=np.int64)
        e[1:] = np.maximum(0, np.maximum(H[i, :-1] + sub, H[i, 1:] + GAP))
        H[i + 1] = np.maximum.accumulate(e - ramp) + ramp
    return H


def bits_for(lo, hi):
    return max(1, math.ceil(math.log2(hi - lo + 1)))


# ------------------------------------------------- 2. attention-guided band

def segment_band(read, ref, d0, w):
    """Per-row band limits that follow per-segment diagonal votes.

    Each segment of the read votes for its own diagonal near d0. A row's band
    spans the diagonals of its segment and both neighbours, plus w either side,
    so the band is narrow where the diagonal is stable and widens only where it
    changes (around an indel).
    """
    eq = read[:, None] == ref[None, :]
    m, n = eq.shape
    order = sorted(range(d0 - RADIUS, d0 + RADIUS + 1), key=lambda d: abs(d - d0))
    seg_diag = []
    for s in range(0, m, SEGMENT):
        rows = np.arange(s, min(s + SEGMENT, m))
        best_d, best_v = None, -1
        for d in order:
            cols = rows + d
            ok = (cols >= 0) & (cols < n)
            valid = int(ok.sum())
            if valid < 10:
                continue
            v = int(eq[rows[ok], cols[ok]].sum())
            if v >= MIN_VOTE_FRAC * valid and v > best_v:
                best_d, best_v = d, v
        seg_diag.append(best_d)

    last = None                                   # fill segments with no clear vote
    for t in range(len(seg_diag)):
        if seg_diag[t] is None:
            seg_diag[t] = last
        last = seg_diag[t]
    nxt = None
    for t in reversed(range(len(seg_diag))):
        if seg_diag[t] is None:
            seg_diag[t] = nxt if nxt is not None else d0
        nxt = seg_diag[t]

    lo = np.empty(m, dtype=np.int64)
    hi = np.empty(m, dtype=np.int64)
    for i in range(m):
        t = i // SEGMENT
        near = seg_diag[max(t - 1, 0):t + 2]
        lo[i] = i + min(near) - w
        hi[i] = i + max(near) + w
    return lo, hi


# ------------------------------------------------- 3. event-driven estimate

def wavefront_edit(a, b):
    """Myers O(ND) alignment. Returns (D, wavefront updates, base compares).

    D counts inserted plus deleted bases, so a substitution costs 2. An update
    is one (D, diagonal) step; compares are the bases checked while sliding
    along runs of matches.
    """
    n, m = len(a), len(b)
    v = {1: 0}
    updates = compares = 0
    for d in range(n + m + 1):
        for k in range(-d, d + 1, 2):
            updates += 1
            if k == -d or (k != d and v.get(k - 1, -1) < v.get(k + 1, -1)):
                x = v.get(k + 1, 0)
            else:
                x = v.get(k - 1, 0) + 1
            y = x - k
            while x < n and y < m:
                compares += 1
                if a[x] != b[y]:
                    break
                x, y = x + 1, y + 1
            v[k] = x
            if x >= n and y >= m:
                return d, updates, compares
    raise AssertionError("unreachable")


# ------------------------------------------------- experiment

def run(n_per_cell, seed):
    rng = random.Random(seed)
    regions = read_fasta(ROOT / "data" / "regions.fasta")
    conditions = [("none", 0)] + [(k, L) for k in ("del", "ins") for L in INDEL_LENGTHS]
    diffs = {"horizontal": [0, 0], "vertical": [0, 0], "diagonal": [0, 0]}
    max_score = 0
    rows = []

    for locus, (centre, half) in LOCI.items():
        region = regions[locus]
        window = region[centre - half:centre + half]
        ref = encode(window)
        for kind, length in conditions:
            for _ in range(n_per_cell):
                read_s, _ = simulate_read(rng, region, centre, half, kind, length)
                read = encode(read_s)

                H = sw_matrix(read, ref)
                full = int(H.max())
                max_score = max(max_score, full)
                for name, d in (("horizontal", H[1:, 1:] - H[1:, :-1]),
                                ("vertical", H[1:, 1:] - H[:-1, 1:]),
                                ("diagonal", H[1:, 1:] - H[:-1, :-1])):
                    diffs[name][0] = min(diffs[name][0], int(d.min()))
                    diffs[name][1] = max(diffs[name][1], int(d.max()))

                d0, _ = diagonal_vote(read, ref)
                full_cells = READ_LEN * len(ref)
                row = {"indel_len": length, "full": full}
                i = np.arange(READ_LEN)
                for w in (4, 8, 16):
                    s, c = sw_score(read, ref, lo=i + d0 - w, hi=i + d0 + w)
                    row[f"fixed{w}"] = s == full
                    row[f"fixed{w}_cells"] = c
                for w in (1, 2, 4):
                    lo, hi = segment_band(read, ref, d0, w)
                    s, c = sw_score(read, ref, lo=lo, hi=hi)
                    row[f"attn{w}"] = s == full
                    row[f"attn{w}_cells"] = c
                row["full_cells"] = full_cells

                _, path = sw_align(read_s, window)
                (i0, j0, _), (i1, j1, _) = path[0], path[-1]
                dist, upd, cmp_ = wavefront_edit(read_s[i0:i1 + 1], window[j0:j1 + 1])
                row.update(edit=dist, wf_updates=upd, wf_compares=cmp_,
                           aligned_rows=i1 - i0 + 1)
                rows.append(row)
    return rows, diffs, max_score


def report(rows, diffs, max_score):
    print("1. Difference encoding")
    print(f"   absolute score: max {max_score} -> {bits_for(0, max_score)} bits unsigned")
    for name, (lo, hi) in diffs.items():
        print(f"   {name:10s} difference: {lo:+d} .. {hi:+d} -> {bits_for(lo, hi)} bits")

    print("\n2. Attention-guided band: fraction of reads equal to full SW")
    keys = ["fixed4", "fixed8", "fixed16", "attn1", "attn2", "attn4"]
    print("   " + "  ".join(f"{h:>7s}" for h in ["indel"] + keys))
    for L in [0] + INDEL_LENGTHS:
        sub = [r for r in rows if r["indel_len"] == L]
        print("   " + "  ".join([f"{L:>7d}"] + [f"{np.mean([r[k] for r in sub]):>7.3f}" for k in keys]))
    print("   mean cells per read")
    print("   " + "  ".join([f"{'':>7s}"] + [f"{np.mean([r[k + '_cells'] for r in rows]):>7.0f}" for k in keys]))
    print(f"   (full matrix: {np.mean([r['full_cells'] for r in rows]):.0f} cells)")

    print("\n3. Event-driven estimate, aligned part of each read (mean per read)")
    print("   " + "  ".join(f"{h:>9s}" for h in ["indel", "D", "updates", "compares", "band8 cells"]))
    for L in [0] + INDEL_LENGTHS:
        sub = [r for r in rows if r["indel_len"] == L]
        band = np.mean([17 * r["aligned_rows"] for r in sub])
        print("   " + "  ".join([f"{L:>9d}"] + [f"{np.mean([r[k] for r in sub]):>9.1f}"
              for k in ("edit", "wf_updates", "wf_compares")] + [f"{band:>9.0f}"]))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=60, help="reads per locus per condition")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    report(*run(args.n, args.seed))
