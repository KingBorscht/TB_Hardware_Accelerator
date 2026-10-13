"""Band-width study for the two-tier alignment proposal.

Question: after a cheap diagonal vote (tier 1) picks the band centre, how wide
does a Smith-Waterman band have to be to return the same score as full SW?
And how often is the ungapped tier-1 score already equal to the full SW score?

Reads are simulated from the H37Rv regions in data/regions.fasta. All reads are
forward-strand (strand handling, SCOPE decision 4, is still open).

Run:  python reference/band_study.py            (from the project root)
Out:  docs/band_study.csv, docs/band_study.png, summary printed to stdout
"""

import argparse
import csv
import random
from pathlib import Path

import numpy as np

from sw_reference import GAP, MATCH, MISMATCH, sw_align

ROOT = Path(__file__).resolve().parent.parent
BASES = "ACGT"

# (hotspot centre index inside the fetched region, window half-width in bp).
# PROVISIONAL: window coordinates belong to the Bio focus (SCOPE decisions 1, 3).
LOCI = {
    "rpoB_RRDR": (240, 90),
    "katG_315": (200, 60),
    "inhA_promoter": (200, 60),
}

READ_LEN = 150
FIXED_WIDTHS = [0, 1, 2, 3, 4, 6, 8, 12, 16]   # band half-widths
ADAPTIVE_WIDTHS = [2, 4, 8]
INDEL_LENGTHS = [1, 2, 3, 5, 8, 12]
SUB_ERROR_RATE = 0.005                          # sequencing substitution errors
STEER_MIN_SCORE = 20                            # adaptive band ignores weaker rows


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = ""
        elif line:
            seqs[name] += line.upper()
    return seqs


def encode(seq):
    return np.frombuffer(seq.encode(), dtype=np.uint8)


# ---------------------------------------------------------------- scoring

def sw_score(read, ref, lo=None, hi=None, adaptive_w=None, d0=0):
    """Score-only SW, one numpy row at a time. Returns (score, cells computed).

    Fixed band: lo[i]..hi[i] are the inclusive reference columns of row i.
    Adaptive band: half-width adaptive_w, centre starts on diagonal d0 and may
    move one extra column left or right per row, towards the row maximum.
    Cells outside the band read as 0, so every in-band score is a real
    alignment whose path stays inside the band.
    """
    m, n = len(read), len(ref)
    cols = np.arange(n + 1)
    ramp = GAP * cols
    prev = np.zeros(n + 1, dtype=np.int64)
    best, cells = 0, 0
    centre = d0

    for i in range(m):
        if adaptive_w is not None:
            a, b = centre - adaptive_w, centre + adaptive_w
        elif lo is not None:
            a, b = lo[i], hi[i]
        else:
            a, b = 0, n - 1
        a, b = max(a, 0), min(b, n - 1)

        sub = np.where(ref == read[i], MATCH, MISMATCH)
        e = np.zeros(n + 1, dtype=np.int64)
        e[1:] = np.maximum(0, np.maximum(prev[:-1] + sub, prev[1:] + GAP))
        if a > b:
            cur = np.zeros(n + 1, dtype=np.int64)
        else:
            e[: a + 1] = 0
            e[b + 2:] = 0
            # H[j] = max(E[j], H[j-1] + GAP), solved with a running maximum
            cur = np.maximum.accumulate(e - ramp) + ramp
            cur[: a + 1] = 0
            cur[b + 2:] = 0
            cells += b - a + 1
        row_max = int(cur.max())
        best = max(best, row_max)

        if adaptive_w is not None:
            step = 0
            if row_max >= STEER_MIN_SCORE:
                jmax = int(cur.argmax()) - 1
                step = (jmax > centre) - (jmax < centre)
            centre += 1 + step
        prev = cur
    return best, cells


def diagonal_vote(read, ref, lam=0):
    """Tier 1: pick the diagonal with the best matches - lam * mismatches.

    lam = 0 is the plain match count. Returns (best diagonal, matches on it).
    Ties go to the smallest diagonal.
    """
    eq = read[:, None] == ref[None, :]
    m, n = eq.shape
    best_d, best_z, best_c = 0, None, 0
    for d in range(-(m - 1), n):
        c = int(np.trace(eq, offset=d))
        compared = min(m, n - d) - max(0, -d)
        z = (1 + lam) * c - lam * compared
        if best_z is None or z > best_z:
            best_d, best_z, best_c = d, z, c
    return best_d, best_c


def ungapped_score(read, ref, d):
    """Tier 1 score: best local score along diagonal d, no gaps (Kadane)."""
    eq = np.diagonal(read[:, None] == ref[None, :], offset=d)
    s = np.where(eq, MATCH, MISMATCH)
    p = np.concatenate(([0], np.cumsum(s)))
    return int((p - np.minimum.accumulate(p)).max())


# ---------------------------------------------------------------- simulation

def simulate_read(rng, region, centre, half, kind, length):
    """One 150 bp read overlapping the window, with SNPs and an optional indel.

    Returns (read, distance from the indel to the nearest end of the part of
    the read that overlaps the window), or None for the distance if no indel.
    """
    start = rng.randint(centre - READ_LEN + 20, centre - 20)
    template = list(region[start:start + READ_LEN + 20])
    ov_a = max(centre - half - start, 0)                 # overlap, read coords
    ov_b = min(centre + half - start, READ_LEN)

    for _ in range(rng.randint(0, 2)):                   # true SNPs
        p = rng.randrange(ov_a, ov_b)
        template[p] = rng.choice([b for b in BASES if b != template[p]])

    edge = None
    if kind != "none":
        p = rng.randrange(ov_a + 1, ov_b - length - 1)
        if kind == "del":
            del template[p:p + length]
            edge = min(p - ov_a, ov_b - length - p)
        else:
            template[p:p] = [rng.choice(BASES) for _ in range(length)]
            edge = min(p - ov_a, ov_b - p - length)

    read = template[:READ_LEN]
    for p in range(READ_LEN):                            # sequencing errors
        if rng.random() < SUB_ERROR_RATE:
            read[p] = rng.choice([b for b in BASES if b != read[p]])
    return "".join(read), edge


# ---------------------------------------------------------------- experiment

def run(n_per_cell, seed):
    rng = random.Random(seed)
    regions = read_fasta(ROOT / "data" / "regions.fasta")
    conditions = [("none", 0)] + [(k, L) for k in ("del", "ins") for L in INDEL_LENGTHS]
    rows = []

    for locus, (centre, half) in LOCI.items():
        region = regions[locus]
        window = region[centre - half:centre + half]
        ref = encode(window)
        n = len(ref)
        for kind, length in conditions:
            for _ in range(n_per_cell):
                read_s, edge = simulate_read(rng, region, centre, half, kind, length)
                read = encode(read_s)
                full, full_cells = sw_score(read, ref)
                d0, votes = diagonal_vote(read, ref)
                row = {
                    "locus": locus, "kind": kind, "indel_len": length,
                    "indel_edge_dist": "" if edge is None else edge,
                    "full_score": full, "d0": d0, "votes": votes,
                    "tier1_score": ungapped_score(read, ref, d0),
                }
                i = np.arange(READ_LEN)
                for w in FIXED_WIDTHS:
                    s, c = sw_score(read, ref, lo=i + d0 - w, hi=i + d0 + w)
                    row[f"fixed_w{w}"] = s
                    row[f"fixed_w{w}_cells"] = c / full_cells
                for w in ADAPTIVE_WIDTHS:
                    s, c = sw_score(read, ref, adaptive_w=w, d0=d0)
                    row[f"adapt_w{w}"] = s
                    row[f"adapt_w{w}_cells"] = c / full_cells
                rows.append(row)
    return rows


def self_check(seed=1, trials=40):
    """The fast scorer must equal the plain-Python reference exactly."""
    rng = random.Random(seed)
    for _ in range(trials):
        a = "".join(rng.choice(BASES) for _ in range(rng.randint(5, 60)))
        b = "".join(rng.choice(BASES) for _ in range(rng.randint(5, 60)))
        if rng.random() < 0.5:                           # related sequences
            p = rng.randrange(len(b))
            a = b[p:p + 30] + a[:5]
        want, _ = sw_align(a, b)
        got, _ = sw_score(encode(a), encode(b))
        assert want == got, (a, b, want, got)
    print(f"self-check: fast scorer == reference on {trials} random pairs")


def frac(rows, pred):
    return sum(1 for r in rows if pred(r)) / len(rows) if rows else float("nan")


def summarize(rows):
    print("\nFraction of reads whose score equals full SW")
    head = ["indel", "n", "tier1"] + [f"w={w}" for w in (1, 2, 4, 8, 16)] + [f"ad{w}" for w in ADAPTIVE_WIDTHS]
    print("  ".join(f"{h:>6s}" for h in head))
    table = {}
    for length in [0] + INDEL_LENGTHS:
        sub = [r for r in rows if r["indel_len"] == length]
        vals = [frac(sub, lambda r: r["tier1_score"] == r["full_score"])]
        vals += [frac(sub, lambda r, w=w: r[f"fixed_w{w}"] == r["full_score"]) for w in (1, 2, 4, 8, 16)]
        vals += [frac(sub, lambda r, w=w: r[f"adapt_w{w}"] == r["full_score"]) for w in ADAPTIVE_WIDTHS]
        table[length] = vals
        print("  ".join([f"{length:>6d}", f"{len(sub):>6d}"] + [f"{v:>6.3f}" for v in vals]))

    print("\nSmallest fixed half-width at which every read equals full SW")
    for length in [0] + INDEL_LENGTHS:
        sub = [r for r in rows if r["indel_len"] == length]
        ok = [w for w in FIXED_WIDTHS if all(r[f"fixed_w{w}"] == r["full_score"] for r in sub)]
        print(f"  indel {length:>2d} bp: {ok[0] if ok else '> ' + str(FIXED_WIDTHS[-1])}")

    print("\nCells computed, as a fraction of the full matrix (mean)")
    for key in [f"fixed_w{w}" for w in (2, 4, 8, 16)] + [f"adapt_w{w}" for w in ADAPTIVE_WIDTHS]:
        print(f"  {key:10s} {np.mean([r[key + '_cells'] for r in rows]):.3f}")

    missed = [r for r in rows if r["kind"] != "none" and r["tier1_score"] == r["full_score"]]
    indel = [r for r in rows if r["kind"] != "none"]
    print(f"\nIndel reads where tier 1 already equals full SW: {len(missed)}/{len(indel)}"
          " (the indel was clipped off or left unaligned)")
    if missed:
        d = [r["indel_edge_dist"] for r in missed]
        print(f"  their indel-to-overlap-end distance: median {np.median(d):.0f} bp, max {max(d)} bp")
    return table


def plot(rows, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    surface, ink, muted, grid = "#fcfcfb", "#0b0b0b", "#898781", "#e6e5e1"
    blues = {2: "#86b6ef", 4: "#5598e7", 8: "#2a78d6", 16: "#104281"}
    orange = "#eb6834"
    lengths = [0] + INDEL_LENGTHS
    x = list(range(len(lengths)))

    def series(key):
        return [frac([r for r in rows if r["indel_len"] == L],
                     lambda r: r[key] == r["full_score"]) for L in lengths]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), facecolor=surface, sharey=True)
    panels = [
        ("Fixed band centred by the diagonal vote",
         [(f"half-width {w}", series(f"fixed_w{w}"), blues[w]) for w in (2, 4, 8, 16)]),
        ("Tier 1 alone, and an adaptive band",
         [("tier 1 (ungapped)", series("tier1_score"), muted),
          ("adaptive, half-width 4", series("adapt_w4"), orange)]),
    ]
    for ax, (title, lines), legend_loc in zip(axes, panels, ("lower left", "upper right")):
        ax.set_facecolor(surface)
        for label, ys, colour in lines:
            ax.plot(x, ys, color=colour, lw=2, marker="o", ms=5, label=label)
        ax.set_title(title, loc="left", fontsize=11, color=ink)
        ax.set_xticks(x, [str(L) for L in lengths])
        ax.set_xlabel("Indel length in the read (bp); 0 = SNPs only", color=muted, fontsize=9)
        ax.set_ylim(-0.03, 1.05)
        ax.grid(axis="y", color=grid, lw=0.8)
        ax.tick_params(colors=muted, length=0, labelsize=9)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color(grid)
        ax.legend(frameon=False, fontsize=9, labelcolor=ink, loc=legend_loc)
    axes[0].set_ylabel("Reads with score equal to full SW", color=muted, fontsize=9)
    fig.suptitle("How wide does the band need to be?", x=0.06, ha="left", fontsize=13, color=ink)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=150, facecolor=surface)
    print(f"\nwrote {path}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=100, help="reads per locus per condition")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    self_check()
    rows = run(args.n, args.seed)
    out = ROOT / "docs" / "band_study.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out} ({len(rows)} reads)")
    summarize(rows)
    plot(rows, ROOT / "docs" / "band_study.png")
