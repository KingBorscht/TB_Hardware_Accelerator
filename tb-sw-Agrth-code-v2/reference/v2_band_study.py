"""How narrow can the v2 band be on long noisy reads, and does attention help?

SCOPE v2 aligns 1-3 kb nanopore reads with a fixed band of 128 diagonals
(+/-64). The band sets the PE count and the traceback memory (3,072 x W x 4
bits), so a narrower band that still finds the full-SW score saves both.

Band placements compared (all deterministic, all decided before or while the
row is computed, so all are possible in hardware):
  fixed    centre = start diagonal (0, known from the primer match), whole read
  tilted   centre on the straight line from the start to the end of the read
           (end diagonal = amplicon length - read length; needs a read that
           runs primer to primer)
  track    starts at 0; after each row the centre moves 1 diagonal toward the
           row's best cell (a running tracker, no look-ahead)
  vote     the attention vote from v0.2, per segment: every SEG read bases,
           Z = matches - mismatches is counted for each diagonal within the
           band, on the next SEG read bases, and the best becomes the centre

Score: v2 affine (v2_spec.py). A read counts as a hit when the banded score
equals the full-SW score. Reads are noisy copies of 2.5 kb H37Rv windows
(data/amplicons_v2.fasta, placeholder coordinates). The error model is an
ASSUMPTION (independent substitutions and geometric-length indels, no
homopolymer effects); replace it with error rates measured from the club's
Dorado reads.

Run:  python reference/v2_band_study.py              (about 10 minutes)
      python reference/v2_band_study.py -n 5         (quick check)
Out:  results/v2_band_<date>_<time>.csv / .md
"""

import argparse
import csv
import random
from datetime import datetime

import numpy as np

from band_study import ROOT, read_fasta
from v2_spec import EXTEND, MATCH, MISMATCH, OPEN, sw_affine

NOISE = [(0.01, 0.01), (0.01, 0.02), (0.02, 0.04), (0.03, 0.06)]   # (substitution, indel) per base
DEL_SHARE = 0.65                                                   # share of indel events that are deletions
MEAN_INDEL = 1.3                                                   # mean indel length, geometric
WIDTHS = [8, 12, 16, 24, 32, 48, 64, 96, 128]                             # band width W (centre +/- W/2)
STRATEGIES = ["fixed", "tilted", "track", "vote"]
SEG = 32
NEG = -(1 << 30)


def noisy_read(rng, ref, sub, indel):
    out, i = [], 0
    while i < len(ref):
        r = rng.random()
        if r < indel * DEL_SHARE:
            i += geometric(rng)
            continue
        if r < indel:
            out.extend(rng.choice("ACGT") for _ in range(geometric(rng)))
        c = ref[i]
        out.append(rng.choice([b for b in "ACGT" if b != c]) if rng.random() < sub else c)
        i += 1
    return "".join(out)


def geometric(rng):
    k = 1
    while rng.random() > 1 / MEAN_INDEL:
        k += 1
    return k


class Aligner:
    """Row-by-row affine SW in numpy. The band for each row comes from a strategy."""

    def __init__(self, read, ref):
        self.read, self.ref = read, ref
        self.m, self.n = len(read), len(ref)
        r = np.frombuffer(ref.encode(), dtype=np.uint8)
        self.row_score = {c: np.where(r == ord(c), MATCH, MISMATCH).astype(np.int64) for c in "ACGT"}
        self.cols = np.arange(self.n + 1, dtype=np.int64)

    def run(self, centre_of, half):
        """centre_of(i, prev_row_best_diag, centre) -> centre for row i. Returns best score."""
        n = self.n
        H_prev = np.zeros(n + 1, dtype=np.int64)
        F_prev = np.full(n + 1, NEG, dtype=np.int64)
        best, centre, row_best_d = 0, 0, None
        for i in range(1, self.m + 1):
            centre = centre_of(i, row_best_d, centre)
            lo, hi = max(1, i + centre - half), min(n, i + centre + half)
            H = np.zeros(n + 1, dtype=np.int64)
            F = np.full(n + 1, NEG, dtype=np.int64)
            if lo > hi:
                H_prev, F_prev, row_best_d = H, F, None
                continue
            D = H_prev[lo - 1:hi] + self.row_score[self.read[i - 1]][lo - 1:hi]
            Fr = np.maximum(H_prev[lo:hi + 1] - (OPEN + EXTEND), F_prev[lo:hi + 1] - EXTEND)
            Hp = np.maximum(np.maximum(D, Fr), 0)
            # E[j] = max over k < j in the band of Hp[k] - OPEN - EXTEND * (j - k)
            A = np.maximum.accumulate(Hp + EXTEND * self.cols[lo:hi + 1])
            E = np.empty_like(Hp)
            E[0] = NEG
            E[1:] = A[:-1] - OPEN - EXTEND * self.cols[lo + 1:hi + 1]
            Hr = np.maximum(Hp, E)
            H[lo:hi + 1], F[lo:hi + 1] = Hr, Fr
            k = int(Hr.argmax())
            if Hr[k] > best:
                best = int(Hr[k])
            row_best_d = (lo + k) - i if Hr[k] > 0 else None
            H_prev, F_prev = H, F
        return best


def strategy(name, read, ref, half):
    m, n = len(read), len(ref)
    if name == "fixed":
        return lambda i, d, c: 0
    if name == "tilted":
        return lambda i, d, c: round(i * (n - m) / m)
    if name == "track":
        return lambda i, d, c: c if d is None else c + (d > c) - (d < c)
    if name == "vote":
        def vote(i, d, c):
            if (i - 1) % SEG:
                return c
            seg = read[i - 1:i - 1 + SEG]
            best_d, best_z = c, None
            for k in sorted(range(c - half, c + half + 1), key=lambda k: (abs(k - c), k)):
                z = 0
                for t, a in enumerate(seg):
                    j = i - 1 + t + k
                    if 0 <= j < n:
                        z += 1 if a == ref[j] else -1
                if best_z is None or z > best_z:
                    best_d, best_z = k, z
            return best_d
        return vote
    raise ValueError(name)


def check_against_spec(rng, trials=60):
    """The numpy scores must equal the spec model's, on random banded cases."""
    for _ in range(trials):
        m, n = rng.randint(5, 60), rng.randint(5, 60)
        read = "".join(rng.choice("ACGT") for _ in range(m))
        ref = "".join(rng.choice("ACGT") for _ in range(n))
        half = rng.randint(1, 12)
        centres = [0] + [rng.randint(-5, 5) for _ in range(m)]
        fast = Aligner(read, ref).run(lambda i, d, c: centres[i], half)
        slow = sw_affine(read, ref, centres, half)[0]
        assert fast == slow, (read, ref, centres, half, fast, slow)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=25, help="reads per amplicon per noise level")
    ap.add_argument("--seed", type=int, default=0)
    opts = ap.parse_args()
    rng = random.Random(opts.seed)
    check_against_spec(rng)

    amplicons = read_fasta(ROOT / "data" / "amplicons_v2.fasta")
    rows = []
    for sub, indel in NOISE:
        for amp, ref in amplicons.items():
            for _ in range(opts.n):
                read = noisy_read(rng, ref, sub, indel)
                al = Aligner(read, ref)
                full = al.run(lambda i, d, c: 0, len(read) + len(ref))
                drift = len(ref) - len(read)
                for name in STRATEGIES:
                    for W in WIDTHS:
                        s = al.run(strategy(name, read, ref, W // 2), W // 2)
                        rows.append({"sub": sub, "indel": indel, "amplicon": amp, "read_len": len(read),
                                     "end_drift": drift, "strategy": name, "width": W,
                                     "score": s, "full_score": full, "hit": s == full})
        print(f"noise sub {sub:.0%} indel {indel:.0%} done")

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out = ROOT / "results"
    out.mkdir(exist_ok=True)
    csv_path, md_path = out / f"v2_band_{stamp}.csv", out / f"v2_band_{stamp}.md"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    reads = len(rows) // (len(STRATEGIES) * len(WIDTHS))
    lines = [f"# v2 band width study ({datetime.now().isoformat(timespec='seconds')})", "",
             f"{reads} simulated reads of 2.5 kb H37Rv windows, v2 affine scoring "
             f"(+{MATCH}/{MISMATCH}, open {OPEN}, extend {EXTEND}), deletions {DEL_SHARE:.0%} of indels, "
             f"mean indel length {MEAN_INDEL}, seed {opts.seed}. Error model is an assumption. "
             "Cells: share of reads whose banded score equals the full-SW score.", ""]
    for sub, indel in NOISE:
        sel = [r for r in rows if r["sub"] == sub and r["indel"] == indel]
        lines += [f"## Substitutions {sub:.0%}, indels {indel:.0%}", "",
                  "| band width W | " + " | ".join(STRATEGIES) + " |",
                  "|---:|" + "---:|" * len(STRATEGIES)]
        for W in WIDTHS:
            cells = []
            for name in STRATEGIES:
                h = [r["hit"] for r in sel if r["strategy"] == name and r["width"] == W]
                cells.append(f"{sum(h) / len(h):.1%}")
            lines.append(f"| {W} | " + " | ".join(cells) + " |")
        lines.append("")
    lines += ["Traceback memory at L_max 3,072 and 4 bits per cell: "
              + ", ".join(f"W={W}: {3072 * W * 4 / 1e6:.2f} Mbit" for W in WIDTHS) + ".",
              "", f"Per-read data: `{csv_path.name}`.", ""]
    md_path.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    print("saved", csv_path, "and", md_path.name)


if __name__ == "__main__":
    main()
