"""How compressible are the v2 traceback pointers? (video-style row prediction)

The traceback memory stores a 4-bit pointer for every cell in the band. Treat
each row of the band as a video frame and ask how many bits per cell a lossless
coder would need, using the entropy of the pointer given what it can be
predicted from:
  raw        4 bits, as SCOPE v2 stores it
  order 0    entropy of the pointer on its own
  left       given the pointer to its left in the same row (intra prediction)
  previous   given the pointer on the same diagonal in the previous row
             (inter prediction; following the diagonal is the motion compensation
             for the moving band)
  both       given both
Entropy is a lower bound for an ideal coder; a simple hardware coder does worse.
'flag code' is one such simple coder: 1 bit "same as the previous row on this
diagonal", plus the 4-bit pointer only when it differs.

Run:  python reference/v2_pointer_entropy.py
"""

import math
import random
from collections import Counter

from band_study import ROOT, read_fasta
from v2_band_study import NOISE, Aligner, noisy_read, strategy
from v2_spec import sw_affine

WIDTH = 32
READS = 4                     # per amplicon per noise level


def entropy(pairs):
    """H(X | context) in bits from a list of (context, x)."""
    joint, ctx = Counter(pairs), Counter(c for c, _ in pairs)
    total = len(pairs)
    return -sum(k / total * math.log2(k / ctx[c]) for (c, _), k in joint.items())


def tracker_centres(read, ref, half):
    """Centres chosen by the 'track' strategy, recorded row by row."""
    centres, step = [0], strategy("track", read, ref, half)

    def record(i, d, c):
        centres.append(step(i, d, c))
        return centres[-1]
    Aligner(read, ref).run(record, half)
    return centres


def main():
    rng = random.Random(0)
    amplicons = read_fasta(ROOT / "data" / "amplicons_v2.fasta")
    half = WIDTH // 2
    print(f"band W = {WIDTH + 1} cells per row, tracker band, {READS} reads per amplicon per noise level")
    print(f"{'sub/indel':>10s} {'raw':>5s} {'order 0':>8s} {'left':>6s} {'previous':>9s} {'both':>6s} {'rows equal to previous':>23s} {'same as previous':>17s} {'flag code':>10s}")
    for sub, indel in NOISE:
        pairs0, pairs_l, pairs_p, pairs_b = [], [], [], []
        same_rows = rows = 0
        for ref in amplicons.values():
            for _ in range(READS):
                read = noisy_read(rng, ref, sub, indel)
                centres = tracker_centres(read, ref, half)
                _, _, ptr, _ = sw_affine(read, ref, centres, half)
                prev_row = None
                for i in range(1, len(read) + 1):
                    row = {j - i: ptr[(i, j)] for j in range(i + centres[i] - half, i + centres[i] + half + 1)
                           if (i, j) in ptr}                      # keyed by diagonal
                    left = None
                    for d in sorted(row):
                        p = row[d]
                        up = prev_row.get(d) if prev_row else None
                        pairs0.append((0, p))
                        pairs_l.append((left, p))
                        pairs_p.append((up, p))
                        pairs_b.append(((left, up), p))
                        left = p
                    rows += 1
                    same_rows += prev_row is not None and all(prev_row.get(d) == p for d, p in row.items())
                    prev_row = row
        hit = sum(c == x for c, x in pairs_p) / len(pairs_p)
        print(f"{sub:>4.0%}/{indel:<4.0%} {4:>5.2f} {entropy(pairs0):>8.2f} {entropy(pairs_l):>6.2f} "
              f"{entropy(pairs_p):>9.2f} {entropy(pairs_b):>6.2f} {same_rows / rows:>23.1%} {hit:>17.1%} {1 + 4 * (1 - hit):>10.2f}")


if __name__ == "__main__":
    main()
