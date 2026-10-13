"""Exact vote margins computed from the reference alone (no simulated reads).

Vote statistic on diagonal k, with mismatch penalty lam >= 0:

    Z(k) = matches(k) - lam * mismatches(k) = (1 + lam) * V(k) - lam * N(k)

V(k) is the number of matching positions and N(k) the number of compared
positions, which depends only on k. lam = 0 is the plain match count.

For every read start p, and every single indel (position q, length g, type),
the error-free model read is built and Z is evaluated on every diagonal. The
margin is

    best Z on an acceptable diagonal - best Z on any other diagonal

A diagonal d is acceptable when the window [d - w, d + w] contains every
diagonal the true alignment visits, that is max(k1, k2) - w <= d <= min(k1, k2) + w.

One substitution changes Z on any diagonal by at most (1 + lam), so a read with
e substitutions still votes for an acceptable diagonal whenever
2 (1 + lam) e < margin. Inserted bases are treated as worst case: matches on
every unacceptable diagonal, mismatches on the acceptable ones.

Run:  python reference/vote_margin.py           (from the project root)
"""

import math

import numpy as np

from band_study import LOCI, READ_LEN, ROOT, read_fasta

HALF_WIDTH = 8                # w, for a W = 17 window
MAX_INDEL = HALF_WIDTH        # longest indel the window can follow
PAD = READ_LEN + 20


def diagonal_sums(region, w0, n):
    """Cumulative match counts along diagonals of (region position, window index)."""
    reg = np.frombuffer(region.encode(), dtype=np.uint8)
    win = reg[w0:w0 + n]
    eq = np.zeros((len(reg), n + 2 * PAD), dtype=np.int32)
    eq[:, PAD:PAD + n] = reg[:, None] == win[None, :]
    cs = np.zeros((len(reg) + 1, n + 2 * PAD + 1), dtype=np.int32)
    cs[1:, 1:] = eq
    for x in range(1, len(reg) + 1):
        cs[x, 1:] += cs[x - 1, :-1]
    return cs


def segment(cs, x0, ks, i0, i1):
    """Matches of read positions [i0, i1), which copy region[x0 + i], on every diagonal in ks."""
    if i1 <= i0:
        return np.zeros(len(ks), dtype=np.int64)
    return (cs[x0 + i1, ks + i1 + PAD] - cs[x0 + i0, ks + i0 + PAD]).astype(np.int64)


def compared(ks, n, m):
    """N(k): read positions i in [0, m) with 0 <= i + k < n."""
    return np.maximum(np.minimum(m, n - ks) - np.maximum(0, -ks), 0)


def margins(locus, region, lam):
    centre, half = LOCI[locus]
    w0, n, m = centre - half, 2 * half, READ_LEN
    cs = diagonal_sums(region, w0, n)
    ks = np.arange(-(m - 1), n)
    penalty = lam * compared(ks, n, m)
    starts = range(centre - m + 20, centre - 20 + 1)      # same read starts as the simulation
    out = {"none": [], "del": {g: [] for g in range(1, MAX_INDEL + 1)},
           "ins": {g: [] for g in range(1, MAX_INDEL + 1)}}

    for p in starts:
        k1 = p - w0
        z = (1 + lam) * segment(cs, p, ks, 0, m) - penalty
        ok = np.abs(ks - k1) <= HALF_WIDTH
        out["none"].append(int(z[ok].max() - z[~ok].max()))

        for g in range(1, MAX_INDEL + 1):
            inside = (ks >= k1 + g - HALF_WIDTH) & (ks <= k1 + HALF_WIDTH)   # deletion: k2 = k1 + g
            for q in range(1, m):
                v = segment(cs, p, ks, 0, q) + segment(cs, p + g, ks, q, m)
                z = (1 + lam) * v - penalty
                out["del"][g].append(int(z[inside].max() - z[~inside].max()))

            inside = (ks >= k1 - HALF_WIDTH) & (ks <= k1 - g + HALF_WIDTH)   # insertion: k2 = k1 - g
            for q in range(1, m - g):
                v = segment(cs, p, ks, 0, q) + segment(cs, p - g, ks, q + g, m)
                lo = np.clip(-ks, q, q + g)                # inserted positions landing in the window
                hi = np.clip(n - ks, q, q + g)
                wild = np.maximum(hi - lo, 0)
                z_true = (1 + lam) * v - penalty
                z_wrong = (1 + lam) * (v + wild) - penalty
                out["ins"][g].append(int(z_true[inside].max() - z_wrong[~inside].max()))
    return out


def tolerance(margin, lam):
    """Largest e with 2 (1 + lam) e < margin (-1 if even e = 0 is not guaranteed)."""
    step = 2 * (1 + lam)
    return math.ceil(margin / step) - 1 if margin > 0 else -1


def binom_tail(n, eps, k):
    """P(Binomial(n, eps) >= k)."""
    return sum(math.comb(n, i) * eps ** i * (1 - eps) ** (n - i) for i in range(max(k, 0), n + 1))


if __name__ == "__main__":
    regions = read_fasta(ROOT / "data" / "regions.fasta")
    for lam, label in ((0, "match count, lam = 0"), (1, "matches minus mismatches, lam = 1")):
        print(f"=== Vote statistic: {label} ===")
        print(f"{'locus':15s}{'case':>14s}{'configs':>10s}{'min margin':>12s}{'tolerated e':>13s}")
        overall, everything = None, []
        for locus in LOCI:
            res = margins(locus, regions[locus], lam)
            rows = [("no indel", res["none"]),
                    (f"del <= {MAX_INDEL} bp", [x for g in res["del"] for x in res["del"][g]]),
                    (f"ins <= {MAX_INDEL} bp", [x for g in res["ins"] for x in res["ins"][g]])]
            for case, vals in rows:
                lo = min(vals)
                overall = lo if overall is None else min(overall, lo)
                everything.extend(vals)
                print(f"{locus:15s}{case:>14s}{len(vals):>10d}{lo:>12d}{tolerance(lo, lam):>13d}")
        e_max = tolerance(overall, lam)
        print(f"overall minimum margin {overall}; guaranteed for up to {e_max} substitutions per read")
        if e_max >= 0:
            for eps in (0.005, 0.01, 0.02):
                print(f"  P(more than {e_max} substitutions in {READ_LEN} bases at rate {eps}) = "
                      f"{binom_tail(READ_LEN, eps, e_max + 1):.1e}")
        tol = np.array([tolerance(x, lam) for x in everything])
        print("  share of configurations guaranteed for at least e substitutions: "
              + ", ".join(f"e={e}: {np.mean(tol >= e):.4f}" for e in (0, 1, 2, 3, 5, 8)))
        print()
