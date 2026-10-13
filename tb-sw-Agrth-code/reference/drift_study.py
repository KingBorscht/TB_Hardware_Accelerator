"""How far does an alignment path drift from its starting diagonal on a long read?

SCOPE v2 uses a fixed band of 128 diagonals (+/-64) on reads of up to 3,072 bp.
The path moves one diagonal for every inserted or deleted base, so over a long
noisy read it performs a random walk. This script simulates only that walk (no
alignment) to see how wide the band has to be.

Two ways to place the band:
  fixed     centred on the starting diagonal for the whole read
  tilted    centred on the straight line from the start to the end of the read,
            which is known in advance when the read spans the whole amplicon
            (final drift = reference length - read length)

The indel rates are ASSUMPTIONS, not measurements. Replace them with rates
measured from the club's own Dorado output before drawing conclusions.

Run:  python reference/drift_study.py
"""

import numpy as np

L = 3072                 # read length, SCOPE v2 L_max
HALF = 64                # SCOPE v2 band half-width
MEAN_LEN = 1.3           # mean indel length (geometric), assumed
N = 20000                # simulated reads per setting
RATES = [0.005, 0.01, 0.02, 0.04]      # indel events per base
DEL_SHARE = [0.5, 0.65, 0.8]           # fraction of indel events that are deletions


def simulate(rng, rate, del_share, n):
    """Return (max |drift| with a fixed band, max |deviation| from the tilted line)."""
    p_del, p_ins = rate * del_share, rate * (1 - del_share)
    fixed, tilted = [], []
    t = np.arange(1, L + 1)
    for start in range(0, n, 2000):
        k = min(2000, n - start)
        step = rng.choice(np.array([-1, 0, 1], dtype=np.int32), p=[p_ins, 1 - rate, p_del], size=(k, L))
        hit = step != 0
        step[hit] *= rng.geometric(1 / MEAN_LEN, size=int(hit.sum())).astype(np.int32)
        d = np.cumsum(step, axis=1)
        fixed.append(np.abs(d).max(axis=1))
        line = d[:, -1:] * t[None, :] / L
        tilted.append(np.abs(d - line).max(axis=1))
    return np.concatenate(fixed), np.concatenate(tilted)


if __name__ == "__main__":
    rng = np.random.default_rng(0)
    print(f"read length {L}, band +/-{HALF}, mean indel length {MEAN_LEN}, {N} simulated reads per row")
    print(f"{'indel rate':>10s} {'deletions':>9s} | {'fixed: leaves +/-64':>20s} {'99.9% needs +/-':>16s} |"
          f" {'tilted: leaves +/-64':>21s} {'99.9% needs +/-':>16s}")
    for rate in RATES:
        for share in DEL_SHARE:
            fx, tl = simulate(rng, rate, share, N)
            print(f"{rate:>10.1%} {share:>9.0%} | {np.mean(fx > HALF):>20.2%} {np.quantile(fx, 0.999):>16.0f} |"
                  f" {np.mean(tl > HALF):>21.2%} {np.quantile(tl, 0.999):>16.0f}")
