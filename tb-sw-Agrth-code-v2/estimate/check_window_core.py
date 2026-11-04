"""Simulate WindowCore and compare its best score with the Python sliding-window model."""

import random
import sys
from pathlib import Path

from amaranth.sim import Simulator

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "reference"))
from area_estimate import WindowCore                                    # noqa: E402
from band_study import LOCI, ROOT, diagonal_vote, encode, read_fasta, simulate_read  # noqa: E402
from window_study import sliding_window, window_offsets                 # noqa: E402

CODE = {"A": 0, "C": 1, "G": 2, "T": 3}


def hardware_best(read, window, d0, W):
    dut = WindowCore(W=W, steps=len(read))
    lo, _ = window_offsets(W)
    n = len(window)
    result = {}

    def base(col):
        ok = 0 <= col < n
        return (CODE[window[col]] if ok else 0), int(ok)

    async def bench(ctx):
        ctx.set(dut.clear, 1)
        await ctx.tick()
        ctx.set(dut.clear, 0)
        ctx.set(dut.preload, 1)
        for k in range(W):                       # fill the window for step 0
            b, v = base(d0 + lo + k)
            ctx.set(dut.ref_in, b)
            ctx.set(dut.ref_in_valid, v)
            await ctx.tick()
        ctx.set(dut.preload, 0)
        ctx.set(dut.step, 1)
        for i, r in enumerate(read):
            b, v = base(i + 1 + d0 + lo + W - 1)
            ctx.set(dut.read_base, CODE[r])
            ctx.set(dut.ref_in, b)
            ctx.set(dut.ref_in_valid, v)
            await ctx.tick()
        ctx.set(dut.step, 0)
        await ctx.tick()
        result["best"] = ctx.get(dut.best)

    sim = Simulator(dut)
    sim.add_clock(1e-6)
    sim.add_testbench(bench)
    sim.run()
    return result["best"]


if __name__ == "__main__":
    rng = random.Random(7)
    regions = read_fasta(ROOT / "data" / "regions.fasta")
    checked = bad = 0
    for locus, (centre, half) in LOCI.items():
        region = regions[locus]
        window = region[centre - half:centre + half]
        for kind, length in [("none", 0), ("del", 2), ("ins", 3), ("del", 8), ("ins", 12)]:
            for W in (5, 17):
                read, _ = simulate_read(rng, region, centre, half, kind, length)
                d0, _ = diagonal_vote(encode(read), encode(window))
                H, _ = sliding_window(encode(read), encode(window), d0, W)
                got = hardware_best(read, window, d0, W)
                checked += 1
                if got != int(H.max()):
                    bad += 1
                    print("MISMATCH", locus, kind, length, W, "model", int(H.max()), "hardware", got)
    print(f"checked {checked} reads, {bad} mismatches")
