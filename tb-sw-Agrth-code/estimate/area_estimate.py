"""Rough area estimate: does the sliding-window design fit a small FPGA?

Two blocks are written in Amaranth and synthesised with Yosys:
  WindowCore    W score cells, one read base per step, pointers out, best tracker
  DiagonalVote  slides the reference past a stored read and counts matches

This is a sizing sketch, not the project RTL. Absolute 9-bit scores, one whole
row evaluated per clock (a long combinational chain), no host interface, no
traceback memory, no place-and-route. Counts are iCE40 LUT4s from the built-in
Yosys mapper (-noabc), which is pessimistic; the ABC mapper would not run in
the YoWASP build on this machine. Use the numbers for "does it fit", not for
speed or power.

Run:  .venv/Scripts/python estimate/area_estimate.py     (from the project root)
Needs: pip install amaranth yowasp-yosys
"""

import re
import subprocess
import sys
import tempfile
from pathlib import Path

from amaranth import Elaboratable, Module, Mux, Signal, signed
from amaranth.back import rtlil

MATCH, MISMATCH, GAP = 2, -1, -2
STOP, DIAG, UP, LEFT = 0, 1, 2, 3


class WindowCore(Elaboratable):
    def __init__(self, W=17, score_bits=9, steps=150):
        self.W, self.score_bits, self.steps = W, score_bits, steps
        self.clear = Signal()
        self.preload = Signal()         # shift the reference window without scoring
        self.step = Signal()
        self.read_base = Signal(2)
        self.ref_in = Signal(2)
        self.ref_in_valid = Signal()
        self.pointers = Signal(2 * W)
        self.best = Signal(score_bits)
        self.best_step = Signal(range(steps))
        self.best_cell = Signal(range(W))

    def ports(self):
        return [self.clear, self.preload, self.step, self.read_base, self.ref_in,
                self.ref_in_valid, self.pointers, self.best, self.best_step, self.best_cell]

    def elaborate(self, platform):
        m = Module()
        W, B = self.W, self.score_bits
        ref = [Signal(2, name=f"ref{k}") for k in range(W)]
        valid = [Signal(name=f"valid{k}") for k in range(W)]
        h = [Signal(B, name=f"h{k}") for k in range(W)]
        new = [Signal(B, name=f"new{k}") for k in range(W)]
        step_count = Signal(range(self.steps))

        for k in range(W):
            diag = Signal(signed(B + 2), name=f"diag{k}")
            up = Signal(signed(B + 2), name=f"up{k}")
            left = Signal(signed(B + 2), name=f"left{k}")
            ptr = self.pointers[2 * k:2 * k + 2]
            m.d.comb += [
                diag.eq(h[k] + Mux(ref[k] == self.read_base, MATCH, MISMATCH)),
                up.eq((h[k + 1] if k + 1 < W else 0) + GAP),
                left.eq((new[k - 1] if k else 0) + GAP),
            ]
            # same tie-breaking as the Python reference: diagonal > up > left
            with m.If(~valid[k]):
                m.d.comb += [new[k].eq(0), ptr.eq(STOP)]
            with m.Elif((diag > 0) & (diag >= up) & (diag >= left)):
                m.d.comb += [new[k].eq(diag), ptr.eq(DIAG)]
            with m.Elif((up > 0) & (up >= left)):
                m.d.comb += [new[k].eq(up), ptr.eq(UP)]
            with m.Elif(left > 0):
                m.d.comb += [new[k].eq(left), ptr.eq(LEFT)]
            with m.Else():
                m.d.comb += [new[k].eq(0), ptr.eq(STOP)]

        row_best = new[0]
        row_cell = Signal(range(W), name="row_cell0")
        for k in range(1, W):
            nb = Signal(B, name=f"row_best{k}")
            nc = Signal(range(W), name=f"row_cell{k}")
            m.d.comb += [nb.eq(Mux(new[k] > row_best, new[k], row_best)),
                         nc.eq(Mux(new[k] > row_best, k, row_cell))]
            row_best, row_cell = nb, nc

        with m.If(self.clear):
            m.d.sync += [s.eq(0) for s in h] + [v.eq(0) for v in valid]
            m.d.sync += [self.best.eq(0), step_count.eq(0)]
        with m.Elif(self.preload | self.step):
            m.d.sync += [ref[k].eq(ref[k + 1]) for k in range(W - 1)] + [ref[W - 1].eq(self.ref_in)]
            m.d.sync += [valid[k].eq(valid[k + 1]) for k in range(W - 1)] + [valid[W - 1].eq(self.ref_in_valid)]
        with m.If(self.step & ~self.clear):
            m.d.sync += [h[k].eq(new[k]) for k in range(W)]
            m.d.sync += step_count.eq(step_count + 1)
            with m.If(row_best > self.best):
                m.d.sync += [self.best.eq(row_best), self.best_step.eq(step_count),
                             self.best_cell.eq(row_cell)]
        return m


class DiagonalVote(Elaboratable):
    def __init__(self, L=150, shifts=330):
        self.L = L
        self.clear = Signal()
        self.load = Signal()            # shift one read base in
        self.read_in = Signal(2)
        self.shift = Signal()           # slide the reference by one base
        self.ref_in = Signal(2)
        self.ref_in_valid = Signal()
        self.best_count = Signal(range(L + 1))
        self.best_shift = Signal(range(shifts))

    def ports(self):
        return [self.clear, self.load, self.read_in, self.shift, self.ref_in,
                self.ref_in_valid, self.best_count, self.best_shift]

    def elaborate(self, platform):
        m = Module()
        L = self.L
        read = [Signal(2, name=f"read{k}") for k in range(L)]
        ref = [Signal(2, name=f"ref{k}") for k in range(L)]
        valid = [Signal(name=f"valid{k}") for k in range(L)]
        shift_count = Signal(self.best_shift.shape())

        terms = [valid[k] & (read[k] == ref[k]) for k in range(L)]
        while len(terms) > 1:                      # balanced adder tree
            terms = [terms[i] + terms[i + 1] if i + 1 < len(terms) else terms[i]
                     for i in range(0, len(terms), 2)]
        count = Signal(range(L + 1))
        m.d.comb += count.eq(terms[0])

        with m.If(self.load):
            m.d.sync += [read[k].eq(read[k + 1]) for k in range(L - 1)] + [read[L - 1].eq(self.read_in)]
        with m.If(self.clear):
            m.d.sync += [v.eq(0) for v in valid] + [self.best_count.eq(0), shift_count.eq(0)]
        with m.Elif(self.shift):
            m.d.sync += [ref[k].eq(ref[k + 1]) for k in range(L - 1)] + [ref[L - 1].eq(self.ref_in)]
            m.d.sync += [valid[k].eq(valid[k + 1]) for k in range(L - 1)] + [valid[L - 1].eq(self.ref_in_valid)]
            m.d.sync += shift_count.eq(shift_count + 1)
            with m.If(count > self.best_count):
                m.d.sync += [self.best_count.eq(count), self.best_shift.eq(shift_count)]
        return m


def synth(design):
    """LUT4 and flip-flop counts from Yosys synth_ice40, built-in LUT mapper."""
    text = rtlil.convert(design, name="top", ports=design.ports())
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / "top.il"
        src.write_bytes(text.encode())      # keep LF line endings on Windows
        cmd = [str(Path(sys.executable).with_name("yowasp-yosys")), "-p",
               f"read_rtlil {src.as_posix()}; synth_ice40 -top top -noabc; stat"]
        out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    stats = out[out.rindex("Printing statistics"):]
    cells = {name: int(n) for n, name in re.findall(r"^\s+(\d+)\s+(SB_\w+)\s*$", stats, re.M)}
    luts = cells.get("SB_LUT4", 0)
    ffs = sum(n for c, n in cells.items() if c.startswith("SB_DFF"))
    return luts, ffs


if __name__ == "__main__":
    blocks = [(f"WindowCore W={W}", lambda W=W: WindowCore(W=W)) for W in (5, 9, 17)]
    blocks.append(("DiagonalVote L=150", lambda: DiagonalVote()))
    print(f"{'block':22s}{'LUT4':>8s}{'FF':>8s}")
    for name, make in blocks:
        luts, ffs = synth(make())
        print(f"{name:22s}{luts:>8d}{ffs:>8d}")
    print()
    print("capacity: iCE40UP5K 5,280 LUT4 | Tang Nano 9K (GW1NR-9) 8,640 LUT4 | Zynq-7020 53,200 LUT6")
