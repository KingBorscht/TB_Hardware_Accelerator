import os

from amaranth import Elaboratable, Module, Signal
from amaranth.back import verilog

from sw_pe import SmithWatermanPE

PE_INPUTS = ("h_diag", "h_up", "f_up", "h_left", "e_left", "match")
PE_OUTPUTS = ("h", "e", "f", "tb")


class PETimingWrapper(Elaboratable):
    def __init__(self):
        self.pe = SmithWatermanPE()
        self.ins = {n: Signal(len(getattr(self.pe, n)), name=n) for n in PE_INPUTS}
        self.outs = {n: Signal(len(getattr(self.pe, n)), name=f"{n}_out") for n in PE_OUTPUTS}

    def elaborate(self, platform):
        m = Module()
        m.submodules.pe = self.pe
        for name, sig in self.ins.items():
            m.d.sync += getattr(self.pe, name).eq(sig)
        for name, sig in self.outs.items():
            m.d.sync += sig.eq(getattr(self.pe, name))
        return m


if __name__ == "__main__":
    top = PETimingWrapper()
    ports = [*top.ins.values(), *top.outs.values()]
    os.makedirs("build", exist_ok=True)
    with open("build/sw_pe_top.v", "w") as f:
        f.write(verilog.convert(top, name="sw_pe_top", ports=ports))