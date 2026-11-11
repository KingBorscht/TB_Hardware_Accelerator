"""AXI4-Stream loopback used to measure PYNQ DMA before the PE array exists."""

from amaranth import Elaboratable, Module, Signal


class AxisPassthrough(Elaboratable):
    """Zero-latency ready/valid AXI4-Stream pass-through.

    Connect ``s_*`` to an AXI DMA MM2S stream and ``m_*`` to its S2MM stream
    in Vivado. This validates the PS/PL DMA path independently of the array.
    """

    def __init__(self, *, data_width=64):
        self.s_valid = Signal()
        self.s_ready = Signal()
        self.s_data = Signal(data_width)
        self.s_last = Signal()

        self.m_valid = Signal()
        self.m_ready = Signal()
        self.m_data = Signal(data_width)
        self.m_last = Signal()

    def elaborate(self, platform):
        m = Module()
        m.d.comb += [
            self.s_ready.eq(self.m_ready),
            self.m_valid.eq(self.s_valid),
            self.m_data.eq(self.s_data),
            self.m_last.eq(self.s_last),
        ]
        return m
