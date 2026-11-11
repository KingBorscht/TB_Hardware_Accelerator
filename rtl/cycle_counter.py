"""Cycle counter core for hardware-only timing in the PL."""

from amaranth import Elaboratable, Module, Signal


class CycleCounter(Elaboratable):
    """Count clocks from ``start`` until ``stop`` and latch the result.

    An AXI-Lite wrapper can expose ``latched_cycles`` to the PYNQ host. The
    wrapper should drive start/stop from the accelerator/DMA transaction, not
    from Python timing, when reporting hardware-only latency.
    """

    def __init__(self, *, width=32):
        self.start = Signal()
        self.stop = Signal()
        self.running = Signal()
        self.latched_cycles = Signal(width)

    def elaborate(self, platform):
        m = Module()
        count = Signal(len(self.latched_cycles))

        with m.If(self.start):
            m.d.sync += [count.eq(0), self.running.eq(1)]
        with m.Elif(self.running & self.stop):
            m.d.sync += [self.latched_cycles.eq(count), self.running.eq(0)]
        with m.Elif(self.running):
            m.d.sync += count.eq(count + 1)

        return m
