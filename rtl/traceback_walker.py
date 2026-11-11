"""One-pointer-per-cycle traceback walker for the banded SW pointer RAM.

Pointer RAM is addressed as ``row * 128 + column``. The RAM is deliberately
external to this module: a BRAM implementation presents ``pointer`` one clock
after ``read_addr`` is registered. After that one-cycle priming delay, this
walker consumes one 4-bit pointer per clock.
"""

from amaranth import Elaboratable, Module, Signal


class TracebackWalker(Elaboratable):
    """Walk 4-bit traceback pointers from a best-score endpoint to STOP.

    Pointer encoding is intentionally part of the FPGA/host ABI:

    * 0: STOP
    * 1: DIAG -- consume one read and one reference base
    * 2: UP   -- consume one read base (insertion relative to reference)
    * 3: LEFT -- consume one reference base (deletion relative to reference)

    Encodings 4--15 fail closed. The affine DP cell can keep its internal
    state however it needs, but it must write one of these final path moves for
    each traceback cell.
    """

    STOP = 0
    DIAG = 1
    UP = 2
    LEFT = 3

    def __init__(self, *, row_width=12, band_width=128):
        if band_width < 2 or band_width & (band_width - 1):
            raise ValueError("band_width must be a power of two of at least 2")
        self.row_width = row_width
        self.band_width = band_width
        self.column_width = (band_width - 1).bit_length()

        # Control and endpoint inputs.
        self.start = Signal()
        self.start_row = Signal(row_width)
        self.start_column = Signal(self.column_width)

        # Synchronous BRAM read interface. ``pointer`` is the response for the
        # address issued in the preceding cycle.
        self.read_addr = Signal(row_width + self.column_width)
        self.pointer = Signal(4)

        # A pulse for every non-STOP pointer consumed.
        self.step_valid = Signal()
        self.step_row = Signal(row_width)
        self.step_column = Signal(self.column_width)
        self.move = Signal(2)

        self.active = Signal()
        self.done = Signal()
        self.failed = Signal()

    def elaborate(self, platform):
        m = Module()

        row = Signal(self.row_width)
        column = Signal(self.column_width)
        priming = Signal()

        # Pulse outputs default low; done and failed remain asserted until a
        # subsequent start so software cannot miss a terminal result.
        m.d.sync += self.step_valid.eq(0)

        with m.If(self.start & ~self.active):
            m.d.sync += [
                row.eq(self.start_row),
                column.eq(self.start_column),
                self.read_addr.eq((self.start_row << self.column_width) | self.start_column),
                priming.eq(1),
                self.active.eq(1),
                self.done.eq(0),
                self.failed.eq(0),
            ]
        with m.Elif(self.active):
            # The RAM response for the start address becomes valid after this
            # first active cycle. Every later cycle both consumes a response
            # and issues the address for the following response.
            with m.If(priming):
                m.d.sync += priming.eq(0)
            with m.Elif(self.pointer == self.STOP):
                m.d.sync += [self.active.eq(0), self.done.eq(1)]
            with m.Elif(self.pointer == self.DIAG):
                with m.If((row == 0) | (column == 0)):
                    m.d.sync += [self.active.eq(0), self.done.eq(1), self.failed.eq(1)]
                with m.Else():
                    m.d.sync += [
                        self.step_valid.eq(1),
                        self.step_row.eq(row),
                        self.step_column.eq(column),
                        self.move.eq(self.DIAG),
                        row.eq(row - 1),
                        column.eq(column - 1),
                        self.read_addr.eq(((row - 1) << self.column_width) | (column - 1)),
                    ]
            with m.Elif(self.pointer == self.UP):
                with m.If(row == 0):
                    m.d.sync += [self.active.eq(0), self.done.eq(1), self.failed.eq(1)]
                with m.Else():
                    m.d.sync += [
                        self.step_valid.eq(1),
                        self.step_row.eq(row),
                        self.step_column.eq(column),
                        self.move.eq(self.UP),
                        row.eq(row - 1),
                        self.read_addr.eq(((row - 1) << self.column_width) | column),
                    ]
            with m.Elif(self.pointer == self.LEFT):
                with m.If(column == 0):
                    m.d.sync += [self.active.eq(0), self.done.eq(1), self.failed.eq(1)]
                with m.Else():
                    m.d.sync += [
                        self.step_valid.eq(1),
                        self.step_row.eq(row),
                        self.step_column.eq(column),
                        self.move.eq(self.LEFT),
                        column.eq(column - 1),
                        self.read_addr.eq((row << self.column_width) | (column - 1)),
                    ]
            with m.Else():
                m.d.sync += [self.active.eq(0), self.done.eq(1), self.failed.eq(1)]

        return m
