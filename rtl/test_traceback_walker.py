"""Simulation checks for the external-BRAM traceback walker."""

from amaranth.sim import Simulator

from traceback_walker import TracebackWalker


def test_walks_one_pointer_per_cycle_after_priming():
    dut = TracebackWalker(row_width=4, band_width=8)
    memory = [TracebackWalker.STOP] * (1 << 7)
    # Endpoint (3, 4) -> DIAG -> (2, 3) -> UP -> (1, 3) -> LEFT -> (1, 2) -> STOP.
    memory[(3 << 3) | 4] = TracebackWalker.DIAG
    memory[(2 << 3) | 3] = TracebackWalker.UP
    memory[(1 << 3) | 3] = TracebackWalker.LEFT

    seen = []

    async def bench(ctx):
        ctx.set(dut.start_row, 3)
        ctx.set(dut.start_column, 4)
        ctx.set(dut.start, 1)
        await ctx.tick()
        ctx.set(dut.start, 0)

        # Feed a synchronous RAM response: the value presented during a clock
        # is the pointer at the address requested by the walker.
        while ctx.get(dut.active):
            ctx.set(dut.pointer, memory[ctx.get(dut.read_addr)])
            await ctx.tick()
            if ctx.get(dut.step_valid):
                seen.append((ctx.get(dut.step_row), ctx.get(dut.step_column), ctx.get(dut.move)))

        assert ctx.get(dut.done)
        assert not ctx.get(dut.failed)

    sim = Simulator(dut)
    sim.add_clock(1e-6)
    sim.add_testbench(bench)
    sim.run()

    assert seen == [
        (3, 4, TracebackWalker.DIAG),
        (2, 3, TracebackWalker.UP),
        (1, 3, TracebackWalker.LEFT),
    ]


def test_invalid_pointer_fails_closed():
    dut = TracebackWalker(row_width=4, band_width=8)

    async def bench(ctx):
        ctx.set(dut.start_row, 2)
        ctx.set(dut.start_column, 2)
        ctx.set(dut.start, 1)
        await ctx.tick()
        ctx.set(dut.start, 0)
        ctx.set(dut.pointer, 15)
        await ctx.tick()  # BRAM priming
        await ctx.tick()  # consume invalid pointer
        assert ctx.get(dut.done)
        assert ctx.get(dut.failed)

    sim = Simulator(dut)
    sim.add_clock(1e-6)
    sim.add_testbench(bench)
    sim.run()
