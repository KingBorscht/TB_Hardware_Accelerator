#!/usr/bin/env python3
"""PYNQ-Z2 HP-port DMA loopback benchmark.

Run this on the PYNQ board after a Vivado overlay connects an AXI DMA MM2S
stream to its S2MM stream through ``rtl/axis_passthrough.py``. It intentionally
does not require the PE array.
"""

import argparse
import asyncio
import json
import time


BYTES_PER_PACKED_READ = 768  # ceil(3,072 bases * 2 bits / 8)


def fill_pattern(buffer, byte_count, seed):
    """Write directly into a PYNQ allocation; no NumPy staging array."""
    for index in range(byte_count):
        buffer[index] = (seed + index) & 0xFF


def verify_pattern(buffer, byte_count, seed):
    for index in range(byte_count):
        expected = (seed + index) & 0xFF
        if int(buffer[index]) != expected:
            raise RuntimeError(f"loopback mismatch at byte {index}: got {int(buffer[index]):#04x}, expected {expected:#04x}")


def simple_dma_addresses(dma):
    """Capture simple-mode addresses as board-side zero-copy evidence."""
    mmio = getattr(dma, "mmio", None)
    if mmio is None:
        return None
    try:
        return {"mm2s_source_low": int(mmio.read(0x18)), "s2mm_destination_low": int(mmio.read(0x48))}
    except Exception:
        return None


async def wait_interrupt(dma):
    await asyncio.gather(dma.sendchannel.wait_async(), dma.recvchannel.wait_async())


def wait_for_dma(dma, completion):
    if completion == "poll":
        dma.sendchannel.wait()
        dma.recvchannel.wait()
    else:
        asyncio.run(wait_interrupt(dma))


def counter_start(counter):
    # Planned AXI-Lite wrapper contract: control bit 0 starts/clears.
    if counter is not None:
        counter.write(0x00, 0x1)


def counter_stop_and_read(counter):
    # Planned wrapper contract: control bit 1 stops/latches; 0x10 holds cycles.
    if counter is None:
        return None
    counter.write(0x00, 0x2)
    return int(counter.read(0x10))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bitstream", help="path to loopback overlay .bit file")
    parser.add_argument("--dma", default="axi_dma_0", help="overlay AXI DMA IP name")
    parser.add_argument("--counter", default=None, help="optional AXI-Lite cycle-counter IP name")
    parser.add_argument("--reads-per-batch", type=int, default=64)
    parser.add_argument("--batches", type=int, default=20)
    parser.add_argument("--completion", choices=("poll", "interrupt"), default="poll")
    parser.add_argument("--no-stale-check", action="store_true", help="skip returned-payload validation")
    args = parser.parse_args()
    if args.reads_per_batch < 1 or args.batches < 1:
        parser.error("reads-per-batch and batches must both be positive")

    try:
        import numpy as np
        from pynq import Overlay, allocate
    except ImportError as error:
        raise SystemExit("Run this script on a PYNQ image with pynq and numpy installed") from error

    overlay = Overlay(args.bitstream)
    dma = getattr(overlay, args.dma)
    counter = getattr(overlay, args.counter) if args.counter else None
    bytes_per_batch = args.reads_per_batch * BYTES_PER_PACKED_READ

    # Two source/destination pairs implement ping-pong bank ownership.
    buffers = [
        (allocate(shape=(bytes_per_batch,), dtype=np.uint8, cacheable=True),
         allocate(shape=(bytes_per_batch,), dtype=np.uint8, cacheable=True))
        for _ in range(2)
    ]

    samples = []
    for batch in range(args.batches):
        tx, rx = buffers[batch & 1]
        seed = 31 * batch + 7
        fill_pattern(tx, bytes_per_batch, seed)
        rx[:bytes_per_batch].fill(0xA5)  # Prime CPU cache with stale-looking data.
        tx.flush()                       # CPU -> non-coherent HP DMA.
        rx.flush()                       # Make stale-data validation deterministic.

        wall_start_ns = time.perf_counter_ns()
        counter_start(counter)
        dma.recvchannel.transfer(rx, nbytes=bytes_per_batch)
        dma.sendchannel.transfer(tx, nbytes=bytes_per_batch)
        programmed_addresses = simple_dma_addresses(dma)
        wait_for_dma(dma, args.completion)
        pl_cycles = counter_stop_and_read(counter)
        wall_end_ns = time.perf_counter_ns()

        rx.invalidate()                  # Non-coherent HP DMA -> CPU.
        if not args.no_stale_check:
            verify_pattern(rx, bytes_per_batch, seed)
        sample = {
            "batch": batch,
            "bytes": bytes_per_batch,
            "end_to_end_ns": wall_end_ns - wall_start_ns,
            "pl_cycles": pl_cycles,
            "tx_physical_address": int(tx.physical_address),
            "rx_physical_address": int(rx.physical_address),
            "tx_contiguous": bool(tx.flags.c_contiguous),
            "rx_contiguous": bool(rx.flags.c_contiguous),
            "dma_register_addresses": programmed_addresses,
            "stale_data_check": not args.no_stale_check,
        }
        samples.append(sample)
        print(json.dumps(sample, sort_keys=True))

    elapsed_ns = sum(item["end_to_end_ns"] for item in samples)
    total_bytes = args.batches * bytes_per_batch
    print(json.dumps({"summary": {
        "completion": args.completion,
        "batches": args.batches,
        "reads_per_batch": args.reads_per_batch,
        "total_bytes": total_bytes,
        "mean_end_to_end_ns": elapsed_ns / args.batches,
        "effective_mib_per_s": (total_bytes / (1024 * 1024)) / (elapsed_ns / 1e9),
        "note": "Do not claim zero-copy until physical-address evidence and payload checks are saved from this board.",
    }}, sort_keys=True))

    for tx, rx in buffers:
        tx.freebuffer()
        rx.freebuffer()


if __name__ == "__main__":
    main()
