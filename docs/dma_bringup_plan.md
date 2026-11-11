# PYNQ-Z2 DMA bring-up plan

## Chosen first path

Use AXI DMA through a Zynq PS HP port for the initial loopback and array
integration. HP is chosen for batched bandwidth. It is not cache coherent, so
PYNQ cacheable allocations always use:

```text
CPU writes TX buffer -> TX.flush() -> DMA reads TX
DMA writes RX buffer -> RX.invalidate() -> CPU reads RX
```

ACP is not the first path. It is coherent but has different bandwidth and
integration trade-offs; measure it later as a separate experiment.

## Loopback before the array

`rtl/axis_passthrough.py` is a ready/valid AXI4-Stream pass-through. In Vivado
connect AXI DMA MM2S to this block and its output to AXI DMA S2MM, exposing the
DMA as `axi_dma_0` in the PYNQ overlay. This measures DDR, interconnect,
streaming, and DMA setup without involving PEs.

Run on the board:

```bash
python3 host/pynq/dma_loopback.py loopback.bit --dma axi_dma_0 --reads-per-batch 64 --batches 20 --completion poll
python3 host/pynq/dma_loopback.py loopback.bit --dma axi_dma_0 --reads-per-batch 64 --batches 20 --completion interrupt
```

The interrupt comparison requires AXI DMA interrupts wired to the PS and
supported by the loaded overlay. It is a measurement, not an assumed win.

## Stale data and zero-copy evidence

The loopback script fills RX with `0xA5`, flushes it, transfers a changing TX
pattern, invalidates RX, and checks every returned byte. This catches missing
cache maintenance. It uses only `pynq.allocate` buffers in the DMA path and
records physical addresses, contiguity, and AXI DMA simple-mode programmed
addresses when available.

Save those JSON records with the benchmark. Do not call the path zero-copy
until the programmed DMA addresses match the PYNQ allocations and profiling
shows no intermediate NumPy copy. Do not add `np.copy`, list conversion, or a
temporary NumPy array to the transfer path.

## Batching, ping-pong, and timing

A maximum packed read is 768 bytes, so one-read DMA transfers are overhead
dominated. The initial 64-read batch is about 48 KiB. Sweep batch size later.
Two independent TX/RX buffer pairs provide ping-pong ownership and prevent the
CPU from overwriting a DMA-owned bank.

End-to-end time comes from `time.perf_counter_ns()` around flush, submission,
completion, invalidation, and validation. Hardware-only time comes from a PL
cycle counter, never a Python clock. `rtl/cycle_counter.py` is the counter
core. Its planned AXI-Lite wrapper exposes latched cycles at `0x10` after
start/stop control at `0x00`; wire the counter to the actual DMA/accelerator
transaction before reporting hardware-only latency.

```text
hardware_latency_seconds = pl_cycles / fabric_clock_hz
```

Keep hardware-only latency, end-to-end latency, throughput, and payload size
as separate result columns.
