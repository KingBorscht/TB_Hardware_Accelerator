# Traceback memory and walker architecture

This is the pre-synthesis BRAM plan for SCOPE v2. The project budget is 140
RAMB36 blocks. Vivado's first synthesis must replace every estimate below.

## BRAM organization

The 128-PE array writes 128 four-bit pointers each cycle, or 512 bits/cycle.
Pack this as eight 72-bit lanes: sixteen pointers per lane. A 72-bit RAMB36
has 512 words of depth. A guarded 3,200-cycle store therefore needs seven
depth banks per lane.

```text
ceil(512 / 72)    = 8 parallel lanes
ceil(3,200 / 512) = 7 depth banks
8 * 7             = 56 RAMB36 blocks per full traceback buffer
```

At exactly 3,072 rows, six banks fit and cost 48 blocks. The guarded 3,200
cycle allocation is the planning case. Raw capacity understates the cost:
`3,072 * 128 * 4 = 1,572,864` bits is about 43 RAMB36 blocks, but it cannot
meet the required 512-bit write width. Two guarded full pointer buffers cost
112 blocks before all other BRAM users.

```text
cycle_index = row
lane        = PE_index / 16
nibble      = PE_index % 16
```

The forward array owns the write port; the walker reads the older read through
the other port. This avoids read/write collision.

## First-synthesis decision rule

Record post-synthesis counts for these quantities:

```text
B_trace_full : one 512-bit/cycle, 3,200-cycle pointer buffer
B_reference  : reference windows and AXI-Lite reload storage
B_read_pp    : both packed read ping-pong buffers
B_dma         : DMA FIFOs, descriptors, and width converters
B_margin      : 14 RAMB36 blocks reserved for implementation variance
```

```text
if 2*B_trace_full + B_reference + B_read_pp + B_dma + B_margin <= 140:
    use two full pointer buffers
elif B_trace_full + B_reference + B_read_pp + B_dma + B_margin <= 140:
    use one full buffer only for functional bring-up;
    use checkpoint-and-recompute for throughput
else:
    switch immediately to checkpoint-and-recompute
```

The preliminary calculation makes double full-buffering unlikely: `2 * 56 =
112`, leaving 28 blocks before margin. References alone use about `32 * 2,500
* 2 = 160 kbit`, or roughly five RAMB36 blocks by capacity. Do not claim that
double full-pointer buffering fits until the synthesis report proves it.

## Walker and verification

`rtl/traceback_walker.py` consumes the pointer RAM with this ABI:

| Value | Meaning | CIGAR operation |
|---|---|---|
| 0 | stop | end |
| 1 | diagonal | `M` |
| 2 | up | `I` |
| 3 | left | `D` |
| 4--15 | invalid | fail read |

BRAM causes one priming cycle. The walker then consumes and issues one pointer
address per clock, so a 3,100-pointer path takes about 3,101 cycles plus
endpoint setup (about 31 microseconds at 100 MHz). It fails closed on illegal
pointers or moves above/left of matrix zero.

`rtl/test_traceback_walker.py` simulates a diagonal/up/left path and an invalid
pointer. Next, generate pointer RAM from the v2 affine-gap specification model
and compare the emitted CIGAR/move stream bit-for-bit for randomized cases.

## Reference and read-buffer ownership

Reference bases are packed in two bits and reloaded through AXI-Lite only when
no affected read is active. A packed 3,072-base read occupies 768 bytes. Start
with 64-read, 48 KiB payload batches in two ping-pong DMA allocations.

```text
PS fills A -> flush -> DMA/forward owns A -> walker owns pointer bank A
PS fills B -> flush -> DMA/forward owns B -> walker owns pointer bank B
walker completes -> invalidate/read result -> bank becomes PS-owned
```

To prove the walker does not stall the forward pass, measure and preserve:

```text
worst_case_walker_cycles + bank_handoff_cycles <= minimum_forward_pass_cycles
```

The current 3,100-cycle walker and 3,200-cycle forward targets leave little
headroom. If this fails, checkpoint/recompute or queueing is mandatory.
