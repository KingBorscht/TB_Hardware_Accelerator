# Host-side mutation-calling design

Status: design proposal for SCOPE v2.0 (2026-10-10).

## Purpose and boundary

The PYNQ-Z2 FPGA performs primer identification, banded affine-gap
Smith-Waterman alignment, and traceback. The host-side mutation caller turns
each completed alignment into evidence at *catalogued* resistance positions,
then combines evidence from many reads into an amplicon-level call.

This is not a clinical reporting system and does not discover novel mutations.
It reports only whether a versioned catalogue mutation is supported by the
aligned reads. Basecalling is upstream (Dorado FASTQ) and is not a host-caller
responsibility.

```
Dorado FASTQ -> FPGA primer match / alignment / traceback
            -> DMA alignment result + CIGAR
            -> C++ evidence extractor
            -> compact per-read evidence records
            -> Python consensus policy
            -> catalogue mutation calls
```

The FPGA does **not** emit finished sequencing data. It emits alignment
metadata. The host retains the normalized basecalled read so it can compare
read bases with H37Rv bases along the traceback path.

## Scope-v2 constraints incorporated here

| Scope decision | Host-software consequence |
|---|---|
| Gene-oriented reference windows | Catalogued coordinates are stored in gene orientation; metadata maps them back to H37Rv genomic coordinates. |
| First 25 bp primer match, <=3 mismatches | Each result carries one amplicon ID and normalized strand. Reverse reads are reverse-complemented before interpretation. |
| Banded affine-gap local alignment with traceback | The FPGA result must include an endpoint and a complete CIGAR/path, not just a score. |
| Six amplicons, parameterized to 32 | All metadata and grouping keys use `amplicon_id`; no locus is hard-coded. |
| Python majority vote; >=50 coverage; 70% resistant; 10-70% mixed | Python owns the final policy decision. C++ creates deterministic per-read evidence only. |
| Versioned WHO catalogue | Every output includes the reference/catalogue version used; catalogue updates must not change historical calls. |

## Division of implementation

SCOPE v2 explicitly locks the final majority vote to Python. The requested C++
component should therefore be a host library/CLI responsible for DMA result
adaptation, traceback validation, and per-read evidence extraction. It should
write line-delimited JSON (or a compact binary equivalent) consumed by a small
Python consensus script. Keeping threshold policy in Python makes the locked
70%/10% values auditable and easy to change without rebuilding PYNQ host code.

```
host/                         # C++17
  include/tb/                 # public types and FPGA adapter API
  src/                        # CIGAR walker, evidence extractor, DMA adapter
  tests/                      # host-only deterministic unit tests
  tools/extract_evidence.cpp  # FASTQ + FPGA-result stream -> evidence JSONL

reference/consensus.py        # Python: group evidence and apply locked policy
bio/catalogue/                # versioned catalogue and amplicon metadata
```

The DMA adapter is deliberately behind an interface so the evidence extractor
can be tested on a workstation from recorded FPGA-result fixtures.

## Required FPGA-to-host result contract

One scalar Smith-Waterman score is insufficient for mutation calling. For each
accepted read, the FPGA/driver boundary must provide:

```cpp
struct FpgaAlignmentResult {
    uint64_t read_id;
    uint8_t  amplicon_id;       // 0..31; assigned by primer match
    bool     reverse_input;     // original FASTQ orientation
    bool     aligned;
    uint16_t best_score;        // 13-bit unsigned score carried in 16 bits
    uint16_t query_begin;       // normalized/reoriented read coordinate
    uint16_t query_end;
    uint16_t reference_begin;   // gene-oriented amplicon-window coordinate
    uint16_t reference_end;
    std::vector<CigarOp> cigar; // complete traceback, not a band-local sketch
};

struct CigarOp {
    char opcode;                // 'M', 'I', or 'D'
    uint16_t length;
};
```

`M` means aligned bases and may be either an equality or mismatch; C++ compares
the retained read/reference bases to distinguish them. `I` consumes read bases;
`D` consumes reference bases. A traceback result with invalid lengths, an
unknown opcode, coordinates outside its selected window, or a path that does
not consume the stated endpoints is rejected as `TRACEBACK_INVALID` and never
contributes to coverage.

The eventual hardware ABI may encode CIGAR operations as 2-bit words plus run
lengths for DMA efficiency, but the logical contract above is the software
boundary and must remain stable.

## Metadata and catalogue contract

Each amplicon needs metadata rather than hard-coded gene assumptions:

```cpp
struct AmpliconMetadata {
    uint8_t amplicon_id;
    std::string name;                 // e.g. "rpoB_RRDR"
    std::string reference_gene_oriented;
    int32_t h37rv_genome_start;       // coordinate of reference[0]
    int8_t genome_step;               // +1 or -1 in gene orientation
    std::string reference_version;    // e.g. NC_000962.3
};

struct CatalogueVariant {
    std::string catalogue_id;
    uint8_t amplicon_id;
    uint16_t reference_offset;        // zero-based gene-oriented coordinate
    std::string reference;
    std::string alternate;
    std::string gene;
    std::string annotation;           // supplied catalogue label, not inferred
    std::string catalogue_version;
};
```

`genome_step = -1` handles genes such as `katG` whose gene orientation is
opposite the genome. The caller should retain both the catalogue's gene-oriented
offset and a derived 1-based H37Rv coordinate, avoiding accidental reversal of
codon numbering.

The catalogue is input data, versioned under `bio/catalogue/`; it is not
embedded in C++ source. Exact amplicon coordinates and the catalogue version
remain required Bio-team inputs before a real call can be produced.

## Per-read evidence extraction

For each catalogued position covered by a valid traceback:

1. Initialize query and reference cursors from `query_begin` and
   `reference_begin`.
2. Walk every CIGAR operation while validating cursor bounds.
3. For `M`, compare each normalized read base with the gene-oriented reference
   base. Emit the observed allele at any matching catalogue coordinate.
4. For `I` or `D`, emit an insertion/deletion allele only when its normalized
   anchoring rule matches a catalogue variant. Record other indels as
   `OTHER_ALLELE` for quality accounting, not as a new mutation discovery.
5. Emit one evidence record per catalogue site per read. A read with an `N` at
   that site is `AMBIGUOUS` and is not counted as reference or alternate.

The extractor output is intentionally small:

```json
{
  "read_id": 1842,
  "amplicon_id": 0,
  "catalogue_id": "WHO-RPOB-S450L",
  "observed_allele": "ALT",
  "h37rv_position": 761155,
  "alignment_score": 4102,
  "reference_version": "NC_000962.3",
  "catalogue_version": "<pinned-version>"
}
```

Allowed `observed_allele` values are `REF`, `ALT`, `OTHER_ALLELE`,
`AMBIGUOUS`, and `NOT_COVERED`. Only `REF` and the catalogue's exact `ALT`
are eligible for the simple allele-fraction numerator/denominator. The raw
counts of all states remain in the result for auditing.

## Consensus policy (Python)

Group evidence by `(amplicon_id, catalogue_id)` after deduplicating by
`read_id`. Let eligible coverage be `REF + ALT`; ambiguous, other, failed, and
not-covered reads are reported separately and do not inflate coverage.

| Condition | Output status |
|---|---|
| eligible coverage < 50 | `INSUFFICIENT_COVERAGE` |
| ALT / eligible coverage >= 0.70 | `RESISTANT_SUPPORTED` |
| 0.10 <= ALT / eligible coverage < 0.70 | `MIXED` |
| ALT / eligible coverage < 0.10 | `NOT_SUPPORTED` |

The target coverage of 100 is a run-quality goal, not a different calling
threshold. The caller reports allele fraction, all read-state counts, metadata
versions, and the threshold values used. It must not label a result as a
clinical diagnosis.

## Test plan before hardware integration

1. **CIGAR-unit tests:** exact match, mismatch, insertion, deletion, soft
   boundary cases, invalid operations, and paths that exceed endpoints.
2. **Orientation tests:** a forward and reverse FASTQ representation of the
   same amplicon yield identical gene-oriented evidence, including `katG`.
3. **Catalogue projection tests:** SNP and indel records map to the intended
   H37Rv coordinate and catalogue ID.
4. **Consensus tests:** 49/50 coverage boundary; 9.9%, 10%, 69.9%, and 70%
   allele-fraction boundaries; duplicate read IDs; ambiguous bases.
5. **Hardware contract tests:** recorded DMA/traceback fixtures reproduce the
   evidence produced by the software specification model bit-for-bit.

## Decisions still needed before implementation

- Exact primer sequences, amplicon windows, and coordinate convention.
- The pinned WHO catalogue release and the representation/normalization rule
  for indels.
- Whether FPGA traceback emits raw pointer memory, a CIGAR stream, or a host
  DMA buffer generated by a traceback walker. This design requires the logical
  CIGAR/path contract, regardless of transport encoding.
- The host/Fabric ownership of primer matching. SCOPE v2 assigns it to the
  FPGA; a software fallback may exist only for test fixtures and must be
  labeled as such.
