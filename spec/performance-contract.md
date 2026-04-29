# Performance Contract

This document defines normative implementation-level constraints for performance-sensitive code paths while preserving result semantics and public contracts.

## Scope

These rules apply to transformations and I/O paths that operate on variant tables/maps and are expected to scale with row count.

## Required Invariants

1. Result semantics are primary.
- Any performance refactor must preserve CLI behavior, file-format contracts, QC/status semantics, and error-message semantics required by tests/spec.

1. Ownership by stage is explicit.
- Input parser/loader owns parse normalization (`PN`) and parse validation (`PV`).
- Transform kernels own semantic validation (`SV`) required for transform correctness.
- Output writers own contract validation (`CV`) required for emitted object guarantees.
- Type conversion policy is part of parser/loader ownership:
  - `PN`: canonicalization conversions required for downstream semantics (for example string normalization before string-domain checks).
  - `PV`: parseability/type checks (for example numeric coercion checks, null/invalid detection after conversion attempts).

1. Redundant `PN`/`PV` is disallowed without boundary change.
- The same `PN`/`PV` operation must not be repeated downstream unless input origin/trust boundary changes or representation changes invalidate earlier guarantees.
- This includes repeated dtype conversion (for example repeated `astype("string")`) in the same lineage without boundary change.

1. Boundary checks remain strict.
- Parse boundaries must validate structural shape/required columns/types.
- Output boundaries must validate schema and object invariants before emission.

## Type Conversion Contract (PN/PV)

1. `PN` type conversion
- Perform conversions needed to establish canonical in-memory representation at parse boundaries.
- Examples: canonical string dtype before string normalization; canonical integer dtype after successful numeric parse.

1. `PV` conversion validation
- Validate conversion outcomes immediately at parse boundary.
- Examples: failed numeric coercion, null-introducing conversions, invalid tokens after canonicalization.

1. Downstream conversion rule
- Downstream stages must not repeat the same conversion unless:
  - data crossed a trust boundary, or
  - representation changed and invalidated prior guarantees.

1. Documentation rule
- Boundary-function comments must state conversion assumptions/ownership in `Assumes/Performs/Guarantees`.

1. Determinism is preserved.
- Ordering-sensitive outputs must use deterministic ordering; when sorting is required, stable sort behavior must be preserved.

## Allowed Redundancy

Redundant normalization/validation is allowed only if at least one holds:

1. Data crossed a trust boundary (external file, subprocess output, connector/plugin output).
1. Representation changed in a way that invalidates prior guarantees.
1. Boundary-local error semantics require a local check.

When retained in hot paths, rationale must be documented in code comments or PR notes.

## Evidence Requirements For Changes

For performance-sensitive changes, contributors must provide:

1. Stage contract statements for touched major functions, expressed as code comments colocated with the implementation:
- `Assumes`
- `Performs` (with `PN`/`PV`/`SV`/`CV` labels)
- `Guarantees`
- Scope: boundary functions only (parsers/loaders, major transform kernels, writers).

1. Justification for retained Python-native loops in hot paths.
- If a loop is retained in a performance-sensitive path, add a short code comment near the loop describing:
  - why the loop is retained
  - why vectorization/batching is not appropriate or not materially better
  - expected cardinality/complexity characteristics (for example small-cardinality control loop)

1. Test evidence for unchanged behavior semantics (targeted tests at minimum; broader suite when warranted).

1. Profiling evidence when a change is motivated by performance.

## Code Comment Conventions

Use concise comments so contracts stay close to code and update with refactors.

Suggested format for boundary functions:

- `Assumes: ...`
- `Performs: PN(...), PV(...), SV(...), CV(...)`
- `Guarantees: ...`

Suggested format for retained loops in hot paths:

- `PERF: loop retained because ...; vectorization not used because ...; expected scale ...`

## Normative Patterns

The following patterns are normative when applicable:

1. First-failing-row behavior in vectorized logic must preserve legacy semantics:
- compute invalid mask
- fail on first invalid row deterministically
- preserve expected error text semantics

1. Multi-status assignment must use explicit precedence:
- assign statuses in ordered precedence with guards that prevent later overwrite unless intended.

1. Stage schema contracts must be explicit:
- major stage boundaries must define required columns and validate them before use.

1. Temporary helper columns should not leak across stage boundaries unless part of the stage contract.

1. Variable-length string assignment must preserve object dtype:
- when assigning potentially long/variable-length string data into DataFrame columns, avoid implicit NumPy fixed-width unicode materialization (for example `<U...>` from `np.array(list_of_strings)`).
- use object-safe assignment paths (for example pandas object Series/arrays or `to_numpy(dtype=object)`) so one extreme string length cannot upcast the full column to huge fixed-width unicode buffers.

1. External sort is allowed for row-ordering tools:
- tools such as `sort_variants.py` may use multi-pass external sorting with bounded chunk size and bounded merge fan-in, provided stable declared-coordinate ordering and all file-format / metadata semantics are preserved.
- external-sort temporary files must be isolated from canonical output artifacts and finalized by temp-then-atomic-rename so failed runs do not leave partial canonical outputs.

1. System temporary directory use must remain small:
- `tempfile` / system temporary-directory use is allowed only for a few small helper objects, approximately 1 MB total per process.
- Large intermediate files must be placed under a user-controlled output or scratch prefix because HPC environments may have small or quota-limited `/tmp` partitions.

1. Streaming set operations are allowed when first-input semantics are preserved:
- `intersect_variants.py` may hold the first input in memory, stream later inputs, and intersect by per-input membership sets.
- streaming must preserve first-input output order, first-input IDs, metadata validation, and exact `chrom:pos:a1:a2` intersection semantics.

1. Streaming source-to-target matching is allowed when target-order semantics are preserved:
- `match_vmap_to_target.py` may hold target rows in memory and stream source `.vmap` rows.
- streaming must preserve target-row output order, target duplicate rejection, source provenance preservation, first-match-wins behavior, and swap-aware allele classification.
