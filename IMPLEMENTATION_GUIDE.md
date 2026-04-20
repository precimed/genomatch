# Implementation Guide

This document is implementation guidance only. Normative constraints are defined in [spec/performance-contract.md](spec/performance-contract.md).

## Performance guidance (non-normative)

### Goals

- Preserve existing behavior contracts while improving scalability.
- Reduce Python-level per-row loops in compute-heavy paths.

### Preferred implementation patterns

1. Keep orchestration and I/O separate from columnar compute kernels.
2. Prefer table-native transformations over side arrays.
3. For expensive per-label transforms, compute on unique values and remap.
4. Use mixed-path execution when data classes differ materially (for example SNV vs indel): split by row properties, run the best kernel per subset, then recombine before shared downstream checks/output.
5. Use mask-based status assignment with explicit precedence order.
6. Preserve deterministic output ordering (stable sorts where applicable).
7. Define required columns at stage boundaries and validate early.
8. Drop temporary helper columns at stage exits unless contractually required.
9. Prefer `pd.read_table` / `pd.read_csv` for tabular inputs when schema is known.
10. Promote reusable table-native helpers (for example sort, duplicate mask, allele-op composition) into shared utilities (`vtable_utils.py`) when used by multiple tools.

### Vectorization intent

- Eliminating `to_rows()` / `from_rows()` calls is not a goal by itself.
- The goal is full table-native/vectorized execution where appropriate for scale.
- Keep non-vectorized paths only when explicitly justified (for example indel-heavy or non-autosome subsets with expected small cardinality), and document the rationale with a local `PERF:` comment.

### First-failing-row pattern

When behavior requires "first invalid row" semantics:

1. Build an invalid boolean mask.
2. Check mask presence (`any()`).
3. Resolve first failing index deterministically.
4. Raise message consistent with existing contract/tests.

### When loops are often acceptable

- File-emission loops (TSV/VCF writing).
- Small-cardinality control loops (for example unique contig labels).
- External-tool orchestration where loop removal does not materially change runtime.

### dtype discipline

- Readers and loaders own dtype canonicalization. They produce `dtype="object"` for all string columns.
- SV/CV kernels must not call `.astype("string")` (pandas StringDtype) on columns that originated from a loader. That dtype was deliberately chosen at load time; escalating it in a kernel silently re-introduces the overhead the loader avoided.
- The only valid place to use `dtype="string"` or `astype("string")` is at a new trust boundary — e.g. reading subprocess VCF output that has not yet passed through a loader.
- Use `astype(str)` (Python `str`, not pandas `"string"`) when you need a plain string representation for display, logging, or passing to a Python function that expects `str`.
- If a `.str.*` accessor fails on a column without a StringDtype cast, that is a signal the column's dtype was set incorrectly at load time. Fix the loader, not the kernel.

### Reusable utility ideas

- required-column checks
- unique-value map/remap helpers
- strict integer parsing to series + validity mask
- first-true-index helper
- batched allele-op composition
