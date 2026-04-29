# Core objects

## `.vtable`

Tab-separated, exactly 5 columns:

1. `chrom`
2. `pos`
3. `id`
4. `a1`
5. `a2`

Row order is meaningful and must be preserved.

`a1` and `a2` are ordered target alleles, not an unordered allele set.

For canonical reference-aware target rows produced by `restrict_build_compatible.py` or `liftover_build.py`, the convention is:

- `a1` = non-reference allele
- `a2` = reference allele

Importer-emitted rows may use source-format-specific allele ordering before reference-aware canonicalization.

Required sidecar metadata:

```json
{
  "object_type": "variant_table",
  "genome_build": "GRCh37",
  "contig_naming": "ncbi"
}
```

`.vtable` is the provenance-free materialized object in v1. It is produced primarily by:

- `convert_vmap_to_target.py`
- `intersect_variants.py`
- `union_variants.py`

It may also be used as the target input to matching operations where source provenance is irrelevant to the matching logic itself.

Canonical `.vtable` and `.vmap` artifacts are always single files. `@` templates are not part of the public object contract for canonical variant objects.

## `.vmap`

Tab-separated, exactly 8 columns:

1. `chrom`
2. `pos`
3. `id`
4. `a1`
5. `a2`
6. `source_shard`
7. `source_index`
8. `allele_op`

The first five columns are the exact target side.

`a1` and `a2` are ordered target alleles, not an unordered allele set.

For canonical reference-aware target rows produced by `restrict_build_compatible.py` or `liftover_build.py`, the convention is:

- `a1` = non-reference allele
- `a2` = reference allele

Importer-emitted rows may use source-format-specific allele ordering before reference-aware canonicalization.

`source_shard` is a provenance locator, not a biological chromosome field:

- for `@`-sharded importer inputs, it stores the exact substitution string that replaced `@` in the discovered filename
- for single-file imports, it is `.`

`source_index` is zero-based within that `source_shard`.

Allowed `allele_op` values in v1:

- `identity`
- `swap`
- `flip`
- `flip_swap`
- `missing`

Required sidecar metadata:

```json
{
  "object_type": "variant_map",
  "target": {
    "genome_build": "GRCh37",
    "contig_naming": "ncbi"
  }
}
```

Importer-emitted `.vmap` is the normal starting artifact for provenance-preserving workflows. Its target side contains retained canonical imported rows; its source side carries exact raw source provenance as `(source_shard, source_index)`.

## Invariants

- Metadata is authoritative for `genome_build` and `contig_naming` once those fields are present.
- Build mismatch must fail clearly unless the user explicitly runs a liftover tool.
- Canonical target-side rows stored in `.vtable` and `.vmap` support biallelic variants with allele strings stored verbatim as non-empty sequences of `A` / `C` / `G` / `T` characters.
- Generic exact-match, swap, intersection, provenance, and apply operations must support both SNP and non-SNP biallelic alleles without rewriting them based on allele length alone.
- Build-guessing is explicit and lives only in `guess_build.py`.
- Reference-dependent tools resolve assets from `match/config.yaml` plus input metadata; there are no CLI reference overrides in v1.
- Current reference-aware operations use UCSC-style internal reference FASTA and top-level chain config entries after strict validation against declared metadata.
- FASTA-consuming tools use bulk reference access by default. The reference access mode is controlled by `MATCH_REFERENCE_ACCESS_MODE`, which accepts `BULK` and `LEGACY` case-insensitively (`BULK` is the default when unset).
- The legacy per-lookup FASTA fetch path is available only in `LEGACY` mode, primarily for tests and profiling.
- `liftover_build.py` is the only cross-build conversion entrypoint.
- `.vmap` target duplicates by `chrom:pos:a1:a2` are not allowed.
- Source duplicates are tolerated; first occurrence wins.
- If `source_index == -1`, then `source_shard` must be `.` and `allele_op` must be `missing`.
- If metadata declares `contig_naming`, every target row must be valid under that convention. This is an object invariant; downstream tools are not required to re-validate it at every boundary.
- `.vmap` metadata is target-side only.
- `.vmap` provenance is stored on disk in each row as `(source_shard, source_index)`.
- `source_shard` must be stored exactly as discovered or assigned by the producing tool. It must not be normalized, canonicalized, or reinterpreted biologically.
- payload-application tools consume that stored provenance exactly; for `apply_vmap_to_bfile.py` and `apply_vmap_to_pfile.py`, `source_shard` lookup is exact and is not chromosome interpretation
- only `import_*` tools originate new provenance in v1
- non-importer tools may preserve existing provenance on `.vmap` input, but they must not originate a new `.vmap` from a standalone `.vtable`
- When any object (`.vmap`, `.vtable`, or any payload format such as sumstats, pfile/bfile/vcf) is represented in-memory as a DataFrame, `source_index` is a semantic provenance value assigned once by the parser/loader and stored in the corresponding table object (for example `SumstatsTable.source_index`, `VMapRowsTable`, `VariantRowsTable`). Downstream code must propagate provenance by reading from that stored field. Deriving provenance from the pandas structural index (`.index`, `RangeIndex`, `.iloc` position) is prohibited even when values are numerically identical — structural index is implementation-internal and carries no provenance guarantee. Resetting the structural DataFrame index (`reset_index`) is safe because it does not affect the provenance data column.

The exact meaning of `source_shard` depends on the first provenance-bearing step:

- importer-emitted `.vmap`: exact raw shard token, or `.` for a single-file import
- `.vmap` produced from an existing `.vmap`: preserve the original upstream `source_shard`

## Row-order contract
- Declared coordinate order means declared contig order, then numeric position, preserving input order for ties with the same contig and position.
- Canonical importers preserve retained discovery order
- For sharded imports this means deterministic shard order plus within-shard row order.
- Single-input filter/repair/reference-compatible restriction tools preserve retained target-row order unless explicitly documented otherwise.
- match_vmap_to_target.py output order follows the target input order, whether the target is a `.vtable` or the target side of a `.vmap`.
- intersect_variants.py output order follows the first input.
- union_variants.py emits declared coordinate order after deduplicating by first occurrence across all inputs.
- sort_variants.py emits declared coordinate order.
- liftover_build.py emits declared coordinate order after liftover.
- Any tool that re-sorts must say so explicitly
