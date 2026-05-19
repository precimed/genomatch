# Mapping

This file defines exact set-operation semantics. Target-side row-transform semantics are defined in [variant-transforms.md](variant-transforms.md).

## Exact set semantics

- Set membership is by exact `chrom:pos:a1:a2`
- rsID is not used for exact set operations
- Prepared target objects are allele-canonical. Filtering and intersection are exact set operations, not harmonization steps.
- Exact set operations use the biallelic target-row allele model defined in [core-objects.md](core-objects.md#invariants).
- For multi-base alleles, `flip` means strand reverse-complement of each allele string: complement each base (`A <-> T`, `C <-> G`) and reverse the resulting string
- For multi-base alleles, `flip_swap` means apply that same reverse-complement operation and then swap `a1` and `a2`
- Build mismatch is an error; no implicit liftover is performed
- `guess_build.py` is the only build-guessing entrypoint
- Exact set tools may use object metadata `variants_count` to choose memory-efficient scan order, but this is an implementation detail and must not change documented output order, ID policy, provenance, or emitted object type.
- If any input to `intersect_variants.py` or `union_variants.py` lacks `variants_count`, the tool must warn and fall back to user-declared input order for its internal scan plan.
- `restrict_vmap.py` may also fall back to user-declared scan order when `variants_count` is absent, but it must not warn solely for missing `variants_count` because the field is only planning metadata for restriction.
- `.vmap` provenance is `source_shard + source_index`, where `source_index` is a non-negative shard-local row index and `source_shard` is stored exactly as emitted by the provenance-bearing source step
- Declared coordinate order is defined in [core-objects.md](core-objects.md) and is reused by both `sort_variants.py` and `liftover_build.py`

## No target-order remapping

- There is no canonical operation that projects a source `.vmap` into another object's row order.
- `restrict_vmap.py` filters a source `.vmap` by exact membership in all supplied restriction inputs; output IDs, row order, provenance, and `allele_op` remain those of the source `.vmap`.
- Users who need a particular `.vmap` order should produce the source `.vmap` in that order, for example with `sort_variants.py` for declared coordinate order.

## Restricting `.vmap`

- `restrict_vmap.py` restricts one source `.vmap` by exact `chrom:pos:a1:a2` membership in all supplied `.vtable` / `.vmap` restriction inputs
- CLI shape is `restrict_vmap.py <source.vmap> <restriction> [<restriction> ...] --output <output.vmap>`
- `restrict_vmap.py` requires the same `genome_build` and the same `contig_naming` across all inputs
- `restrict_vmap.py` performs no implicit normalization
- mismatched build or contig naming fails clearly
- the source input must be `.vmap`
- restriction inputs may be `.vtable` or `.vmap`
- a source row is retained only if its exact `chrom:pos:a1:a2` key appears in every restriction input; multi-restriction semantics are intersection, not union
- restriction input `.vmap` provenance is ignored; only target-side rows are used as membership filters
- output IDs, row order, provenance, and `allele_op` come from the source `.vmap`
- output is mapped-only `.vmap`
- newly written output metadata must include `variants_count`
- `restrict_vmap.py` does not compute, infer, or change `allele_op`
- `restrict_vmap.py` does not provide a sort mode; users who need sorted output should sort the source `.vmap` before restriction
- for memory efficiency, `restrict_vmap.py` may intersect smaller restriction inputs before scanning the source `.vmap`; if the source `.vmap` is not the smallest object, the implementation may continue reducing candidate membership until it reaches the source, then introduce the source `.vmap` IDs, order, provenance, and `allele_op`, and continue filtering against remaining restriction inputs
- if the source `.vmap` is the smallest object, the implementation may read or stream the source first and then stream restriction inputs as membership filters while preserving source `.vmap` output semantics

## Assigning `.vmap` IDs

- `assign_vmap_ids.py` assigns IDs to one source `.vmap` from an ID source `.vmap` / `.vtable`
- CLI shape is `assign_vmap_ids.py --vmap <source.vmap> --id-source <source.vmap|source.vtable> --output <output.vmap>`
- optional `--unmatched-id-policy {drop,variant-key,missing}` controls source `.vmap` rows with no exact `--id-source` match; default is `drop`
- `assign_vmap_ids.py` requires the same `genome_build` and the same `contig_naming` across inputs
- `assign_vmap_ids.py` performs no implicit normalization, allele flipping, liftover, or ID-based matching
- matching is by exact `chrom:pos:a1:a2`
- the input named by `--vmap` must be a `.vmap`
- the input named by `--id-source` may be `.vmap` or `.vtable`; when it is `.vmap`, only target-side rows are used and its provenance is ignored
- duplicate `chrom:pos:a1:a2` keys in `--id-source` may be ignored when that key is absent from the source `.vmap`
- if `--id-source` contains duplicate `chrom:pos:a1:a2` keys for any key present in the source `.vmap`, fail clearly
- retained output row order follows the source `.vmap`
- retained output provenance (`source_shard`, `source_index`) and `allele_op` come from the source `.vmap`
- retained output `chrom`, `pos`, `a1`, and `a2` come from the source `.vmap`
- retained output `id` is copied from the matching `--id-source` row
- with `--unmatched-id-policy drop`, rows with no exact `chrom:pos:a1:a2` match in `--id-source` are dropped and audited in `<output>.qc.tsv` with status `id_not_found`
- with `--unmatched-id-policy variant-key`, rows with no exact `chrom:pos:a1:a2` match in `--id-source` are retained with `id=chrom:pos:a1:a2` from the source `.vmap` row and are not audited
- with `--unmatched-id-policy missing`, rows with no exact `chrom:pos:a1:a2` match in `--id-source` are retained with `id=.` and are not audited
- rows whose matching `--id-source` ID is empty or `.` are dropped and audited in `<output>.qc.tsv` with status `missing_id`
- `<output>.qc.tsv` contains dropped rows only; rows are keyed by source `.vmap` provenance and contain `source_shard`, `source_index`, `source_id`, and `status`
- if no rows are dropped, no QC sidecar is emitted and any stale `<output>.qc.tsv` is removed
- output is `.vmap`
- newly written output metadata must include `variants_count`

## Intersections

- `intersect_variants.py` intersects exact `chrom:pos:a1:a2`
- CLI shape is `intersect_variants.py <input> <input> [<input> ...] --output <output.vtable>`
- `intersect_variants.py` is the symmetric provenance-free set-intersection tool
- `intersect_variants.py` requires the same `genome_build` and the same `contig_naming` across all inputs
- `intersect_variants.py` performs no implicit normalization
- mismatched build or contig naming fails clearly
- output IDs are target-derived `chrom:pos:a1:a2`, not copied from any particular input row
- `intersect_variants.py` sorts the intersection by declared coordinate order using the same sorting contract as `sort_variants.py`, then emits one row per exact target key
- `intersect_variants.py` always emits `.vtable`
- input `.vmap` provenance and `allele_op` are ignored and dropped
- this is an accepted breaking change from older first-input ID/order behavior, chosen to keep provenance-free set materialization consistent with `union_variants.py`
- `intersect_variants.py` does not provide a sort mode because sorted target-derived output is the only output contract
- for memory efficiency, `intersect_variants.py` may start candidate-key intersection from the smallest input; when `variants_count` is absent, fallback to user-declared input order affects only internal candidate aggregation, not observable output order or IDs

## Unions

- `union_variants.py` unions exact `chrom:pos:a1:a2`
- CLI shape is `union_variants.py <input> <input> [<input> ...] --output <output.vtable>`
- `union_variants.py` requires at least two inputs
- `union_variants.py` requires the same `genome_build` and the same `contig_naming` across all inputs
- `union_variants.py` performs no implicit normalization
- mismatched build or contig naming fails clearly
- `union_variants.py` computes the provenance-free set union of exact target keys across all inputs
- output IDs are target-derived `chrom:pos:a1:a2`, not copied from any particular input row
- `union_variants.py` sorts the union by declared coordinate order using the same sorting contract as `sort_variants.py`, then emits one row per exact target key
- `union_variants.py` emits `.vtable`
- this is an accepted breaking change from older first-retained-occurrence ID behavior, chosen so provenance-free set materialization is source-agnostic
- for memory efficiency, `union_variants.py` may use smaller inputs first for internal candidate aggregation; when `variants_count` is absent, fallback to user-declared input order affects only internal candidate aggregation, not observable output order or IDs
