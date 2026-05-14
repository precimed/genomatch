# Mapping

This file defines exact set-operation semantics. Target-side row-transform semantics are defined in [variant-transforms.md](variant-transforms.md).

## Exact set semantics

- Set membership is by exact `chrom:pos:a1:a2`
- rsID is not used for exact set operations
- Prepared target objects are allele-canonical. Filtering and intersection are exact set operations, not harmonization steps.
- Exact set operations must support both SNP and non-SNP biallelic alleles without rewriting them based on allele length alone.
- For multi-base alleles, `flip` means strand reverse-complement of each allele string: complement each base (`A <-> T`, `C <-> G`) and reverse the resulting string
- For multi-base alleles, `flip_swap` means apply that same reverse-complement operation and then swap `a1` and `a2`
- Build mismatch is an error; no implicit liftover is performed
- `guess_build.py` is the only build-guessing entrypoint
- Exact set tools may use object metadata `variants_count` to choose memory-efficient scan order, but this is an implementation detail and must not change documented output order, ID policy, provenance, or emitted object type.
- If any input to `intersect_variants.py`, `restrict_vmap.py`, or `union_variants.py` lacks `variants_count`, the tool must warn and fall back to user-declared input order for its internal scan plan.
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
- `restrict_vmap.py` does not compute, infer, or change `allele_op`
- `restrict_vmap.py` does not provide a sort mode; users who need sorted output should sort the source `.vmap` before restriction
- for memory efficiency, `restrict_vmap.py` may intersect smaller restriction inputs before scanning the source `.vmap`; if the source `.vmap` is not the smallest object, the implementation may continue reducing candidate membership until it reaches the source, then introduce the source `.vmap` IDs, order, provenance, and `allele_op`, and continue filtering against remaining restriction inputs

## Intersections

- `intersect_variants.py` intersects exact `chrom:pos:a1:a2`
- CLI shape is `intersect_variants.py <input> <input> [<input> ...] --output <output.vtable>`
- `intersect_variants.py` is the symmetric provenance-free set-intersection tool
- `intersect_variants.py` requires the same `genome_build` and the same `contig_naming` across all inputs
- `intersect_variants.py` performs no implicit normalization
- mismatched build or contig naming fails clearly
- output IDs come from the first input
- output order follows the first input and is not re-sorted by declared coordinate order
- `intersect_variants.py` always emits `.vtable`
- input `.vmap` provenance and `allele_op` are ignored and dropped
- `intersect_variants.py` does not provide a sort mode; users who need sorted intersection output should sort the desired first input before intersection
- for memory efficiency, `intersect_variants.py` may start candidate-key intersection from the smallest input, but output IDs and row order still come from the first CLI input

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
- for memory efficiency, `union_variants.py` may use smaller inputs first for internal candidate aggregation; output remains the same source-agnostic target-key union in declared coordinate order
