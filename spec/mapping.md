# Mapping

This file defines exact matching and set-operation semantics. Target-side row-transform semantics are defined in [variant-transforms.md](variant-transforms.md).

## Matching semantics

- Matching is by `chrom:pos:a1:a2`
- rsID does not drive matching
- Target row order defines output row order
- Same-build exact matching emits only `identity`, `swap` or `missing`
- Same-build exact matching must support both SNP and non-SNP biallelic alleles; `swap` is defined by allele ordering, not by allele length
- For multi-base alleles, `flip` means strand reverse-complement of each allele string: complement each base (`A <-> T`, `C <-> G`) and reverse the resulting string
- For multi-base alleles, `flip_swap` means apply that same reverse-complement operation and then swap `a1` and `a2`
- Missing target rows remain in `.vmap` with `source_shard=.` `source_index=-1` and `allele_op=missing`
- Build mismatch is an error; no implicit liftover is performed
- `guess_build.py` is the only build-guessing entrypoint
- `.vmap` provenance is `source_shard + source_index`, where `source_index` is shard-local and `source_shard` is stored exactly as emitted by the provenance-bearing source step
- Declared coordinate order is defined in [core-objects.md](core-objects.md) and is reused by both `sort_variants.py` and `liftover_build.py`

## Mapping

- `match_vmap_to_target.py` emits only `identity`, `swap`, or `missing`
- `match_vmap_to_target.py` requires the same `genome_build` and the same `contig_naming` across all inputs; wrong-strand resolution is not part of this tool
- `match_vmap_to_target.py` requires the source input to be a `.vmap`
- `match_vmap_to_target.py` accepts the target input as a `.vtable` or a `.vmap`
- matching operates only on the target side of the source `.vmap`
- if the target input is a `.vmap`, matching also operates only on its target side and ignores its provenance entirely
- if the target input is a `.vmap`, the tool must emit a warning that target provenance is ignored
- users who want a provenance-free target and no warning should materialize it first with `convert_vmap_to_target.py`
- the source provenance carried by the input `.vmap` is irrelevant to the matching logic itself but must still be preserved in the emitted `.vmap`
- output order follows the target input order, whether that order comes from a `.vtable` or from the target side of a `.vmap`
- `match_vmap_to_target.py` emits `.vmap`

## Intersections

- `intersect_variants.py` intersects exact `chrom:pos:a1:a2`
- `intersect_variants.py` requires the same `genome_build` and the same `contig_naming` across all inputs
- `intersect_variants.py` performs no implicit normalization
- mismatched build or contig naming fails clearly
- output IDs come from the first input
- output order follows the first input and is not re-sorted by declared coordinate order
- `intersect_variants.py` does not provide a sort mode; users who need sorted intersection output should sort the desired first input before intersection
- `intersect_variants.py` emits `.vtable`

## Unions

- `union_variants.py` unions exact `chrom:pos:a1:a2`
- `union_variants.py` requires at least two inputs
- `union_variants.py` requires the same `genome_build` and the same `contig_naming` across all inputs
- `union_variants.py` performs no implicit normalization
- mismatched build or contig naming fails clearly
- duplicate exact rows are deduplicated by first occurrence across the full input stream, scanning inputs in CLI order and rows in file order
- output IDs come from that first retained occurrence
- after deduplication, `union_variants.py` emits declared coordinate order using the same sorting contract as `sort_variants.py`
- stable ordering for ties is inherited from that declared-coordinate sort, so rows with the same declared contig and numeric position retain first-occurrence order
- `union_variants.py` emits `.vtable`
