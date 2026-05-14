# Concepts

## Core mental model

Conceptually, the toolkit treats raw input as three pieces that can be separated and recombined: the variant table itself (`.vtable`), the attached payload, and the provenance linking variants back to the raw source:

- `raw data = .vtable + payload + provenance`
- `.vtable` is the ordered variant table, one row per biallelic variant
- `payload` is the external format-specific data attached to those variants, such as summary-stat columns, PLINK 1 BED/BIM/FAM genotype matrices, or PLINK 2 PGEN/PVAR/PSAM genotype matrices; payloads can be split across multiple files (shards)
- `provenance` is the link back from variant to its raw payload, represented as `(source_shard, source_index)` plus `allele_op` to capture `a1 <-> a2` allele swaps and strand flips
- `.vmap = .vtable + provenance`; `.vmap` is therefore the object emitted by import operations, and used for most downstream transforms: it keeps transformed table of mapped variants together with enough provenance to rewrite the original payload later, but it does not contain the payload itself. A canonical `.vmap` is mapped-only; rows missing from a source payload are absent rather than represented by sentinel rows.

The toolkit provides a set of tools for transforming tables of variants: importing raw inputs, normalizing representation, restricting rows, liftover between builds, performing intersections and unions of variant tables. These transformations can carry provenance linking each transformed variant row back to the original payload. That is what `.vmap` stores: a transformed variant table together with provenance, but without the payload itself. We refer to that transformed variant table as the *target side* of `.vmap`, while its *source side* is the original payload that the stored provenance links back to. Target side of `.vmap` can be seen as a `.vtable` object that only stores the transformed variant table, without provenance. `apply_vmap_*` reconnects the transformed variant rows in `.vmap` to the original payload and writes a rewritten payload in that transformed row order.

- In the following diagram `.vmap` is shown on the right, and its provenance points to source payload on the left:
  ```text
                         +---- .vmap -----+
       payload           |   .vtable +    |
       (source) <--------|---provenance   |
                         +----------------+
  ```
- Transformations will update `.vmap` while still keeping its provenance pointing to the same source payload

Intersection of variants is performed on exact `chr:bp:a1:a2`, ignoring variant `id`. A `prepare_variants.py` pipeline combines operations that standardize `chr:bp:a1:a2`, making it into a unique variant identifier.
Intersection supports biallelic SNPs and non-SNPs, while reference-aware restriction and liftover are narrower and are specified in the spec docs.

For exact schema and edge-case rules, see [SPEC.md](../SPEC.md), [core-objects.md](../spec/core-objects.md), [contigs-and-metadata.md](../spec/contigs-and-metadata.md), [shard-discovery.md](../spec/shard-discovery.md), [importers.md](../spec/importers.md), [variant-transforms.md](../spec/variant-transforms.md), [mapping.md](../spec/mapping.md), [payload-application.md](../spec/payload-application.md), and [workflow.md](../spec/workflow.md). End-user installation is described in [install.md](install.md), and reference-aware assets/config in [downloads.md](downloads.md).

## Core objects and metadata

Each `.vtable` or `.vmap` object is a single file with JSON sidecars at `<object>.meta.json`.

| Object | Purpose | Required metadata |
| --- | --- | --- |
| `.vmap` | Mapped-only provenance-carrying object for import, cleanup, reference-compatible restriction, liftover, `restrict_vmap.py`, and application | `object_type=variant_map`, `target.genome_build`, `target.contig_naming`, `variants_count` |
| `.vtable` | Provenance-free table, produced mainly by `convert_vmap_to_target.py`, `intersect_variants.py`, and `union_variants.py` | `object_type=variant_table`, `genome_build`, `contig_naming`, `variants_count` |

Downstream transforms other than `normalize_contigs.py` require declared `contig_naming`. If an imported object is unresolved or mixed, clean it with `normalize_contigs.py` before continuing. See [core-objects.md](../spec/core-objects.md) for the on-disk column contracts.

Newly written `.vmap` and `.vtable` sidecars include top-level `variants_count` row-count metadata. Readers accept missing `variants_count` for backward compatibility. Set-operation tools may use it to choose memory-efficient scan order, but not to change documented output order, ID policy, provenance, or emitted object type.

`a1` and `a2` are ordered alleles, not an unordered allele set. Importers preserve source-format allele meaning at import time, while reference-aware tools such as `restrict_build_compatible.py` and `liftover_build.py` emit reference-aligned variants with `a1=non-reference` and `a2=reference`. This keeps canonical reference-aware rows aligned with ALT/effect/counting-style semantics while still making the reference allele explicit in `a2`.

## Using .vmap vs .vtable as input across tools

- `import_*` tools are the only tools that originate new provenance, so they emit `.vmap`.
- For tools that accept both `.vmap` and `.vtable`, the default rule is to act on target rows and emit the same artifact type as the input, unless the tool explicitly defines a provenance-free `.vtable` output.
- `intersect_variants.py` always emits provenance-free `.vtable`.
- `restrict_vmap.py` is the provenance-preserving exact restriction tool for `.vmap`; restriction inputs are membership filters only and multiple restriction inputs are combined by intersection.
- `apply_vmap_*` consumes `.vmap` as a mapping object, not as a table. The payload must still match the `.vmap` source provenance exactly.
