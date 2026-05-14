from __future__ import annotations

from typing import Dict, Iterable, Set


def build_needed_source_indices(vmap_rows: Iterable[object]) -> Dict[str, Set[int]]:
    needed: Dict[str, Set[int]] = {}
    for row in vmap_rows:
        source_index = int(row.source_index)
        needed.setdefault(str(row.source_shard), set()).add(source_index)
    return needed


def output_variant_id(vmap_row: object, *, retain_snp_id: bool) -> str:
    if retain_snp_id:
        return str(vmap_row.id)
    return f"{vmap_row.chrom}:{vmap_row.pos}:{vmap_row.a1}:{vmap_row.a2}"
