import importlib
import json
from pathlib import Path
import sys

import pytest

from genomatch import vtable_utils
from utils import read_tsv, write_object_metadata, write_lines


def write_object(path, rows, *, include_count):
    target = {"genome_build": "GRCh37", "contig_naming": "ncbi"}
    if path.suffix == ".vmap":
        metadata = {"object_type": "variant_map", "target": target}
    else:
        metadata = {"object_type": "variant_table", **target}
    if include_count:
        metadata["variants_count"] = len(rows)
    write_lines(path, rows)
    write_object_metadata(path.with_name(path.name + ".meta.json"), metadata)


# "declared" has no counts and repeats reach the chunked reader; "membership-first"
# makes the repeated input the driver. restrict_vmap has a separate source-first branch.
TOOLS = ["intersect_variants", "union_variants", "restrict_vmap"]
REPEATED_KEY_CASES = [
    *[(tool, "vmap", plan) for tool in TOOLS for plan in ("declared", "membership-first")],
    ("restrict_vmap", "vmap", "source-first"),
    *[(tool, "vtable", "membership-first") for tool in TOOLS],
]


@pytest.mark.parametrize("tool,object_type,scan_plan", REPEATED_KEY_CASES)
def test_repeated_membership_keys_preserve_set_results(tmp_path, monkeypatch, tool, object_type, scan_plan):
    module = importlib.import_module(f"genomatch.{tool}")
    monkeypatch.setattr(vtable_utils, "WRITE_CHUNK_ROWS", 2)
    source = tmp_path / ("source.vmap" if tool == "restrict_vmap" else f"source.{object_type}")
    membership = tmp_path / f"membership.{object_type}"
    output = tmp_path / ("output.vmap" if tool == "restrict_vmap" else "output.vtable")
    source_rows = [f"1\t{pos}\tsrc{pos}\tAC\tA" for pos in [600, 100, 500, 200, 400, 300]]
    if source.suffix == ".vmap":
        source_rows = [f"{row}\tsource\t{idx}\tswap" for idx, row in enumerate(source_rows)]
    # First two occurrences share a chunk; later copies cross chunk boundaries.
    # The absent key exercises repeats that never match any driver row.
    positions = [100, 100, 700, 200, 100]
    if scan_plan == "source-first":
        positions += [700, 200]
    membership_rows = [f"1\t{pos}\tmember{idx}\tAC\tA" for idx, pos in enumerate(positions)]
    if object_type == "vmap":
        membership_rows = [f"{row}\tignored\t{idx}\tflip" for idx, row in enumerate(membership_rows)]
    for path, rows in [(source, source_rows), (membership, membership_rows)]:
        write_object(path, rows, include_count=scan_plan != "declared")

    monkeypatch.setattr(sys, "argv", [tool, str(source), str(membership), "--output", str(output)])
    assert module.main() == 0

    if tool == "restrict_vmap":
        expected = [source_rows[1].split("\t"), source_rows[3].split("\t")]
    else:
        positions = range(100, 701, 100) if tool == "union_variants" else [100, 200]
        expected = [["1", str(pos), f"1:{pos}:A:AC", "AC", "A"] for pos in positions]
    assert read_tsv(output) == expected
    metadata = json.loads(output.with_name(output.name + ".meta.json").read_text())
    assert metadata["variants_count"] == len(expected)
    assert not output.with_name(output.name + ".qc.tsv").exists()


@pytest.mark.parametrize("reader", ["whole", "chunked"])
@pytest.mark.parametrize(
    "row,error",
    [
        ("1\t100\tbad\tAC\tA\t.\t2\tinvalid", "invalid allele_op"),
        ("1\t0\tbad\tAC\tA\t.\t2\tidentity", "invalid pos"),
    ],
)
def test_membership_readers_still_validate_rows(tmp_path, reader, row, error):
    from genomatch.exact_set_utils import read_target_table, restrict_frame_to_membership_chunks

    path = tmp_path / "membership.vmap"
    write_lines(path, ["1\t100\tfirst\tAC\tA\t.\t0\tidentity", "1\t100\tsecond\tAC\tA\t.\t1\tidentity", row])
    with pytest.raises(ValueError, match=error):
        if reader == "whole":
            read_target_table(path)
        else:
            import pandas as pd

            driver = pd.DataFrame([["1", "100", "AC", "A"]], columns=["chrom", "pos", "a1", "a2"])
            restrict_frame_to_membership_chunks(driver, path)


@pytest.mark.parametrize("builds", [("unknown", "unknown"), ("unknown", "GRCh37")], ids=["all-unknown", "mixed"])
def test_exact_matching_rejects_unknown_builds(builds):
    from genomatch.exact_set_utils import TargetObjectInfo, require_shared_target_metadata

    infos = [
        TargetObjectInfo(
            path=Path(f"input{idx}.vmap"),
            object_type="variant_map",
            target_metadata={"genome_build": build, "contig_naming": "ncbi"},
            raw_metadata={},
            variants_count=None,
        )
        for idx, build in enumerate(builds)
    ]
    with pytest.raises(ValueError, match="known genome_build; run guess_build.py --write"):
        require_shared_target_metadata(infos)
