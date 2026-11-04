import importlib
import json
import sys

import pytest

from genomatch import vtable_utils
from utils import read_tsv, run_py, write_json, write_lines


def write_object(path, rows, *, include_count):
    target = {"genome_build": "GRCh37", "contig_naming": "ncbi"}
    if path.suffix == ".vmap":
        metadata = {"object_type": "variant_map", "target": target}
    else:
        metadata = {"object_type": "variant_table", **target}
    if include_count:
        metadata["variants_count"] = len(rows)
    write_lines(path, rows)
    write_json(path.with_name(path.name + ".meta.json"), metadata)


@pytest.mark.parametrize("tool", ["intersect_variants", "union_variants", "restrict_vmap"])
@pytest.mark.parametrize("object_type", ["vmap", "vtable"])
@pytest.mark.parametrize("scan_plan", ["declared", "source-first", "membership-first"])
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
        ("1\t100\tbad\tAC\tA\t.\t-1\tidentity", "source_index out of range"),
        ("1\t100\tbad\tAC\tA\t.\t2\tinvalid", "invalid allele_op"),
        ("1\t0\tbad\tAC\tA\t.\t2\tidentity", "invalid pos"),
        ("1\t100\tbad\tN\tA\t.\t2\tidentity", "allele"),
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


@pytest.mark.parametrize('tool', ['intersect_variants', 'union_variants', 'restrict_vmap', 'assign_vmap_ids'])
@pytest.mark.parametrize('membership_type', ['vmap', 'vtable'])
@pytest.mark.parametrize('builds', [('unknown', 'unknown'), ('unknown', 'GRCh37'), ('GRCh37', 'unknown')])
def test_exact_matching_rejects_unknown_builds_before_reading_rows(tmp_path, tool, membership_type, builds):
    source = tmp_path / 'source.vmap'
    membership = tmp_path / f'membership.{membership_type}'
    output = tmp_path / ('output.vmap' if tool in {'restrict_vmap', 'assign_vmap_ids'} else 'output.vtable')
    for path, build in zip((source, membership), builds):
        target = {'genome_build': build, 'contig_naming': 'ncbi'}
        metadata = ({'object_type': 'variant_map', 'target': target} if path.suffix == '.vmap'
                    else {'object_type': 'variant_table', **target})
        # Metadata rejection must happen before malformed payloads are parsed.
        write_lines(path, ['malformed'])
        write_json(path.with_name(path.name + '.meta.json'), metadata)
    args = (['--vmap', source, '--id-source', membership] if tool == 'assign_vmap_ids'
            else [source, membership])
    result = run_py(f'{tool}.py', *args, '--output', output)
    assert result.returncode != 0
    assert 'all inputs must have a known genome_build' in result.stderr
    assert 'guess_build.py' in result.stderr
    assert f'{source}: genome_build={builds[0]}' in result.stderr
    assert f'{membership}: genome_build={builds[1]}' in result.stderr
    assert not output.exists()
    assert not output.with_name(output.name + '.meta.json').exists()
