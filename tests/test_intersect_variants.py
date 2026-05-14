import json

import pandas as pd

import genomatch.intersect_variants as intersect_module
from genomatch.exact_set_utils import TargetObjectInfo
from genomatch.tabular_rows import VariantRowsTable
from utils import read_tsv, run_py, write_json, write_lines


def test_intersect_variants_assigns_target_ids_and_sorts(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    meta = {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"}
    write_lines(first, ["1\t200\tid2\tC\tT", "1\t100\tid1\tA\tG", "1\t300\tid3\tA\tC"])
    write_lines(second, ["1\t100\tx\tA\tG", "1\t200\ty\tC\tT"])
    write_json(first.with_name(first.name + ".meta.json"), meta)
    write_json(second.with_name(second.name + ".meta.json"), meta)
    result = run_py("intersect_variants.py", first, second, "--output", out)
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "100", "1:100:A:G", "A", "G"], ["1", "200", "1:200:C:T", "C", "T"]]
    assert f"loaded 3 variants from {first}" in result.stderr
    assert f"loaded 2 variants from {second}" in result.stderr
    assert f"after intersecting {second}, 2 variants remain" in result.stderr
    assert "variants_count missing from metadata" in result.stderr
    assert "using CLI input order" in result.stderr


def test_intersect_variants_uses_count_order_and_emits_sorted_target_ids(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    write_lines(first, ["1\t300\tlast\tA\tC", "1\t100\tkeep\tA\tG", "1\t200\tdrop\tC\tT"])
    write_lines(second, ["1\t100\tother\tA\tG"])
    write_json(
        first.with_name(first.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi", "variants_count": 3},
    )
    write_json(
        second.with_name(second.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi", "variants_count": 1},
    )

    result = run_py("intersect_variants.py", first, second, "--output", out)

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "100", "1:100:A:G", "A", "G"]]
    assert result.stderr.index(f"loaded 1 variants from {second}") < result.stderr.index(f"loaded 3 variants from {first}")
    assert "variants_count missing from metadata" not in result.stderr


def test_intersect_count_order_does_not_reread_first_input_when_it_is_driver(tmp_path, monkeypatch):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    first_frame = pd.DataFrame([["1", "100", "keep", "A", "G"]], columns=["chrom", "pos", "id", "a1", "a2"])
    second_frame = pd.DataFrame(
        [["1", "100", "other", "A", "G"], ["1", "200", "other2", "C", "T"]],
        columns=["chrom", "pos", "id", "a1", "a2"],
    )
    frames = {first: first_frame, second: second_frame}
    read_counts = {first: 0, second: 0}
    infos = [
        TargetObjectInfo(
            path=first,
            object_type="variant_table",
            target_metadata={"genome_build": "GRCh37", "contig_naming": "ncbi"},
            raw_metadata={
                "object_type": "variant_table",
                "genome_build": "GRCh37",
                "contig_naming": "ncbi",
                "variants_count": 1,
            },
            variants_count=1,
        ),
        TargetObjectInfo(
            path=second,
            object_type="variant_table",
            target_metadata={"genome_build": "GRCh37", "contig_naming": "ncbi"},
            raw_metadata={
                "object_type": "variant_table",
                "genome_build": "GRCh37",
                "contig_naming": "ncbi",
                "variants_count": 2,
            },
            variants_count=2,
        ),
    ]

    def fake_read_target_table(path):
        read_counts[path] += 1
        return VariantRowsTable.from_frame(frames[path].copy(), copy=False)

    def fake_restrict_frame_to_membership_chunks(frame, path):
        return frame, len(frames[path])

    monkeypatch.setattr(intersect_module, "read_target_table", fake_read_target_table)
    monkeypatch.setattr(intersect_module, "restrict_frame_to_membership_chunks", fake_restrict_frame_to_membership_chunks)

    result = intersect_module.intersect_count_order([first, second], infos)

    assert read_counts[first] == 1
    assert result.reset_index(drop=True).equals(first_frame)


def test_intersect_variants_streams_later_input_by_full_input_membership(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    meta = {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"}
    write_lines(first, ["1\t1\tfirst\tA\tG", "1\t100001\tlast\tC\tT"])
    with open(second, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("1\t1\tx\tA\tG\n")
        for idx in range(2, 100001):
            handle.write(f"1\t{idx}\tfiller{idx}\tA\tC\n")
        handle.write("1\t100001\ty\tC\tT\n")
    write_json(first.with_name(first.name + ".meta.json"), meta)
    write_json(second.with_name(second.name + ".meta.json"), meta)

    result = run_py("intersect_variants.py", first, second, "--output", out)

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "1", "1:1:A:G", "A", "G"], ["1", "100001", "1:100001:C:T", "C", "T"]]
    assert f"loaded 100001 variants from {second}" in result.stderr


def test_intersect_variants_checks_metadata_before_loading_rows(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    write_lines(first, ["1\t100\tid1\tA\tG"])
    write_lines(second, ["not\ta\tvalid\tvtable\trow"])
    write_json(first.with_name(first.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"})
    write_json(second.with_name(second.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh38", "contig_naming": "ncbi"})

    result = run_py("intersect_variants.py", first, second, "--output", out)

    assert result.returncode != 0
    assert "same genome_build" in result.stderr
    assert f"{first}: genome_build=GRCh37, contig_naming=ncbi" in result.stderr
    assert f"{second}: genome_build=GRCh38, contig_naming=ncbi" in result.stderr
    assert "invalid vtable row" not in result.stderr


def test_intersect_variants_validates_later_vtable_like_first_input(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    meta = {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"}
    write_lines(first, ["1\t100\tid1\tA\tG"])
    write_lines(second, ["1\t100\t\tA\tG"])
    write_json(first.with_name(first.name + ".meta.json"), meta)
    write_json(second.with_name(second.name + ".meta.json"), meta)

    result = run_py("intersect_variants.py", first, second, "--output", out)

    assert result.returncode != 0
    assert "row is missing id" in result.stderr


def test_intersect_variants_does_not_canonicalize_contig_spellings(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    meta = {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"}
    write_lines(first, ["chr1\t100\tid1\tA\tG"])
    write_lines(second, ["1\t100\tid2\tA\tG"])
    write_json(first.with_name(first.name + ".meta.json"), meta)
    write_json(second.with_name(second.name + ".meta.json"), meta)

    result = run_py("intersect_variants.py", first, second, "--output", out)

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == []


def test_intersect_variants_does_not_accept_sort_flag(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    meta = {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"}
    write_lines(first, ["1\t100\tid1\tA\tG"])
    write_lines(second, ["1\t100\tid2\tA\tG"])
    write_json(first.with_name(first.name + ".meta.json"), meta)
    write_json(second.with_name(second.name + ".meta.json"), meta)

    result = run_py("intersect_variants.py", first, second, "--output", out, "--sort")

    assert result.returncode != 0
    assert "unrecognized arguments: --sort" in result.stderr


def test_intersect_variants_checks_contig_naming_before_loading_rows(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    write_lines(first, ["1\t100\tid1\tA\tG"])
    write_lines(second, ["not\ta\tvalid\tvtable\trow"])
    write_json(first.with_name(first.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"})
    write_json(second.with_name(second.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ucsc"})

    result = run_py("intersect_variants.py", first, second, "--output", out)

    assert result.returncode != 0
    assert "same contig_naming" in result.stderr
    assert f"{first}: genome_build=GRCh37, contig_naming=ncbi" in result.stderr
    assert f"{second}: genome_build=GRCh37, contig_naming=ucsc" in result.stderr
    assert "invalid vtable row" not in result.stderr


def test_intersect_variants_fails_cleanly_on_missing_contig_naming(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    write_lines(first, ["1\t100\tid1\tA\tG"])
    write_lines(second, ["not\ta\tvalid\tvtable\trow"])
    write_json(first.with_name(first.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37"})
    write_json(second.with_name(second.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"})

    result = run_py("intersect_variants.py", first, second, "--output", out)

    assert result.returncode != 0
    assert "contig_naming" in result.stderr
    assert "normalize_contigs.py" in result.stderr
    assert f"{first}: genome_build=GRCh37, contig_naming=<missing>" in result.stderr
    assert "invalid vtable row" not in result.stderr


def test_intersect_variants_empty_later_input_yields_empty_result(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    write_lines(first, ["1\t100\tid1\tA\tG"])
    write_lines(second, [])
    write_json(first.with_name(first.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"})
    write_json(second.with_name(second.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"})

    result = run_py("intersect_variants.py", first, second, "--output", out)

    assert result.returncode == 0, result.stderr
    assert out.exists()
    from utils import read_tsv
    assert read_tsv(out) == []


def test_intersect_variants_rejects_over_wide_later_vtable_row(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    write_lines(first, ["1\t100\tid1\tA\tG"])
    write_lines(second, ["1\t100\tid1\tA\tG\textra_column"])
    write_json(first.with_name(first.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"})
    write_json(second.with_name(second.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"})

    result = run_py("intersect_variants.py", first, second, "--output", out)

    assert result.returncode != 0
    assert "invalid variant table row" in result.stderr


def test_intersect_variants_accepts_vmap_and_drops_provenance(tmp_path):
    first = tmp_path / "a.vmap"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    write_lines(first, ["1\t100\tfrom_vmap\tA\tG\t.\t0\tidentity", "1\t200\tdrop\tC\tT\t.\t1\tswap"])
    write_lines(second, ["1\t100\tfrom_vtable\tA\tG"])
    write_json(
        first.with_name(first.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )
    write_json(second.with_name(second.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"})

    result = run_py("intersect_variants.py", first, second, "--output", out)

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "100", "1:100:A:G", "A", "G"]]


def test_intersect_variants_ignores_ids_and_vmap_provenance_in_membership(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vmap"
    out = tmp_path / "out.vtable"
    write_lines(
        first,
        [
            "1\t100\tfirst_id\tA\tG",
            "1\t200\tfirst_only\tC\tT",
        ],
    )
    write_lines(second, ["1\t100\tdifferent_id\tA\tG\tignored_shard\t99\tflip_swap"])
    write_json(first.with_name(first.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"})
    write_json(
        second.with_name(second.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )

    result = run_py("intersect_variants.py", first, second, "--output", out)

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "100", "1:100:A:G", "A", "G"]]


def test_intersect_variants_multiple_inputs_are_intersection_not_union(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    third = tmp_path / "c.vtable"
    out = tmp_path / "out.vtable"
    meta = {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"}
    write_lines(
        first,
        [
            "1\t100\tboth\tA\tG",
            "1\t200\tsecond_only\tC\tT",
            "1\t300\tthird_only\tG\tT",
        ],
    )
    write_lines(second, ["1\t100\tx\tA\tG", "1\t200\ty\tC\tT"])
    write_lines(third, ["1\t100\tx\tA\tG", "1\t300\tz\tG\tT"])
    for path in (first, second, third):
        write_json(path.with_name(path.name + ".meta.json"), meta)

    result = run_py("intersect_variants.py", first, second, third, "--output", out)

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "100", "1:100:A:G", "A", "G"]]


def test_intersect_variants_supports_non_snp_biallelic_exact_intersection(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    meta = {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"}
    write_lines(first, ["1\t100\tindel\tAC\tA", "1\t200\tsnp\tG\tT"])
    write_lines(second, ["1\t100\tother\tAC\tA"])
    write_json(first.with_name(first.name + ".meta.json"), meta)
    write_json(second.with_name(second.name + ".meta.json"), meta)

    result = run_py("intersect_variants.py", first, second, "--output", out)

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "100", "1:100:AC:A", "AC", "A"]]


def test_intersect_variants_first_vmap_output_metadata_is_vtable(tmp_path):
    first = tmp_path / "a.vmap"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    write_lines(first, ["1\t100\tfrom_vmap\tA\tG\t.\t0\tidentity"])
    write_lines(second, ["1\t100\tfrom_vtable\tA\tG"])
    write_json(
        first.with_name(first.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )
    write_json(second.with_name(second.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"})

    result = run_py("intersect_variants.py", first, second, "--output", out)

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "100", "1:100:A:G", "A", "G"]]
    metadata = json.loads(out.with_name(out.name + ".meta.json").read_text(encoding="utf-8"))
    assert metadata == {
        "object_type": "variant_table",
        "genome_build": "GRCh37",
        "contig_naming": "ncbi",
        "variants_count": 1,
    }


def test_intersect_variants_rejects_duplicate_first_vmap_target_keys(tmp_path):
    first = tmp_path / "a.vmap"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    write_lines(first, ["1\t100\tfirst\tA\tG\t.\t0\tidentity", "1\t100\tsecond\tA\tG\t.\t1\tidentity"])
    write_lines(second, ["1\t100\tid1\tA\tG"])
    write_json(
        first.with_name(first.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )
    write_json(second.with_name(second.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"})

    result = run_py("intersect_variants.py", first, second, "--output", out)

    assert result.returncode != 0
    assert "duplicate chrom:pos:a1:a2 in vmap target rows" in result.stderr


def test_intersect_variants_rejects_duplicate_vmap_target_keys(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vmap"
    out = tmp_path / "out.vtable"
    write_lines(first, ["1\t100\tid1\tA\tG"])
    write_lines(second, ["1\t100\tfirst\tA\tG\t.\t0\tidentity", "1\t100\tsecond\tA\tG\t.\t1\tidentity"])
    write_json(first.with_name(first.name + ".meta.json"), {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"})
    write_json(
        second.with_name(second.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )

    result = run_py("intersect_variants.py", first, second, "--output", out)

    assert result.returncode != 0
    assert "duplicate chrom:pos:a1:a2 in vmap target rows" in result.stderr


def test_intersect_variants_rejects_template_paths(tmp_path):
    first = tmp_path / "a.vtable"
    second = tmp_path / "b.vtable"
    out = tmp_path / "out.vtable"
    meta = {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"}
    write_lines(first, ["1\t100\tid1\tA\tG"])
    write_lines(second, ["1\t100\tid2\tA\tG"])
    write_json(first.with_name(first.name + ".meta.json"), meta)
    write_json(second.with_name(second.name + ".meta.json"), meta)

    result = run_py("intersect_variants.py", first, "b.@.vtable", "--output", out)

    assert result.returncode != 0
    assert "does not accept '@' paths" in result.stderr
