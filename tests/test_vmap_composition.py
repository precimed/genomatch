from utils import read_tsv, run_py, write_json, write_lines


def test_vmap_match_target_composition_preserves_missing_rows(tmp_path):
    source = tmp_path / "src.vmap"
    target = tmp_path / "tgt.vtable"
    out = tmp_path / "map.vmap"
    write_lines(source, ["1\t100\tmid1\tA\tG\t.\t-1\tmissing", "1\t200\tmid2\tT\tC\tshard1\t4\tidentity"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {
            "object_type": "variant_map",
            "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"},
        },
    )
    write_lines(target, ["1\t100\tfinal1\tA\tG", "1\t200\tfinal2\tT\tC"])
    write_json(
        target.with_name(target.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    result = run_py("match_vmap_to_target.py", "--source", source, "--target", target, "--output", out)
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["1", "100", "final1", "A", "G", ".", "-1", "missing"],
        ["1", "200", "final2", "T", "C", "shard1", "4", "identity"],
    ]
