from utils import read_tsv, run_py, write_json, write_lines
import json


def test_convert_vmap_to_target_materializes_target_side(tmp_path):
    source = tmp_path / "map.vmap"
    out = tmp_path / "out.vtable"
    write_lines(source, ["1\t100\tt1\tA\tG\t1\t0\tidentity", "1\t200\tt2\tC\tT\t.\t-1\tmissing"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {
            "object_type": "variant_map",
            "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"},
        },
    )
    result = run_py("convert_vmap_to_target.py", "--source", source, "--output", out)
    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "100", "t1", "A", "G"], ["1", "200", "t2", "C", "T"]]
    meta = json.loads(out.with_name(out.name + ".meta.json").read_text(encoding="utf-8"))
    assert meta == {
        "object_type": "variant_table",
        "genome_build": "GRCh37",
        "contig_naming": "ncbi",
        "created_by": "convert_vmap_to_target.py",
        "derived_from": str(source),
    }
