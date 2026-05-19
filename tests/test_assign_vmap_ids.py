import json
from pathlib import Path

import pytest

from genomatch import vtable_utils
from genomatch.assign_vmap_ids import assign_ids_from_id_source
from genomatch.exact_set_utils import load_target_object_info
from genomatch.vtable_utils import read_vmap_table
from utils import read_tsv, run_py, write_json, write_lines


def write_vmap(path: Path, lines: list[str], *, genome_build: str = "GRCh38", contig_naming: str = "ncbi") -> None:
    write_lines(path, lines)
    write_json(
        path.with_name(path.name + ".meta.json"),
        {
            "object_type": "variant_map",
            "target": {"genome_build": genome_build, "contig_naming": contig_naming},
            "variants_count": len(lines),
            "created_by": "test",
        },
    )


def write_vtable(path: Path, lines: list[str], *, genome_build: str = "GRCh38", contig_naming: str = "ncbi") -> None:
    write_lines(path, lines)
    write_json(
        path.with_name(path.name + ".meta.json"),
        {
            "object_type": "variant_table",
            "genome_build": genome_build,
            "contig_naming": contig_naming,
            "variants_count": len(lines),
        },
    )


def test_assign_vmap_ids_copies_ids_and_drops_unassigned_rows_to_qc(tmp_path):
    source = tmp_path / "source.vmap"
    id_source = tmp_path / "ids.vtable"
    out = tmp_path / "out.vmap"
    write_vmap(
        source,
        [
            "1\t10\told1\tA\tG\t.\t0\tidentity",
            "1\t20\told2\tC\tT\t.\t1\tswap",
            "1\t25\told25\tT\tC\t.\t2\tidentity",
            "1\t30\told3\tG\tA\t.\t3\tidentity",
        ],
    )
    write_vtable(
        id_source,
        [
            "1\t10\trs10\tA\tG",
            "1\t25\t.\tT\tC",
        ],
    )

    result = run_py("assign_vmap_ids.py", "--vmap", source, "--id-source", id_source, "--output", out)

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "10", "rs10", "A", "G", ".", "0", "identity"]]
    assert read_tsv(out.with_name(out.name + ".qc.tsv")) == [
        ["source_shard", "source_index", "source_id", "status"],
        [".", "1", "old2", "id_not_found"],
        [".", "2", "old25", "missing_id"],
        [".", "3", "old3", "id_not_found"],
    ]
    metadata = json.loads(out.with_name(out.name + ".meta.json").read_text(encoding="utf-8"))
    assert metadata["object_type"] == "variant_map"
    assert metadata["target"] == {"genome_build": "GRCh38", "contig_naming": "ncbi"}
    assert metadata["variants_count"] == 1
    assert metadata["created_by"] == "test"


def test_assign_vmap_ids_variant_key_policy_retains_unmatched_rows_without_qc(tmp_path):
    source = tmp_path / "source.vmap"
    id_source = tmp_path / "ids.vtable"
    out = tmp_path / "out.vmap"
    write_vmap(
        source,
        [
            "1\t10\told1\tA\tG\t.\t0\tidentity",
            "1\t20\told2\tC\tT\t.\t1\tidentity",
            "1\t25\told25\tT\tC\t.\t2\tidentity",
        ],
    )
    write_vtable(
        id_source,
        [
            "1\t10\trs10\tA\tG",
            "1\t25\t.\tT\tC",
        ],
    )

    result = run_py(
        "assign_vmap_ids.py",
        "--vmap",
        source,
        "--id-source",
        id_source,
        "--output",
        out,
        "--unmatched-id-policy",
        "variant-key",
    )

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["1", "10", "rs10", "A", "G", ".", "0", "identity"],
        ["1", "20", "1:20:C:T", "C", "T", ".", "1", "identity"],
    ]
    assert read_tsv(out.with_name(out.name + ".qc.tsv")) == [
        ["source_shard", "source_index", "source_id", "status"],
        [".", "2", "old25", "missing_id"],
    ]


def test_assign_vmap_ids_missing_policy_retains_unmatched_rows_with_missing_id_without_qc(tmp_path):
    source = tmp_path / "source.vmap"
    id_source = tmp_path / "ids.vtable"
    out = tmp_path / "out.vmap"
    write_vmap(
        source,
        [
            "1\t10\told1\tA\tG\t.\t0\tidentity",
            "1\t20\told2\tC\tT\t.\t1\tidentity",
        ],
    )
    write_vtable(id_source, ["1\t10\trs10\tA\tG"])

    result = run_py(
        "assign_vmap_ids.py",
        "--vmap",
        source,
        "--id-source",
        id_source,
        "--output",
        out,
        "--unmatched-id-policy",
        "missing",
    )

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["1", "10", "rs10", "A", "G", ".", "0", "identity"],
        ["1", "20", ".", "C", "T", ".", "1", "identity"],
    ]
    assert not out.with_name(out.name + ".qc.tsv").exists()


def test_assign_vmap_ids_accepts_vmap_id_source_and_ignores_its_provenance(tmp_path):
    source = tmp_path / "source.vmap"
    id_source = tmp_path / "ids.vmap"
    out = tmp_path / "out.vmap"
    write_vmap(
        source,
        [
            "1\t20\told2\tC\tT\tpayload\t5\tswap",
            "1\t10\told1\tA\tG\tpayload\t4\tidentity",
        ],
    )
    write_vmap(
        id_source,
        [
            "1\t10\trs10\tA\tG\tlookup\t100\tidentity",
            "1\t20\trs20\tC\tT\tlookup\t101\tidentity",
        ],
    )

    result = run_py("assign_vmap_ids.py", "--vmap", source, "--id-source", id_source, "--output", out)

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [
        ["1", "20", "rs20", "C", "T", "payload", "5", "swap"],
        ["1", "10", "rs10", "A", "G", "payload", "4", "identity"],
    ]
    assert not out.with_name(out.name + ".qc.tsv").exists()


def test_assign_vmap_ids_removes_stale_qc_when_no_rows_are_dropped(tmp_path):
    source = tmp_path / "source.vmap"
    id_source = tmp_path / "ids.vtable"
    out = tmp_path / "out.vmap"
    write_vmap(source, ["1\t10\told1\tA\tG\t.\t0\tidentity"])
    write_vtable(id_source, ["1\t10\trs10\tA\tG"])
    write_lines(out.with_name(out.name + ".qc.tsv"), ["stale"])

    result = run_py("assign_vmap_ids.py", "--vmap", source, "--id-source", id_source, "--output", out)

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "10", "rs10", "A", "G", ".", "0", "identity"]]
    assert not out.with_name(out.name + ".qc.tsv").exists()


def test_assign_vmap_ids_rejects_duplicate_used_id_source_target_keys(tmp_path):
    source = tmp_path / "source.vmap"
    id_source = tmp_path / "ids.vtable"
    out = tmp_path / "out.vmap"
    write_vmap(source, ["1\t10\told1\tA\tG\t.\t0\tidentity"])
    write_vtable(
        id_source,
        [
            "1\t10\trs10a\tA\tG",
            "1\t10\trs10b\tA\tG",
        ],
    )

    result = run_py("assign_vmap_ids.py", "--vmap", source, "--id-source", id_source, "--output", out)

    assert result.returncode != 0
    assert "duplicate used chrom:pos:a1:a2 key: 1:10:A:G" in result.stderr
    assert not out.exists()


def test_assign_vmap_ids_rejects_duplicate_used_key_across_id_source_chunks(tmp_path, monkeypatch):
    source = tmp_path / "source.vmap"
    id_source = tmp_path / "ids.vtable"
    write_vmap(source, ["1\t10\told1\tA\tG\t.\t0\tidentity"])
    write_vtable(
        id_source,
        [
            "1\t10\trs10a\tA\tG",
            "1\t10\trs10b\tA\tG",
        ],
    )
    monkeypatch.setattr(vtable_utils, "WRITE_CHUNK_ROWS", 1)
    source_frame = read_vmap_table(source).to_frame(copy=False)
    id_source_info = load_target_object_info(id_source)

    with pytest.raises(ValueError, match="duplicate used chrom:pos:a1:a2 key: 1:10:A:G"):
        assign_ids_from_id_source(
            source_frame,
            id_source,
            id_source_info,
            unmatched_id_policy="drop",
        )


def test_assign_vmap_ids_ignores_duplicate_unused_id_source_target_keys(tmp_path):
    source = tmp_path / "source.vmap"
    id_source = tmp_path / "ids.vtable"
    out = tmp_path / "out.vmap"
    write_vmap(source, ["1\t10\told1\tA\tG\t.\t0\tidentity"])
    write_vtable(
        id_source,
        [
            "1\t20\trs20a\tC\tT",
            "1\t20\trs20b\tC\tT",
            "1\t10\trs10\tA\tG",
        ],
    )

    result = run_py("assign_vmap_ids.py", "--vmap", source, "--id-source", id_source, "--output", out)

    assert result.returncode == 0, result.stderr
    assert read_tsv(out) == [["1", "10", "rs10", "A", "G", ".", "0", "identity"]]


def test_assign_vmap_ids_rejects_metadata_mismatch(tmp_path):
    source = tmp_path / "source.vmap"
    id_source = tmp_path / "ids.vtable"
    out = tmp_path / "out.vmap"
    write_vmap(source, ["1\t10\told1\tA\tG\t.\t0\tidentity"], genome_build="GRCh37")
    write_vtable(id_source, ["1\t10\trs10\tA\tG"], genome_build="GRCh38")

    result = run_py("assign_vmap_ids.py", "--vmap", source, "--id-source", id_source, "--output", out)

    assert result.returncode != 0
    assert "all inputs must have the same genome_build" in result.stderr


def test_assign_vmap_ids_rejects_contig_naming_mismatch(tmp_path):
    source = tmp_path / "source.vmap"
    id_source = tmp_path / "ids.vtable"
    out = tmp_path / "out.vmap"
    write_vmap(source, ["1\t10\told1\tA\tG\t.\t0\tidentity"], contig_naming="ncbi")
    write_vtable(id_source, ["1\t10\trs10\tA\tG"], contig_naming="ucsc")

    result = run_py("assign_vmap_ids.py", "--vmap", source, "--id-source", id_source, "--output", out)

    assert result.returncode != 0
    assert "all inputs must have the same contig_naming" in result.stderr


def test_assign_vmap_ids_rejects_non_vmap_source(tmp_path):
    source = tmp_path / "source.vtable"
    id_source = tmp_path / "ids.vtable"
    out = tmp_path / "out.vmap"
    write_vtable(source, ["1\t10\told1\tA\tG"])
    write_vtable(id_source, ["1\t10\trs10\tA\tG"])

    result = run_py("assign_vmap_ids.py", "--vmap", source, "--id-source", id_source, "--output", out)

    assert result.returncode != 0
    assert "requires --vmap to be a .vmap" in result.stderr
