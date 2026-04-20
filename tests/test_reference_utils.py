from pathlib import Path

import pytest

from genomatch import reference_utils
from utils import write_fasta, write_lines


def test_resolve_config_path_uses_match_config_override(tmp_path, monkeypatch):
    config = tmp_path / "custom-config.yaml"
    config.write_text("references: {}\nchain: {}\n", encoding="utf-8")

    monkeypatch.setenv("MATCH_CONFIG", str(config))
    monkeypatch.setattr(reference_utils, "DEFAULT_MATCH_CONFIG", tmp_path / "missing-default.yaml")

    assert reference_utils.resolve_config_path() == config


def test_resolve_config_path_uses_default_ref_config(monkeypatch, tmp_path):
    ref = tmp_path / "ref"
    ref.mkdir()
    config = ref / "config.yaml"
    config.write_text("references: {}\nchain: {}\n", encoding="utf-8")

    monkeypatch.delenv("MATCH_CONFIG", raising=False)
    monkeypatch.setattr(reference_utils, "DEFAULT_MATCH_CONFIG", config)

    assert reference_utils.resolve_config_path() == config


def test_resolve_config_path_requires_override_or_default(monkeypatch, tmp_path):
    monkeypatch.delenv("MATCH_CONFIG", raising=False)
    monkeypatch.setattr(reference_utils, "DEFAULT_MATCH_CONFIG", tmp_path / "missing-config.yaml")

    with pytest.raises(ValueError, match="MATCH_CONFIG is required unless /ref/config.yaml exists"):
        reference_utils.resolve_config_path()


def test_relative_paths_resolve_from_ref_config_directory(monkeypatch, tmp_path):
    ref = tmp_path / "ref"
    config = ref / "config.yaml"
    grch37 = ref / "ucsc" / "GRCh37" / "hg19.fa"
    grch38 = ref / "ucsc" / "GRCh38" / "hg38.fa"
    chain37to38 = ref / "chain" / "hg19ToHg38.over.chain.gz"
    chain38to37 = ref / "chain" / "hg38ToHg19.over.chain.gz"

    grch37.parent.mkdir(parents=True)
    grch38.parent.mkdir(parents=True)
    chain37to38.parent.mkdir(parents=True)

    write_fasta(grch37, {"chr1": "A"})
    write_fasta(grch38, {"chr1": "A"})
    write_lines(chain37to38, ["dummy"])
    write_lines(chain38to37, ["dummy"])
    config.write_text(
        "\n".join(
            [
                "references:",
                "  ucsc:",
                "    GRCh37:",
                "      fasta: ucsc/GRCh37/hg19.fa",
                "    GRCh38:",
                "      fasta: ucsc/GRCh38/hg38.fa",
                "chain:",
                "  hg19ToHg38: chain/hg19ToHg38.over.chain.gz",
                "  hg38ToHg19: chain/hg38ToHg19.over.chain.gz",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    monkeypatch.delenv("MATCH_CONFIG", raising=False)
    monkeypatch.setattr(reference_utils, "DEFAULT_MATCH_CONFIG", config)

    assert reference_utils.resolve_internal_reference_fasta("GRCh37") == grch37
    assert reference_utils.resolve_internal_reference_fasta("GRCh38") == grch38
    assert reference_utils.resolve_liftover_chain("GRCh37", "GRCh38") == chain37to38
    assert reference_utils.resolve_liftover_chain("GRCh38", "GRCh37") == chain38to37


def test_resolve_reference_access_mode_defaults_to_bulk(monkeypatch):
    monkeypatch.delenv("MATCH_REFERENCE_ACCESS_MODE", raising=False)
    assert reference_utils.resolve_reference_access_mode() == "BULK"


@pytest.mark.parametrize(
    "raw_value,expected",
    [
        ("bulk", "BULK"),
        ("BULK", "BULK"),
        ("legacy", "LEGACY"),
        ("LeGaCy", "LEGACY"),
    ],
)
def test_resolve_reference_access_mode_accepts_case_insensitive_values(monkeypatch, raw_value, expected):
    monkeypatch.setenv("MATCH_REFERENCE_ACCESS_MODE", raw_value)
    assert reference_utils.resolve_reference_access_mode() == expected


def test_resolve_reference_access_mode_rejects_unknown_value(monkeypatch):
    monkeypatch.setenv("MATCH_REFERENCE_ACCESS_MODE", "fast")
    with pytest.raises(ValueError, match="unsupported MATCH_REFERENCE_ACCESS_MODE"):
        reference_utils.resolve_reference_access_mode()
