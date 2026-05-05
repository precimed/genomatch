from pathlib import Path

import pytest

from genomatch import reference_utils
from utils import write_fasta, write_lines, write_match_config, write_primary_ucsc_fasta


def test_resolve_config_path_uses_match_config_override(tmp_path, monkeypatch):
    config = tmp_path / "custom-config.yaml"
    config.write_text("builds: {}\nliftover: []\n", encoding="utf-8")

    monkeypatch.setenv("MATCH_CONFIG", str(config))
    monkeypatch.setattr(reference_utils, "DEFAULT_MATCH_CONFIG", tmp_path / "missing-default.yaml")

    assert reference_utils.resolve_config_path() == config


def test_resolve_config_path_uses_default_ref_config(monkeypatch, tmp_path):
    ref = tmp_path / "ref"
    ref.mkdir()
    config = ref / "config.yaml"
    config.write_text("builds: {}\nliftover: []\n", encoding="utf-8")

    monkeypatch.delenv("MATCH_CONFIG", raising=False)
    monkeypatch.setattr(reference_utils, "DEFAULT_MATCH_CONFIG", config)

    assert reference_utils.resolve_config_path() == config


def test_resolve_config_path_requires_override_or_default(monkeypatch, tmp_path):
    monkeypatch.delenv("MATCH_CONFIG", raising=False)
    monkeypatch.setattr(reference_utils, "DEFAULT_MATCH_CONFIG", tmp_path / "missing-config.yaml")

    with pytest.raises(ValueError, match="MATCH_CONFIG is required unless /ref/config.yaml exists"):
        reference_utils.resolve_config_path()


def test_load_match_config_caches_by_config_path(tmp_path, monkeypatch):
    config = tmp_path / "config.yaml"
    config.write_text("builds: {}\nliftover: []\n", encoding="utf-8")
    calls = []

    def fake_load_yaml(path):
        calls.append(path)
        return {"builds": {}, "liftover": []}

    monkeypatch.setenv("MATCH_CONFIG", str(config))
    monkeypatch.setattr(reference_utils, "_load_yaml", fake_load_yaml)

    assert reference_utils.load_match_config() == {"builds": {}, "liftover": []}
    assert reference_utils.load_match_config() == {"builds": {}, "liftover": []}
    assert calls == [config]


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

    write_primary_ucsc_fasta(grch37, {"chr1": "A"})
    write_primary_ucsc_fasta(grch38, {"chr1": "A"})
    write_lines(chain37to38, ["dummy"])
    write_lines(chain38to37, ["dummy"])
    write_match_config(
        config,
        grch37_ucsc_fasta=grch37,
        grch38_ucsc_fasta=grch38,
        chain37to38=chain37to38,
        chain38to37=chain38to37,
    )

    monkeypatch.delenv("MATCH_CONFIG", raising=False)
    monkeypatch.setattr(reference_utils, "DEFAULT_MATCH_CONFIG", config)

    assert reference_utils.resolve_internal_reference_fasta("GRCh37") == grch37
    assert reference_utils.resolve_internal_reference_fasta("GRCh38") == grch38
    assert reference_utils.resolve_liftover_chain("GRCh37", "GRCh38") == chain37to38
    assert reference_utils.resolve_liftover_chain("GRCh38", "GRCh37") == chain38to37


def test_resolve_internal_reference_fasta_accepts_t2t_build(monkeypatch, tmp_path):
    grch37 = tmp_path / "hg19.fa"
    grch38 = tmp_path / "hg38.fa"
    t2t = tmp_path / "chm13v2.fa"
    config = tmp_path / "config.yaml"
    write_primary_ucsc_fasta(grch37, {"chr1": "A"})
    write_primary_ucsc_fasta(grch38, {"chr1": "A"})
    write_primary_ucsc_fasta(t2t, {"chr1": "A"})
    write_match_config(
        config,
        grch37_ucsc_fasta=grch37,
        grch38_ucsc_fasta=grch38,
        t2t_chm13v2_ucsc_fasta=t2t,
    )

    monkeypatch.setenv("MATCH_CONFIG", str(config))

    assert reference_utils.resolve_internal_reference_fasta("T2T-CHM13v2.0") == t2t


def test_legacy_config_shape_is_rejected(monkeypatch, tmp_path):
    config = tmp_path / "config.yaml"
    config.write_text("references: {}\nchain: {}\n", encoding="utf-8")
    monkeypatch.setenv("MATCH_CONFIG", str(config))

    with pytest.raises(ValueError, match="legacy match config shape"):
        reference_utils.resolve_internal_reference_fasta("GRCh37")


def test_duplicate_liftover_edges_fail_clearly(monkeypatch, tmp_path):
    grch37 = tmp_path / "hg19.fa"
    grch38 = tmp_path / "hg38.fa"
    chain_a = tmp_path / "a.chain"
    chain_b = tmp_path / "b.chain"
    config = tmp_path / "config.yaml"
    write_primary_ucsc_fasta(grch37, {"chr1": "A"})
    write_primary_ucsc_fasta(grch38, {"chr1": "A"})
    write_lines(chain_a, ["dummy"])
    write_lines(chain_b, ["dummy"])
    write_match_config(config, grch37_ucsc_fasta=grch37, grch38_ucsc_fasta=grch38)
    config.write_text(
        "\n".join(
            [
                "builds:",
                "  GRCh37:",
                f"    ucsc_fasta: {grch37}",
                "  GRCh38:",
                f"    ucsc_fasta: {grch38}",
                "liftover:",
                "  - source: GRCh37",
                "    target: GRCh38",
                f"    chain: {chain_a}",
                "  - source: GRCh37",
                "    target: GRCh38",
                f"    chain: {chain_b}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    monkeypatch.setenv("MATCH_CONFIG", str(config))

    with pytest.raises(ValueError, match="duplicate liftover edges"):
        reference_utils.resolve_liftover_chain("GRCh37", "GRCh38")


def test_primary_contig_validation_reports_missing_contigs(monkeypatch, tmp_path):
    grch37 = tmp_path / "hg19.fa"
    config = tmp_path / "config.yaml"
    write_fasta(grch37, {"chr1": "A"})
    config.write_text(
        "\n".join(
            [
                "builds:",
                "  GRCh37:",
                f"    ucsc_fasta: {grch37}",
                "liftover: []",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    monkeypatch.setenv("MATCH_CONFIG", str(config))

    with pytest.raises(ValueError, match="missing required primary contigs"):
        reference_utils.resolve_internal_reference_fasta("GRCh37")


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
