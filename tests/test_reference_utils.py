from pathlib import Path

import numpy as np
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


@pytest.mark.parametrize("mode", [None, "bulk", "LeGaCy"])
def test_reference_alleles_match_full_alleles_and_bounds(tmp_path, monkeypatch, mode):
    fasta = tmp_path / "reference.fa"
    write_fasta(fasta, {"chr1": "agcNryta", "chr2": "t"})
    if mode is None:
        monkeypatch.delenv("MATCH_REFERENCE_ACCESS_MODE", raising=False)
    else:
        monkeypatch.setenv("MATCH_REFERENCE_ACCESS_MODE", mode)
    queries = [
        ("chr1", 1, "A", True),
        ("chr1", 1, "AG", True),
        ("chr1", 1, "AGC", True),
        ("chr1", 1, "AT", False),  # First base alone agrees.
        ("chr1", 2, "GC", True),
        ("chr1", 3, "CAAA", False),  # Reference ambiguity characters.
        ("chr1", 7, "TA", True),
        ("chr1", 8, "A", True),
        ("chr1", 8, "AG", False),
        ("chr1", 9, "A", False),
        ("chr1", 1, "AGCAAAAAA", False),
        ("chr1", 0, "A", False),  # Negative startswith offset would match the end.
        ("chr1", -1, "TA", False),
        ("chr1", 1, "", False),
        ("chr1", 9, "", False),
        ("chr2", 1, "T", True),
        ("chr2", 1, "TA", False),
        ("missing", 1, "A", False),
        ("missing", 1, "", False),
        ("chr1", 1, "AG", True),  # Preserve duplicates and original order.
    ]
    contigs, positions, alleles, expected = zip(*queries)
    result = reference_utils.reference_alleles_match(fasta, contigs, positions, alleles)
    assert result.dtype == np.dtype(bool)
    np.testing.assert_array_equal(result, expected)
    # Also accept arrays directly, as future table-based callers will supply them.
    np.testing.assert_array_equal(
        reference_utils.reference_alleles_match(
            fasta, np.array(contigs, dtype=object), np.array(positions), np.array(alleles, dtype=object)
        ),
        expected,
    )
    # Existing one-base APIs retain their results, including ambiguity characters.
    base_queries = [("chr1", 1), ("chr1", 4), ("chr1", 9), ("missing", 1)]
    bases = reference_utils.fetch_reference_bases(fasta, base_queries)
    assert bases == {("chr1", 1): "A", ("chr1", 4): "N", ("chr1", 9): "", ("missing", 1): ""}
    for (contig, pos), base in bases.items():
        assert reference_utils.fetch_reference_base(fasta, contig, pos) == base


@pytest.mark.parametrize("mode", ["BULK", "LEGACY"])
def test_reference_alleles_match_empty_batch_does_not_open_fasta(tmp_path, monkeypatch, mode):
    monkeypatch.setenv("MATCH_REFERENCE_ACCESS_MODE", mode)
    result = reference_utils.reference_alleles_match(tmp_path / "absent.fa", [], [], [])
    assert result.shape == (0,)
    assert result.dtype == np.dtype(bool)


@pytest.mark.parametrize("contigs,positions,alleles", [(["chr1"], [], ["A"]), ([], [], ["A"])])
def test_reference_alleles_match_rejects_misaligned_queries(tmp_path, contigs, positions, alleles):
    with pytest.raises(ValueError, match="must have equal lengths"):
        reference_utils.reference_alleles_match(tmp_path / "absent.fa", contigs, positions, alleles)


def test_reference_alleles_match_bulk_batches_distinct_queries_and_isolates_paths(tmp_path, monkeypatch):
    monkeypatch.setenv("MATCH_REFERENCE_ACCESS_MODE", "BULK")
    fasta_a = tmp_path / "a.fa"
    fasta_b = tmp_path / "b.fa"
    write_fasta(fasta_a, {"chr1": "AGC" * 1000, "chr2": "TAC" * 1000})
    write_fasta(fasta_b, {"chr1": "TTC" * 1000, "chr2": "GAC" * 1000})
    open_fasta = reference_utils.open_indexed_fasta
    fetches = []

    class TrackedFasta:
        def __init__(self, path):
            self.path = path
            self.fasta = open_fasta(path)
            self.references = self.fasta.references

        def fetch(self, contig, start=None, end=None):
            fetches.append((self.path, contig, start, end))
            return self.fasta.fetch(contig, start, end)

        def __getattr__(self, name):
            return getattr(self.fasta, name)

    monkeypatch.setattr(reference_utils, "open_indexed_fasta", TrackedFasta)
    # Interleave contigs to exercise positional alignment, using distinct positions
    # so a per-unique-variant indexed fetch implementation cannot pass this test.
    contigs = ["chr1", "chr2"] * 1000
    positions = [3 * (index // 2) + 1 for index in range(len(contigs))]
    alleles = ["AG", "TAC"] * 1000
    result_a = reference_utils.reference_alleles_match(fasta_a, contigs, positions, alleles)
    result_b = reference_utils.reference_alleles_match(fasta_b, contigs, positions, alleles)
    assert result_a.all()
    assert not result_b.any()
    np.testing.assert_array_equal(
        reference_utils.reference_alleles_match(fasta_a, contigs, positions, alleles), result_a
    )
    assert reference_utils.fetch_reference_bases(fasta_a, [("chr1", 1)]) == {("chr1", 1): "A"}
    assert fetches
    assert len(fetches) < len(positions) // 10


@pytest.mark.parametrize("mode", ["BULK", "LEGACY"])
def test_reference_alleles_match_missing_index_fails_clearly(tmp_path, monkeypatch, mode):
    monkeypatch.setenv("MATCH_REFERENCE_ACCESS_MODE", mode)
    fasta = tmp_path / "unindexed.fa"
    fasta.write_text(">chr1\nAGC\n", encoding="utf-8")
    with pytest.raises(ValueError, match="FASTA index not found"):
        reference_utils.reference_alleles_match(fasta, ["chr1"], [1], ["AG"])


@pytest.mark.parametrize("mode", ["BULK", "LEGACY"])
def test_reference_alleles_match_fetch_failure_is_not_a_mismatch(tmp_path, monkeypatch, mode):
    monkeypatch.setenv("MATCH_REFERENCE_ACCESS_MODE", mode)

    class BrokenFasta:
        references = ("chr1",)

        def get_reference_length(self, contig):
            return 3

        def fetch(self, *args):
            raise OSError("unreadable reference")

    monkeypatch.setattr(reference_utils, "open_indexed_fasta", lambda path: BrokenFasta())
    with pytest.raises(ValueError, match="unreadable reference"):
        reference_utils.reference_alleles_match(tmp_path / "broken.fa", ["chr1"], [1], ["AG"])
