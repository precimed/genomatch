"""Regenerate Tutorial 1 independently and verify its documented command path."""
import importlib.util
import json
import re
import shlex
import sys

import pytest
import yaml

from real_liftover_helpers import (
    CHAIN_37_TO_38, CHAIN_38_TO_37, GRCH37_UCSC_SOURCE_FASTA,
    GRCH38_UCSC_SOURCE_FASTA, ensure_real_liftover_assets,
    resolve_bcftools_with_liftover,
)
from utils import REPO_ROOT, write_primary_ucsc_fasta


spec = importlib.util.spec_from_file_location("tutorial_fixture", REPO_ROOT / "scripts/generate_tutorial_1_fixture.py")
tutorial = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = tutorial
spec.loader.exec_module(tutorial)


def test_tutorial_expectations_use_full_ref_and_preserve_supplied_order(tmp_path):
    resolve_bcftools_with_liftover()
    fasta = tmp_path / "reference.fa"
    write_primary_ucsc_fasta(fasta, {"chr1": "AGCAGCAGC"})
    scenarios = [
        tutorial.Scenario("deletion", "ref_cohort", tutorial.Candidate("1", 1, "AG", "A", 1, "A", "C"),
                          (), True, False, True, variant_kind="deletion"),
        tutorial.Scenario("opposite", "ref_cohort", tutorial.Candidate("1", 4, "AG", "A", 4, "A", "C"),
                          (), True, False, True, variant_kind="deletion", allele_mode_cohort="flip"),
        tutorial.Scenario("both_match", "ref_cohort", tutorial.Candidate("1", 7, "AG", "AGC", 7, "A", "C"),
                          (), True, False, True, variant_kind="insertion", allele_mode_cohort="swap"),
    ]
    expected = tutorial.expected_preparation(scenarios, {"grch38_fasta": fasta}, tmp_path)
    assert expected["cohort"]["deletion"] == ("1", "1", "A", "AG")
    assert "opposite" not in expected["cohort"]
    assert expected["cohort"]["both_match"] == ("1", "8", "G", "GC")
    assert expected["reference"]["both_match"] != expected["cohort"]["both_match"]


def test_tutorial_rejects_first_base_only_fixture_ref(tmp_path):
    fasta = tmp_path / "reference.fa"
    write_primary_ucsc_fasta(fasta, {"chr1": "ATC"})
    scenario = tutorial.Scenario("bad_ref", "ref_cohort", tutorial.Candidate("1", 1, "AG", "A", 1, "A", "C"),
                                 (), True, False, True, variant_kind="deletion")
    with pytest.raises(RuntimeError, match="full FASTA"):
        tutorial.expected_preparation([scenario], {"grch38_fasta": fasta}, tmp_path)


def test_tutorial_regeneration_and_documented_roundtrip(tmp_path, monkeypatch):
    ensure_real_liftover_assets()
    resolve_bcftools_with_liftover()
    assets = {"grch37_fasta": GRCH37_UCSC_SOURCE_FASTA, "grch38_fasta": GRCH38_UCSC_SOURCE_FASTA,
              "chain37to38": CHAIN_37_TO_38, "chain38to37": CHAIN_38_TO_37}
    config = tmp_path / "config.yaml"
    config.write_text(yaml.safe_dump({
        "builds": {"GRCh37": {"ucsc_fasta": str(assets["grch37_fasta"])},
                   "GRCh38": {"ucsc_fasta": str(assets["grch38_fasta"])}},
        "liftover": [{"source": source, "target": target, "chain": str(assets[chain])}
                     for source, target, chain in (("GRCh37", "GRCh38", "chain37to38"),
                                                   ("GRCh38", "GRCh37", "chain38to37"))],
    }))
    monkeypatch.setenv("MATCH_CONFIG", str(config))
    monkeypatch.setattr(tutorial, "FIXTURE_DIR", tmp_path / "fixture")
    scenarios = tutorial.build_scenarios(assets)
    expected = tutorial.expected_preparation(scenarios, assets, tmp_path)
    assert {role: len(keys) for role, keys in expected.items()} == {"reference": 200, "study": 94, "cohort": 149}
    assert {role: sum(len(key[2]) > 1 or len(key[3]) > 1 for key in keys.values())
            for role, keys in expected.items()} == {"reference": 6, "study": 0, "cohort": 5}
    fixture = tmp_path / "fixture"
    fixture.mkdir()
    tutorial.write_reference_vcf(scenarios, fixture / "reference.grch38.vcf")
    tutorial.write_sumstats(scenarios, fixture / "study.grch37.tsv")
    tutorial.write_sumstats_metadata(fixture / "study.grch37.meta.yaml")
    tutorial.write_cohort_bfiles(scenarios, fixture)
    tutorial.write_manifest(scenarios, expected, fixture / "manifest.tsv")
    tutorial.write_fixture_readme(tutorial.venn_counts(expected), fixture / "README.md")
    checked_in = REPO_ROOT / "examples/tutorial-1"
    assert {path.name for path in fixture.iterdir()} == {path.name for path in checked_in.iterdir()}
    for path in fixture.iterdir():
        assert path.read_bytes() == (checked_in / path.name).read_bytes(), path.name

    commands = []
    original_run = tutorial.run_cli

    def capture(cmd, **kwargs):
        commands.append(cmd)
        original_run(cmd, **kwargs)

    monkeypatch.setattr(tutorial, "run_cli", capture)
    work = tmp_path / "roundtrip"
    tutorial.validate_roundtrip(fixture, work, scenarios, expected, assets)

    # Run the documented command arguments, including flags: doc drift must fail.
    blocks = re.findall(r"```bash\n(.*?)```", (REPO_ROOT / "docs/tutorial-1.md").read_text(), re.S)
    documented = []
    for block in blocks:
        for line in block.replace("\\\n", "").splitlines():
            if line.startswith(("prepare_variants", "restrict_vmap.py", "project_payload.py")):
                tokens = shlex.split(line)
                documented.append([token.replace("examples/tutorial-1", str(fixture))
                                   .replace("work/tutorial-1", str(work / "work"))
                                   .replace("out/tutorial-1", str(work / "out")) for token in tokens])
    assert commands[:len(documented)] == documented
    assert len(documented) == 7

    # Successful exit codes cannot hide wrong keys, locators, metadata or recodes.
    corruptions = [
        (work / "work/cohort.vmap", lambda s: s.replace("\tidentity\n", "\tswap\n", 1), "retained sort stage|operation"),
        (work / "work/reference.vmap.meta.json", lambda s: json.dumps({**json.loads(s), "variants_count": 0}), "metadata"),
        (work / "out/cohort.21.bed", None, "genotypes"),
        (work / "out/study.reference.tsv", lambda s: s.replace("\n21\t", "\n22\t", 1), "variant/order"),
    ]
    for path, change, message in corruptions:
        original = path.read_bytes()
        if change:
            altered = change(original.decode()).encode()
        else:
            altered = original[:3] + bytes([original[3] ^ 3]) + original[4:]
        assert altered != original
        path.write_bytes(altered)
        try:
            with pytest.raises(RuntimeError, match=message):
                tutorial.validate_outputs(fixture, work, scenarios, expected, assets)
        finally:
            path.write_bytes(original)
