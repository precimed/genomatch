from plink_vcf_helpers import build_vmap
from utils import read_tsv, write_json, write_lines


def test_build_vmap_composes_swap_with_existing_source_allele_op(tmp_path):
    source = tmp_path / "source.vmap"
    target = tmp_path / "target.vtable"
    output = tmp_path / "output.vmap"
    write_lines(source, ["1\t100\tsource\tA\tG\t.\t0\tswap"])
    write_json(
        source.with_name(source.name + ".meta.json"),
        {"object_type": "variant_map", "target": {"genome_build": "GRCh37", "contig_naming": "ncbi"}},
    )
    write_lines(target, ["1\t100\ttarget\tG\tA"])
    write_json(
        target.with_name(target.name + ".meta.json"),
        {"object_type": "variant_table", "genome_build": "GRCh37", "contig_naming": "ncbi"},
    )

    build_vmap(source, target, output)

    assert read_tsv(output) == [["1", "100", "target", "G", "A", ".", "0", "identity"]]
