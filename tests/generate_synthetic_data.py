#!/usr/bin/env python3
import argparse
import os
import random
from typing import List

BASES = ["A", "C", "G", "T"]


def sample_genotype(p: float, rng: random.Random) -> int:
    g = 0
    if rng.random() < p:
        g += 1
    if rng.random() < p:
        g += 1
    return g


def generate_vcf(
    out_path: str,
    n_samples: int,
    n_snps: int,
    seed: int,
    missing_rate: float,
    pos_start: int,
    pos_step: int,
) -> List[str]:
    rng = random.Random(seed)
    samples = [f"ID{idx:04d}" for idx in range(1, n_samples + 1)]

    with open(out_path, "w", newline="\n") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)
        handle.write(header + "\n")

        for snp_idx in range(n_snps):
            pos = pos_start + snp_idx * pos_step
            ref = rng.choice(BASES)
            alt = rng.choice([b for b in BASES if b != ref])
            maf = rng.uniform(0.05, 0.5)

            row = ["1", str(pos), f"rs{snp_idx + 1}", ref, alt, ".", "PASS", ".", "GT"]
            for _ in range(n_samples):
                if rng.random() < missing_rate:
                    gt = "./."
                else:
                    g = sample_genotype(maf, rng)
                    if g == 0:
                        gt = "0/0"
                    elif g == 1:
                        gt = "0/1"
                    else:
                        gt = "1/1"
                row.append(gt)
            handle.write("\t".join(row) + "\n")

    return samples


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate synthetic VCF data for match tests.")
    parser.add_argument("--out-dir", default=os.path.join(os.path.dirname(__file__), "synthetic"))
    parser.add_argument("--vcf-name", default="base.vcf")
    parser.add_argument("--n-samples", type=int, default=50)
    parser.add_argument("--n-snps", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20240219)
    parser.add_argument("--missing-rate", type=float, default=0.01)
    parser.add_argument("--pos-start", type=int, default=100000)
    parser.add_argument("--pos-step", type=int, default=10)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    vcf_path = os.path.join(args.out_dir, args.vcf_name)
    samples = generate_vcf(
        vcf_path,
        args.n_samples,
        args.n_snps,
        args.seed,
        args.missing_rate,
        args.pos_start,
        args.pos_step,
    )

    with open(os.path.join(args.out_dir, "samples.txt"), "w", newline="\n") as handle:
        for sample in samples:
            handle.write(sample + "\n")


if __name__ == "__main__":
    main()
