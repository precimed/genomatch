# Variant Preparation Examples

Assumptions:

- run from the repository root
- activate the `match-liftover` conda environment first
- reference assets are configured as described in `match/INSTALL.md`
- examples write outputs under `variant_lists/prepared/`
- all examples assume the default destination build setting, `--dst-build GRCh38`

## Data source

The raw input files used in the examples below come from these source collections:

- `CMIG`: `eadb_core_chrALL.bim`, `DEMGENEv230501_all_link1.bim`, `huntstudy_HNT_CHR1_22_PID107822.pvar`, `tromsostudy_TU_imputed_8030_00676.pvar`, `husk_all.autosomes.v240619.plink1.bim`
- `TSD`: `HGDP_1KG_v260302_all.pvar`
- `CMIG` shard-based collections: UKB genotype `.bim` shards and UKB HRC-imputed QCed `.bim` shards

Keep environment-specific source paths out of this file. If you need the original local-path note, check [`variant_lists/README.md`](/home/oleksanf/github/precimed/genotools/variant_lists/README.md).

## Recommended workflow

Use `match/prepare_variants.py` as the default entrypoint. Typical pattern:

```bash
python3 match/prepare_variants.py --input-format <bim|pvar|vcf|sumstats> \
  --input <raw-input> \
  --output variant_lists/prepared/<name>
```

`--output` is a stem. The final prepared artifact is written as `<output>.vmap`, and retained stage outputs default to the same stem unless you set `--prefix` explicitly.
If you are resuming a previously canceled run, add ``--resume``.
For force delete outputs and re-run, use ``--force``.

For intersection of two or more of prepared artifacts, use:

```bash
python3 match/intersect_variants.py --inputs \
    variant_lists/prepared/<first>.vmap \
    variant_lists/prepared/<second>.vmap \
  --output variant_lists/shared/<shared-name>.vtable
```

For projecting a payload into an existing shared target universe, use:

```bash
python3 match/project_payload.py --input-format <bfile|pfile|sumstats|sumstats-clean> \
  --input <raw-input> \
  --sumstats-metadata <sumstats-metadata-if-needed> \
  --source-vmap variant_lists/prepared/<source-name>.vmap \
  --target variant_lists/shared/<shared-name>.vtable \
  --output <projected-output>
```

## Workflow steps

### Reference and cohort inputs

```bash
python3 match/prepare_variants.py --input-format bim \
  --input variant_lists/eadb_core_chrALL.bim \
  --output variant_lists/prepared/eadb_core_chrALL

python3 match/prepare_variants.py --input-format vcf --resume \
  --input variant_lists/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz \
  --output variant_lists/prepared/HRC.r1-1.GRCh37.wgs.mac5.sites

python3 match/prepare_variants.py --input-format bim --resume \
  --input variant_lists/DEMGENEv230501_all_link1.bim \
  --output variant_lists/prepared/DEMGENEv230501_all_link1

python3 match/prepare_variants.py --input-format bim --resume \
  --input variant_lists/husk_all.autosomes.v240619.plink1.bim \
  --output variant_lists/prepared/husk_all.autosomes.v240619.plink1

python3 match/prepare_variants.py --input-format pvar --resume \
  --input variant_lists/huntstudy_HNT_CHR1_22_PID107822.pvar \
  --output variant_lists/prepared/huntstudy_HNT_CHR1_22_PID107822

python3 match/prepare_variants.py --input-format pvar --resume \
  --input variant_lists/tromsostudy_TU_imputed_8030_00676.pvar \
  --output variant_lists/prepared/tromsostudy_TU_imputed_8030_00676

python3 match/prepare_variants.py --input-format pvar --resume \
  --input variant_lists/HGDP_1KG_v260302_all.pvar \
  --output variant_lists/prepared/HGDP_1KG_v260302_all
```

Notes:

- `HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz` is VCF-like enough for the current importer: it starts with `#CHROM` and its first five columns are `CHROM POS ID REF ALT`.
- `import_vcf.py` does not require a `.vcf` extension. It reads gzipped text and uses the first five tab-separated columns.
- rows with multiallelic `ALT` or non-`A/C/G/T` alleles will be dropped to import QC
- `HGDP_1KG_v260302_all.pvar` may contain `PAR1` and `PAR2`; the current normalization path aliases those labels to `X`

### Summary-stat inputs

#### UKB height

```bash
python3 match/prepare_variants.py --input-format sumstats --resume \
  --input variant_lists/sumstats/50_irnt.gwas.imputed_v3.both_sexes.tsv.gz \
  --sumstats-metadata variant_lists/sumstats/UKB_HEIGHT_2018_irnt.variant_only.yaml \
  --output variant_lists/prepared/UKB_HEIGHT_2018_irnt
```

Notes:
- this helper metadata imports `variant` as a joined `CHR:POS:allele1:allele2` field

#### GLGC HDL without UKB

```bash
python3 match/prepare_variants.py --input-format sumstats --resume \
  --input variant_lists/sumstats/without_UKB_HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.indels.gz \
  --sumstats-metadata variant_lists/sumstats/GLGC_LIPIDS_2021_HDL_EUR_noUKB.yaml \
  --output variant_lists/prepared/GLGC_LIPIDS_2021_HDL_EUR_noUKB
```

### Generated summary-stat YAML for `phs_models/`

Generated metadata files:

- `variant_lists/phs_models/eadb_ALZ_phs_apoe_genome_2026_01_29_model.yaml`
- `variant_lists/phs_models/DesikanModel_noAPOE.yaml`
- `variant_lists/phs_models/01_geneticScore.finalModel.snplist.GRCh38.yaml`

Commands:

```bash
python3 match/prepare_variants.py --input-format sumstats --resume \
  --input variant_lists/phs_models/eadb_ALZ_phs_apoe_genome_2026_01_29_model.txt \
  --sumstats-metadata variant_lists/phs_models/eadb_ALZ_phs_apoe_genome_2026_01_29_model.yaml \
  --output variant_lists/prepared/eadb_ALZ_phs_apoe_genome_2026_01_29_model

python3 match/prepare_variants.py --input-format sumstats --resume \
  --input variant_lists/phs_models/DesikanModel_noAPOE.txt \
  --sumstats-metadata variant_lists/phs_models/DesikanModel_noAPOE.yaml \
  --output variant_lists/prepared/DesikanModel_noAPOE

python3 match/prepare_variants.py --input-format sumstats --resume \
  --input variant_lists/phs_models/01_geneticScore.finalModel.snplist.csv \
  --sumstats-metadata variant_lists/phs_models/01_geneticScore.finalModel.snplist.GRCh38.yaml \
  --output variant_lists/prepared/01_geneticScore.finalModel.snplist.GRCh38
```

Notes:

- `DesikanModel_noAPOE.txt` begins with many `#` comment lines; the summary-stat importer skips comment lines before reading the header row
- `eadb_ALZ_phs_apoe_genome_2026_01_29_model.txt` has both `POSITION` and `BP`; the generated metadata uses `BP` because `POSITION` is `0` in sampled rows
- `01_geneticScore.finalModel.snplist.csv` contains both `Position (GRCh37)` and `Position (GRCh38)`; the generated metadata uses `Position (GRCh38)` to align with `--dst-build GRCh38`


## Sharded raw imports with `@`

You can still use raw-input shard discovery with `prepare_variants.py`:

```bash
python3 match/prepare_variants.py --input-format bim --resume \
  --input variant_lists/ukb_genotype/ukb_snp_@_v2.bim \
  --output variant_lists/prepared/ukb_genotype

python3 match/prepare_variants.py --input-format bim --resume \
  --input variant_lists/UKB_imputed_QCed/ukb27412_imp_chr@_v3.bim \
  --output variant_lists/prepared/UKB_imputed_QCed

python3 match/prepare_variants.py --input-format bim --resume \
  --input variant_lists/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
  --output variant_lists/prepared/1000G_EUR_Phase3_plink
```

Shard-discovery notes:

- `variant_lists/ukb_genotype/ukb_snp_@_v2.bim` captures `ukb_snp_chr1_v2.bim` through `ukb_snp_chr22_v2.bim`, plus `chrX`, `chrY`, `chrXY`, and `chrMT`
- `variant_lists/UKB_imputed_QCed/ukb27412_imp_chr@_v3.bim` captures `ukb27412_imp_chr23_v3.bim`, because token `23` is substituted into the `@` slot after the literal `chr` prefix
- `variant_lists/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim` captures the numbered shard files; the literal placeholder file `1000G.EUR.QC.@.bim` is not used by shard discovery


### Intersect lifted 1000G EUR Phase3, UKB imputed QCed, and UKB height

Use the prepared artifacts to define one shared `GRCh38` target universe:

```bash
python3 match/intersect_variants.py --inputs \
    variant_lists/prepared/1000G_EUR_Phase3_plink.vmap \
    variant_lists/prepared/ukb_genotype.vmap \
  --output variant_lists/shared/grch38_shared.vtable
```

Then use `project_payload.py` for the payloads you want to rewrite into that shared universe.

For the UKB height summary-stat payload:

```bash
python3 match/project_payload.py --input-format sumstats-clean \
  --input variant_lists/sumstats/50_irnt.gwas.imputed_v3.both_sexes.tsv.gz \
  --sumstats-metadata variant_lists/sumstats/UKB_HEIGHT_2018_irnt.variant_only.yaml \
  --source-vmap variant_lists/prepared/UKB_HEIGHT_2018_irnt.vmap \
  --target variant_lists/shared/grch38_shared.vtable \
  --output variant_lists/shared/UKB_HEIGHT_2018_irnt.tsv.gz
```

For the GLGC HDL summary-stat payload:

```bash
python3 match/project_payload.py --input-format sumstats-clean \
  --input variant_lists/sumstats/without_UKB_HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.indels.gz \
  --sumstats-metadata variant_lists/sumstats/GLGC_LIPIDS_2021_HDL_EUR_noUKB.yaml \
  --source-vmap variant_lists/prepared/GLGC_LIPIDS_2021_HDL_EUR_noUKB.vmap \
  --target variant_lists/shared/grch38_shared.vtable \
  --output variant_lists/shared/GLGC_LIPIDS_2021_HDL_EUR_noUKB.tsv.gz
```

For the 1000G EUR Phase3 PLINK payload:

```bash
python3 match/project_payload.py --input-format bfile \
  --input variant_lists/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
  --source-vmap variant_lists/prepared/1000G_EUR_Phase3_plink.vmap \
  --target variant_lists/shared/grch38_shared.vtable \
  --output variant_lists/shared/1000G_EUR_Phase3_plink/1000G.EUR.QC.@
```

For the UKB genotypes PLINK payload (can be executed on a machine where full `.bed/.bim/.fam` inputs are present ):

```bash
python3 match/project_payload.py --input-format bfile \
  --input variant_lists/ukb_genotype/ukb_snp_@_v2.bim \
  --source-vmap variant_lists/prepared/ukb_genotype.vmap \
  --target variant_lists/shared/grch38_shared.vtable \
  --output variant_lists/shared/ukb_genotype/ukb_snp_@_v2
```
