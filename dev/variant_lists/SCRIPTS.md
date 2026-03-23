# Variant Preparation Examples

Assumptions:

- run from the local `variant_lists/` directory
- `MATCH_CONFIG` is already set
- the local `match-liftover` development environment is already active and provides the `genomatch` CLI tools
- examples write outputs under `prepared/`
- all examples assume the default destination build setting, `--dst-build GRCh38`

## Data source

The raw input files used in the examples below come from these source collections:

- `CMIG`: `eadb_core_chrALL.bim`, `DEMGENEv230501_all_link1.bim`, `huntstudy_HNT_CHR1_22_PID107822.pvar`, `tromsostudy_TU_imputed_8030_00676.pvar`, `husk_all.autosomes.v240619.plink1.bim`
- `TSD`: `HGDP_1KG_v260302_all.pvar`
- `CMIG` shard-based collections: UKB genotype `.bim` shards and UKB HRC-imputed QCed `.bim` shards

## Recommended workflow

Use `prepare_variants.py` as the default entrypoint. Typical pattern:

```bash
prepare_variants.py --input-format <bim|pvar|vcf|sumstats> \
  --input <raw-input> \
  --output prepared/<name>
```

`--output` is a stem. The final prepared artifact is written as `<output>.vmap`, and retained stage outputs default to the same stem unless you set `--prefix` explicitly.
If you are resuming a previously canceled run, add ``--resume``.
For force delete outputs and re-run, use ``--force``.

For intersection of two or more of prepared artifacts, use:

```bash
intersect_variants.py --inputs \
    prepared/<first>.vmap \
    prepared/<second>.vmap \
  --output shared/<shared-name>.vtable
```

For projecting a payload into an existing shared target universe, use:

```bash
project_payload.py --input-format <bfile|pfile|sumstats|sumstats-clean> \
  --input <raw-input> \
  --sumstats-metadata <sumstats-metadata-if-needed> \
  --source-vmap prepared/<source-name>.vmap \
  --target shared/<shared-name>.vtable \
  --output <projected-output>
```

## Workflow steps

### Reference and cohort inputs

```bash
prepare_variants.py --input-format bim \
  --input eadb_core_chrALL.bim \
  --output prepared/eadb_core_chrALL

prepare_variants.py --input-format vcf --resume \
  --input HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz \
  --output prepared/HRC.r1-1.GRCh37.wgs.mac5.sites

prepare_variants.py --input-format bim --resume \
  --input DEMGENEv230501_all_link1.bim \
  --output prepared/DEMGENEv230501_all_link1

prepare_variants.py --input-format bim --resume \
  --input husk_all.autosomes.v240619.plink1.bim \
  --output prepared/husk_all.autosomes.v240619.plink1

prepare_variants.py --input-format pvar --resume \
  --input huntstudy_HNT_CHR1_22_PID107822.pvar \
  --output prepared/huntstudy_HNT_CHR1_22_PID107822

prepare_variants.py --input-format pvar --resume \
  --input tromsostudy_TU_imputed_8030_00676.pvar \
  --output prepared/tromsostudy_TU_imputed_8030_00676

prepare_variants.py --input-format pvar --resume \
  --input HGDP_1KG_v260302_all.pvar \
  --output prepared/HGDP_1KG_v260302_all
```

Notes:

- `HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz` is VCF-like enough for the current importer: it starts with `#CHROM` and its first five columns are `CHROM POS ID REF ALT`.
- `import_vcf.py` does not require a `.vcf` extension. It reads gzipped text and uses the first five tab-separated columns.
- rows with multiallelic `ALT` or non-`A/C/G/T` alleles will be dropped to import QC
- `HGDP_1KG_v260302_all.pvar` may contain `PAR1` and `PAR2`; the current normalization path aliases those labels to `X`

### Summary-stat inputs

#### UKB height

```bash
prepare_variants.py --input-format sumstats --resume \
  --input sumstats/50_irnt.gwas.imputed_v3.both_sexes.tsv.gz \
  --sumstats-metadata sumstats/UKB_HEIGHT_2018_irnt.variant_only.yaml \
  --output prepared/UKB_HEIGHT_2018_irnt
```

Notes:
- this helper metadata imports `variant` as a joined `CHR:POS:allele1:allele2` field

#### GLGC HDL without UKB

```bash
prepare_variants.py --input-format sumstats --resume \
  --input sumstats/without_UKB_HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.indels.gz \
  --sumstats-metadata sumstats/GLGC_LIPIDS_2021_HDL_EUR_noUKB.yaml \
  --output prepared/GLGC_LIPIDS_2021_HDL_EUR_noUKB
```

### Generated summary-stat YAML for `phs_models/`

Generated metadata files:

- `phs_models/eadb_ALZ_phs_apoe_genome_2026_01_29_model.yaml`
- `phs_models/DesikanModel_noAPOE.yaml`
- `phs_models/01_geneticScore.finalModel.snplist.GRCh38.yaml`

Commands:

```bash
prepare_variants.py --input-format sumstats --resume \
  --input phs_models/eadb_ALZ_phs_apoe_genome_2026_01_29_model.txt \
  --sumstats-metadata phs_models/eadb_ALZ_phs_apoe_genome_2026_01_29_model.yaml \
  --output prepared/eadb_ALZ_phs_apoe_genome_2026_01_29_model

prepare_variants.py --input-format sumstats --resume \
  --input phs_models/DesikanModel_noAPOE.txt \
  --sumstats-metadata phs_models/DesikanModel_noAPOE.yaml \
  --output prepared/DesikanModel_noAPOE

prepare_variants.py --input-format sumstats --resume \
  --input phs_models/01_geneticScore.finalModel.snplist.csv \
  --sumstats-metadata phs_models/01_geneticScore.finalModel.snplist.GRCh38.yaml \
  --output prepared/01_geneticScore.finalModel.snplist.GRCh38
```

Notes:

- `DesikanModel_noAPOE.txt` begins with many `#` comment lines; the summary-stat importer skips comment lines before reading the header row
- `eadb_ALZ_phs_apoe_genome_2026_01_29_model.txt` has both `POSITION` and `BP`; the generated metadata uses `BP` because `POSITION` is `0` in sampled rows
- `01_geneticScore.finalModel.snplist.csv` contains both `Position (GRCh37)` and `Position (GRCh38)`; the generated metadata uses `Position (GRCh38)` to align with `--dst-build GRCh38`


## Sharded raw imports with `@`

You can still use raw-input shard discovery with `prepare_variants.py`:

```bash
prepare_variants.py --input-format bim --resume \
  --input ukb_genotype/ukb_snp_@_v2.bim \
  --output prepared/ukb_genotype

prepare_variants.py --input-format bim --resume \
  --input UKB_imputed_QCed/ukb27412_imp_chr@_v3.bim \
  --output prepared/UKB_imputed_QCed

prepare_variants.py --input-format bim --resume \
  --input 1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
  --output prepared/1000G_EUR_Phase3_plink
```

Shard-discovery notes:

- `ukb_genotype/ukb_snp_@_v2.bim` captures `ukb_snp_chr1_v2.bim` through `ukb_snp_chr22_v2.bim`, plus `chrX`, `chrY`, `chrXY`, and `chrMT`
- `UKB_imputed_QCed/ukb27412_imp_chr@_v3.bim` captures `ukb27412_imp_chr23_v3.bim`, because token `23` is substituted into the `@` slot after the literal `chr` prefix
- `1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim` captures the numbered shard files; the literal placeholder file `1000G.EUR.QC.@.bim` is not used by shard discovery


### Intersect lifted 1000G EUR Phase3, UKB imputed QCed, and UKB height

Use the prepared artifacts to define one shared `GRCh38` target universe:

```bash
intersect_variants.py --inputs \
    prepared/1000G_EUR_Phase3_plink.vmap \
    prepared/ukb_genotype.vmap \
  --output shared/grch38_shared.vtable
```

Then use `project_payload.py` for the payloads you want to rewrite into that shared universe.

For the UKB height summary-stat payload:

```bash
project_payload.py --input-format sumstats-clean \
  --input sumstats/50_irnt.gwas.imputed_v3.both_sexes.tsv.gz \
  --sumstats-metadata sumstats/UKB_HEIGHT_2018_irnt.variant_only.yaml \
  --source-vmap prepared/UKB_HEIGHT_2018_irnt.vmap \
  --target shared/grch38_shared.vtable \
  --output shared/UKB_HEIGHT_2018_irnt.tsv.gz
```

For the GLGC HDL summary-stat payload:

```bash
project_payload.py --input-format sumstats-clean \
  --input sumstats/without_UKB_HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.indels.gz \
  --sumstats-metadata sumstats/GLGC_LIPIDS_2021_HDL_EUR_noUKB.yaml \
  --source-vmap prepared/GLGC_LIPIDS_2021_HDL_EUR_noUKB.vmap \
  --target shared/grch38_shared.vtable \
  --output shared/GLGC_LIPIDS_2021_HDL_EUR_noUKB.tsv.gz
```

For the 1000G EUR Phase3 PLINK payload:

```bash
project_payload.py --input-format bfile \
  --input 1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
  --source-vmap prepared/1000G_EUR_Phase3_plink.vmap \
  --target shared/grch38_shared.vtable \
  --output shared/1000G_EUR_Phase3_plink/1000G.EUR.QC.@
```

For the UKB genotypes PLINK payload (can be executed on a machine where full `.bed/.bim/.fam` inputs are present ):

```bash
project_payload.py --input-format bfile \
  --input ukb_genotype/ukb_snp_@_v2.bim \
  --source-vmap prepared/ukb_genotype.vmap \
  --target shared/grch38_shared.vtable \
  --output shared/ukb_genotype/ukb_snp_@_v2
```
