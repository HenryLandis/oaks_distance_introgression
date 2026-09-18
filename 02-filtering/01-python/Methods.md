Input is VCF files that have been genotyped with GATK pipeline.

Split variant and invariant sites.
```bash
sbatch 00-variant_split.srun
```

Scripts to remove corrupted lines from VCF files.
```bash
sbatch 01a-variant_remove_null.srun
sbatch 01b-variant_split_remove_null.srun
sbatch 01c-invariant_split_remove_null.srun
```

Script to prepare and run biallelic-only files through GATK variant filtration. Parameters control for read/mapping quality, strand bias, and excess heterozygosity. Applies PASS or other markers to reads.
```bash
sbatch 02a-variant_filtering.srun # variants only
sbatch 02b-variant_filtering.srun # all sites
```

Further processing the variants and allsites files.
```bash
sbatch 03a-vcftools_variants.srun
sbatch 03b-vcftools_allsites.srun
```

Rename columns for downstream analysis, and remove columns with too much missing data.
```bash
sbatch 04-rename_columns.srun
sbatch 05-remove_columns.srun
```