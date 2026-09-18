Symbolic link to VCF file.
```bash
ln -s ../../03-popstructure/01-pca/vcf_with_outgroups/genotypes_all_sites_variant_retained_final2.vcf.gz input.vcf.gz
```

Species and populations assignment text files comprise the sample list and species/population assignments, with a tab between columns.

Import Newick tree output from IQ-TREE into FigTree. Root on the reference and export.

Run Dsuite and F-branch stats.
```bash
sbatch 01-dsuite_inds.srun
sbatch 02-dsuite_corrections.srun
sbatch 04-fbranch.srun
```