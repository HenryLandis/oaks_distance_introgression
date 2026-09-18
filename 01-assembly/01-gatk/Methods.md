Initial input is read files aligned with bwa, then marked for duplicates and combined with picard. GATK is used to call variants. **Do not use default GATK (v4.2.6.1): its output is not compliant with VCF format specifications used by other software.** All scripts specify GATK v4.2.2.0.

Run HaploCaller on the marked individual reads:
```bash
sbatch 01-haplo_caller.srun
```

Run following script twice, first to convert the combined BAM file to BED format, and the second to merge intervals within the BED file. Merging the intervals speeds the creation of genomic databases later.
```bash
sbatch 02-bedtools_intervals.srun
```

Before creating the database, build namelist.sample_map file. It consists of two columns separated by a tab: the first column contains the sample names, and the other contains the filepath to the associated gzipped VCF file. For example:
HL-058_S64_R1_001	vcf_files/HL-058_S64_R1_001.g.vcf.gz

This is the point to determine which samples proceed in the analysis (specific populations, species, etc).

Creating the database requires a large amount of temporary writing space. Create a temporary folder in Kamiak's scratch space:
```bash
cd /scratch/user/henry.landis
export myScratch = "$(mkworkspace )" 
```

Add the scratch folder to the following script and run:
```bash
sbatch 03-genomics_db_import.srun
```

Run script to jointly call SNPs and indels. By default, calls variant sites only. For specific downstream analyses, use -all-sites tag to additionally include invariant sites.
```bash
sbatch 04-genotype_gvcfs.srun
```