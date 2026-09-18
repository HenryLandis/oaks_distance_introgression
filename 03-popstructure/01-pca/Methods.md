Filtered VCF file is input to PLINK for calculating population statistics.

Run a PCA, and also generate output files for further analysis. 
Take .eigenvec and .eigenval files to R for visualizing.
```bash
sbatch 01a-plink_pca_all.srun
```

Version of PCA runs that consider only two of the three species at a time. First requires filtering out columns of the unused species.
```bash
sbatch 01b-plink_pca_pairwise.swun
```

To get a phylogenetic tree, the VCF input can be converted to phylip format. The version of PLINK on Kamiak doesn't support output in phylip format. Instead, use vcf2phylip.
```bash
git clone https://github.com/edgardomortiz/vcf2phylip.git
```

Make versions of the VCF files with a pseudo-outgroup. Script retrived from:
https://www.biostars.org/p/9541926/#:~:text=what%20should%20be%20the%20genotype,not%20heterozygous%20at%20any%20site
Assmes the loci are diploid/unphased and the pseudo-outgroup (reference) is not heterozygous at any site.
```bash
module load htslib
bgzip -d ../../02-filtering/01-python/genotypes_all_sites_variant_retained_final2.vcf.gz

mkdir vcf_with_outgroups
cat ../../02-filtering/01-python/genotypes_all_sites_variant_retained_final2.vcf | ./add_reference.awk > vcf_with_outgroups/genotypes_all_sites_variant_retained_final2.vcf

bgzip ../../02-filtering/01-python/genotypes_all_sites_variant_retained_final2.vcf
bgzip vcf_with_outgroups/genotypes_all_sites_variant_retained_final2.vcf
```

Script is modified to add the following tag at every locus:
0/0:20,0:20:PASS:99:0,0,255
[GT:AD:DP:FT:GQ:PL]
[genotype, read depth per allele, read depth, filter, genotype quality, genotype likelihood]

Convert files with pseudo-outgroups to phylip.
```bash
sbatch 02-vcf2phylip.srun
```