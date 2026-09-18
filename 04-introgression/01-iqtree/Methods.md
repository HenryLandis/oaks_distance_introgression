Create symbolic link to phylip file.
```bash
ln -s ../../03-popstructure/01-pca/phylip_files/genotypes_all_sites_variant_retained_final2.min4.phy data.phy
```

Run IQ-TREE.
```bash
sbatch iqtree.srun
```