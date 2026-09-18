The D3 implementation by Hibbins & Hahn (2019) starts from a phylip file. Create a symbolic link:
```bash
ln -s ../../03-popstructure/01-pca/phylip_files/genotypes_variants_only_maxmissing_0_75_outgroup.min4.phy data.phy
```

The script assumes exactly four lines. Subset the file to the outgroup and each relevant trio to test.