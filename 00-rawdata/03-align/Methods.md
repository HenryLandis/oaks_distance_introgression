Single-end RAD-seq FASTQ files from University of Minnesota have been previously demultiplexed, as well as trimmed with fastp to remove adapters and restriction overhangs.

Reads are aligned to the *Quercus lobata* genome. Place in separate "genome" directory.
```bash
wget -O Qlobata_genome.fasta.gz https://valleyoak.ucla.edu/genomic-resources/Qlobata.v3.0.RptMsk4.0.6.on-RptMdl1.0.8.softmasked.fasta.gz
```

The genome must be indexed prior to aligning reads.
```bash
bwa index ../../../genome/Qlobata_genome.fasta
```

Run script to align trimmed reads to genome, then take the output SAM files and convert them to BAM files, followed by sorting.
```bash
sbatch 01-align.srun
```

Mark and remove duplicate reads, then sort/index marked files.
```bash
sbatch 02-picard_dupes.srun
```

Merge the marked BAM files into one file, then sort/index that file.
```bash
sbatch 03-merge.srun
```

Output in picard_files should be retained for further processing. For disk space considerations, intermediate SAM/BAM files can be deleted if necessary.
