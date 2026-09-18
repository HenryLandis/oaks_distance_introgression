Install Admixture.
```bash
wget https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz
gunzip admixture_linux-1.3.0.tar.gz
tar xvf admixture_linux-1.3.0.tar
```

Make folders to contain output of Admixture. Then, run Admixture for K=1-12.
```bash
mkdir admixk1 admixk2 admixk3 admixk4 admixk5 admixk6 admixk7 admixk8 admixk9 admixk10 admixk11 admixk12
sbatch admixturek1.srun
sbatch admixturek2.srun
sbatch admixturek3.srun
sbatch admixturek4.srun
sbatch admixturek5.srun
sbatch admixturek6.srun
sbatch admixturek7.srun
sbatch admixturek8.srun
sbatch admixturek9.srun
sbatch admixturek10.srun
sbatch admixturek11.srun
sbatch admixturek12.srun
```

Average CV values per k.
```bash
sbatch error.srun
```