# Software required
## Sequence downloading
```sh
fastq-dump
```

## QC

### `illumina-utils`

Link: https://github.com/merenlab/illumina-utils

```sh
conda create -n illumina-utils python=3
source activate illumina-utils
pip install illumina-utils
source deactivate
```

```sh
bowtie2 --version | head -n 1 | awk '{print $3}'
2.3.3.1

source activate illumina-utils
iu-filter-quality-minoche -v
Illumina-utils v2.0.2
source deactivate

samtools --version
samtools 1.5
Using htslib 1.6

megahit -v
MEGAHIT v1.1.2

source activate anvio3
anvi-init-bam -v
Anvi'o version ...............................: 3
Profile DB version ...........................: 20
Contigs DB version ...........................: 9
Pan DB version ...............................: 5
Samples information DB version ...............: 2
Genome data storage version ..................: 1
Auxiliary data storage version ...............: 4
Anvi'server users data storage version .......: 1
source deactivate
```

### MaxBin2
Link: https://downloads.jbei.org/data/microbial_communities/MaxBin/MaxBin.html
Link readme: https://downloads.jbei.org/data/microbial_communities/MaxBin/README.txt
```sh
run_MaxBin.pl -v
MaxBin 2.2.4
MaxBin 2.2.4
```

### CheckM

Install CheckM

```sh
conda create -n CheckM python=2.7
source activate CheckM
conda install -n CheckM numpy scipy matplotlib pysam dendropy ScreamingBackpack
pip install checkm-genome
checkm data setRoot <data_directory>
checkm data update
checkm
source deactivate
```

Test
```sh
checkm test ~/checkm_test_results
```

### PhyloPhlAn

Link: https://huttenhower.sph.harvard.edu/phylophlan

Install

```bash
conda create -n phylophlan python=2.7
source activate phylophlan
conda install -n phylophlan numpy scipy matplotlib
pip install biopython

# download phylophlan

cd programs/phylophlan/
# add phylophlan to $PATH
./phylophlan.py -u example_corynebacteria --nproc 32
./phylophlan.py -i -t example_insertion --nproc 32
source deactivate
```
