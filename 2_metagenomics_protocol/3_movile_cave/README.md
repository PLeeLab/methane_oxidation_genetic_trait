# Metagenomics from Movile Cave in Mangalia, Romania

<!-- TOC -->

- [Metagenomics from Movile Cave in Mangalia, Romania](#metagenomics-from-movile-cave-in-mangalia-romania)
    - [[Stage 0] Preliminary](#stage-0-preliminary)
    - [[Stage 1] Quality Control (QC)](#stage-1-quality-control-qc)
    - [[Stage 2] Co-assembly](#stage-2-co-assembly)
    - [[Stage 3] Binning](#stage-3-binning)
    - [[Stage 4] Refine bins](#stage-4-refine-bins)
    - [[Stage 5] Functional annotation - `prokka`](#stage-5-functional-annotation---prokka)
    - [[Stage 6] Taxonomy classification of bins](#stage-6-taxonomy-classification-of-bins)
    - [[Stage 7]  Refine functional annotation - `EggNOG-mapper`](#stage-7--refine-functional-annotation---eggnog-mapper)

<!-- /TOC -->

## [Stage 0] Preliminary
### Create project directory
```sh
jv~$ mkdir metagenomics-movile_cave
jv~$ cd metagenomics-movile_cave/
```

### Retrieve sequences from ENA system

```sh
jv~$ mkdir 00_RAW
jv~$ cd 00_RAW/
```

> This dataset was published in:
> Kumaresan, D., Stephenson, J., Doxey, A. C., Bandukwala, H., Brooks, E., Hillebrand-Voiculescu, A., ... & Murrell, J. C. (2018). Aerobic proteobacterial methylotrophs in Movile Cave: genomic and metagenomic analyses. Microbiome, 6(1), 1.
> https://doi.org/10.1186/s40168-017-0383-2

### Download metagenomics data

>"The shotgun metagenome sequences have been deposited at the European Nucleotide Archive under the project accession PRJEB12283."

Microbial mat:
- ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR119/001/ERR1198911/ERR1198911_1.fastq.gz
- ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR119/001/ERR1198911/ERR1198911_2.fastq.gz
- ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR119/002/ERR1198912/ERR1198912_1.fastq.gz
- ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR119/002/ERR1198912/ERR1198912_2.fastq.gz

Sediment:
- ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR119/003/ERR1198913/ERR1198913_1.fastq.gz
- ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR119/003/ERR1198913/ERR1198913_2.fastq.gz
- ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR119/004/ERR1198914/ERR1198914_1.fastq.gz
- ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR119/004/ERR1198914/ERR1198914_2.fastq.gz

```sh
jv~$ ls -1

ERR1198911_1.fastq.gz
ERR1198911_2.fastq.gz
ERR1198912_1.fastq.gz
ERR1198912_2.fastq.gz
ERR1198913_1.fastq.gz
ERR1198913_2.fastq.gz
ERR1198914_1.fastq.gz
ERR1198914_2.fastq.gz
```

Extract `fastq` files

```sh
jv~$ gzip -d ERR1198911_1.fastq.gz
jv~$ gzip -d ERR1198911_2.fastq.gz
jv~$ gzip -d ERR1198912_1.fastq.gz
jv~$ gzip -d ERR1198912_2.fastq.gz
jv~$ gzip -d ERR1198913_1.fastq.gz
jv~$ gzip -d ERR1198913_2.fastq.gz
jv~$ gzip -d ERR1198914_1.fastq.gz
jv~$ gzip -d ERR1198914_2.fastq.gz
```

### Check number of reads

#### Sample 1
```sh
jv~$ head -1 ERR1198911_1.fastq
# check the header "@ERR1198911.1"

jv~$ grep -c "@ERR1198911" ERR1198911_1.fastq
1408891

jv~$ grep -c "@ERR1198911" ERR1198911_2.fastq
1408891
```

#### Sample 2
```sh
jv~$ head -1 ERR1198912_1.fastq
# @ERR1198912.1

jv~$ grep -c "@ERR1198912" ERR1198912_1.fastq
597510
jv~$ grep -c "@ERR1198912" ERR1198912_2.fastq
597510
```

#### Sample 3
```sh
jv~$ grep -c "@ERR1198913" ERR1198913_1.fastq
700716
jv~$ grep -c "@ERR1198913" ERR1198913_2.fastq
700716
```

#### Sample 4
```sh
jv~$ grep -c "@ERR1198914" ERR1198914_1.fastq
577022
jv~$ grep -c "@ERR1198914" ERR1198914_2.fastq
577022
```


#### Summary of reads
Generate some stats of the samples in `R`

```sh
jv~$ R
```

```R
x <- c(1408891, 597510, 700716, 577022)
mean(x)
[1] 821034.8

median(x)
[1] 649113

sd(x)
[1] 395624.9

quit()
```

|#|Sample    |reads  |
|-|----------|-------|
|1|ERR1198911|1408891|
|2|ERR1198912| 597510|
|3|ERR1198913| 700716|
|4|ERR1198914| 577022|
|**Average**|**821034.8**|**~0.82 M reads**|
|Std.Dev|395624.9|~0.4 M reads|

> Low depth of reads (current datasets are ~ 20-90 M reads)

## [Stage 1] Quality Control (QC)

### Quality Filtering

> This protocol section is partially based on: http://merenlab.org/tutorials/assembly-based-metagenomics/

Generate a TAB-delimited `samples.txt` file to point out where are your raw `R1` and `R2` files for each sample.

|sample    |r1                         |r2                          |
|----------|---------------------------|----------------------------|
|Sample_01 | 00_RAW/ERR1198911_1.fastq | 00_RAW/ERR1198911_2.fastq  |
|Sample_02 | 00_RAW/ERR1198912_1.fastq | 00_RAW/ERR1198912_2.fastq  |
|Sample_03 | 00_RAW/ERR1198913_1.fastq | 00_RAW/ERR1198913_2.fastq  |
|Sample_04 | 00_RAW/ERR1198914_1.fastq | 00_RAW/ERR1198914_2.fastq  |


Create a directory for quality-filtered `R1` and `R2`
```sh
jv~$ mkdir 01_QC

jv~$ source activate illumina-utils

jv~$ iu-gen-configs samples.txt -o 01_QC

Report .......................................: Read for 4 samples is read
Output directory set in configs ..............: /disk/rdisk09/jvilladaa2/metagenomics-movile_cave/01_QC
Prefix for R1 ................................: None
Prefix for R2 ................................: None

jv~$ ls -1 01_QC/
Sample_01.ini
Sample_02.ini
Sample_03.ini
Sample_04.ini
```

Run quality filtering for all your samples at once:

> Note that we had to include the `parameter` `--ignore-deflines`
```sh
jv~$ for ini in 01_QC/*.ini; do iu-filter-quality-minoche --ignore-deflines $ini; done

 110% -- (num pairs processed: 1,408,000)
Read ID tracker dict is being stored ...
 108% -- (num pairs processed: 597,000)
Read ID tracker dict is being stored ...
 106% -- (num pairs processed: 700,000)
Read ID tracker dict is being stored ...
 104% -- (num pairs processed: 577,000)
Read ID tracker dict is being stored ...
```

Ref to the method used:

> Minoche, A. E., Dohm, J. C., & Himmelbauer, H. (2011). Evaluation of genomic high-throughput sequencing data generated on Illumina HiSeq and genome analyzer systems. *Genome biology*, 12(11), R112.

The contents of the `01_QC/` directory should look like this:
```sh
jv~$ cd ..
jv~$ ls 01_QC/
Sample_01.ini
Sample_01-QUALITY_PASSED_R1.fastq
Sample_01-QUALITY_PASSED_R2.fastq
Sample_01-READ_IDs.cPickle.z
Sample_01-STATS.txt
Sample_02.ini
Sample_02-QUALITY_PASSED_R1.fastq
Sample_02-QUALITY_PASSED_R2.fastq
Sample_02-READ_IDs.cPickle.z
Sample_02-STATS.txt
Sample_03.ini
Sample_03-QUALITY_PASSED_R1.fastq
Sample_03-QUALITY_PASSED_R2.fastq
Sample_03-READ_IDs.cPickle.z
Sample_03-STATS.txt
Sample_04.ini
Sample_04-QUALITY_PASSED_R1.fastq
Sample_04-QUALITY_PASSED_R2.fastq
Sample_04-READ_IDs.cPickle.z
Sample_04-STATS.txt
```

#### Checking QC results stats in `*_STATS.txt` files:

##### Sample_01-STATS
```sh
jv~$ cat 01_QC/Sample_01-STATS.txt

number of pairs analyzed      : 1408891
total pairs passed            : 1306751 (%92.75 of all pairs)
  total pair_1 trimmed        : 0 (%0.00 of all passed pairs)
  total pair_2 trimmed        : 0 (%0.00 of all passed pairs)
total pairs failed            : 102140 (%7.25 of all pairs)
  pairs failed due to pair_1  : 13387 (%13.11 of all failed pairs)
  pairs failed due to pair_2  : 79382 (%77.72 of all failed pairs)
  pairs failed due to both    : 9371 (%9.17 of all failed pairs)
  FAILED_REASON_P             : 897 (%0.88 of all failed pairs)
  FAILED_REASON_C33           : 101242 (%99.12 of all failed pairs)
  FAILED_REASON_N             : 1 (%0.00 of all failed pairs)
```

##### Sample_02-STATS
```sh
jv~$ cat 01_QC/Sample_02-STATS.txt

number of pairs analyzed      : 597510
total pairs passed            : 537816 (%90.01 of all pairs)
  total pair_1 trimmed        : 0 (%0.00 of all passed pairs)
  total pair_2 trimmed        : 0 (%0.00 of all passed pairs)
total pairs failed            : 59694 (%9.99 of all pairs)
  pairs failed due to pair_1  : 5939 (%9.95 of all failed pairs)
  pairs failed due to pair_2  : 48536 (%81.31 of all failed pairs)
  pairs failed due to both    : 5219 (%8.74 of all failed pairs)
  FAILED_REASON_P             : 648 (%1.09 of all failed pairs)
  FAILED_REASON_N             : 1 (%0.00 of all failed pairs)
  FAILED_REASON_C33           : 59045 (%98.91 of all failed pairs)
```

##### Sample_03-STATS
```sh
jv~$ cat 01_QC/Sample_03-STATS.txt

number of pairs analyzed      : 700716
total pairs passed            : 622470 (%88.83 of all pairs)
  total pair_1 trimmed        : 0 (%0.00 of all passed pairs)
  total pair_2 trimmed        : 0 (%0.00 of all passed pairs)
total pairs failed            : 78246 (%11.17 of all pairs)
  pairs failed due to pair_1  : 5800 (%7.41 of all failed pairs)
  pairs failed due to pair_2  : 66422 (%84.89 of all failed pairs)
  pairs failed due to both    : 6024 (%7.70 of all failed pairs)
  FAILED_REASON_C33           : 78119 (%99.84 of all failed pairs)
  FAILED_REASON_P             : 127 (%0.16 of all failed pairs)
```

##### Sample_04-STATS
```sh
jv~$ cat 01_QC/Sample_04-STATS.txt

number of pairs analyzed      : 577022
total pairs passed            : 502078 (%87.01 of all pairs)
  total pair_1 trimmed        : 0 (%0.00 of all passed pairs)
  total pair_2 trimmed        : 0 (%0.00 of all passed pairs)
total pairs failed            : 74944 (%12.99 of all pairs)
  pairs failed due to pair_1  : 4864 (%6.49 of all failed pairs)
  pairs failed due to pair_2  : 64403 (%85.93 of all failed pairs)
  pairs failed due to both    : 5677 (%7.57 of all failed pairs)
  FAILED_REASON_C33           : 74800 (%99.81 of all failed pairs)
  FAILED_REASON_P             : 144 (%0.19 of all failed pairs)
```

##### Deactivate conda environment

```sh
source deactivate
```

##### Summary Stage 1 - QC Filtering

```sh
jv~$ grep 'total pairs passed' 01_QC/*STATS.txt
```

|#|Sample    |number of pairs analyzed   | total pairs passed  |
|-|----------|---------------------------|---------------------|
|1|ERR1198911|1408891                    |1306751 (%92.75)     |
|2|ERR1198912|597510                     | 537816 (%90.01)     |
|3|ERR1198913|700716                     | 622470 (%88.83)     |
|4|ERR1198914|577022                     | 502078 (%87.01)     |

## [Stage 2] Co-assembly

Generating `contigs.fasta`

> This protocol section is partially based on: http://merenlab.org/tutorials/assembly-based-metagenomics/

Note that no conda environment is activated at this stage, and the `python` version in which these scripts are working is `Python 2.7.14`.

List your filtered reads
```sh
jv~$ ls -1 01_QC/*fastq

01_QC/Sample_01-QUALITY_PASSED_R1.fastq
01_QC/Sample_01-QUALITY_PASSED_R2.fastq
01_QC/Sample_02-QUALITY_PASSED_R1.fastq
01_QC/Sample_02-QUALITY_PASSED_R2.fastq
01_QC/Sample_03-QUALITY_PASSED_R1.fastq
01_QC/Sample_03-QUALITY_PASSED_R2.fastq
01_QC/Sample_04-QUALITY_PASSED_R1.fastq
01_QC/Sample_04-QUALITY_PASSED_R2.fastq
```

Create two environment variables:
```sh
jv~$ R1s=`ls 01_QC/*QUALITY_PASSED_R1* | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`

jv~$ R2s=`ls 01_QC/*QUALITY_PASSED_R2* | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`
```

The environment variables should look like this:

```sh
jv~$ echo $R1s

01_QC/Sample_01-QUALITY_PASSED_R1.fastq,01_QC/Sample_02-QUALITY_PASSED_R1.fastq,01_QC/Sample_03-QUALITY_PASSED_R1.fastq,01_QC/Sample_04-QUALITY_PASSED_R1.fastq

jv~$ echo $R2s

01_QC/Sample_01-QUALITY_PASSED_R2.fastq,01_QC/Sample_02-QUALITY_PASSED_R2.fastq,01_QC/Sample_03-QUALITY_PASSED_R2.fastq,01_QC/Sample_04-QUALITY_PASSED_R2.fastq
```

### Run MEGAHIT

The code has the form:

```sh
jv~$ NUM_THREADS=32
jv~$ echo $NUM_THREADS
32

jv~$ MIN_CONTIG_SIZE=1000
jv~$ echo $MIN_CONTIG_SIZE
1000

jv~$ megahit -1 $R1s -2 $R2s --min-contig-len $MIN_CONTIG_SIZE -m 0.95 -o 02_ASSEMBLY/ -t $NUM_THREADS

251.0Gb memory in total.
Using: 239.057Gb.
MEGAHIT v1.1.2
.
.
.

--- [STAT] 16848 contigs, total 32383346 bp, min 1000 bp, max 95199 bp, avg 1922 bp, N50 1845 bp
--- ALL DONE. Time elapsed: 5555.232126 seconds ---
```
_**Time elapsed: 5555.232126 seconds = ~1.5 hours**_

### Refining contigs with `anvio`

```sh
jv~$ mkdir 03_CONTIGS
jv~$ source activate anvio3
jv~$ anvi-script-reformat-fasta 02_ASSEMBLY/final.contigs.fa -o 03_CONTIGS/contigs.fa --min-len 1000 --simplify-names --report name_conversions_1000.txt

Input ........................................: 02_ASSEMBLY/final.contigs.fa
Output .......................................: 03_CONTIGS/contigs.fa
Minimum length ...............................: 1,000
Total num contigs ............................: 16,848
Total num nucleotides ........................: 32,383,346
Contigs removed ..............................: 0 (0.00% of all)
Nucleotides removed ..........................: 0 (0.00% of all)
Deflines simplified ..........................: True
```

Testing for `--min-len 1500`

```sh
jv~$ anvi-script-reformat-fasta 02_ASSEMBLY/final.contigs.fa -o 03_CONTIGS/contigs_1500.fa --min-len 1500 --simplify-names --report name_conversions_1500.txt
Input ........................................: 02_ASSEMBLY/final.contigs.fa
Output .......................................: 03_CONTIGS/contigs_1500.fa
Minimum length ...............................: 1,500
Total num contigs ............................: 16,848
Total num nucleotides ........................: 32,383,346
Contigs removed ..............................: 10162 (60.32% of all)
Nucleotides removed ..........................: 12134435 (37.47% of all)
Deflines simplified ..........................: True
```

Testing for `--min-len 2500`
```sh
jv~$ anvi-script-reformat-fasta 02_ASSEMBLY/final.contigs.fa -o 03_CONTIGS/contigs_2500.fa --min-len 2500 --simplify-names --report name_conversions_2000.txt
Input ........................................: 02_ASSEMBLY/final.contigs.fa
Output .......................................: 03_CONTIGS/contigs_2500.fa
Minimum length ...............................: 2,500
Total num contigs ............................: 16,848
Total num nucleotides ........................: 32,383,346
Contigs removed ..............................: 14702 (87.26% of all)
Nucleotides removed ..........................: 20596179 (63.60% of all)
Deflines simplified ..........................: True
```

We will proceed with `03_CONTIGS/contigs.fa` which has **`--min-len` = 1000**

##### Deactivate the `anvio` conda environment

```sh
jv~$ source deactivate
```

## [Stage 3] Binning

### Binning with `MaxBin v2.2.4`

```sh
jv~$ mkdir 04_BINNING

```

Test `Maxbin2` installation

```sh
jv~$ run_MaxBin.pl -v
MaxBin 2.2.4
```

Create `reads_list.txt` to be passed to `MaxBin2`:

```sh
jv~$ ls 01_QC/*.fastq > 04_BINNING/reads_list.txt
```

Corroborate the `reads_list.txt` content:

```sh
jv~$ cat 04_BINNING/reads_list.txt

01_QC/Sample_01-QUALITY_PASSED_R1.fastq
01_QC/Sample_01-QUALITY_PASSED_R2.fastq
01_QC/Sample_02-QUALITY_PASSED_R1.fastq
01_QC/Sample_02-QUALITY_PASSED_R2.fastq
01_QC/Sample_03-QUALITY_PASSED_R1.fastq
01_QC/Sample_03-QUALITY_PASSED_R2.fastq
01_QC/Sample_04-QUALITY_PASSED_R1.fastq
01_QC/Sample_04-QUALITY_PASSED_R2.fastq
```

Run MaxBin binning on the assembled contigs:

> Note that here we are indicating `-thread 32`

```sh
jv~$ run_MaxBin.pl -contig 03_CONTIGS/contigs.fa -reads_list 04_BINNING/reads_list.txt -out 04_BINNING/MOVILE_CAVE -thread 32
========== Job finished ==========
Yielded 11 bins for contig (scaffold) file 03_CONTIGS/contigs.fa
...
Summary file: 04_BINNING/MOVILE_CAVE.summary
Genome abundance info file: 04_BINNING/MOVILE_CAVE.abundance
Marker counts: 04_BINNING/MOVILE_CAVE.marker
Marker genes for each bin: 04_BINNING/MOVILE_CAVE.marker_of_each_gene.tar.gz
Bin files: 04_BINNING/MOVILE_CAVE.001.fasta - 04_BINNING/MOVILE_CAVE.011.fasta
Unbinned sequences: 04_BINNING/MOVILE_CAVE.noclass
...
========== Elapsed Time ==========
0 hours 10 minutes and 6 seconds.
```
**Elapsed time: ~11 minutes**

## [Stage 4] Refine bins

CheckM

```sh
jv~$ source activate CheckM
jv~$ checkm lineage_wf -t 32 -x fasta 04_BINNING/ 05_CHECKM

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Bin Id                     Marker lineage            # genomes   # markers   # marker sets    0     1    2    3    4   5+   Completeness   Contamination   Strain heterogeneity
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  MOVILE_CAVE.009   c__Betaproteobacteria (UID3959)       235         418           210         7    321   87   3    0   0       97.95           24.20              50.00
  MOVILE_CAVE.003   c__Betaproteobacteria (UID3888)       323         386           233         91   190   80   22   2   1       80.29           30.63              10.12
  MOVILE_CAVE.001   c__Alphaproteobacteria (UID3305)      564         349           230         75   256   18   0    0   0       80.00            5.10               0.00
  MOVILE_CAVE.010         k__Bacteria (UID203)            5449        104            58         13    36   43   10   2   0       78.45           40.99              52.94
  MOVILE_CAVE.004         k__Bacteria (UID203)            5449        104            58         20    40   32   10   2   0       70.53           34.50               2.70
  MOVILE_CAVE.006        k__Bacteria (UID2565)            2921        152            93         73    54   22   1    2   0       53.73           27.71               8.11
  MOVILE_CAVE.007         k__Bacteria (UID203)            5449        104            58         38    34   25   5    1   1       51.61           22.52               3.57
  MOVILE_CAVE.011         k__Bacteria (UID203)            5449        104            58         75    21   7    1    0   0       45.69           15.52              40.00
  MOVILE_CAVE.002   c__Gammaproteobacteria (UID4267)      119         544           284        337   185   18   4    0   0       37.99            3.74               3.33
  MOVILE_CAVE.005   c__Betaproteobacteria (UID3888)       323         387           234        234   128   21   3    1   0       34.57            3.80               2.78
  MOVILE_CAVE.008        k__Bacteria (UID3187)            2258        181           110        129    50   2    0    0   0       24.38            1.36              50.00
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source deactivate
```
**> Time: ~9 minutes**

One high quality bin (complete ≥ 70% and contam ≤10%) was retrieved: `MOVILE_CAVE.001`. We will also include `MOVILE_CAVE.009`, `MOVILE_CAVE.011`, `MOVILE_CAVE.002` and `MOVILE_CAVE.005` in further analyses.

## [Stage 5] Functional annotation - `prokka`

```sh
jv~$ mkdir 06_ANNOTATION
jv~$ prokka --outdir 06_ANNOTATION/MOVILE_CAVE.001 --prefix MOVILE_CAVE.001 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/MOVILE_CAVE.001.fasta
jv~$ prokka --outdir 06_ANNOTATION/MOVILE_CAVE.009 --prefix MOVILE_CAVE.009 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/MOVILE_CAVE.009.fasta
jv~$ prokka --outdir 06_ANNOTATION/MOVILE_CAVE.011 --prefix MOVILE_CAVE.011 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/MOVILE_CAVE.011.fasta
jv~$ prokka --outdir 06_ANNOTATION/MOVILE_CAVE.002 --prefix MOVILE_CAVE.002 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/MOVILE_CAVE.002.fasta
jv~$ prokka --outdir 06_ANNOTATION/MOVILE_CAVE.005 --prefix MOVILE_CAVE.005 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/MOVILE_CAVE.005.fasta
```
## [Stage 6] Taxonomy classification of bins

**Input of `phylophlan` should be bins in `.faa` format**

Create a directory to save the results
```sh
jv~$ mkdir 07_TAXONOMY/
```

Copy the bins in amoni acids `.faa` format from the `prokka` output:
```sh
jv~$ cp 06_ANNOTATION/MOVILE_CAVE.001/MOVILE_CAVE.001.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/MOVILE_CAVE.009/MOVILE_CAVE.009.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/MOVILE_CAVE.011/MOVILE_CAVE.011.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/MOVILE_CAVE.002/MOVILE_CAVE.002.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/MOVILE_CAVE.005/MOVILE_CAVE.005.faa 07_TAXONOMY/

jv~$ ls 07_TAXONOMY/

MOVILE_CAVE.001.faa  MOVILE_CAVE.005.faa  MOVILE_CAVE.011.faa  MOVILE_CAVE.002.faa  MOVILE_CAVE.009.faa
```

Go to the `phylophlan` directory (`$DIR_PHYLOPHLAN`) and create a new folder inside the `input` folder containing the bins in `.faa` format:

```sh
jv~$ cd $DIR_PHYLOPHLAN
jv~$ mkdir input/07_TAXONOMY-MOVILE_CAVE/
jv~$ cp $DIR_PROJECT/07_TAXONOMY/* input/07_TAXONOMY-MOVILE_CAVE/
jv~$ ls input/07_TAXONOMY-MOVILE_CAVE/
MOVILE_CAVE.001.faa  MOVILE_CAVE.002.faa  MOVILE_CAVE.005.faa  MOVILE_CAVE.009.faa  MOVILE_CAVE.011.faa
```

```sh
jv~$ source activate phylophlan
jv~$ ./phylophlan.py -i -t 07_TAXONOMY-MOVILE_CAVE --nproc 32
jv~$ source deactivate
jv~$ ls -1 input/07_TAXONOMY-MOVILE_CAVE/
```

Results in `imputed_conf_high-conf.txt`
```sh
jv~$ cat output/07_TAXONOMY-MOVILE_CAVE/imputed_conf_high-conf.txt
MOVILE_CAVE.009	d__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Gallionellales.f__Gallionellaceae.g__?.s__?.t__?
MOVILE_CAVE.005	d__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Methylophilales.f__Methylophilaceae.g__?.s__?.t__?
MOVILE_CAVE.011	d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__?.f__?.g__?.s__?.t__?
```

Results in `imputed_conf_low-conf.txt`
```sh
jv~$ cat output/07_TAXONOMY-MOVILE_CAVE/imputed_conf_low_conf.txt
MOVILE_CAVE.002	d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__?.f__?.g__?.s__?.t__?
MOVILE_CAVE.001	d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__?.f__?.g__?.s__?.t__?
```
**Total time: ~ 1 h**

Copy the `phylophlan` results from `$DIR_PHYLOPHLAN` to the `$DIR_PROJECT`

```sh 
jv~$ cp output/07_TAXONOMY-MOVILE_CAVE/* ~/data/metagenomics-movile_cave/07_TAXONOMY/
```

## [Stage 7]  Refine functional annotation - `EggNOG-mapper`
`EggNOG-mapper v4.5.1` works on `python2.7`

`$DIR_EGGNOG` is the directory to `eggnog-mapper/emapper.py` script.

```sh
jv~$ mkdir 08_ANNOTATION_REFINE

jv~$ cd 08_ANNOTATION_REFINE

jv~$ python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/MOVILE_CAVE.009.faa --output MOVILE_CAVE.009 -m diamond --cpu 32
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/MOVILE_CAVE.005.faa --output MOVILE_CAVE.005 -m diamond --cpu 32
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/MOVILE_CAVE.011.faa --output MOVILE_CAVE.011 -m diamond --cpu 32
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/MOVILE_CAVE.002.faa --output MOVILE_CAVE.002 -m diamond --cpu 32
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/MOVILE_CAVE.001.faa --output MOVILE_CAVE.001 -m diamond --cpu 32
```