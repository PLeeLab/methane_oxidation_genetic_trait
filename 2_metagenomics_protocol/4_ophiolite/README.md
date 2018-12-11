# Metagenomics from Santa Elena Ophiolite alkaline spring, Costa Rica

<!-- TOC -->

- [Metagenomics from Santa Elena Ophiolite alkaline spring, Costa Rica](#metagenomics-from-santa-elena-ophiolite-alkaline-spring-costa-rica)
    - [[Stage 0] Preliminary](#stage-0-preliminary)
    - [[Stage 1] Quality Control (QC)](#stage-1-quality-control-qc)
    - [[Stage 2] Co-assembly (pooled-assembly)](#stage-2-co-assembly-pooled-assembly)
    - [[Stage 3] Binning](#stage-3-binning)
    - [[Stage 4] Refine bins](#stage-4-refine-bins)
    - [[Stage 5] Functional annotation - `prokka`](#stage-5-functional-annotation---prokka)
    - [[Stage 6] Taxonomy classification of bins](#stage-6-taxonomy-classification-of-bins)
    - [[Stage 7]  Refine functional annotation of methanotrophic bacteria using `EggNOG-mapper`](#stage-7--refine-functional-annotation-of-methanotrophic-bacteria-using-eggnog-mapper)

<!-- /TOC -->

## [Stage 0] Preliminary
### Create project directory
```sh
jv~$ mkdir metagenomics-ophiolite/
jv~$ cd metagenomics-ophiolite/
```

### Retrieve sequences from SRA NCBI system

> This dataset was published in:
> Crespo-Medina, M., Twing, K. I., Sánchez-Murillo, R., Brazelton, W. J., McCollom, T. M., & Schrenk, M. O. (2017). Methane Dynamics in a Tropical Serpentinizing Environment: The Santa Elena Ophiolite, Costa Rica. Frontiers in microbiology, 8, 916.
> https://doi.org/10.3389/fmicb.2017.00916

### Download metagenomics data

The data of the project are in the NCBI at https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA340462

This code will retrieve the reads in separate ``FASTQ`` files:

```sh
jv~$ mkdir 00_RAW
jv~$ cd 00_RAW/
jv~$ fastq-dump -I --split-files SRR4101184 SRR4101185
jv~$ ls -1
```

### Check number of reads

#### Sample 1
```sh
jv~$ head -1 SRR4101184_1.fastq
jv~$ grep -c "@SRR4101184" SRR4101184_1.fastq
45476317
jv~$ grep -c "@SRR4101184" SRR4101184_2.fastq
45476317
```

#### Sample 2
```sh
jv~$ grep -c "@SRR4101185" SRR4101185_1.fastq
60780480
jv~$ grep -c "@SRR4101185" SRR4101185_2.fastq
60780480
```

```sh
jv~$ R
```

```R
x <- c(45476317, 60780480)
mean(x)
[1] 53128398

sd(x)
[1] 10821677

quit()
```

|#|Sample    |reads  |
|-|----------|-------|
|1|SRR4101184|53128398|
|2|SRR4101185|10821677|
|**Average**|**53128398**|**~53 M reads**|
|Std.Dev|10821677|~10 M reads|

> Sufficient number of reads (cat feces 20-90 M reads)

## [Stage 1] Quality Control (QC)

### Quality Filtering

> Paired-end  sequencing  with a 100 cycle Illumina HiSeq run generated partial ∼30 bp overlaps. 
> This protocol section is partially based on: http://merenlab.org/tutorials/assembly-based-metagenomics/

Generate a TAB-delimited `samples.txt` file to point out where are your raw `R1` and `R2` files for each sample.

|sample    |r1                         |r2                          |
|----------|---------------------------|----------------------------|
|Sample_01 | 00_RAW/SRR4101184_1.fastq | 00_RAW/SRR4101184_2.fastq  |
|Sample_02 | 00_RAW/SRR4101185_1.fastq | 00_RAW/SRR4101185_2.fastq  |


Create a directory for quality-filtered `R1` and `R2`
```sh
jv~$ cd ..
jv~$ mkdir 01_QC/
jv~$ source activate illumina-utils
jv~$ iu-gen-configs samples.txt -o 01_QC
Report .......................................: Read for 2 samples is read
Output directory set in configs ..............: /disk/rdisk09/jvilladaa2/metagenomics-ophiolite/01_QC
Prefix for R1 ................................: None
Prefix for R2 ................................: Nonee
jv~$ ls 01_QC/
Sample_01.ini  Sample_02.ini
```

Run quality filtering for all your samples at once:

Note that we are applying here the method by Minoche et al. and we also include the `parameter` `--ignore-deflines`

```sh
jv~$ for ini in 01_QC/*.ini; do iu-filter-quality-minoche --ignore-deflines $ini; done
 96% -- (num pairs processed: 45,476,000)
Read ID tracker dict is being stored ...
 96% -- (num pairs processed: 60,780,000)
Read ID tracker dict is being stored ...
```

Ref to the method used:

> Minoche, A. E., Dohm, J. C., & Himmelbauer, H. (2011). Evaluation of genomic high-throughput sequencing data generated on Illumina HiSeq and genome analyzer systems. *Genome biology*, 12(11), R112.

The contents of the `01_QC/` directory should look like this:
```sh
jv~$ ls 01_QC/
Sample_01.ini                      Sample_02.ini
Sample_01-QUALITY_PASSED_R1.fastq  Sample_02-QUALITY_PASSED_R1.fastq
Sample_01-QUALITY_PASSED_R2.fastq  Sample_02-QUALITY_PASSED_R2.fastq
Sample_01-READ_IDs.cPickle.z       Sample_02-READ_IDs.cPickle.z
Sample_01-STATS.txt                Sample_02-STATS.txt
```

#### Checking QC results stats in `*_STATS.txt` files:

##### Sample_01-STATS
```sh
jv~$ cat 01_QC/Sample_01-STATS.txt
number of pairs analyzed      : 45476317
total pairs passed            : 40308912 (%88.64 of all pairs)
  total pair_1 trimmed        : 7573622 (%18.79 of all passed pairs)
  total pair_2 trimmed        : 7370258 (%18.28 of all passed pairs)
total pairs failed            : 5167405 (%11.36 of all pairs)
  pairs failed due to pair_1  : 1567375 (%30.33 of all failed pairs)
  pairs failed due to pair_2  : 2177160 (%42.13 of all failed pairs)
  pairs failed due to both    : 1422870 (%27.54 of all failed pairs)
  FAILED_REASON_N             : 145516 (%2.82 of all failed pairs)
  FAILED_REASON_C33           : 1724261 (%33.37 of all failed pairs)
  FAILED_REASON_P             : 3297628 (%63.82 of all failed pairs)
```

##### Sample_02-STATS
```sh
jv~$ cat 01_QC/Sample_02-STATS.txt
number of pairs analyzed      : 60780480
total pairs passed            : 51813498 (%85.25 of all pairs)
  total pair_1 trimmed        : 10549832 (%20.36 of all passed pairs)
  total pair_2 trimmed        : 10897854 (%21.03 of all passed pairs)
total pairs failed            : 8966982 (%14.75 of all pairs)
  pairs failed due to pair_1  : 3080989 (%34.36 of all failed pairs)
  pairs failed due to pair_2  : 2960630 (%33.02 of all failed pairs)
  pairs failed due to both    : 2925363 (%32.62 of all failed pairs)
  FAILED_REASON_N             : 175857 (%1.96 of all failed pairs)
  FAILED_REASON_P             : 6059759 (%67.58 of all failed pairs)
  FAILED_REASON_C33           : 2731366 (%30.46 of all failed pairs)
```

Deactivate conda environment

```sh
source deactivate
```

##### Summary Stage 1 - QC Filtering

```sh
jv~$ grep 'total pairs passed' 01_QC/*STATS.txt
```

|#|Sample    |number of pairs analyzed   | total pairs passed  |
|-|----------|---------------------------|---------------------|
|1|SRR4101184|45476317                   |40308912 (%88.64)    |
|2|SRR4101185|60780480                   |51813498 (%85.25)    |

## [Stage 2] Co-assembly (pooled-assembly)

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
```

Create two environment variables:
```sh
jv~$ R1s=`ls 01_QC/*QUALITY_PASSED_R1* | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`

jv~$ R2s=`ls 01_QC/*QUALITY_PASSED_R2* | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`
```

The environment variables should look like this:

```sh
jv~$ echo $R1s
01_QC/Sample_01-QUALITY_PASSED_R1.fastq,01_QC/Sample_02-QUALITY_PASSED_R1.fastq

jv~$ echo $R2s
01_QC/Sample_01-QUALITY_PASSED_R2.fastq,01_QC/Sample_02-QUALITY_PASSED_R2.fastq
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
--- [STAT] 76785 contigs, total 223941517 bp, min 1000 bp, max 405231 bp, avg 2916 bp, N50 3916 bp
--- [Wed Feb 14 20:08:35 2018] ALL DONE. Time elapsed: 16063.311761 seconds ---
```
**Time elapsed: ~5 hours**

### Refining contigs with `anvio`

```sh
jv~$ mkdir 03_CONTIGS/
jv~$ source activate anvio3
jv~$ anvi-script-reformat-fasta 02_ASSEMBLY/final.contigs.fa -o 03_CONTIGS/contigs.fa --min-len 1000 --simplify-names --report name_conversions_1000.txt
Input ........................................: 02_ASSEMBLY/final.contigs.fa
Output .......................................: 03_CONTIGS/contigs.fa
Minimum length ...............................: 1,000
Total num contigs ............................: 76,785
Total num nucleotides ........................: 223,941,517
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
Total num contigs ............................: 76,785
Total num nucleotides ........................: 223,941,517
Contigs removed ..............................: 35003 (45.59% of all)
Nucleotides removed ..........................: 42234776 (18.86% of all)
Deflines simplified ..........................: True
```

Testing for `--min-len 2500`
```sh
jv~$ anvi-script-reformat-fasta 02_ASSEMBLY/final.contigs.fa -o 03_CONTIGS/contigs_2500.fa --min-len 2500 --simplify-names --report name_conversions_2500.txt

Input ........................................: 02_ASSEMBLY/final.contigs.fa
Output .......................................: 03_CONTIGS/contigs_2500.fa
Minimum length ...............................: 2,500
Total num contigs ............................: 94,436
Total num nucleotides ........................: 256,767,472
Contigs removed ..............................: 70325 (74.47% of all)
Nucleotides removed ..........................: 103004121 (40.12% of all)
Deflines simplified ..........................: True
```

We will proceed with `03_CONTIGS/contigs.fa` which has **`--min-len` = 1000**

Deactivate the `anvio` conda environment
```sh
jv~$ source deactivate
```

## [Stage 3] Binning

### Binning with `MaxBin v2.2.4`

```sh
jv~$ mkdir 04_BINNING/

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
```
Run MaxBin binning on the assembled contigs:

> Note that here we are indicating `-thread 32`

```sh
jv~$ run_MaxBin.pl -contig 03_CONTIGS/contigs.fa -reads_list 04_BINNING/reads_list.txt -out 04_BINNING/OPHIOLITE -thread 32
========== Job finished ==========
Yielded 94 bins for contig (scaffold) file 03_CONTIGS/contigs.fa

Here are the output files for this run.
Please refer to the README file for further details.

Summary file: 04_BINNING/OPHIOLITE.summary
Genome abundance info file: 04_BINNING/OPHIOLITE.abundance
Marker counts: 04_BINNING/OPHIOLITE.marker
Marker genes for each bin: 04_BINNING/OPHIOLITE.marker_of_each_gene.tar.gz
Bin files: 04_BINNING/OPHIOLITE.001.fasta - 04_BINNING/OPHIOLITE.094.fasta
Unbinned sequences: 04_BINNING/OPHIOLITE.noclass
.
.
.
========== Elapsed Time ==========
2 hours 7 minutes and 37 seconds.
```
**Elapsed time: ~2 hours**

## [Stage 4] Refine bins

CheckM

```sh
jv~$ source activate CheckM
jv~$ checkm lineage_wf -t 32 -x fasta 04_BINNING/ 05_CHECKM
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Bin Id                   Marker lineage            # genomes   # markers   # marker sets    0     1     2    3    4    5+   Completeness   Contamination   Strain heterogeneity
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  OPHIOLITE.040        k__Bacteria (UID2569)            434         278           186         2    273    3    0    0    0       99.19            1.61               0.00
  OPHIOLITE.085         k__Bacteria (UID203)            5449        104            58         1     45    49   9    0    0       99.14           35.31               0.00
  OPHIOLITE.059        k__Bacteria (UID2569)            434         278           186         2    266    10   0    0    0       98.92            4.03               0.00
  OPHIOLITE.010      p__Bacteroidetes (UID2591)         364         303           203         4    293    6    0    0    0       98.52            2.22              16.67
  OPHIOLITE.091      p__Cyanobacteria (UID2182)          89         544           424         7    535    2    0    0    0       98.35            0.08               0.00
  OPHIOLITE.084         k__Bacteria (UID203)            5449        103            58         1     79    19   4    0    0       98.28           28.10               0.00
  OPHIOLITE.076         k__Bacteria (UID203)            5449        103            57         1     64    28   8    1    1       98.25           36.44               0.00
  OPHIOLITE.050         k__Bacteria (UID203)            5449        104            58         2     48    26   25   3    0       97.41           37.07               9.24
  OPHIOLITE.029     o__Burkholderiales (UID4000)        193         427           214         20   386    20   1    0    0       97.16            7.29              78.26
  OPHIOLITE.013     o__Actinomycetales (UID1572)        580         286           171         5    273    8    0    0    0       97.08            1.87              62.50
  OPHIOLITE.086     p__Proteobacteria (UID3880)         1495        261           164         14   229    18   0    0    0       96.99            8.91              11.11
  OPHIOLITE.045      p__Bacteroidetes (UID2605)         350         316           210         7    289    20   0    0    0       96.90            6.17              45.00
  OPHIOLITE.011         k__Bacteria (UID203)            5449        104            58         7     59    34   2    2    0       96.63           25.49              53.85
  OPHIOLITE.049         k__Bacteria (UID203)            5449        104            58         2     78    15   9    0    0       96.55           20.98               2.38
  OPHIOLITE.015         k__Bacteria (UID203)            5449        103            57         2     73    26   2    0    0       96.49           13.16              34.38
  OPHIOLITE.047         k__Bacteria (UID203)            5449        104            58         3     54    45   2    0    0       94.83           18.81               3.92
  OPHIOLITE.094        k__Bacteria (UID2570)            433         273           183         14   255    4    0    0    0       94.54            1.14              25.00
  OPHIOLITE.005     o__Burkholderiales (UID4000)        193         427           214         38   339    47   3    0    0       93.89           12.19              69.64
  OPHIOLITE.009     o__Burkholderiales (UID4001)        108         570           250         36   367   107   48   11   1       93.87           42.38              47.40
  OPHIOLITE.008     o__Burkholderiales (UID4000)        193         427           214         37   301    84   4    1    0       93.78           22.08              41.18
  OPHIOLITE.062         k__Bacteria (UID203)            5449        104            58         5     72    27   0    0    0       93.10           23.75               3.70
  OPHIOLITE.014         k__Bacteria (UID203)            5449        104            58         9     38    48   8    1    0       92.79           71.65              70.51
  OPHIOLITE.080        k__Bacteria (UID2495)            2993        147            91         13   128    5    1    0    0       91.44            2.11              25.00
  OPHIOLITE.033        p__Firmicutes (UID241)           930         213           118         20   163    28   2    0    0       90.47           14.41              52.94
  OPHIOLITE.038         k__Bacteria (UID203)            5449        104            58         6     39    56   3    0    0       89.66           42.24               4.62
  OPHIOLITE.087         k__Bacteria (UID203)            5449        104            58         7     70    24   3    0    0       89.31           22.34               6.06
  OPHIOLITE.043        k__Bacteria (UID2329)            174         149            89         18   119    9    2    1    0       88.99           11.77               4.76
  OPHIOLITE.004     o__Burkholderiales (UID4000)        193         427           214         53   348    25   1    0    0       88.31            3.35              53.57
  OPHIOLITE.006         k__Bacteria (UID203)            5449        104            58         8     67    23   6    0    0       87.93           19.44              60.98
  OPHIOLITE.089         k__Bacteria (UID203)            5449        104            58         14    42    47   1    0    0       87.70           59.90              12.00
  OPHIOLITE.069     f__Spirochaetaceae (UID2535)         33         249           143         33   195    21   0    0    0       87.49            6.23              23.81
  OPHIOLITE.041     o__Burkholderiales (UID4000)        193         426           214         67   340    16   3    0    0       87.40            4.88              20.00
  OPHIOLITE.022     o__Actinomycetales (UID1572)        580         286           171         48   233    5    0    0    0       86.65            1.95              20.00
  OPHIOLITE.063   c__Gammaproteobacteria (UID4274)      112         580           289        106   421    49   4    0    0       86.29            9.82               8.20
  OPHIOLITE.056         k__Bacteria (UID203)            5449        104            58         11    47    37   7    2    0       86.21           45.86              52.86
  OPHIOLITE.075        k__Bacteria (UID2569)            434         278           186         53   140    81   4    0    0       86.17           34.63               5.38
  OPHIOLITE.083   c__Alphaproteobacteria (UID3422)       26         529           308         94   253   126   40   15   1       85.66           55.07               0.00
  OPHIOLITE.001     o__Burkholderiales (UID4000)        193         427           214         65   355    7    0    0    0       85.60            0.92              42.86
  OPHIOLITE.054       p__Euryarchaeota (UID49)           95         228           153         26   193    9    0    0    0       84.97            4.95              22.22
  OPHIOLITE.081         k__Bacteria (UID203)            5449        104            58         46    56    2    0    0    0       84.48            3.45               0.00
  OPHIOLITE.042         k__Bacteria (UID203)            5449        104            58         20    30    30   15   5    4       83.78           82.26               2.07
  OPHIOLITE.077         k__Bacteria (UID203)            5449        104            58         11    43    39   11   0    0       83.46           50.47               0.00
  OPHIOLITE.046         k__Bacteria (UID203)            5449        103            58         17    61    22   3    0    0       83.45           14.92               3.23
  OPHIOLITE.078        p__Firmicutes (UID241)           930         213           118         27   161    24   1    0    0       83.39           12.20              18.52
  OPHIOLITE.051         k__Bacteria (UID203)            5449        104            58         13    27    32   24   8    0       82.18           76.27               0.00
  OPHIOLITE.082         k__Bacteria (UID203)            5449        104            58         13    91    0    0    0    0       82.03            0.00               0.00
  OPHIOLITE.065         k__Bacteria (UID203)            5449        103            58         18    33    24   11   6    11      79.23           68.18               0.82
  OPHIOLITE.036         k__Bacteria (UID203)            5449        104            58         16    70    18   0    0    0       78.13           11.77               0.00
  OPHIOLITE.079         k__Bacteria (UID203)            5449        104            58         51    53    0    0    0    0       77.59            0.00               0.00
  OPHIOLITE.092      p__Bacteroidetes (UID2605)         350         316           210         79   220    14   3    0    0       76.74            4.48               0.00
  OPHIOLITE.012   c__Betaproteobacteria (UID3971)       223         425           211        120   277    26   2    0    0       76.06            7.69              65.62
  OPHIOLITE.037          k__Archaea (UID2)              207         149           107         28    60    38   18   4    1       75.93           59.38               2.38
  OPHIOLITE.068   c__Deltaproteobacteria (UID3216)       83         247           155         48   178    20   1    0    0       75.39            8.33               4.35
  OPHIOLITE.023    f__Rhodobacteraceae (UID3340)         84         568           330        133   392    42   1    0    0       73.72            7.64              13.33
  OPHIOLITE.007     o__Burkholderiales (UID4001)        108         570           250        157   337    73   3    0    0       73.10           13.08              37.80
  OPHIOLITE.052          k__Archaea (UID2)              207         149           107         47    81    16   4    1    0       69.61           15.89               0.00
  OPHIOLITE.053        p__Firmicutes (UID241)           930         213           118         56   135    22   0    0    0       69.34           11.11              77.27
  OPHIOLITE.028        k__Bacteria (UID2495)            2993        143            89         34    92    17   0    0    0       68.56            9.70              82.35
  OPHIOLITE.061         k__Bacteria (UID203)            5449        104            58         30    35    21   15   3    0       68.04           60.03               3.57
  OPHIOLITE.071        p__Firmicutes (UID239)           1324        175           101         58    80    33   4    0    0       67.67           23.72               6.67
  OPHIOLITE.066         k__Bacteria (UID203)            5449        104            58         33    44    13   10   4    0       63.89           23.98               0.00
  OPHIOLITE.093        k__Bacteria (UID2495)            2993        140            85         42    86    12   0    0    0       62.47            6.28              33.33
  OPHIOLITE.058         k__Bacteria (UID203)            5449        104            58         32    30    25   10   5    2       62.35           42.12               0.00
  OPHIOLITE.026     o__Burkholderiales (UID4000)        193         427           214        176   200    42   7    2    0       60.78           14.24               9.33
  OPHIOLITE.032         k__Bacteria (UID203)            5449        100            55         31    38    22   9    0    0       59.01           38.55               0.00
  OPHIOLITE.048        k__Bacteria (UID1452)            924         151           101         50    84    15   2    0    0       58.80           13.20               9.52
  OPHIOLITE.073       c__Clostridia (UID1085)            35         420           196        171   206    42   1    0    0       58.80           11.06              11.11
  OPHIOLITE.057         k__Bacteria (UID203)            5449        104            58         30    42    21   9    1    1       57.99           34.11               1.56
  OPHIOLITE.016         k__Bacteria (UID203)            5449        104            58         62    25    11   5    1    0       57.87           36.21               3.12
  OPHIOLITE.088        k__Bacteria (UID2495)            2993        143            89         66    61    14   2    0    0       56.72           17.42               0.00
  OPHIOLITE.031         k__Bacteria (UID203)            5449        103            57         39    32    26   6    0    0       54.68           29.73              11.36
  OPHIOLITE.074        k__Bacteria (UID2495)            2993        147            91         63    55    19   7    2    1       53.95           22.83               1.61
  OPHIOLITE.064         k__Bacteria (UID203)            5449        103            57         52    34    17   0    0    0       52.34           11.24               0.00
  OPHIOLITE.018   c__Betaproteobacteria (UID3971)       223         425           211        221   187    17   0    0    0       50.99            4.87              29.41
  OPHIOLITE.021     o__Burkholderiales (UID4000)        193         427           214        215   152    45   10   5    0       50.55           17.20               7.62
  OPHIOLITE.020         k__Bacteria (UID203)            5449        104            58         40    37    17   8    2    0       50.47           28.86              24.53
  OPHIOLITE.055        k__Bacteria (UID1452)            924         151           101         68    73    10   0    0    0       49.41            8.53              10.00
  OPHIOLITE.044        k__Bacteria (UID2329)            174         149            89         67    46    19   10   4    3       49.21           48.71               5.56
  OPHIOLITE.024         k__Bacteria (UID203)            5449        100            55         43    31    16   10   0    0       47.58           27.15               0.00
  OPHIOLITE.067     p__Actinobacteria (UID1454)         732         206           120         89    82    28   6    1    0       45.23           12.46               0.00
  OPHIOLITE.035        k__Bacteria (UID1452)            924         151           101         86    44    19   2    0    0       41.84           16.53              24.00
  OPHIOLITE.027   c__Gammaproteobacteria (UID4274)      112         581           290        340   195    41   5    0    0       40.33            6.61               5.36
  OPHIOLITE.034         k__Bacteria (UID203)            5449        104            58         55    30    15   4    0    0       36.44           12.93               0.00
  OPHIOLITE.090        k__Bacteria (UID2495)            2993        147            91        106    40    1    0    0    0       35.06            0.18               0.00
  OPHIOLITE.072        k__Bacteria (UID2569)            434         278           186        183    73    21   1    0    0       31.88            7.22               0.00
  OPHIOLITE.039         k__Bacteria (UID203)            5449        104            58         65    37    1    1    0    0       28.70            5.17              25.00
  OPHIOLITE.025   c__Gammaproteobacteria (UID4274)      112         580           289        393   133    44   9    1    0       28.59            8.45               1.30
  OPHIOLITE.060         k__Bacteria (UID203)            5449        104            58         52    34    9    7    2    0       27.48            9.09               4.76
  OPHIOLITE.019         k__Bacteria (UID203)            5449        103            57         83    20    0    0    0    0       20.33            0.00               0.00
  OPHIOLITE.002     o__Burkholderiales (UID4000)        193         427           214        348    44    35   0    0    0       17.45            7.91              62.86
  OPHIOLITE.070   c__Deltaproteobacteria (UID3217)       62         280           168        209    55    15   1    0    0       17.24            3.53               0.00
  OPHIOLITE.017         k__Bacteria (UID203)            5449        104            58         85    19    0    0    0    0       15.03            0.00               0.00
  OPHIOLITE.030         k__Bacteria (UID203)            5449        104            58         94    10    0    0    0    0       10.86            0.00               0.00
  OPHIOLITE.003             root (UID1)                 5656         56            24         56    0     0    0    0    0        0.00            0.00               0.00
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source deactivate
```
**> Time: ~25 minutes**

We filtered out other than high-quality bins. We considered high-quality bins those meeting the criteria: **Completeness ≥ 70 && Contamination ≤ 10**. 24 bins are of high quality:

```sh
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Bin Id                   Marker lineage            # genomes   # markers   # marker sets    0     1     2    3    4    5+   Completeness   Contamination   Strain heterogeneity
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  OPHIOLITE.040        k__Bacteria (UID2569)            434         278           186         2    273    3    0    0    0       99.19            1.61               0.00
  OPHIOLITE.059        k__Bacteria (UID2569)            434         278           186         2    266    10   0    0    0       98.92            4.03               0.00
  OPHIOLITE.010      p__Bacteroidetes (UID2591)         364         303           203         4    293    6    0    0    0       98.52            2.22              16.67
  OPHIOLITE.091      p__Cyanobacteria (UID2182)          89         544           424         7    535    2    0    0    0       98.35            0.08               0.00
  OPHIOLITE.029     o__Burkholderiales (UID4000)        193         427           214         20   386    20   1    0    0       97.16            7.29              78.26
  OPHIOLITE.013     o__Actinomycetales (UID1572)        580         286           171         5    273    8    0    0    0       97.08            1.87              62.50
  OPHIOLITE.086     p__Proteobacteria (UID3880)         1495        261           164         14   229    18   0    0    0       96.99            8.91              11.11
  OPHIOLITE.045      p__Bacteroidetes (UID2605)         350         316           210         7    289    20   0    0    0       96.90            6.17              45.00
  OPHIOLITE.094        k__Bacteria (UID2570)            433         273           183         14   255    4    0    0    0       94.54            1.14              25.00
  OPHIOLITE.080        k__Bacteria (UID2495)            2993        147            91         13   128    5    1    0    0       91.44            2.11              25.00
  OPHIOLITE.004     o__Burkholderiales (UID4000)        193         427           214         53   348    25   1    0    0       88.31            3.35              53.57
  OPHIOLITE.069     f__Spirochaetaceae (UID2535)         33         249           143         33   195    21   0    0    0       87.49            6.23              23.81
  OPHIOLITE.041     o__Burkholderiales (UID4000)        193         426           214         67   340    16   3    0    0       87.40            4.88              20.00
  OPHIOLITE.022     o__Actinomycetales (UID1572)        580         286           171         48   233    5    0    0    0       86.65            1.95              20.00
  OPHIOLITE.063   c__Gammaproteobacteria (UID4274)      112         580           289        106   421    49   4    0    0       86.29            9.82               8.20
  OPHIOLITE.001     o__Burkholderiales (UID4000)        193         427           214         65   355    7    0    0    0       85.60            0.92              42.86
  OPHIOLITE.054       p__Euryarchaeota (UID49)           95         228           153         26   193    9    0    0    0       84.97            4.95              22.22
  OPHIOLITE.081         k__Bacteria (UID203)            5449        104            58         46    56    2    0    0    0       84.48            3.45               0.00
  OPHIOLITE.082         k__Bacteria (UID203)            5449        104            58         13    91    0    0    0    0       82.03            0.00               0.00
  OPHIOLITE.079         k__Bacteria (UID203)            5449        104            58         51    53    0    0    0    0       77.59            0.00               0.00
  OPHIOLITE.092      p__Bacteroidetes (UID2605)         350         316           210         79   220    14   3    0    0       76.74            4.48               0.00
  OPHIOLITE.012   c__Betaproteobacteria (UID3971)       223         425           211        120   277    26   2    0    0       76.06            7.69              65.62
  OPHIOLITE.068   c__Deltaproteobacteria (UID3216)       83         247           155         48   178    20   1    0    0       75.39            8.33               4.35
  OPHIOLITE.023    f__Rhodobacteraceae (UID3340)         84         568           330        133   392    42   1    0    0       73.72            7.64              13.33
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
```


## [Stage 5] Functional annotation - `prokka`

Annotating high-quality bins:

```sh
jv~$ mkdir 06_ANNOTATION/
jv~$ prokka --outdir 06_ANNOTATION/OPHIOLITE.040 --prefix OPHIOLITE.040 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/OPHIOLITE.040.fasta
jv~$ prokka --outdir 06_ANNOTATION/OPHIOLITE.059 --prefix OPHIOLITE.059 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/OPHIOLITE.059.fasta
jv~$ prokka --outdir 06_ANNOTATION/OPHIOLITE.010 --prefix OPHIOLITE.010 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/OPHIOLITE.010.fasta
jv~$ prokka --outdir 06_ANNOTATION/OPHIOLITE.091 --prefix OPHIOLITE.091 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/OPHIOLITE.091.fasta
jv~$ prokka --outdir 06_ANNOTATION/OPHIOLITE.029 --prefix OPHIOLITE.029 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/OPHIOLITE.029.fasta
jv~$ prokka --outdir 06_ANNOTATION/OPHIOLITE.013 --prefix OPHIOLITE.013 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/OPHIOLITE.013.fasta
jv~$ prokka --outdir 06_ANNOTATION/OPHIOLITE.086 --prefix OPHIOLITE.086 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/OPHIOLITE.086.fasta
jv~$ prokka --outdir 06_ANNOTATION/OPHIOLITE.045 --prefix OPHIOLITE.045 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/OPHIOLITE.045.fasta
jv~$ prokka --outdir 06_ANNOTATION/OPHIOLITE.094 --prefix OPHIOLITE.094 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/OPHIOLITE.094.fasta
jv~$ prokka --outdir 06_ANNOTATION/OPHIOLITE.080 --prefix OPHIOLITE.080 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/OPHIOLITE.080.fasta
jv~$ prokka --outdir 06_ANNOTATION/OPHIOLITE.004 --prefix OPHIOLITE.004 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/OPHIOLITE.004.fasta
jv~$ prokka --outdir 06_ANNOTATION/OPHIOLITE.069 --prefix OPHIOLITE.069 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/OPHIOLITE.069.fasta
jv~$ prokka --outdir 06_ANNOTATION/OPHIOLITE.041 --prefix OPHIOLITE.041 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/OPHIOLITE.041.fasta
jv~$ prokka --outdir 06_ANNOTATION/OPHIOLITE.022 --prefix OPHIOLITE.022 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/OPHIOLITE.022.fasta
jv~$ prokka --outdir 06_ANNOTATION/OPHIOLITE.063 --prefix OPHIOLITE.063 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/OPHIOLITE.063.fasta
jv~$ prokka --outdir 06_ANNOTATION/OPHIOLITE.001 --prefix OPHIOLITE.001 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/OPHIOLITE.001.fasta
jv~$ prokka --outdir 06_ANNOTATION/OPHIOLITE.054 --prefix OPHIOLITE.054 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/OPHIOLITE.054.fasta
jv~$ prokka --outdir 06_ANNOTATION/OPHIOLITE.081 --prefix OPHIOLITE.081 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/OPHIOLITE.081.fasta
jv~$ prokka --outdir 06_ANNOTATION/OPHIOLITE.082 --prefix OPHIOLITE.082 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/OPHIOLITE.082.fasta
jv~$ prokka --outdir 06_ANNOTATION/OPHIOLITE.079 --prefix OPHIOLITE.079 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/OPHIOLITE.079.fasta
jv~$ prokka --outdir 06_ANNOTATION/OPHIOLITE.092 --prefix OPHIOLITE.092 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/OPHIOLITE.092.fasta
jv~$ prokka --outdir 06_ANNOTATION/OPHIOLITE.012 --prefix OPHIOLITE.012 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/OPHIOLITE.012.fasta
jv~$ prokka --outdir 06_ANNOTATION/OPHIOLITE.068 --prefix OPHIOLITE.068 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/OPHIOLITE.068.fasta
jv~$ prokka --outdir 06_ANNOTATION/OPHIOLITE.023 --prefix OPHIOLITE.023 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/OPHIOLITE.023.fasta
```

## [Stage 6] Taxonomy classification of bins

**> Input for `` should be bins in `.faa` format**

Create a directory to save the results
```sh
jv~$ mkdir 07_TAXONOMY/
```

Copy the bins in amoni acids `.faa` format from the `prokka` output:
```sh
jv~$ cp 06_ANNOTATION/OPHIOLITE.040/OPHIOLITE.040.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/OPHIOLITE.059/OPHIOLITE.059.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/OPHIOLITE.010/OPHIOLITE.010.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/OPHIOLITE.091/OPHIOLITE.091.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/OPHIOLITE.029/OPHIOLITE.029.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/OPHIOLITE.013/OPHIOLITE.013.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/OPHIOLITE.086/OPHIOLITE.086.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/OPHIOLITE.045/OPHIOLITE.045.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/OPHIOLITE.094/OPHIOLITE.094.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/OPHIOLITE.080/OPHIOLITE.080.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/OPHIOLITE.004/OPHIOLITE.004.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/OPHIOLITE.069/OPHIOLITE.069.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/OPHIOLITE.041/OPHIOLITE.041.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/OPHIOLITE.022/OPHIOLITE.022.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/OPHIOLITE.063/OPHIOLITE.063.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/OPHIOLITE.001/OPHIOLITE.001.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/OPHIOLITE.054/OPHIOLITE.054.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/OPHIOLITE.081/OPHIOLITE.081.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/OPHIOLITE.082/OPHIOLITE.082.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/OPHIOLITE.079/OPHIOLITE.079.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/OPHIOLITE.092/OPHIOLITE.092.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/OPHIOLITE.012/OPHIOLITE.012.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/OPHIOLITE.068/OPHIOLITE.068.faa 07_TAXONOMY/
jv~$ cp 06_ANNOTATION/OPHIOLITE.023/OPHIOLITE.023.faa 07_TAXONOMY/
jv~$ ls -1 07_TAXONOMY/
```

Go to the `phylophlan` directory (`$DIR_PHYLOPHLAN`), then create a new folder inside the `input` folder containing the bins in `.faa` format:

```sh
jv~$ cd ~/programs/phylophlan/
jv~$ source activate phylophlan
jv~$ mkdir input/07_TAXONOMY_OPHIOLITE/
jv~$ cp ~/data/metagenomics-ophiolite/07_TAXONOMY/* input/07_TAXONOMY_OPHIOLITE/
jv~$ ls -1 input/07_TAXONOMY_OPHIOLITE/
OPHIOLITE.001.faa
OPHIOLITE.004.faa
OPHIOLITE.010.faa
OPHIOLITE.012.faa
OPHIOLITE.013.faa
OPHIOLITE.022.faa
OPHIOLITE.023.faa
OPHIOLITE.029.faa
OPHIOLITE.040.faa
OPHIOLITE.041.faa
OPHIOLITE.045.faa
OPHIOLITE.054.faa
OPHIOLITE.059.faa
OPHIOLITE.063.faa
OPHIOLITE.068.faa
OPHIOLITE.069.faa
OPHIOLITE.079.faa
OPHIOLITE.080.faa
OPHIOLITE.081.faa
OPHIOLITE.082.faa
OPHIOLITE.086.faa
OPHIOLITE.091.faa
OPHIOLITE.092.faa
OPHIOLITE.094.faa
```

```sh
jv~$ ./phylophlan.py -i -t 07_TAXONOMY_OPHIOLITE --nproc 32
Total time: 2317.37 seconds Unique: 3163/3195 Bad splits: 6/3160 Worst delta-LogLk 0.379
Tree built! The output newick file is in output/07_TAXONOMY_OPHIOLITE/07_TAXONOMY_OPHIOLITE.tree.int.nwk
Tree building finished in 2487 secs (3564 total time).
Trying to impute taxonomic labels for taxa newly integrated into the tree... Done!
Writing taxonomic imputation outputs ... Done!
jv~$ source deactivate
```

Results with high confidence are in `imputed_conf_high-conf.txt`
```sh
jv~$ cat output/07_TAXONOMY_OPHIOLITE/imputed_conf_high-conf.txt
OPHIOLITE.022	d__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Actinomycetales.f__Microbacteriaceae.g__?.s__?.t__?
OPHIOLITE.023	d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhodobacterales.f__Rhodobacteraceae.g__Rhodobacter.s__sphaeroides.t__?
OPHIOLITE.013	d__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Actinomycetales.f__Microbacteriaceae.g__?.s__?.t__?
OPHIOLITE.012	d__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Rhodocyclales.f__Rhodocyclaceae.g__?.s__?.t__?
OPHIOLITE.063	d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Methylococcales.f__Methylococcaceae.g__?.s__?.t__?
OPHIOLITE.068	d__Bacteria.p__Proteobacteria.c__Deltaproteobacteria.o__?.f__?.g__?.s__?.t__?
OPHIOLITE.040	d__Bacteria.p__Bacteroidetes.c__?.o__?.f__?.g__?.s__?.t__?
OPHIOLITE.086	d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Thiotrichales.f__Thiotrichaceae.g__Beggiatoa.s__?.t__?
```
Results with medium confidence are in `imputed_conf_medium-conf.txt`

```sh
jv~$ cat output/07_TAXONOMY_OPHIOLITE/imputed_conf_medium-conf.txt
OPHIOLITE.004	d__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__?.f__?.g__?.s__?.t__?
OPHIOLITE.001	d__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__?.f__?.g__?.s__?.t__?
OPHIOLITE.079	d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhodospirillales.f__?.g__?.s__?.t__?
OPHIOLITE.054	d__Archaea.p__Euryarchaeota.c__Methanomicrobia.o__?.f__?.g__?.s__?.t__?
OPHIOLITE.041	d__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__?.f__?.g__?.s__?.t__?
OPHIOLITE.010	d__Bacteria.p__Bacteroidetes.c__?.o__?.f__?.g__?.s__?.t__?
OPHIOLITE.094	d__Bacteria.p__Chlorobi.c__Chlorobia.o__Chlorobiales.f__Chlorobiaceae.g__?.s__?.t__?
OPHIOLITE.081	d__Bacteria.p__Chlorobi.c__Chlorobia.o__Chlorobiales.f__Chlorobiaceae.g__?.s__?.t__?
OPHIOLITE.069	d__Bacteria.p__Spirochaetes.c__Spirochaetes.o__Spirochaetales.f__Treponemaceae.g__Treponema.s__?.t__?
```
**> Total time: ~ 1 h**

Copy the `phylophlan` results from `$DIR_PHYLOPHLAN` to the `$DIR_PROJECT`

```sh 
jv~$ cp output/07_TAXONOMY_OPHIOLITE/* ~/data/metagenomics-ophiolite/07_TAXONOMY/
```

## [Stage 7]  Refine functional annotation of methanotrophic bacteria using `EggNOG-mapper`
`EggNOG-mapper v4.5.1` works on `python2.7`

`$DIR_EGGNOG` is the directory to `eggnog-mapper/emapper.py` script.

```sh
jv~$ cd ~/data/metagenomics-ophiolite/
jv~$ mkdir 08_ANNOTATION_REFINE/
jv~$ cd 08_ANNOTATION_REFINE/
jv~$ python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/OPHIOLITE.063.faa --output OPHIOLITE.063 -m diamond --cpu 32
jv~$ python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/OPHIOLITE.022.faa --output OPHIOLITE.022 -m diamond --cpu 32
jv~$ python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/OPHIOLITE.023.faa --output OPHIOLITE.023 -m diamond --cpu 32
jv~$ python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/OPHIOLITE.013.faa --output OPHIOLITE.013 -m diamond --cpu 32
jv~$ python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/OPHIOLITE.012.faa --output OPHIOLITE.012 -m diamond --cpu 32
jv~$ python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/OPHIOLITE.068.faa --output OPHIOLITE.068 -m diamond --cpu 32
jv~$ python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/OPHIOLITE.040.faa --output OPHIOLITE.040 -m diamond --cpu 32
jv~$ python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/OPHIOLITE.086.faa --output OPHIOLITE.086 -m diamond --cpu 32
jv~$ python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/OPHIOLITE.004.faa --output OPHIOLITE.004 -m diamond --cpu 32
jv~$ python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/OPHIOLITE.001.faa --output OPHIOLITE.001 -m diamond --cpu 32
jv~$ python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/OPHIOLITE.079.faa --output OPHIOLITE.079 -m diamond --cpu 32
jv~$ python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/OPHIOLITE.054.faa --output OPHIOLITE.054 -m diamond --cpu 32
jv~$ python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/OPHIOLITE.041.faa --output OPHIOLITE.041 -m diamond --cpu 32
jv~$ python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/OPHIOLITE.010.faa --output OPHIOLITE.010 -m diamond --cpu 32
jv~$ python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/OPHIOLITE.094.faa --output OPHIOLITE.094 -m diamond --cpu 32
jv~$ python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/OPHIOLITE.081.faa --output OPHIOLITE.081 -m diamond --cpu 32
jv~$ python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/OPHIOLITE.069.faa --output OPHIOLITE.069 -m diamond --cpu 32
jv~$ python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/OPHIOLITE.029.faa --output OPHIOLITE.029 -m diamond --cpu 32
jv~$ python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/OPHIOLITE.045.faa --output OPHIOLITE.045 -m diamond --cpu 32
jv~$ python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/OPHIOLITE.059.faa --output OPHIOLITE.059 -m diamond --cpu 32
jv~$ python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/OPHIOLITE.080.faa --output OPHIOLITE.080 -m diamond --cpu 32
jv~$ python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/OPHIOLITE.082.faa --output OPHIOLITE.082 -m diamond --cpu 32
jv~$ python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/OPHIOLITE.091.faa --output OPHIOLITE.091 -m diamond --cpu 32
jv~$ python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/OPHIOLITE.092.faa --output OPHIOLITE.092 -m diamond --cpu 32
```
