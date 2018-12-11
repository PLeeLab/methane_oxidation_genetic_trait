# Metagenomics from Serpentinite Springs of the Voltri Massif, Italy

<!-- TOC -->

- [Metagenomics from Serpentinite Springs of the Voltri Massif, Italy](#metagenomics-from-serpentinite-springs-of-the-voltri-massif-italy)
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
jv~$ mkdir metagenomics-serpentinite-springs
jv~$ cd metagenomics-serpentinite-springs/
```

### Retrieve sequences from SRA NCBI system

> This dataset was published in:
> Brazelton, W. J., Thornton, C. N., Hyer, A., Twing, K. I., Longino, A. A., Lang, S. Q., ... & Schrenk, M. O. (2017). Metagenomic identification of active methanogens and methanotrophs in serpentinite springs of the Voltri Massif, Italy. _PeerJ_, 5, e2945.
> https://peerj.com/articles/2945/

### Download metagenomics data

>Open the `SraAccList.txt` file. Replace `\n` by ` ` . Then paste the resulting line in the field below.

This code will retrieve the reads in separate ``FASTQ`` files:

```sh
jv~$ mkdir 00_RAW
jv~$ cd 00_RAW/
jv~$ fastq-dump -I --split-files SRR1636510 SRR1636512 SRR1636513 SRR1636514 SRR2058405 SRR2058406 SRR2058407 SRR2058408
jv~$ ls -1
```

### Check number of reads

#### Sample 1
```sh
jv~$ head -1 SRR1636510_1.fastq
# check the header "@SRR1636510.1.1"

jv~$ grep -c "@SRR1636510" SRR1636510_1.fastq
10501350

jv~$ grep -c "@SRR1636510" SRR1636510_2.fastq
10501350
```

#### Sample 2
```sh
jv~$ head -1 SRR1636512_1.fastq
# @SRR1636512.1.1

jv~$ grep -c "@SRR1636512" SRR1636512_1.fastq
22635241
jv~$ grep -c "@SRR1636512" SRR1636512_2.fastq
22635241
```

#### Sample 3
```sh
jv~$ grep -c "@SRR1636513" SRR1636513_1.fastq
30005532
jv~$ grep -c "@SRR1636513" SRR1636513_2.fastq
30005532
```

#### Sample 4
```sh
jv~$ grep -c "@SRR1636514" SRR1636514_1.fastq
38071931
jv~$ grep -c "@SRR1636514" SRR1636514_2.fastq
38071931
```

#### Sample 5
```sh
jv~$ grep -c "@SRR2058405" SRR2058405_1.fastq
15452030
jv~$ grep -c "@SRR2058405" SRR2058405_2.fastq
15452030
```

#### Sample 6
```sh
jv~$ grep -c "@SRR2058406" SRR2058406_1.fastq
19439690
jv~$ grep -c "@SRR2058406" SRR2058406_2.fastq
19439690
```

#### Sample 7
```sh
jv~$ grep -c "@SRR2058407" SRR2058407_1.fastq
20259258
jv~$ grep -c "@SRR2058407" SRR2058407_2.fastq
20259258
```

#### Sample 8
```sh
jv~$ grep -c "@SRR2058408" SRR2058408_1.fastq
23890776
jv~$ grep -c "@SRR2058408" SRR2058408_2.fastq
23890776
```

#### Summary of reads
Generate some stats of the samples in `R`

```sh
jv~$ R
```

```R
x <- c(10501350, 22635241, 30005532, 38071931, 15452030, 19439690, 20259258, 23890776)
mean(x)
[1] 22531976

median(x)
21447250

sd(x)
[1] 8525511

quit()
```

|#|Sample    |reads  |
|-|----------|-------|
|1|SRR1636510|10501350|
|2|SRR1636512|22635241|
|3|SRR1636513|30005532|
|4|SRR1636514|38071931|
|5|SRR2058405|15452030|
|6|SRR2058406|19439690|
|7|SRR2058407|20259258|
|8|SRR2058408|23890776|
|**Average**|**22531976**|**~22 M reads**|
|Std.Dev|8525511|~8 M reads|

> Sufficient number of reads (cat feces 20-90 M reads)

## [Stage 1] Quality Control (QC)

### Quality Filtering - overlapped reads

The authors of the paper say:
>"Paired-end sequencing with a 100 cycle Illumina HiSeq run generated partial ∼30 bp overlaps, and six libraries were multiplexed per lane."

> This protocol section is partially based on: http://merenlab.org/tutorials/assembly-based-metagenomics/

Generate a TAB-delimited `samples.txt` file to point out where are your raw `R1` and `R2` files for each sample.

|sample    |r1                         |r2                          |
|----------|---------------------------|----------------------------|
|Sample_01 | 00_RAW/SRR1636510_1.fastq | 00_RAW/SRR1636510_2.fastq  |
|Sample_02 | 00_RAW/SRR1636512_1.fastq | 00_RAW/SRR1636512_2.fastq  |
|Sample_03 | 00_RAW/SRR1636513_1.fastq | 00_RAW/SRR1636513_2.fastq  |
|Sample_04 | 00_RAW/SRR1636514_1.fastq | 00_RAW/SRR1636514_2.fastq  |
|Sample_05 | 00_RAW/SRR2058405_1.fastq | 00_RAW/SRR2058405_2.fastq  |
|Sample_06 | 00_RAW/SRR2058406_1.fastq | 00_RAW/SRR2058406_2.fastq  |
|Sample_07 | 00_RAW/SRR2058407_1.fastq | 00_RAW/SRR2058407_2.fastq  |
|Sample_08 | 00_RAW/SRR2058408_1.fastq | 00_RAW/SRR2058408_2.fastq  |

Create a directory for quality-filtered `R1` and `R2`
```sh
jv~$ mkdir 01_QC_overlap/

jv~$ source activate illumina-utils

jv~$ iu-gen-configs samples.txt -o 01_QC_overlap

Report .......................................: Read for 8 samples is read
Output directory set in configs ..............: /metagenomics-serpentinite-springs/01_QC_overlap
Prefix for R1 ................................: None
Prefix for R2 ................................: None

jv~$ ls 01_QC_overlap/
Sample_01.ini  Sample_03.ini  Sample_05.ini  Sample_07.ini
Sample_02.ini  Sample_04.ini  Sample_06.ini  Sample_08.ini
```

Run quality filtering for all the samples at once:

The authors of the paper say:
>"Paired-end sequencing with a 100 cycle Illumina HiSeq run generated partial ∼30 bp overlaps, and six libraries were multiplexed per lane."

Therefore, we will use `iu-merge-pairs` for partially overlapping reads (instead of `iu-filter-quality-minoche`) in constrast to the method that has been used in  other protocols within this repository.

Note also that we include the `parameter` `--ignore-deflines`

```sh
jv~$ for ini in 01_QC_overlap/*.ini; do iu-merge-pairs --ignore-deflines $ini; done
jv~$ source deactivate
```

The contents of the `01_QC_overlap/` directory should look like this:
```sh
jv~$ ls 01_QC_overlap/
Sample_01_FAILED                Sample_03_MISMATCHES_BREAKDOWN  Sample_06.ini
Sample_01_FAILED_WITH_Ns        Sample_03_STATS                 Sample_06_MERGED
Sample_01.ini                   Sample_04_FAILED                Sample_06_MISMATCHES_BREAKDOWN
Sample_01_MERGED                Sample_04_FAILED_WITH_Ns        Sample_06_STATS
Sample_01_MISMATCHES_BREAKDOWN  Sample_04.ini                   Sample_07_FAILED
Sample_01_STATS                 Sample_04_MERGED                Sample_07_FAILED_WITH_Ns
Sample_02_FAILED                Sample_04_MISMATCHES_BREAKDOWN  Sample_07.ini
Sample_02_FAILED_WITH_Ns        Sample_04_STATS                 Sample_07_MERGED
Sample_02.ini                   Sample_05_FAILED                Sample_07_MISMATCHES_BREAKDOWN
Sample_02_MERGED                Sample_05_FAILED_WITH_Ns        Sample_07_STATS
Sample_02_MISMATCHES_BREAKDOWN  Sample_05.ini                   Sample_08_FAILED
Sample_02_STATS                 Sample_05_MERGED                Sample_08_FAILED_WITH_Ns
Sample_03_FAILED                Sample_05_MISMATCHES_BREAKDOWN  Sample_08.ini
Sample_03_FAILED_WITH_Ns        Sample_05_STATS                 Sample_08_MERGED
Sample_03.ini                   Sample_06_FAILED                Sample_08_MISMATCHES_BREAKDOWN
Sample_03_MERGED                Sample_06_FAILED_WITH_Ns        Sample_08_STATS
```

#### Checking QC results stats in `*_STATS.txt` files:

##### Sample_01-STATS
```sh
jv~$ cat 01_QC_overlap/Sample_01_STATS
Number of pairs analyzed ...............................	10501350
Prefix failed in read 1 ................................	0
Prefix failed in read 2 ................................	0
Prefix failed in both ..................................	0
Passed prefix total ....................................	10501350
Failed prefix total ....................................	0
Merged total ...........................................	7780452
```

##### Sample_02-STATS
```sh
jv~$ cat 01_QC_overlap/Sample_02_STATS
Number of pairs analyzed ...............................	22635241
Prefix failed in read 1 ................................	0
Prefix failed in read 2 ................................	0
Prefix failed in both ..................................	0
Passed prefix total ....................................	22635241
Failed prefix total ....................................	0
Merged total ...........................................	19387360
```

##### Sample_03-STATS
```sh
jv~$ head -7 01_QC_overlap/Sample_03_STATS
Number of pairs analyzed ...............................	30005532
Prefix failed in read 1 ................................	0
Prefix failed in read 2 ................................	0
Prefix failed in both ..................................	0
Passed prefix total ....................................	30005532
Failed prefix total ....................................	0
Merged total ...........................................	25469249
```

##### Sample_04-STATS
```sh
jv~$ head -7 01_QC_overlap/Sample_04_STATS
Number of pairs analyzed ...............................	38071931
Prefix failed in read 1 ................................	0
Prefix failed in read 2 ................................	0
Prefix failed in both ..................................	0
Passed prefix total ....................................	38071931
Failed prefix total ....................................	0
Merged total ...........................................	32960531
```

##### Sample_05-STATS
```sh
jv~$ head -7 01_QC_overlap/Sample_05_STATS
Number of pairs analyzed ...............................	15452030
Prefix failed in read 1 ................................	0
Prefix failed in read 2 ................................	0
Prefix failed in both ..................................	0
Passed prefix total ....................................	15452030
Failed prefix total ....................................	0
Merged total ...........................................	15062989
```

##### Sample_06-STATS
```sh
jv~$ head -7 01_QC_overlap/Sample_06_STATS
Number of pairs analyzed ...............................	19439690
Prefix failed in read 1 ................................	0
Prefix failed in read 2 ................................	0
Prefix failed in both ..................................	0
Passed prefix total ....................................	19439690
Failed prefix total ....................................	0
Merged total ...........................................	18685905
```

##### Sample_07-STATS
```sh
jv~$ head -7 01_QC_overlap/Sample_07_STATS
Number of pairs analyzed ...............................	20259258
Prefix failed in read 1 ................................	0
Prefix failed in read 2 ................................	0
Prefix failed in both ..................................	0
Passed prefix total ....................................	20259258
Failed prefix total ....................................	0
Merged total ...........................................	19286425
```

##### Sample_08-STATS
```sh
jv~$ head -7 01_QC_overlap/Sample_08_STATS
Number of pairs analyzed ...............................	23890776
Prefix failed in read 1 ................................	0
Prefix failed in read 2 ................................	0
Prefix failed in both ..................................	0
Passed prefix total ....................................	23890776
Failed prefix total ....................................	0
Merged total ...........................................	23125398
```

##### Summary Stage 1 - QC Filtering
```sh
jv~$ grep 'total pairs passed' 01_QC_overlap/*STATS.txt
```

|#|Sample    |number of pairs analyzed   | Merged total  |
|-|----------|---------------------------|---------------|
|1|SRR1636510|10501350                   | 7780452|
|2|SRR1636512|22635241                   |19387360|
|3|SRR1636513|30005532                   |25469249|
|4|SRR1636514|38071931                   |32960531|
|5|SRR2058405|15452030                   |15062989|
|6|SRR2058406|19439690                   |18685905|
|7|SRR2058407|2652260                    |19286425|
|8|SRR2058408|23890776                   |23125398|


## [Stage 2] Co-assembly (pooled-assembly)

Generating `contigs.fasta`

### Run MEGAHIT

Create an environment variable:
```sh
jv~$ R12s=`ls 01_QC_overlap/*MERGED* | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`
```

```sh
jv~$ echo $R12s
01_QC_overlap/Sample_01_MERGED,01_QC_overlap/Sample_02_MERGED,01_QC_overlap/Sample_03_MERGED,01_QC_overlap/Sample_04_MERGED,01_QC_overlap/Sample_05_MERGED,01_QC_overlap/Sample_06_MERGED,01_QC_overlap/Sample_07_MERGED,01_QC_overlap/Sample_08_MERGED
```

Note that as the pairs were merged, we will use the parameter `-r in this assembly setup.

```sh
jv~$ megahit -r $R12s --min-contig-len 1000 -m 0.99 -o 02_ASSEMBLY_overlap/ -t 32
251.0Gb memory in total.
Using: 249.122Gb.
MEGAHIT v1.1.2
.
.
.
--- [STAT] 111139 contigs, total 263683255 bp, min 1000 bp, max 169126 bp, avg 2373 bp, N50 2653 bp
--- [Wed Feb 14 00:46:27 2018] ALL DONE. Time elapsed: 26600.867049 seconds ---
```

_**Time elapsed: 26600.867049 seconds = ~8 hours**_

### Refining contigs with `anvio`

```sh
jv~$ mkdir 03_CONTIGS_overlap/
jv~$ source activate anvio3
jv~$ anvi-script-reformat-fasta 02_ASSEMBLY_overlap/final.contigs.fa -o 03_CONTIGS_overlap/contigs.fa --min-len 1000 --simplify-names --report name_conversions_1000.txt
Input ........................................: 02_ASSEMBLY_overlap/final.contigs.fa
Output .......................................: 03_CONTIGS_overlap/contigs.fa
Minimum length ...............................: 1,000
Total num contigs ............................: 111,139
Total num nucleotides ........................: 263,683,255
Contigs removed ..............................: 0 (0.00% of all)
Nucleotides removed ..........................: 0 (0.00% of all)
Deflines simplified ..........................: True
```

Testing for `--min-len 1500`

```sh
jv~$ anvi-script-reformat-fasta 02_ASSEMBLY_overlap/final.contigs.fa -o 03_CONTIGS_overlap/contigs_1500.fa --min-len 1500 --simplify-names --report name_conversions_1500.txt
Input ........................................: 02_ASSEMBLY_overlap/final.contigs.fa
Output .......................................: 03_CONTIGS_overlap/contigs_1500.fa
Minimum length ...............................: 1,500
Total num contigs ............................: 111,139
Total num nucleotides ........................: 263,683,255
Contigs removed ..............................: 54465 (49.01% of all)
Nucleotides removed ..........................: 65534078 (24.85% of all)
Deflines simplified ..........................: True
```

Testing for `--min-len 2500`
```sh
jv~$ anvi-script-reformat-fasta 02_ASSEMBLY_overlap/final.contigs.fa -o 03_CONTIGS_overlap/contigs_2500.fa --min-len 2500 --simplify-names --report name_conversions_2500.txt
Input ........................................: 02_ASSEMBLY_overlap/final.contigs.fa
Output .......................................: 03_CONTIGS_overlap/contigs_2500.fa
Minimum length ...............................: 2,500
Total num contigs ............................: 111,139
Total num nucleotides ........................: 263,683,255
Contigs removed ..............................: 86451 (77.79% of all)
Nucleotides removed ..........................: 125926650 (47.76% of all)
Deflines simplified ..........................: True
```

We will proceed with `03_CONTIGS_overlap/contigs.fa` which has **`--min-len` = 1000**

Deactivate the `anvio` conda environment

```sh
jv~$ source deactivate
```

## [Stage 3] Binning

### Binning with `MaxBin v2.2.4`

```sh
jv~$ mkdir 04_BINNING_overlap/

```

Test `Maxbin2` installation

```sh
jv~$ run_MaxBin.pl -v
MaxBin 2.2.4
```

Create `reads_list.txt` to be passed to `MaxBin2`:

```sh
jv~$ ls 01_QC_overlap/*MERGED > 04_BINNING_overlap/reads_list.txt
```

Corroborate the `reads_list.txt` content:

```sh
jv~$ cat 04_BINNING_overlap/reads_list.txt
01_QC_overlap/Sample_01_MERGED
01_QC_overlap/Sample_02_MERGED
01_QC_overlap/Sample_03_MERGED
01_QC_overlap/Sample_04_MERGED
01_QC_overlap/Sample_05_MERGED
01_QC_overlap/Sample_06_MERGED
01_QC_overlap/Sample_07_MERGED
01_QC_overlap/Sample_08_MERGED
```

Run MaxBin binning on the assembled contigs:

> Note that here we are indicating `-thread 32`

```sh
jv~$ run_MaxBin.pl -contig 03_CONTIGS_overlap/contigs.fa -reads_list 04_BINNING_overlap/reads_list.txt -out 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap -thread 32
========== Job finished ==========
Yielded 107 bins for contig (scaffold) file 03_CONTIGS_overlap/contigs.fa
Here are the output files for this run.
Please refer to the README file for further details.
Summary file: 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.summary
Genome abundance info file: 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.abundance
Marker counts: 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.marker
Marker genes for each bin: 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.marker_of_each_gene.tar.gz
Bin files: 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.001.fasta - 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.107.fasta
Unbinned sequences: 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.noclass
.
.
.
========== Elapsed Time ==========
2 hours 43 minutes and 59 seconds.
```
**Elapsed time: ~3 hours**

We got 107 bins at this stage. Using the `minoche` method in the other protocol we retrieved 106 bins.

## [Stage 4] Refine bins
CheckM

```sh
jv~$ source activate CheckM
jv~$ checkm lineage_wf -t 32 -x fasta 04_BINNING_overlap/ 05_CHECKM_overlap/
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Bin Id                                      Marker lineage            # genomes   # markers   # marker sets    0     1     2     3    4    5+   Completeness   Contamination   Strain heterogeneity
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  SERPENTINITE_SPRINGS_overlap.044         k__Bacteria (UID203)            5449        104            58         0     60    39    5    0    0       100.00          40.53               1.85
  SERPENTINITE_SPRINGS_overlap.027   c__Gammaproteobacteria (UID4274)      112         580           289         26   539    15    0    0    0       98.10            2.60              60.00
  SERPENTINITE_SPRINGS_overlap.052      p__Bacteroidetes (UID2605)         350         316           210         4    277    28    7    0    0       98.10           16.67               4.08
  SERPENTINITE_SPRINGS_overlap.004   c__Deltaproteobacteria (UID3217)       62         280           168         6    272    2     0    0    0       97.00            0.89               0.00
  SERPENTINITE_SPRINGS_overlap.091    o__Sphingomonadales (UID3310)         26         569           293         18   529    22    0    0    0       96.75            2.80              31.82
  SERPENTINITE_SPRINGS_overlap.065   c__Gammaproteobacteria (UID4274)      112         580           289         53   517    10    0    0    0       96.50            1.78              60.00
  SERPENTINITE_SPRINGS_overlap.006        k__Bacteria (UID3187)            2258        181           110         5    167    8     1    0    0       96.31            3.43               9.09
  SERPENTINITE_SPRINGS_overlap.032         k__Bacteria (UID203)            5449        103            58         5     67    31    0    0    0       96.21           16.99               0.00
  SERPENTINITE_SPRINGS_overlap.103      p__Bacteroidetes (UID2591)         364         303           203         16   276    11    0    0    0       95.21            2.19              72.73
  SERPENTINITE_SPRINGS_overlap.005         k__Bacteria (UID203)            5449        104            58         3     53    46    2    0    0       94.83           44.00              94.23
  SERPENTINITE_SPRINGS_overlap.098    f__Flavobacteriaceae (UID2817)        81         511           283         25   427    56    3    0    0       94.59           12.06              67.69
  SERPENTINITE_SPRINGS_overlap.096     o__Burkholderiales (UID4001)        108         570           250         33   518    18    1    0    0       93.42            3.57              42.86
  SERPENTINITE_SPRINGS_overlap.056      p__Bacteroidetes (UID2605)         350         316           210         18   261    35    2    0    0       93.00            9.42              75.61
  SERPENTINITE_SPRINGS_overlap.070     o__Burkholderiales (UID4000)        193         427           214         28   208   140    48   3    0       91.95           60.58               7.95
  SERPENTINITE_SPRINGS_overlap.077    o__Sphingomonadales (UID3310)         26         569           293         58   195   120   105   64   27      91.67           119.92             23.36
  SERPENTINITE_SPRINGS_overlap.018     p__Actinobacteria (UID1454)         732         206           120         13   185    8     0    0    0       91.67            4.06              62.50
  SERPENTINITE_SPRINGS_overlap.058        k__Bacteria (UID2329)            174         149            89         9    127    8     2    3    0       91.57            8.46              12.50
  SERPENTINITE_SPRINGS_overlap.072        k__Bacteria (UID2570)            433         274           183         24   240    10    0    0    0       91.38            2.47              50.00
  SERPENTINITE_SPRINGS_overlap.068      o__Clostridiales (UID1120)         304         247           141         16   201    27    3    0    0       91.13           13.62               5.56
  SERPENTINITE_SPRINGS_overlap.094      p__Bacteroidetes (UID2591)         364         303           203         62   227    14    0    0    0       91.11            3.97              71.43
  SERPENTINITE_SPRINGS_overlap.078     o__Burkholderiales (UID4000)        193         427           214         71   303    46    7    0    0       90.79           16.05              35.82
  SERPENTINITE_SPRINGS_overlap.063         k__Bacteria (UID203)            5449        104            58         9     34    23    28   6    4       90.05           93.45               2.13
  SERPENTINITE_SPRINGS_overlap.013         k__Bacteria (UID209)            5443        105            59         6     88    11    0    0    0       89.83            9.94               0.00
  SERPENTINITE_SPRINGS_overlap.050    f__Rhodobacteraceae (UID3340)         84         568           330         74   189   253    42   8    2       89.17           63.67               2.91
  SERPENTINITE_SPRINGS_overlap.041    f__Xanthomonadaceae (UID4214)         55         659           290         68   545    45    1    0    0       88.41            7.44              58.33
  SERPENTINITE_SPRINGS_overlap.024         k__Bacteria (UID203)            5449        104            58         8     46    44    6    0    0       88.22           39.89               6.45
  SERPENTINITE_SPRINGS_overlap.022         k__Bacteria (UID203)            5449        104            58         20    45    38    1    0    0       86.65           37.52               2.44
  SERPENTINITE_SPRINGS_overlap.102     o__Actinomycetales (UID1663)        488         310           185         55    96    83    42   25   9       86.40           103.22             18.81
  SERPENTINITE_SPRINGS_overlap.073        k__Bacteria (UID2372)            131         177           106         46   122    9     0    0    0       86.32            4.75              11.11
  SERPENTINITE_SPRINGS_overlap.047         k__Bacteria (UID203)            5449        104            58         10    44    30    18   2    0       85.63           67.65              26.04
  SERPENTINITE_SPRINGS_overlap.017          k__Archaea (UID2)              207         149           107         25    94    28    2    0    0       85.62           21.02              55.88
  SERPENTINITE_SPRINGS_overlap.088        k__Bacteria (UID2569)            434         275           184         31   231    11    2    0    0       85.56            4.37              17.65
  SERPENTINITE_SPRINGS_overlap.092     o__Actinomycetales (UID1572)        580         286           171         57   117   102    9    1    0       84.98           47.56              48.15
  SERPENTINITE_SPRINGS_overlap.064        p__Firmicutes (UID241)           930         213           118         43   154    16    0    0    0       83.69            8.86              50.00
  SERPENTINITE_SPRINGS_overlap.076     o__Burkholderiales (UID4000)        193         427           214         70   295    55    7    0    0       83.43           13.11              46.05
  SERPENTINITE_SPRINGS_overlap.046        p__Firmicutes (UID241)           930         213           118         25   174    13    1    0    0       83.05            6.84              25.00
  SERPENTINITE_SPRINGS_overlap.045         k__Bacteria (UID203)            5449        103            58         12    67    16    7    1    0       82.59           16.70              11.63
  SERPENTINITE_SPRINGS_overlap.082     o__Actinomycetales (UID1572)        580         286           171         37   205    43    1    0    0       82.50           14.05              52.17
  SERPENTINITE_SPRINGS_overlap.087     o__Burkholderiales (UID4001)        108         549           242        124   359    60    4    0    2       82.47           13.68              18.62
  SERPENTINITE_SPRINGS_overlap.010        k__Bacteria (UID1452)            924         151           101         20   127    4     0    0    0       82.34            2.07              50.00
  SERPENTINITE_SPRINGS_overlap.055         k__Bacteria (UID203)            5449        104            58         13    53    32    5    1    0       81.87           46.37               1.89
  SERPENTINITE_SPRINGS_overlap.033         k__Bacteria (UID203)            5449        100            55         13    48    31    4    4    0       81.79           43.22              19.40
  SERPENTINITE_SPRINGS_overlap.086    f__Flavobacteriaceae (UID2817)        81         511           283        118   326    51    16   0    0       80.85           15.79               4.04
  SERPENTINITE_SPRINGS_overlap.095     o__Burkholderiales (UID4000)        193         427           214        101   241    57    22   4    2       80.79           23.10              27.54
  SERPENTINITE_SPRINGS_overlap.100     o__Burkholderiales (UID4001)        108         570           250        114   341    94    18   2    1       80.45           28.04              17.65
  SERPENTINITE_SPRINGS_overlap.057         k__Bacteria (UID203)            5449        103            58         20    55    20    8    0    0       79.08           23.16               0.00
  SERPENTINITE_SPRINGS_overlap.104     o__Actinomycetales (UID1696)        455         315           190         84    84    99    42   5    1       78.82           71.25              36.98
  SERPENTINITE_SPRINGS_overlap.054         k__Bacteria (UID203)            5449        104            58         14    80    10    0    0    0       78.74           11.38              30.00
  SERPENTINITE_SPRINGS_overlap.090    f__Flavobacteriaceae (UID2817)        81         511           283        153   321    32    5    0    0       78.25            8.55              36.17
  SERPENTINITE_SPRINGS_overlap.036        k__Bacteria (UID2565)            2921        143            88         32    56    46    7    2    0       77.24           52.86               3.80
  SERPENTINITE_SPRINGS_overlap.085      p__Bacteroidetes (UID2591)         364         303           203         89   131    63    19   1    0       75.91           39.08               2.38
  SERPENTINITE_SPRINGS_overlap.020         k__Bacteria (UID203)            5449        104            58         16    63    21    4    0    0       75.71           19.05              12.12
  SERPENTINITE_SPRINGS_overlap.049         k__Bacteria (UID203)            5449        103            58         26    45    29    3    0    0       75.70           29.55               5.26
  SERPENTINITE_SPRINGS_overlap.059        p__Firmicutes (UID241)           930         213           118         55   142    15    1    0    0       75.19            6.64              38.89
  SERPENTINITE_SPRINGS_overlap.066       c__Clostridia (UID1118)           387         223           124         83   125    15    0    0    0       75.17            6.72               0.00
  SERPENTINITE_SPRINGS_overlap.071         k__Bacteria (UID203)            5449        104            58         46    30    24    4    0    0       73.84           50.34              25.00
  SERPENTINITE_SPRINGS_overlap.093     o__Burkholderiales (UID4001)        108         570           250        139   385    46    0    0    0       73.12            7.93              78.26
  SERPENTINITE_SPRINGS_overlap.081         k__Bacteria (UID203)            5449        104            58         19    45    29    6    5    0       72.26           23.00               1.30
  SERPENTINITE_SPRINGS_overlap.043       c__Clostridia (UID1118)           387         223           124         63   139    19    2    0    0       72.07            7.65              60.00
  SERPENTINITE_SPRINGS_overlap.025        p__Firmicutes (UID239)           1324        175           101         43   115    17    0    0    0       71.72           12.23               5.88
  SERPENTINITE_SPRINGS_overlap.030         k__Bacteria (UID203)            5449        104            58         21    41    35    6    1    0       70.38           42.79               0.00
  SERPENTINITE_SPRINGS_overlap.028        k__Bacteria (UID2570)            433         273           183         72   137    56    7    1    0       69.97           25.21               2.41
  SERPENTINITE_SPRINGS_overlap.001         k__Bacteria (UID203)            5449        104            58         21    73    9     1    0    0       68.94            8.54              16.67
  SERPENTINITE_SPRINGS_overlap.060         k__Bacteria (UID209)            5443        105            59         23    78    4     0    0    0       66.92            2.85              25.00
  SERPENTINITE_SPRINGS_overlap.069         k__Bacteria (UID203)            5449        103            57         25    40    23    7    5    3       66.91           35.16               0.00
  SERPENTINITE_SPRINGS_overlap.038         k__Bacteria (UID203)            5449        104            58         26    44    27    4    3    0       66.77           23.07               7.02
  SERPENTINITE_SPRINGS_overlap.014         k__Bacteria (UID203)            5449        104            58         35    45    22    2    0    0       66.19           17.60              21.43
  SERPENTINITE_SPRINGS_overlap.067     o__Burkholderiales (UID4000)        193         427           214        172   228    27    0    0    0       65.27            7.26              48.15
  SERPENTINITE_SPRINGS_overlap.002         k__Bacteria (UID203)            5449        104            58         25    76    3     0    0    0       64.63            3.61               0.00
  SERPENTINITE_SPRINGS_overlap.009         k__Bacteria (UID203)            5449        104            58         60    42    2     0    0    0       63.66            3.45              50.00
  SERPENTINITE_SPRINGS_overlap.089        k__Bacteria (UID2982)             88         229           147         83   112    31    3    0    0       63.45           12.50               0.00
  SERPENTINITE_SPRINGS_overlap.074          k__Archaea (UID2)              207         145           103         51    68    21    3    1    1       63.44           17.65              10.87
  SERPENTINITE_SPRINGS_overlap.053     o__Burkholderiales (UID4000)        193         426           214        157   255    14    0    0    0       62.75            3.43              28.57
  SERPENTINITE_SPRINGS_overlap.037        p__Firmicutes (UID241)           930         213           118         91    84    30    4    4    0       61.80           26.38              21.21
  SERPENTINITE_SPRINGS_overlap.016     o__Burkholderiales (UID4000)        193         427           214        150   259    17    1    0    0       60.85            4.40              20.00
  SERPENTINITE_SPRINGS_overlap.079         k__Bacteria (UID203)            5449        104            58         33    36    22    10   2    1       60.76           43.77               0.00
  SERPENTINITE_SPRINGS_overlap.051        k__Bacteria (UID2495)            2993        143            89         54    54    30    3    0    2       60.41           32.80               1.69
  SERPENTINITE_SPRINGS_overlap.023         k__Bacteria (UID203)            5449        103            57         32    45    19    5    2    0       59.14           20.88               4.35
  SERPENTINITE_SPRINGS_overlap.021         k__Bacteria (UID203)            5449        104            58         31    48    25    0    0    0       58.54           21.64              20.00
  SERPENTINITE_SPRINGS_overlap.099      p__Bacteroidetes (UID2605)         350         316           210        150   150    15    1    0    0       55.88            6.83               0.00
  SERPENTINITE_SPRINGS_overlap.062        k__Bacteria (UID1452)            924         151           101         68    78    5     0    0    0       52.97            4.46               0.00
  SERPENTINITE_SPRINGS_overlap.042        k__Bacteria (UID2329)            174         149            89         60    69    20    0    0    0       52.08           12.20              40.00
  SERPENTINITE_SPRINGS_overlap.040         k__Bacteria (UID203)            5449        103            58         33    37    20    13   0    0       51.72           15.02               1.69
  SERPENTINITE_SPRINGS_overlap.039         k__Bacteria (UID203)            5449        100            55         40    32    14    12   2    0       46.02           33.64              41.94
  SERPENTINITE_SPRINGS_overlap.061         k__Bacteria (UID203)            5449        104            58         51    33    15    5    0    0       44.69           17.43               0.00
  SERPENTINITE_SPRINGS_overlap.015         k__Bacteria (UID203)            5449        104            58         49    51    3     1    0    0       43.05            6.03              16.67
  SERPENTINITE_SPRINGS_overlap.029         k__Bacteria (UID203)            5449        104            58         73    25    4     2    0    0       38.20            5.96               0.00
  SERPENTINITE_SPRINGS_overlap.031        k__Bacteria (UID2570)            433         274           183        156    89    23    5    1    0       35.70            6.08               2.27
  SERPENTINITE_SPRINGS_overlap.083         c__Bacilli (UID259)             750         273           152        158   107    8     0    0    0       35.00            0.97              25.00
  SERPENTINITE_SPRINGS_overlap.035        k__Bacteria (UID2495)            2993        143            89         91    35    14    3    0    0       33.74           12.88              34.78
  SERPENTINITE_SPRINGS_overlap.084      o__Cytophagales (UID2938)           28         439           266        275   148    13    2    1    0       31.85            4.43               4.00
  SERPENTINITE_SPRINGS_overlap.075        k__Bacteria (UID2329)            174         149            89         93    37    18    1    0    0       31.48            8.92               9.52
  SERPENTINITE_SPRINGS_overlap.007         k__Bacteria (UID203)            5449        104            58         85    18    1     0    0    0       21.69            1.72              100.00
  SERPENTINITE_SPRINGS_overlap.003         k__Bacteria (UID203)            5449        104            58         79    25    0     0    0    0       21.54            0.00               0.00
  SERPENTINITE_SPRINGS_overlap.026             root (UID1)                 5656         56            24         51    5     0     0    0    0       20.83            0.00               0.00
  SERPENTINITE_SPRINGS_overlap.008         k__Bacteria (UID203)            5449        104            58         90    14    0     0    0    0       15.82            0.00               0.00
  SERPENTINITE_SPRINGS_overlap.097         k__Bacteria (UID203)            5449        104            58         91    13    0     0    0    0       14.66            0.00               0.00
  SERPENTINITE_SPRINGS_overlap.080             root (UID1)                 5656         56            24         51    3     2     0    0    0       14.58            8.33               0.00
  SERPENTINITE_SPRINGS_overlap.106         k__Bacteria (UID203)            5449        104            58         95    7     2     0    0    0       13.79            3.45              100.00
  SERPENTINITE_SPRINGS_overlap.011         k__Bacteria (UID203)            5449        104            58         90    12    1     1    0    0        7.73            0.47               0.00
  SERPENTINITE_SPRINGS_overlap.012         k__Bacteria (UID203)            5449        104            58         98    6     0     0    0    0        7.21            0.00               0.00
  SERPENTINITE_SPRINGS_overlap.107         k__Bacteria (UID203)            5449        103            58         94    6     3     0    0    0        6.90            1.29              100.00
  SERPENTINITE_SPRINGS_overlap.101         k__Bacteria (UID203)            5449        104            58         91    12    1     0    0    0        6.74            1.72               0.00
  SERPENTINITE_SPRINGS_overlap.034   c__Gammaproteobacteria (UID4274)      112         580           289        516    60    4     0    0    0        5.10            0.22              25.00
  SERPENTINITE_SPRINGS_overlap.019         k__Bacteria (UID203)            5449        104            58        100    4     0     0    0    0        0.63            0.00               0.00
  SERPENTINITE_SPRINGS_overlap.105             root (UID1)                 5656         56            24         56    0     0     0    0    0        0.00            0.00               0.00
  SERPENTINITE_SPRINGS_overlap.048             root (UID1)                 5656         56            24         56    0     0     0    0    0        0.00            0.00               0.00
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source deactivate
```
**> Time: ~30 minutes**

We filtered out other than high-quality bins. We considered high-quality bins those meeting the criteria: **Completeness ≥ 70 && Contamination ≤ 10**. The 24 (2 bins less than in the `minoche` protocol) selected bins were:

```sh
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Bin Id                                      Marker lineage            # genomes   # markers   # marker sets    0     1     2     3    4    5+   Completeness   Contamination   Strain heterogeneity
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  SERPENTINITE_SPRINGS_overlap.027   c__Gammaproteobacteria (UID4274)      112         580           289         26   539    15    0    0    0       98.10            2.60              60.00
  SERPENTINITE_SPRINGS_overlap.004   c__Deltaproteobacteria (UID3217)       62         280           168         6    272    2     0    0    0       97.00            0.89               0.00
  SERPENTINITE_SPRINGS_overlap.091    o__Sphingomonadales (UID3310)         26         569           293         18   529    22    0    0    0       96.75            2.80              31.82
  SERPENTINITE_SPRINGS_overlap.065   c__Gammaproteobacteria (UID4274)      112         580           289         53   517    10    0    0    0       96.50            1.78              60.00
  SERPENTINITE_SPRINGS_overlap.006        k__Bacteria (UID3187)            2258        181           110         5    167    8     1    0    0       96.31            3.43               9.09
  SERPENTINITE_SPRINGS_overlap.103      p__Bacteroidetes (UID2591)         364         303           203         16   276    11    0    0    0       95.21            2.19              72.73
  SERPENTINITE_SPRINGS_overlap.096     o__Burkholderiales (UID4001)        108         570           250         33   518    18    1    0    0       93.42            3.57              42.86
  SERPENTINITE_SPRINGS_overlap.056      p__Bacteroidetes (UID2605)         350         316           210         18   261    35    2    0    0       93.00            9.42              75.61
  SERPENTINITE_SPRINGS_overlap.018     p__Actinobacteria (UID1454)         732         206           120         13   185    8     0    0    0       91.67            4.06              62.50
  SERPENTINITE_SPRINGS_overlap.058        k__Bacteria (UID2329)            174         149            89         9    127    8     2    3    0       91.57            8.46              12.50
  SERPENTINITE_SPRINGS_overlap.072        k__Bacteria (UID2570)            433         274           183         24   240    10    0    0    0       91.38            2.47              50.00
  SERPENTINITE_SPRINGS_overlap.094      p__Bacteroidetes (UID2591)         364         303           203         62   227    14    0    0    0       91.11            3.97              71.43
  SERPENTINITE_SPRINGS_overlap.013         k__Bacteria (UID209)            5443        105            59         6     88    11    0    0    0       89.83            9.94               0.00
  SERPENTINITE_SPRINGS_overlap.041    f__Xanthomonadaceae (UID4214)         55         659           290         68   545    45    1    0    0       88.41            7.44              58.33
  SERPENTINITE_SPRINGS_overlap.073        k__Bacteria (UID2372)            131         177           106         46   122    9     0    0    0       86.32            4.75              11.11
  SERPENTINITE_SPRINGS_overlap.088        k__Bacteria (UID2569)            434         275           184         31   231    11    2    0    0       85.56            4.37              17.65
  SERPENTINITE_SPRINGS_overlap.064        p__Firmicutes (UID241)           930         213           118         43   154    16    0    0    0       83.69            8.86              50.00
  SERPENTINITE_SPRINGS_overlap.046        p__Firmicutes (UID241)           930         213           118         25   174    13    1    0    0       83.05            6.84              25.00
  SERPENTINITE_SPRINGS_overlap.010        k__Bacteria (UID1452)            924         151           101         20   127    4     0    0    0       82.34            2.07              50.00
  SERPENTINITE_SPRINGS_overlap.090    f__Flavobacteriaceae (UID2817)        81         511           283        153   321    32    5    0    0       78.25            8.55              36.17
  SERPENTINITE_SPRINGS_overlap.059        p__Firmicutes (UID241)           930         213           118         55   142    15    1    0    0       75.19            6.64              38.89
  SERPENTINITE_SPRINGS_overlap.066       c__Clostridia (UID1118)           387         223           124         83   125    15    0    0    0       75.17            6.72               0.00
  SERPENTINITE_SPRINGS_overlap.093     o__Burkholderiales (UID4001)        108         570           250        139   385    46    0    0    0       73.12            7.93              78.26
  SERPENTINITE_SPRINGS_overlap.043       c__Clostridia (UID1118)           387         223           124         63   139    19    2    0    0       72.07            7.65              60.00
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
```

## [Stage 5] Functional annotation - `prokka`

Annotating high-quality bins:

```sh
jv~$ mkdir 06_ANNOTATION_overlap/
jv~$ prokka --outdir 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.027 --prefix SERPENTINITE_SPRINGS_overlap.027 --kingdom Bacteria --metagenome --cpus 32 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.027.fasta
jv~$ prokka --outdir 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.004 --prefix SERPENTINITE_SPRINGS_overlap.004 --kingdom Bacteria --metagenome --cpus 32 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.004.fasta
jv~$ prokka --outdir 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.091 --prefix SERPENTINITE_SPRINGS_overlap.091 --kingdom Bacteria --metagenome --cpus 32 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.091.fasta
jv~$ prokka --outdir 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.065 --prefix SERPENTINITE_SPRINGS_overlap.065 --kingdom Bacteria --metagenome --cpus 32 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.065.fasta
jv~$ prokka --outdir 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.006 --prefix SERPENTINITE_SPRINGS_overlap.006 --kingdom Bacteria --metagenome --cpus 32 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.006.fasta
jv~$ prokka --outdir 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.103 --prefix SERPENTINITE_SPRINGS_overlap.103 --kingdom Bacteria --metagenome --cpus 32 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.103.fasta
jv~$ prokka --outdir 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.096 --prefix SERPENTINITE_SPRINGS_overlap.096 --kingdom Bacteria --metagenome --cpus 32 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.096.fasta
jv~$ prokka --outdir 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.072 --prefix SERPENTINITE_SPRINGS_overlap.072 --kingdom Bacteria --metagenome --cpus 32 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.072.fasta
jv~$ prokka --outdir 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.056 --prefix SERPENTINITE_SPRINGS_overlap.056 --kingdom Bacteria --metagenome --cpus 32 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.056.fasta
jv~$ prokka --outdir 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.018 --prefix SERPENTINITE_SPRINGS_overlap.018 --kingdom Bacteria --metagenome --cpus 32 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.018.fasta
jv~$ prokka --outdir 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.058 --prefix SERPENTINITE_SPRINGS_overlap.058 --kingdom Bacteria --metagenome --cpus 32 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.058.fasta
jv~$ prokka --outdir 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.094 --prefix SERPENTINITE_SPRINGS_overlap.094 --kingdom Bacteria --metagenome --cpus 32 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.094.fasta
jv~$ prokka --outdir 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.013 --prefix SERPENTINITE_SPRINGS_overlap.013 --kingdom Bacteria --metagenome --cpus 32 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.013.fasta
jv~$ prokka --outdir 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.041 --prefix SERPENTINITE_SPRINGS_overlap.041 --kingdom Bacteria --metagenome --cpus 32 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.041.fasta
jv~$ prokka --outdir 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.073 --prefix SERPENTINITE_SPRINGS_overlap.073 --kingdom Bacteria --metagenome --cpus 32 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.073.fasta
jv~$ prokka --outdir 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.088 --prefix SERPENTINITE_SPRINGS_overlap.088 --kingdom Bacteria --metagenome --cpus 32 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.088.fasta
jv~$ prokka --outdir 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.064 --prefix SERPENTINITE_SPRINGS_overlap.064 --kingdom Bacteria --metagenome --cpus 32 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.064.fasta
jv~$ prokka --outdir 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.046 --prefix SERPENTINITE_SPRINGS_overlap.046 --kingdom Bacteria --metagenome --cpus 32 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.046.fasta
jv~$ prokka --outdir 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.010 --prefix SERPENTINITE_SPRINGS_overlap.010 --kingdom Bacteria --metagenome --cpus 32 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.010.fasta
jv~$ prokka --outdir 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.090 --prefix SERPENTINITE_SPRINGS_overlap.090 --kingdom Bacteria --metagenome --cpus 32 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.090.fasta
jv~$ prokka --outdir 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.059 --prefix SERPENTINITE_SPRINGS_overlap.059 --kingdom Bacteria --metagenome --cpus 32 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.059.fasta
jv~$ prokka --outdir 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.066 --prefix SERPENTINITE_SPRINGS_overlap.066 --kingdom Bacteria --metagenome --cpus 32 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.066.fasta
jv~$ prokka --outdir 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.093 --prefix SERPENTINITE_SPRINGS_overlap.093 --kingdom Bacteria --metagenome --cpus 32 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.093.fasta
jv~$ prokka --outdir 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.043 --prefix SERPENTINITE_SPRINGS_overlap.043 --kingdom Bacteria --metagenome --cpus 32 04_BINNING_overlap/SERPENTINITE_SPRINGS_overlap.043.fasta
```

## [Stage 6] Taxonomy classification of bins

Input for `phylophlan` should be bins in `.faa` format

Create a directory to save the results
```sh
jv~$ mkdir 07_TAXONOMY_overlap/
```

Copy the bins in amino acids `.faa` format from the `prokka` output:
```sh
jv~$ cp 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.027/SERPENTINITE_SPRINGS_overlap.027.faa 07_TAXONOMY_overlap/
jv~$ cp 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.004/SERPENTINITE_SPRINGS_overlap.004.faa 07_TAXONOMY_overlap/
jv~$ cp 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.091/SERPENTINITE_SPRINGS_overlap.091.faa 07_TAXONOMY_overlap/
jv~$ cp 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.065/SERPENTINITE_SPRINGS_overlap.065.faa 07_TAXONOMY_overlap/
jv~$ cp 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.006/SERPENTINITE_SPRINGS_overlap.006.faa 07_TAXONOMY_overlap/
jv~$ cp 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.103/SERPENTINITE_SPRINGS_overlap.103.faa 07_TAXONOMY_overlap/
jv~$ cp 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.096/SERPENTINITE_SPRINGS_overlap.096.faa 07_TAXONOMY_overlap/
jv~$ cp 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.072/SERPENTINITE_SPRINGS_overlap.072.faa 07_TAXONOMY_overlap/
jv~$ cp 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.056/SERPENTINITE_SPRINGS_overlap.056.faa 07_TAXONOMY_overlap/
jv~$ cp 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.018/SERPENTINITE_SPRINGS_overlap.018.faa 07_TAXONOMY_overlap/
jv~$ cp 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.058/SERPENTINITE_SPRINGS_overlap.058.faa 07_TAXONOMY_overlap/
jv~$ cp 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.094/SERPENTINITE_SPRINGS_overlap.094.faa 07_TAXONOMY_overlap/
jv~$ cp 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.013/SERPENTINITE_SPRINGS_overlap.013.faa 07_TAXONOMY_overlap/
jv~$ cp 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.041/SERPENTINITE_SPRINGS_overlap.041.faa 07_TAXONOMY_overlap/
jv~$ cp 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.073/SERPENTINITE_SPRINGS_overlap.073.faa 07_TAXONOMY_overlap/
jv~$ cp 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.088/SERPENTINITE_SPRINGS_overlap.088.faa 07_TAXONOMY_overlap/
jv~$ cp 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.064/SERPENTINITE_SPRINGS_overlap.064.faa 07_TAXONOMY_overlap/
jv~$ cp 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.046/SERPENTINITE_SPRINGS_overlap.046.faa 07_TAXONOMY_overlap/
jv~$ cp 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.010/SERPENTINITE_SPRINGS_overlap.010.faa 07_TAXONOMY_overlap/
jv~$ cp 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.090/SERPENTINITE_SPRINGS_overlap.090.faa 07_TAXONOMY_overlap/
jv~$ cp 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.059/SERPENTINITE_SPRINGS_overlap.059.faa 07_TAXONOMY_overlap/
jv~$ cp 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.066/SERPENTINITE_SPRINGS_overlap.066.faa 07_TAXONOMY_overlap/
jv~$ cp 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.093/SERPENTINITE_SPRINGS_overlap.093.faa 07_TAXONOMY_overlap/
jv~$ cp 06_ANNOTATION_overlap/SERPENTINITE_SPRINGS_overlap.043/SERPENTINITE_SPRINGS_overlap.043.faa 07_TAXONOMY_overlap/
jv~$ ls -1 07_TAXONOMY_overlap/
jv~$ ls -1 07_TAXONOMY_overlap/ | wc -l
24
```

Go to the `phylophlan` directory, then create a new folder inside the `input` folder containing the bins in `.faa` format:

```sh
jv~$ cd ~/programs/phylophlan/
jv~$ source activate phylophlan
jv~$ mkdir input/07_TAXONOMY_overlap-serpentinite-springs
jv~$ cp ~/data/metagenomics-serpentinite-springs/07_TAXONOMY_overlap/* input/07_TAXONOMY_overlap-serpentinite-springs/
jv~$ ls -1 input/07_TAXONOMY_overlap-serpentinite-springs/ | wc -l
24
```

```sh
jv~$ ./phylophlan.py -i -t 07_TAXONOMY_overlap-serpentinite-springs --nproc 32
Total time: 2862.76 seconds Unique: 3163/3195 Bad splits: 5/3160 Worst delta-LogLk 0.385
Tree built! The output newick file is in output/07_TAXONOMY_overlap-serpentinite-springs/07_TAXONOMY_overlap-serpentinite-springs.tree.int.nwk
Tree building finished in 3234 secs (4367 total time).

Trying to impute taxonomic labels for taxa newly integrated into the tree... Done!
Writing taxonomic imputation outputs ... Done!
jv~$ source deactivate
```

Results with high confidence are in `imputed_conf_high-conf.txt`
```sh
jv~$ cat output/07_TAXONOMY_overlap-serpentinite-springs/imputed_conf_high-conf.txt
SERPENTINITE_SPRINGS_overlap.027	d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Methylococcales.f__Methylococcaceae.g__?.s__?.t__?
SERPENTINITE_SPRINGS_overlap.010	d__Bacteria.p__Chloroflexi.c__Dehalococcoidetes.o__Dehalococcoidales.f__?.g__?.s__?.t__?
SERPENTINITE_SPRINGS_overlap.004	d__Bacteria.p__Proteobacteria.c__Deltaproteobacteria.o__?.f__?.g__?.s__?.t__?
SERPENTINITE_SPRINGS_overlap.006	d__Bacteria.p__Nitrospirae.c__Nitrospira.o__Nitrospirales.f__Nitrospiraceae.g__?.s__?.t__?
SERPENTINITE_SPRINGS_overlap.103	d__Bacteria.p__Bacteroidetes.c__Sphingobacteria.o__Sphingobacteriales.f__?.g__?.s__?.t__?
SERPENTINITE_SPRINGS_overlap.093	d__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Burkholderiales.f__Burkholderiaceae.g__Polynucleobacter.s__?.t__?
SERPENTINITE_SPRINGS_overlap.090	d__Bacteria.p__Bacteroidetes.c__Flavobacteria.o__Flavobacteriales.f__Flavobacteriaceae.g__Flavobacterium.s__?.t__?
SERPENTINITE_SPRINGS_overlap.096	d__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Burkholderiales.f__Burkholderiaceae.g__Polynucleobacter.s__?.t__?
SERPENTINITE_SPRINGS_overlap.094	d__Bacteria.p__Bacteroidetes.c__Sphingobacteria.o__Sphingobacteriales.f__?.g__?.s__?.t__?
SERPENTINITE_SPRINGS_overlap.072	d__Bacteria.p__Bacteroidetes.c__Cytophagia.o__Cytophagales.f__?.g__?.s__?.t__?
SERPENTINITE_SPRINGS_overlap.073	d__Bacteria.p__Firmicutes.c__Erysipelotrichi.o__Erysipelotrichales.f__Erysipelotrichaceae.g__Erysipelothrix.s__rhusiopathiae.t__?
SERPENTINITE_SPRINGS_overlap.065	d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Methylococcales.f__Methylococcaceae.g__?.s__?.t__?
```

At this stage we could observe that the main difference in the taxonomy classfication with high confidence is that with the `minoche` method were detected 3 other c__Clostridia bins:
```sh
SERPENTINITE_SPRINGS.060	d__Bacteria.p__Firmicutes.c__Clostridia.o__?.f__?.g__?.s__?.t__?
SERPENTINITE_SPRINGS.071	d__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Clostridiaceae.g__Alkaliphilus.s__?.t__?
SERPENTINITE_SPRINGS.043	d__Bacteria.p__Firmicutes.c__Clostridia.o__?.f__?.g__?.s__?.t__?
```

While with the `merge` (overlap) method were classified 3 p__Bacteroidetes with more detail:
```sh
SERPENTINITE_SPRINGS_overlap.103	d__Bacteria.p__Bacteroidetes.c__Sphingobacteria.o__Sphingobacteriales.f__?.g__?.s__?.t__?
SERPENTINITE_SPRINGS_overlap.090	d__Bacteria.p__Bacteroidetes.c__Flavobacteria.o__Flavobacteriales.f__Flavobacteriaceae.g__Flavobacterium.s__?.t__?
SERPENTINITE_SPRINGS_overlap.094	d__Bacteria.p__Bacteroidetes.c__Sphingobacteria.o__Sphingobacteriales.f__?.g__?.s__?.t__?
```

Results with medium confidence are in `imputed_conf_medium-conf.txt`

```sh
jv~$ cat output/07_TAXONOMY_overlap-serpentinite-springs/imputed_conf_medium-conf.txt
SERPENTINITE_SPRINGS_overlap.043	d__Bacteria.p__Firmicutes.c__Clostridia.o__Thermoanaerobacterales.f__?.g__?.s__?.t__?
SERPENTINITE_SPRINGS_overlap.018	d__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Coriobacteriales.f__Coriobacteriaceae.g__?.s__?.t__?
SERPENTINITE_SPRINGS_overlap.091	d__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Sphingomonadales.f__?.g__?.s__?.t__?
SERPENTINITE_SPRINGS_overlap.066	d__Bacteria.p__Firmicutes.c__Clostridia.o__Thermoanaerobacterales.f__?.g__?.s__?.t__?
SERPENTINITE_SPRINGS_overlap.058	d__Bacteria.p__Firmicutes.c__Erysipelotrichi.o__Erysipelotrichales.f__Erysipelotrichaceae.g__?.s__?.t__?
```

The main difference in the taxonomy classfication with medium confidence is that with the `minoche` method were detected 2 other bins:
```sh
SERPENTINITE_SPRINGS.087	d__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__?.f__?.g__?.s__?.t__?
SERPENTINITE_SPRINGS.033	d__Bacteria.p__Thermi.c__?.o__?.f__?.g__?.s__?.t__?
```

While with the `merge` (overlap) method were classified 2 o__Thermoanaerobacterales bins with high detail:
```sh
SERPENTINITE_SPRINGS_overlap.043	d__Bacteria.p__Firmicutes.c__Clostridia.o__Thermoanaerobacterales.f__?.g__?.s__?.t__?
SERPENTINITE_SPRINGS_overlap.066	d__Bacteria.p__Firmicutes.c__Clostridia.o__Thermoanaerobacterales.f__?.g__?.s__?.t__?
```
**> Total time: ~ 1.5 h**

Altogether, based on the `phylophlan` classification, it seems the `merge` method performed better.

Copy the `phylophlan` results from `$DIR_PHYLOPHLAN` to the `$DIR_PROJECT`
```sh 
jv~$ cp output/07_TAXONOMY_overlap-serpentinite-springs/* ~/data/metagenomics-serpentinite-springs/07_TAXONOMY_overlap/
```

It seems that Bins 027 and 065 correspond to methanotrophic bacteria. Both bins were classified with high cofidence by `phylophlan`, and their recovered genomes are highly complete with low contamination according to `CheckM`:

**SERPENTINITE_SPRINGS_overlap.027**:
```
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Bin Id                              Marker lineage            # genomes   # markers   # marker sets    0     1     2    3    4   5+   Completeness   Contamination   Strain heterogeneity
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  SERPENTINITE_SPRINGS_overlap.027   c__Gammaproteobacteria (UID4274)      112         580           289         26   539    15    0    0    0       98.10            2.60              60.00

SERPENTINITE_SPRINGS_overlap.027	d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Methylococcales.f__Methylococcaceae.g__?.s__?.t__?
```

**SERPENTINITE_SPRINGS_overlap.065**:
```sh
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Bin Id                              Marker lineage            # genomes   # markers   # marker sets    0     1     2    3    4   5+   Completeness   Contamination   Strain heterogeneity
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  SERPENTINITE_SPRINGS_overlap.065   c__Gammaproteobacteria (UID4274)      112         580           289         53   517    10    0    0    0       96.50            1.78              60.00

SERPENTINITE_SPRINGS_overlap.065	d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Methylococcales.f__Methylococcaceae.g__?.s__?.t__?
```

## [Stage 7]  Refine functional annotation of methanotrophic bacteria using `EggNOG-mapper`
`EggNOG-mapper v4.5.1` works on `python2.7`

`$DIR_EGGNOG` is the directory to `eggnog-mapper/emapper.py` script.

```sh
jv~$ cd ~/data/metagenomics-serpentinite-springs/
jv~$ mkdir 08_ANNOTATION_REFINE_overlap/
jv~$ cd 08_ANNOTATION_REFINE_overlap/
# f__Methylococcaceae
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY_overlap/SERPENTINITE_SPRINGS_overlap.027.faa --output SERPENTINITE_SPRINGS_overlap.027 -m diamond --cpu 32
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY_overlap/SERPENTINITE_SPRINGS_overlap.065.faa --output SERPENTINITE_SPRINGS_overlap.065 -m diamond --cpu 32
# other bins
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY_overlap/SERPENTINITE_SPRINGS_overlap.006.faa --output SERPENTINITE_SPRINGS_overlap.006 -m diamond --cpu 32
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY_overlap/SERPENTINITE_SPRINGS_overlap.004.faa --output SERPENTINITE_SPRINGS_overlap.004 -m diamond --cpu 32
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY_overlap/SERPENTINITE_SPRINGS_overlap.091.faa --output SERPENTINITE_SPRINGS_overlap.091 -m diamond --cpu 32
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY_overlap/SERPENTINITE_SPRINGS_overlap.103.faa --output SERPENTINITE_SPRINGS_overlap.103 -m diamond --cpu 32
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY_overlap/SERPENTINITE_SPRINGS_overlap.096.faa --output SERPENTINITE_SPRINGS_overlap.096 -m diamond --cpu 32
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY_overlap/SERPENTINITE_SPRINGS_overlap.056.faa --output SERPENTINITE_SPRINGS_overlap.056 -m diamond --cpu 32
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY_overlap/SERPENTINITE_SPRINGS_overlap.018.faa --output SERPENTINITE_SPRINGS_overlap.018 -m diamond --cpu 32
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY_overlap/SERPENTINITE_SPRINGS_overlap.058.faa --output SERPENTINITE_SPRINGS_overlap.058 -m diamond --cpu 32
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY_overlap/SERPENTINITE_SPRINGS_overlap.072.faa --output SERPENTINITE_SPRINGS_overlap.072 -m diamond --cpu 32
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY_overlap/SERPENTINITE_SPRINGS_overlap.094.faa --output SERPENTINITE_SPRINGS_overlap.094 -m diamond --cpu 32
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY_overlap/SERPENTINITE_SPRINGS_overlap.013.faa --output SERPENTINITE_SPRINGS_overlap.013 -m diamond --cpu 32
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY_overlap/SERPENTINITE_SPRINGS_overlap.041.faa --output SERPENTINITE_SPRINGS_overlap.041 -m diamond --cpu 32
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY_overlap/SERPENTINITE_SPRINGS_overlap.073.faa --output SERPENTINITE_SPRINGS_overlap.073 -m diamond --cpu 32
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY_overlap/SERPENTINITE_SPRINGS_overlap.088.faa --output SERPENTINITE_SPRINGS_overlap.088 -m diamond --cpu 32
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY_overlap/SERPENTINITE_SPRINGS_overlap.064.faa --output SERPENTINITE_SPRINGS_overlap.064 -m diamond --cpu 32
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY_overlap/SERPENTINITE_SPRINGS_overlap.046.faa --output SERPENTINITE_SPRINGS_overlap.046 -m diamond --cpu 32
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY_overlap/SERPENTINITE_SPRINGS_overlap.010.faa --output SERPENTINITE_SPRINGS_overlap.010 -m diamond --cpu 32
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY_overlap/SERPENTINITE_SPRINGS_overlap.090.faa --output SERPENTINITE_SPRINGS_overlap.090 -m diamond --cpu 32
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY_overlap/SERPENTINITE_SPRINGS_overlap.059.faa --output SERPENTINITE_SPRINGS_overlap.059 -m diamond --cpu 32
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY_overlap/SERPENTINITE_SPRINGS_overlap.066.faa --output SERPENTINITE_SPRINGS_overlap.066 -m diamond --cpu 32
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY_overlap/SERPENTINITE_SPRINGS_overlap.093.faa --output SERPENTINITE_SPRINGS_overlap.093 -m diamond --cpu 32
python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY_overlap/SERPENTINITE_SPRINGS_overlap.043.faa --output SERPENTINITE_SPRINGS_overlap.043 -m diamond --cpu 32
```
