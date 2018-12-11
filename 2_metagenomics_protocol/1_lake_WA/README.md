# Metagenomics from Lake Washington

<!-- TOC -->

- [Metagenomics from Lake Washington](#metagenomics-from-lake-washington)
    - [[Stage 0] Preliminary](#stage-0-preliminary)
    - [[Stage 1] Quality Control (QC)](#stage-1-quality-control-qc)
    - [[Stage 2] Co-assembly](#stage-2-co-assembly)
    - [[Stage 3] Binning](#stage-3-binning)
    - [[Stage 4] Refine bins](#stage-4-refine-bins)
    - [[Stage 5] Functional annotation - `prokka`](#stage-5-functional-annotation---prokka)
    - [[Stage 6] Taxonomy classification of bins](#stage-6-taxonomy-classification-of-bins)
    - [[Stage 7] Refine functional annotation - `EggNOG-mapper`](#stage-7-refine-functional-annotation---eggnog-mapper)
- [Results](#results)

<!-- /TOC -->

## [Stage 0] Preliminary
### Create project directory
```sh
jv~$ mkdir metagenomics_lake_WA
jv~$ cd metagenomics_lake_WA/
```

### Retrieve sequences from SRA NCBI system

> This dataset was published in:
> Oshkin, I. Y., Beck, D. A., Lamb, A. E., Tchesnokova, V., Benuska, G., McTaggart, T. L., ... & Chistoserdova, L. (2015). Methane-fed microbial microcosms show differential community dynamics and pinpoint taxa involved in communal response. _The ISME journal_, 9(5), 1119-1129.
> https://www.nature.com/articles/ismej2014203

### Download metagenomics data

This code will retrieve the reads in separate ``FASTQ`` files:

```sh
jv~$ mkdir 00_RAW
jv~$ cd 00_RAW/
jv~$ fastq-dump -I --split-files SRR1505166 SRR1505167 SRR1505168 SRR1505169 SRR1505170 SRR1505171 SRR1505172 SRR1505173
jv~$ ls 00_RAW/
```

### Check number of reads

#### Sample 1
```sh
jv~$ head SRR1505166_1.fastq
# check the header "@SRR1505166.1.1"

jv~$ grep -c "@SRR1505166" SRR1505166_1.fastq
1827066

jv~$ grep -c "@SRR1505166" SRR1505166_2.fastq
1827066
```

#### Sample 2
```sh
jv~$ head  SRR1505167_1.fastq
# @SRR1505167

jv~$ grep -c "@SRR1505167" SRR1505167_1.fastq
1712886
jv~$ grep -c "@SRR1505167" SRR1505167_2.fastq
1712886
```

#### Sample 3
```sh
jv~$ grep -c "@SRR1505168" SRR1505168_1.fastq
1685155
jv~$ grep -c "@SRR1505168" SRR1505168_2.fastq
1685155
```

#### Sample 4
```sh
jv~$ grep -c "@SRR1505169" SRR1505169_1.fastq
2360078
jv~$ grep -c "@SRR1505169" SRR1505169_2.fastq
2360078
```

#### Sample 5
```sh
jv~$ grep -c "@SRR1505170" SRR1505170_1.fastq
1351287
jv~$ grep -c "@SRR1505170" SRR1505170_2.fastq
1351287
```

#### Sample 6
```sh
jv~$ grep -c "@SRR1505171" SRR1505171_1.fastq
2864984
jv~$ grep -c "@SRR1505171" SRR1505171_2.fastq
2864984
```

#### Sample 7
```sh
jv~$ grep -c "@SRR1505172" SRR1505172_1.fastq
2652260
jv~$ grep -c "@SRR1505172" SRR1505172_2.fastq
2652260
```

#### Sample 8
```sh
jv~$ grep -c "@SRR1505173" SRR1505173_1.fastq
2101967
jv~$ grep -c "@SRR1505173" SRR1505173_2.fastq
2101967
```

#### Summary of reads
Generate some stats of the samples in `R`

```sh
jv~$ R
```

```R
x <- c(1827066, 1712886, 1685155, 2360078, 1351287, 2864984, 2652260, 2101967)
mean(x)
[1] 2069460
sd(x)
[1] 522000.6
```

|#|Sample    |reads  |
|-|----------|-------|
|1|SRR1505166|1827066|
|2|SRR1505167|1712886|
|3|SRR1505168|1685155|
|4|SRR1505170|2360078|
|5|SRR1505169|1351287|
|6|SRR1505171|2864984|
|7|SRR1505172|2652260|
|8|SRR1505173|2101967|
|**Average**|**2069460**|**~2M reads**|
|Std.Dev|522000.6|~0.5M reads|

> Find a recent dataset with more reads! (cats feces 20-90M reads)

## [Stage 1] Quality Control (QC)

### Quality Filtering

> This protocol section is partially based on: http://merenlab.org/tutorials/assembly-based-metagenomics/

Generate a TAB-delimited `samples.txt` file to point out where are your raw `R1` and `R2` files for each sample.

|sample    |r1                         |r2                          |
|----------|---------------------------|----------------------------|
|Sample_01 | 00_RAW/SRR1505166_1.fastq | 00_RAW/SRR1505166_2.fastq  |
|Sample_02 | 00_RAW/SRR1505167_1.fastq | 00_RAW/SRR1505167_2.fastq  |
|Sample_03 | 00_RAW/SRR1505168_1.fastq | 00_RAW/SRR1505168_2.fastq  |
|Sample_04 | 00_RAW/SRR1505169_1.fastq | 00_RAW/SRR1505169_2.fastq  |
|Sample_05 | 00_RAW/SRR1505170_1.fastq | 00_RAW/SRR1505170_2.fastq  |
|Sample_06 | 00_RAW/SRR1505171_1.fastq | 00_RAW/SRR1505171_2.fastq  |
|Sample_07 | 00_RAW/SRR1505172_1.fastq | 00_RAW/SRR1505172_2.fastq  |
|Sample_08 | 00_RAW/SRR1505173_1.fastq | 00_RAW/SRR1505173_2.fastq  |

Create a directory for quality-filtered `R1` and `R2`
```sh
jv~$ mkdir 01_QC

jv~$ source activate illumina-utils

jv~$ iu-gen-configs samples.txt -o 01_QC
Report .......................................: Read for 8 samples is read
Output directory set in configs ..............: /metagenomics_lake_WA/01_QC
Prefix for R1 ................................: None
Prefix for R2 ................................: None

jv~$ ls 01_QC/
Sample_01.ini  Sample_02.ini  Sample_03.ini  Sample_04.ini  Sample_05.ini  Sample_06.ini  Sample_07.ini  Sample_08.ini
```

Run quality filtering for all your samples at once:

> Note that we had to include the `parameter` `--ignore-deflines`
```sh
jv~$ for ini in 01_QC/*.ini; do iu-filter-quality-minoche --ignore-deflines $ini; done
91% -- (num pairs processed: 1,827,000)
Read ID tracker dict is being stored ...
 91% -- (num pairs processed: 1,712,000)
Read ID tracker dict is being stored ...
 91% -- (num pairs processed: 1,685,000)
Read ID tracker dict is being stored ...
 91% -- (num pairs processed: 2,360,000)
Read ID tracker dict is being stored ...
 90% -- (num pairs processed: 1,351,000)
Read ID tracker dict is being stored ...
 91% -- (num pairs processed: 2,864,000)
Read ID tracker dict is being stored ...
 91% -- (num pairs processed: 2,652,000)
Read ID tracker dict is being stored ...
 90% -- (num pairs processed: 2,101,000)
Read ID tracker dict is being stored ...
```
Ref to the method used:
> Minoche, A. E., Dohm, J. C., & Himmelbauer, H. (2011). Evaluation of genomic high-throughput sequencing data generated on Illumina HiSeq and genome analyzer systems. *Genome biology*, 12(11), R112.

The contents of the `01_QC/` directory should look like this:
```sh
jv~$ ls 01_QC/
```

#### Checking QC results stats in `*_STATS.txt` files:

##### Sample_01-STATS
```sh
jv~$ nano 01_QC/Sample_01-STATS.txt

number of pairs analyzed      : 1827066
total pairs passed            : 1801218 (%98.59 of all pairs)
  total pair_1 trimmed        : 0 (%0.00 of all passed pairs)
  total pair_2 trimmed        : 0 (%0.00 of all passed pairs)
total pairs failed            : 25848 (%1.41 of all pairs)
  pairs failed due to pair_1  : 8763 (%33.90 of all failed pairs)
  pairs failed due to pair_2  : 12893 (%49.88 of all failed pairs)
  pairs failed due to both    : 4192 (%16.22 of all failed pairs)
  FAILED_REASON_C33           : 23019 (%89.06 of all failed pairs)
  FAILED_REASON_P             : 2829 (%10.94 of all failed pairs)
```

##### Sample_02-STATS
```sh
jv~$ nano 01_QC/Sample_02-STATS.txt

number of pairs analyzed      : 1712886
total pairs passed            : 1628445 (%95.07 of all pairs)
  total pair_1 trimmed        : 0 (%0.00 of all passed pairs)
  total pair_2 trimmed        : 0 (%0.00 of all passed pairs)
total pairs failed            : 84441 (%4.93 of all pairs)
  pairs failed due to pair_1  : 1577 (%1.87 of all failed pairs)
  pairs failed due to pair_2  : 78032 (%92.41 of all failed pairs)
  pairs failed due to both    : 4832 (%5.72 of all failed pairs)
  FAILED_REASON_P             : 1055 (%1.25 of all failed pairs)
  FAILED_REASON_C33           : 83384 (%98.75 of all failed pairs)
  FAILED_REASON_N             : 2 (%0.00 of all failed pairs)
```

##### Sample_03-STATS
```sh
jv~$ 

number of pairs analyzed      : 1685155
total pairs passed            : 1595021 (%94.65 of all pairs)
  total pair_1 trimmed        : 0 (%0.00 of all passed pairs)
  total pair_2 trimmed        : 0 (%0.00 of all passed pairs)
total pairs failed            : 90134 (%5.35 of all pairs)
  pairs failed due to pair_1  : 1382 (%1.53 of all failed pairs)
  pairs failed due to pair_2  : 81934 (%90.90 of all failed pairs)
  pairs failed due to both    : 6818 (%7.56 of all failed pairs)
  FAILED_REASON_N             : 3 (%0.00 of all failed pairs)
  FAILED_REASON_P             : 3155 (%3.50 of all failed pairs)
  FAILED_REASON_C33           : 86976 (%96.50 of all failed pairs)
```

##### Sample_04-STATS
```sh
jv~$ nano 01_QC/Sample_04-STATS.txt

number of pairs analyzed      : 2360078
total pairs passed            : 2243502 (%95.06 of all pairs)
  total pair_1 trimmed        : 0 (%0.00 of all passed pairs)
  total pair_2 trimmed        : 0 (%0.00 of all passed pairs)
total pairs failed            : 116576 (%4.94 of all pairs)
  pairs failed due to pair_1  : 3167 (%2.72 of all failed pairs)
  pairs failed due to pair_2  : 104984 (%90.06 of all failed pairs)
  pairs failed due to both    : 8425 (%7.23 of all failed pairs)
  FAILED_REASON_P             : 2754 (%2.36 of all failed pairs)
  FAILED_REASON_N             : 1 (%0.00 of all failed pairs)
  FAILED_REASON_C33           : 113821 (%97.64 of all failed pairs)
```

##### Sample_05-STATS
```sh
jv~$ nano 01_QC/Sample_05-STATS.txt

number of pairs analyzed      : 1351287
total pairs passed            : 1328989 (%98.35 of all pairs)
  total pair_1 trimmed        : 0 (%0.00 of all passed pairs)
  total pair_2 trimmed        : 0 (%0.00 of all passed pairs)
total pairs failed            : 22298 (%1.65 of all pairs)
  pairs failed due to pair_1  : 4929 (%22.11 of all failed pairs)
  pairs failed due to pair_2  : 15614 (%70.02 of all failed pairs)
  pairs failed due to both    : 1755 (%7.87 of all failed pairs)
  FAILED_REASON_C33           : 21933 (%98.36 of all failed pairs)
  FAILED_REASON_P             : 365 (%1.64 of all failed pairs)
```

##### Sample_06-STATS
```sh
jv~$ nano 01_QC/Sample_06-STATS.txt

number of pairs analyzed      : 2864984
total pairs passed            : 2769733 (%96.68 of all pairs)
  total pair_1 trimmed        : 0 (%0.00 of all passed pairs)
  total pair_2 trimmed        : 0 (%0.00 of all passed pairs)
total pairs failed            : 95251 (%3.32 of all pairs)
  pairs failed due to pair_1  : 3643 (%3.82 of all failed pairs)
  pairs failed due to pair_2  : 80829 (%84.86 of all failed pairs)
  pairs failed due to both    : 10779 (%11.32 of all failed pairs)
  FAILED_REASON_P             : 5277 (%5.54 of all failed pairs)
  FAILED_REASON_N             : 1 (%0.00 of all failed pairs)
  FAILED_REASON_C33           : 89973 (%94.46 of all failed pairs)
```

##### Sample_07-STATS
```sh
jv~$ nano 01_QC/Sample_07-STATS.txt

number of pairs analyzed      : 2652260
total pairs passed            : 2552257 (%96.23 of all pairs)
  total pair_1 trimmed        : 0 (%0.00 of all passed pairs)
  total pair_2 trimmed        : 0 (%0.00 of all passed pairs)
total pairs failed            : 100003 (%3.77 of all pairs)
  pairs failed due to pair_1  : 3232 (%3.23 of all failed pairs)
  pairs failed due to pair_2  : 89397 (%89.39 of all failed pairs)
  pairs failed due to both    : 7374 (%7.37 of all failed pairs)
  FAILED_REASON_P             : 2249 (%2.25 of all failed pairs)
  FAILED_REASON_C33           : 97754 (%97.75 of all failed pairs)
```

##### Sample_08-STATS
```sh
jv~$ nano 01_QC/Sample_08-STATS.txt

number of pairs analyzed      : 2101967
total pairs passed            : 1978716 (%94.14 of all pairs)
  total pair_1 trimmed        : 0 (%0.00 of all passed pairs)
  total pair_2 trimmed        : 0 (%0.00 of all passed pairs)
total pairs failed            : 123251 (%5.86 of all pairs)
  pairs failed due to pair_1  : 2285 (%1.85 of all failed pairs)
  pairs failed due to pair_2  : 113716 (%92.26 of all failed pairs)
  pairs failed due to both    : 7250 (%5.88 of all failed pairs)
  FAILED_REASON_N             : 2 (%0.00 of all failed pairs)
  FAILED_REASON_P             : 2003 (%1.63 of all failed pairs)
  FAILED_REASON_C33           : 121246 (%98.37 of all failed pairs)
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
|1|SRR1505166|1827066                    |1801218 (%98.59)     |
|2|SRR1505166|1712886                    |1628445 (%95.07)     |
|3|SRR1505166|1685155                    |1595021 (%94.65)     |
|4|SRR1505166|2360078                    |2243502 (%95.06)     |
|5|SRR1505166|1351287                    |1328989 (%98.35)     |
|6|SRR1505166|2864984                    |2769733 (%96.68)     |
|7|SRR1505166|2652260                    |2552257 (%96.23)     |
|8|SRR1505166|2101967                    |1978716 (%94.14)     |


## [Stage 2] Co-assembly

Generating `contigs.fasta`

> This protocol section is partially based on: http://merenlab.org/tutorials/assembly-based-metagenomics/

Note that no conda environment is activated at this stage, and the `python` version in where these scripts are working is in `Python 2.7.14`.

List your filtered reads
```sh
jv~$ ls 01_QC/*fastq

01_QC/Sample_01-QUALITY_PASSED_R1.fastq
01_QC/Sample_01-QUALITY_PASSED_R2.fastq
01_QC/Sample_02-QUALITY_PASSED_R1.fastq
01_QC/Sample_02-QUALITY_PASSED_R2.fastq
01_QC/Sample_03-QUALITY_PASSED_R1.fastq
01_QC/Sample_03-QUALITY_PASSED_R2.fastq
01_QC/Sample_04-QUALITY_PASSED_R1.fastq
01_QC/Sample_04-QUALITY_PASSED_R2.fastq
01_QC/Sample_05-QUALITY_PASSED_R1.fastq
01_QC/Sample_05-QUALITY_PASSED_R2.fastq
01_QC/Sample_06-QUALITY_PASSED_R1.fastq
01_QC/Sample_06-QUALITY_PASSED_R2.fastq
01_QC/Sample_07-QUALITY_PASSED_R1.fastq
01_QC/Sample_07-QUALITY_PASSED_R2.fastq
01_QC/Sample_08-QUALITY_PASSED_R1.fastq
01_QC/Sample_08-QUALITY_PASSED_R2.fastq
```

Create two environment variables:
```sh
jv~$ R1s=`ls 01_QC/*QUALITY_PASSED_R1* | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`
jv~$ R2s=`ls 01_QC/*QUALITY_PASSED_R2* | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`
```

The environment variables should look like this:

```sh
jv~$ echo $R1s

01_QC/Sample_01-QUALITY_PASSED_R1.fastq,01_QC/Sample_02-QUALITY_PASSED_R1.fastq,01_QC/Sample_03-QUALITY_PASSED_R1.fastq,01_QC/Sample_04-QUALITY_PASSED_R1.fastq,01_QC/Sample_05-QUALITY_PASSED_R1.fastq,01_QC/Sample_06-QUALITY_PASSED_R1.fastq,01_QC/Sample_07-QUALITY_PASSED_R1.fastq,01_QC/Sample_08-QUALITY_PASSED_R1.fastq

jv~$ echo $R2s

01_QC/Sample_01-QUALITY_PASSED_R2.fastq,01_QC/Sample_02-QUALITY_PASSED_R2.fastq,01_QC/Sample_03-QUALITY_PASSED_R2.fastq,01_QC/Sample_04-QUALITY_PASSED_R2.fastq,01_QC/Sample_05-QUALITY_PASSED_R2.fastq,01_QC/Sample_06-QUALITY_PASSED_R2.fastq,01_QC/Sample_07-QUALITY_PASSED_R2.fastq,01_QC/Sample_08-QUALITY_PASSED_R2.fastq
```

### Run MEGAHIT

The code has the form:
```sh
jv~$ megahit -1 $R1s -2 $R2s --min-contig-len $MIN_CONTIG_SIZE -m 0.85 -o 02_ASSEMBLY/ -t $NUM_THREADS
```

We are goint to set the parameters `MIN_CONTIG_SIZE = 1000` and `NUM_THREADS = 32`.

```sh
jv~$ megahit -1 $R1s -2 $R2s --min-contig-len 1000 -m 0.85 -o 02_ASSEMBLY/ -t 32

251.0Gb memory in total.
Using: 213.893Gb.
MEGAHIT v1.1.2
.
.
.
--- [STAT] 11599 contigs, total 27001230 bp, min 1000 bp, max 65169 bp, avg 2328 bp, N50 2482 bp
--- [Tue Nov 21 11:37:39 2017] ALL DONE. Time elapsed: 1784.559302 seconds ---
```

_**Time elapsed: 1784.559302 seconds = ~30 minutes**_

### Refining contigs with `anvio`

```sh
jv~$ mkdir 03_CONTIGS
jv~$ source activate anvio3
jv~$ anvi-script-reformat-fasta 02_ASSEMBLY/final.contigs.fa -o 03_CONTIGS/contigs.fa --min-len 1000 --simplify-names --report name_conversions.txt

Input ........................................: 02_ASSEMBLY/final.contigs.fa
Output .......................................: 03_CONTIGS/contigs.fa
Minimum length ...............................: 1,000
Total num contigs ............................: 11,599
Total num nucleotides ........................: 27,001,230
Contigs removed ..............................: 0 (0.00% of all)
Nucleotides removed ..........................: 0 (0.00% of all)
Deflines simplified ..........................: True
```

Testing for `--min-len 1500`

```sh
jv~$ anvi-script-reformat-fasta 02_ASSEMBLY/final.contigs.fa -o 03_CONTIGS/contigs_1500.fa --min-len 1500 --simplify-names --report name_conversions.txt
Input ........................................: 02_ASSEMBLY/final.contigs.fa
Output .......................................: 03_CONTIGS/contigs_1500.fa
Minimum length ...............................: 1,500
Total num contigs ............................: 11,599
Total num nucleotides ........................: 27,001,230
Contigs removed ..............................: 5796 (49.97% of all)
Nucleotides removed ..........................: 6991926 (25.89% of all)
Deflines simplified ..........................: True
```

Testing for `--min-len 2500`
```sh
jv~$ anvi-script-reformat-fasta 02_ASSEMBLY/final.contigs.fa -o 03_CONTIGS/contigs_2500.fa --min-len 2500 --simplify-names --report name_conversions.txt

Input ........................................: 02_ASSEMBLY/final.contigs.fa
Output .......................................: 03_CONTIGS/contigs_2500.fa
Minimum length ...............................: 2,500
Total num contigs ............................: 11,599
Total num nucleotides ........................: 27,001,230
Contigs removed ..............................: 9279 (80.00% of all)
Nucleotides removed ..........................: 13555420 (50.20% of all)
Deflines simplified ..........................: True
```

We will proceed with `03_CONTIGS/contigs.fa` which has **`--min-len` = 1000**

##### Deactivate the `anvio` conda environment

```sh
jv~$ source deactivate
```


## [Stage 4] Binning

### Binning with `MaxBin v2.2.4`

```sh
jv~$ mkdir 05_BINNING

```

Test `Maxbin2` installation

```sh
jv~$ run_MaxBin.pl -v
MaxBin 2.2.4
```

Create `reads_list.txt` to be passed to `MaxBin2`:

```sh
jv~$ ls 01_QC/*.fastq > 05_BINNING/reads_list.txt
```

Corroborate the `reads_list.txt` content:

```sh
jv~$ nano 05_BINNING/reads_list.txt

01_QC/Sample_01-QUALITY_PASSED_R1.fastq
01_QC/Sample_01-QUALITY_PASSED_R2.fastq
01_QC/Sample_02-QUALITY_PASSED_R1.fastq
01_QC/Sample_02-QUALITY_PASSED_R2.fastq
01_QC/Sample_03-QUALITY_PASSED_R1.fastq
01_QC/Sample_03-QUALITY_PASSED_R2.fastq
01_QC/Sample_04-QUALITY_PASSED_R1.fastq
01_QC/Sample_04-QUALITY_PASSED_R2.fastq
01_QC/Sample_05-QUALITY_PASSED_R1.fastq
01_QC/Sample_05-QUALITY_PASSED_R2.fastq
01_QC/Sample_06-QUALITY_PASSED_R1.fastq
01_QC/Sample_06-QUALITY_PASSED_R2.fastq
01_QC/Sample_07-QUALITY_PASSED_R1.fastq
01_QC/Sample_07-QUALITY_PASSED_R2.fastq
01_QC/Sample_08-QUALITY_PASSED_R1.fastq
01_QC/Sample_08-QUALITY_PASSED_R2.fastq
```


Run MaxBin binning on the assembled contigs:

> Note that here we are indicating `-thread 32`

```sh
jv~$ run_MaxBin.pl -contig 03_CONTIGS/contigs.fa -reads_list 05_BINNING/reads_list.txt -out 05_BINNING/LAKE_WA -thread 32

========== Job finished ==========
Yielded 9 bins for contig (scaffold) file 03_CONTIGS/contigs.fa
...
Summary file: 05_BINNING/LAKE_WA.summary
Genome abundance info file: 05_BINNING/LAKE_WA.abundance
Marker counts: 05_BINNING/LAKE_WA.marker
Marker genes for each bin: 05_BINNING/LAKE_WA.marker_of_each_gene.tar.gz
Bin files: 05_BINNING/LAKE_WA.001.fasta - 05_BINNING/LAKE_WA.009.fasta
Unbinned sequences: 05_BINNING/LAKE_WA.noclass
...
Store abundance information of reads file [01_QC/Sample_01-QUALITY_PASSED_R1.fastq] in [05_BINNING/LAKE_WA.abund1].
Store abundance information of reads file [01_QC/Sample_01-QUALITY_PASSED_R2.fastq] in [05_BINNING/LAKE_WA.abund2].
...
========== Elapsed Time ==========
0 hours 23 minutes and 46 seconds.
```
**Elapsed time: ~25 minutes**

The `manual mode` of the code above can be as follows:
**> DO NOT do this if the previous code worked well**
```sh
jv~$ run_MaxBin.pl -contig 03_CONTIGS/contigs.fa -reads 01_QC/Sample_01-QUALITY_PASSED_R1.fastq -reads2 01_QC/Sample_01-QUALITY_PASSED_R2.fastq -reads3 01_QC/Sample_02-QUALITY_PASSED_R1.fastq -reads4 01_QC/Sample_02-QUALITY_PASSED_R2.fastq -reads5 01_QC/Sample_03-QUALITY_PASSED_R1.fastq -reads6 01_QC/Sample_03-QUALITY_PASSED_R2.fastq -reads7 01_QC/Sample_04-QUALITY_PASSED_R1.fastq -reads8 01_QC/Sample_04-QUALITY_PASSED_R2.fastq -reads9 01_QC/Sample_05-QUALITY_PASSED_R1.fastq -reads10 01_QC/Sample_05-QUALITY_PASSED_R2.fastq -reads11 01_QC/Sample_06-QUALITY_PASSED_R1.fastq -reads12 01_QC/Sample_06-QUALITY_PASSED_R2.fastq -reads13 01_QC/Sample_07-QUALITY_PASSED_R1.fastq -reads14 01_QC/Sample_07-QUALITY_PASSED_R2.fastq -reads15 01_QC/Sample_08-QUALITY_PASSED_R1.fastq -reads16 01_QC/Sample_08-QUALITY_PASSED_R2.fastq -out 05_BINNING/LAKE_WA -thread 32
```


## [Stage 4] Refine bins

CheckM

```sh
jv~$ source activate CheckM
jv~$ checkm lineage_wf -t 32 -x fasta 05_BINNING/ 06_CHECKM

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Bin Id                 Marker lineage            # genomes   # markers   # marker sets    0     1    2    3   4   5+   Completeness   Contamination   Strain heterogeneity
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  LAKE_WA.006   c__Betaproteobacteria (UID3888)       323         387           234         4    294   81   7   1   0       99.09           30.05               8.33
  LAKE_WA.008   c__Betaproteobacteria (UID3888)       323         387           234         14   353   20   0   0   0       95.26            6.39              45.00
  LAKE_WA.001   c__Betaproteobacteria (UID3888)       323         387           234         89   284   14   0   0   0       83.34            5.77              57.14
  LAKE_WA.002         k__Bacteria (UID203)            5449        104            58         15    71   14   4   0   0       79.00            8.86              92.31
  LAKE_WA.004   c__Gammaproteobacteria (UID4274)      112         580           289        145   406   28   1   0   0       70.79            4.61              70.97
  LAKE_WA.007     o__Burkholderiales (UID4105)         54         553           264        228   245   73   6   1   0       56.07           13.32              22.68
  LAKE_WA.009       g__Pseudomonas (UID4576)           50         1090          364        444   590   56   0   0   0       52.97            3.06              32.14
  LAKE_WA.003   c__Betaproteobacteria (UID3888)       323         387           234        225   127   33   2   0   0       37.93            8.28              35.90
  LAKE_WA.005         k__Bacteria (UID203)            5449        103            57         98    5    0    0   0   0        6.30            0.00               0.00
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source deactivate
```

**> Time: ~10 minutes**

## [Stage 5] Functional annotation - `prokka`

```sh
jv~$ mkdir 07_ANNOTATION
jv~$ prokka --outdir 07_ANNOTATION/LAKE_WA.001 --prefix LAKE_WA.001 --kingdom Bacteria --metagenome --cpus 32 05_BINNING/LAKE_WA.001.fasta
jv~$ prokka --outdir 07_ANNOTATION/LAKE_WA.002 --prefix LAKE_WA.002 --kingdom Bacteria --metagenome --cpus 32 05_BINNING/LAKE_WA.002.fasta
jv~$ prokka --outdir 07_ANNOTATION/LAKE_WA.003 --prefix LAKE_WA.003 --kingdom Bacteria --metagenome --cpus 32 05_BINNING/LAKE_WA.003.fasta
jv~$ prokka --outdir 07_ANNOTATION/LAKE_WA.004 --prefix LAKE_WA.004 --kingdom Bacteria --metagenome --cpus 32 05_BINNING/LAKE_WA.004.fasta
jv~$ prokka --outdir 07_ANNOTATION/LAKE_WA.005 --prefix LAKE_WA.005 --kingdom Bacteria --metagenome --cpus 32 05_BINNING/LAKE_WA.005.fasta
jv~$ prokka --outdir 07_ANNOTATION/LAKE_WA.006 --prefix LAKE_WA.006 --kingdom Bacteria --metagenome --cpus 32 05_BINNING/LAKE_WA.006.fasta
jv~$ prokka --outdir 07_ANNOTATION/LAKE_WA.007 --prefix LAKE_WA.007 --kingdom Bacteria --metagenome --cpus 32 05_BINNING/LAKE_WA.007.fasta
jv~$ prokka --outdir 07_ANNOTATION/LAKE_WA.008 --prefix LAKE_WA.008 --kingdom Bacteria --metagenome --cpus 32 05_BINNING/LAKE_WA.008.fasta
jv~$ prokka --outdir 07_ANNOTATION/LAKE_WA.009 --prefix LAKE_WA.009 --kingdom Bacteria --metagenome --cpus 32 05_BINNING/LAKE_WA.009.fasta
```

| Bin         | Walltime used (minustes) |
|-------------|--------------------------|
| LAKE_WA.001 | 1.07                     |
| LAKE_WA.002 | 0.60                     |
| LAKE_WA.003 | 0.63                     |
| LAKE_WA.004 | 0.90                     |
| LAKE_WA.005 | 0.42                     |
| LAKE_WA.006 | 1.50                     |
| LAKE_WA.007 | 1.70                     |
| LAKE_WA.008 | 1.32                     |
| LAKE_WA.009 | 1.03                     |


## [Stage 6] Taxonomy classification of bins

**> Input for `` should be bins in `.faa` format**

Create a directory to save the results
```sh
jv~$ mkdir 08_TAXONOMY
```

Copy the bins in amino acids `.faa` format from the `prokka` output:
```sh
jv~$ cp 07_ANNOTATION/LAKE_WA.001/LAKE_WA.001.faa 08_TAXONOMY
jv~$ cp 07_ANNOTATION/LAKE_WA.001/LAKE_WA.001.faa 08_TAXONOMY
jv~$ cp 07_ANNOTATION/LAKE_WA.002/LAKE_WA.002.faa 08_TAXONOMY
jv~$ cp 07_ANNOTATION/LAKE_WA.003/LAKE_WA.003.faa 08_TAXONOMY
jv~$ cp 07_ANNOTATION/LAKE_WA.004/LAKE_WA.004.faa 08_TAXONOMY
jv~$ cp 07_ANNOTATION/LAKE_WA.005/LAKE_WA.005.faa 08_TAXONOMY
jv~$ cp 07_ANNOTATION/LAKE_WA.006/LAKE_WA.006.faa 08_TAXONOMY
jv~$ cp 07_ANNOTATION/LAKE_WA.007/LAKE_WA.007.faa 08_TAXONOMY
jv~$ cp 07_ANNOTATION/LAKE_WA.008/LAKE_WA.008.faa 08_TAXONOMY
jv~$ cp 07_ANNOTATION/LAKE_WA.009/LAKE_WA.009.faa 08_TAXONOMY

jv~$ ls 08_TAXONOMY/

LAKE_WA.001.faa  LAKE_WA.002.faa  LAKE_WA.003.faa  LAKE_WA.004.faa  LAKE_WA.005.faa  LAKE_WA.006.faa  LAKE_WA.007.faa  LAKE_WA.008.faa  LAKE_WA.009.faa
```

Go to the `phylophlan` directory (`$DIR_PHYLOPHLAN`) and create a new folder inside the `input` folder containing the bins in `.faa` format:

```sh
jv~$ cd $DIR_PHYLOPHLAN
jv~$ mkdir input/08_TAXONOMY
jv~$ cp $DIR_PROJECT/08_TAXONOMY/* input/08_TAXONOMY
jv~$ ls input/08_TAXONOMY/

LAKE_WA.001.faa  LAKE_WA.002.faa  LAKE_WA.003.faa  LAKE_WA.004.faa  LAKE_WA.005.faa  LAKE_WA.006.faa  LAKE_WA.007.faa  LAKE_WA.008.faa  LAKE_WA.009.faa
```

```sh
jv~$ source activate phylophlan
jv~$ ./phylophlan.py -i -t 08_TAXONOMY --nproc 32
jv~$ source deactivate
jv~$ ls -lh input/08_TAXONOMY
```

Results in `imputed_conf_high-conf.txt`
```sh
jv~$ cd output/08_TAXONOMY/
jv~$ nano imputed_conf_high-conf.txt

LAKE_WA.006     d__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Methylophilales.f__Methylophilaceae.g__Methylotenera.s__?.t__?
LAKE_WA.008     d__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Methylophilales.f__Methylophilaceae.g__Methylotenera.s__?.t__?
LAKE_WA.001     d__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Methylophilales.f__Methylophilaceae.g__Methylotenera.s__?.t__?
LAKE_WA.003     d__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Methylophilales.f__Methylophilaceae.g__Methylotenera.s__?.t__?
LAKE_WA.002     d__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Methylophilales.f__Methylophilaceae.g__Methylotenera.s__?.t__?
LAKE_WA.005     d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Methylococcales.f__Methylococcaceae.g__?.s__?.t__?
LAKE_WA.004     d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Methylococcales.f__Methylococcaceae.g__?.s__?.t__?
LAKE_WA.009     d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pseudomonadales.f__Pseudomonadaceae.g__Pseudomonas.s__fluorescens.t__?

```

Results in `imputed_conf_medium-conf.txt`
```sh
jv~$ nano imputed_conf_medium-conf.txt

LAKE_WA.007     d__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__?.f__?.g__?.s__?.t__?
```

**> Time:**
> Aligning finished in 1538 secs (1561 total time)
> Tree building finished in 3096 secs (4657 total time).
**> Total: ~ 1h 20min**

Copy the `phylophlan` results from `$DIR_PHYLOPHLAN` to the `$DIR_PROJECT`

```sh 
jv~$ cp output/08_TAXONOMY/* $PROJECT_DIR/08_TAXONOMY/
```


## [Stage 7] Refine functional annotation - `EggNOG-mapper`
`EggNOG-mapper v4.5.1` works on `python2.7`

`$DIR_EGGNOG` is the directory to `eggnog-mapper/emapper.py` script.

```sh
jv~$ mkdir 09_ANNOTATION_REFINE

jv~$ cd 09_ANNOTATION_REFINE

jv~$ python2.7 $DIR_EGGNOG/emapper.py -i ../07_ANNOTATION/LAKE_WA.001/LAKE_WA.001.faa --output LAKE_WA.001  -m diamond --cpu 32
jv~$ python2.7 $DIR_EGGNOG/emapper.py -i ../07_ANNOTATION/LAKE_WA.002/LAKE_WA.002.faa --output LAKE_WA.002  -m diamond --cpu 32
jv~$ python2.7 $DIR_EGGNOG/emapper.py -i ../07_ANNOTATION/LAKE_WA.003/LAKE_WA.003.faa --output LAKE_WA.003  -m diamond --cpu 32
jv~$ python2.7 $DIR_EGGNOG/emapper.py -i ../07_ANNOTATION/LAKE_WA.004/LAKE_WA.004.faa --output LAKE_WA.004  -m diamond --cpu 32
jv~$ python2.7 $DIR_EGGNOG/emapper.py -i ../07_ANNOTATION/LAKE_WA.005/LAKE_WA.005.faa --output LAKE_WA.005  -m diamond --cpu 32
jv~$ python2.7 $DIR_EGGNOG/emapper.py -i ../07_ANNOTATION/LAKE_WA.006/LAKE_WA.006.faa --output LAKE_WA.006  -m diamond --cpu 32
jv~$ python2.7 $DIR_EGGNOG/emapper.py -i ../07_ANNOTATION/LAKE_WA.007/LAKE_WA.007.faa --output LAKE_WA.007  -m diamond --cpu 32
jv~$ python2.7 $DIR_EGGNOG/emapper.py -i ../07_ANNOTATION/LAKE_WA.008/LAKE_WA.008.faa --output LAKE_WA.008  -m diamond --cpu 32
jv~$ python2.7 $DIR_EGGNOG/emapper.py -i ../07_ANNOTATION/LAKE_WA.009/LAKE_WA.009.faa --output LAKE_WA.009  -m diamond --cpu 32
```