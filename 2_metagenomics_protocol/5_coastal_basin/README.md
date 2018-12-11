# The northern head of the coastal basin of GD (Golfo Dulce), Costa Rica

- https://doi.org/10.3389/fmars.2017.00023
- Assembly: https://www.ncbi.nlm.nih.gov/assembly?LinkName=genome_assembly&from_uid=50906
- https://www.ncbi.nlm.nih.gov/assembly/GCA_001901525.1
- Taxonomy: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1920888
- cellular organisms; Bacteria; Proteobacteria; Gammaproteobacteria; Methylococcales; unclassified Methylococcales
- Methylococcales bacterium OPU3_GD_OMZ (g-proteobacteria)

>Padilla, C. C., Bertagnolli, A. D., Bristow, L. A., Sarode, N., Glass, J. B., Thamdrup, B., & Stewart, F. J. (2017). Metagenomic binning recovers a transcriptionally active Gammaproteobacterium linking methanotrophy to partial denitrification in an anoxic oxygen minimum zone. Frontiers in Marine Science, 4, 23.


## [Stage 0] Download the available binned genome

```s
mkdir 04_BINNING
# download the available genome here
```

## [Stage 4] Refine bins
```sh
jv~$ source activate CheckM
jv~$ checkm lineage_wf -t 32 -x fasta 04_BINNING/ 05_CHECKM/
source deactivate
```

```s
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Bin Id                             Marker lineage            # genomes   # markers   # marker sets   0     1    2    3   4   5+   Completeness   Contamination   Strain heterogeneity
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Golfo_Dulce_Bin010_OPU3   c__Gammaproteobacteria (UID4274)      112         581           290        24   522   34   1   0   0       95.27            5.46              21.62
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
```

## [Stage 5] Functional annotation - `prokka`

Annotating the bin:

```sh
jv~$ mkdir 06_ANNOTATION/
jv~$ prokka --outdir 06_ANNOTATION/Golfo_Dulce_Bin010_OPU3 --prefix Golfo_Dulce_Bin010_OPU3 --kingdom Bacteria --metagenome --cpus 32 04_BINNING/Golfo_Dulce_Bin010_OPU3.fasta
jv~$ prokka --outdir 06_ANNOTATION_genome//Golfo_Dulce_Bin010_OPU3 --prefix Golfo_Dulce_Bin010_OPU3 --kingdom Bacteria --cpus 32 04_BINNING/Golfo_Dulce_Bin010_OPU3.fasta
```

## [Stage 6] Taxonomy classification of the bin

Create a directory to save the results
```sh
jv~$ mkdir 07_TAXONOMY/
```

Copy the bins in amino acids `.faa` format from the `prokka` output:
```sh
jv~$ cp 06_ANNOTATION/Golfo_Dulce_Bin010_OPU3/Golfo_Dulce_Bin010_OPU3.faa 07_TAXONOMY/
jv~$ ls -1 07_TAXONOMY/
Golfo_Dulce_Bin010_OPU3.faa
```

Go to the `phylophlan` directory (`$DIR_PHYLOPHLAN`), then create a new folder inside the `input` folder containing the bins in `.faa` format:

```sh
jv~$ cd ~/programs/phylophlan/
jv~$ source activate phylophlan
jv~$ mkdir input/07_TAXONOMY_Golfo_Dulce_/
jv~$ cp ~/data/metagenomics-coastal_basin/07_TAXONOMY/Golfo_Dulce_Bin010_OPU3.faa input/07_TAXONOMY_Golfo_Dulce_/
jv~$ ls -1 input/07_TAXONOMY_OPHIOLITE_overlap/
Golfo_Dulce_Bin010_OPU3.faa
```

```sh
jv~$ ./phylophlan.py -i -t 07_TAXONOMY_Golfo_Dulce_ --nproc 32
Total time: 3293.99 seconds Unique: 3140/3172 Bad splits: 5/3137 Worst delta-LogLk 0.796
Tree built! The output newick file is in output/07_TAXONOMY_Golfo_Dulce_/07_TAXONOMY_Golfo_Dulce_.tree.int.nwk
Tree building finished in 3492 secs (4688 total time).
Trying to impute taxonomic labels for taxa newly integrated into the tree... Done!
Writing taxonomic imputation outputs ... Done!
jv~$ source deactivate
```


```sh
jv~$ cat output/07_TAXONOMY_Golfo_Dulce_/imputed_conf_high-conf.txt
```

```s
Golfo_Dulce_Bin010_OPU3	d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Methylococcales.f__Methylococcaceae.g__?.s__?.t__?
```

Copy the `phylophlan` results from `$DIR_PHYLOPHLAN` to the `$DIR_PROJECT`

```sh 
jv~$ cp output/07_TAXONOMY_Golfo_Dulce_/* ~/data/metagenomics-coastal_basin/07_TAXONOMY/
```

## [Stage 7]  Refine functional annotation of methanotrophic bacteria using `EggNOG-mapper`
`EggNOG-mapper v4.5.1` works on `python2.7`

`$DIR_EGGNOG` is the directory to `eggnog-mapper/emapper.py` script.

```sh
jv~$ cd ~/data/metagenomics-coastal_basin/
jv~$ mkdir 08_ANNOTATION
jv~$ cd 08_ANNOTATION
jv~$ python2.7 ~/programs/eggnog-mapper/emapper.py -i ../07_TAXONOMY/Golfo_Dulce_Bin010_OPU3.faa --output Golfo_Dulce_Bin010_OPU3 -m diamond --cpu 32
```