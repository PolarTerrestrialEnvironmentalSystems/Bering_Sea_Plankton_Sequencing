# Code for Diatom Amplicon sequencing

This directory contains the code for filtering, processing and analyzing diatom amplicon sequencing data from a marine sediment core from the Bering Sea for the manuscript:
Plankton community change during the last 124 000 years in the subarctic Bering Sea derived from sedimentary ancient DNA (add link after publication) <br>
It contains all necessary files and scripts for reproducing the results presented in the manuscript.

## Raw sequencing data
The raw sequencing data are stored in the European Nucleotide Archive in bioproject XXX (add when raw data are uploaded).

## Taxonomic classification with ObiTools
We used the Python package [OBITools 3.0.1](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12428) to analyze the raw paired-end sequencing data of the amplicon-sequencing of the rbcL gene of diatoms. It merges paired-end reads, demultiplexes samples, and performs the taxonomic classification based on sequence similarity to the customized rbcL-EMBL nucleotide reference database `rbcl_embl143_db.fasta`.
Before running the ObiTools script `Obitools_APMG_41_embl143_rbcL.sh`, the follwoing adjustments have to be made to the  `rbcl_embl143_db.fasta` file:

```
## define path for output.dms
DATAPATH=
 
obi import rbcl_embl143_db.fasta $DATAPATH/rbcl_embl143_ref
obi import --taxdump taxdump.tar.gz $DATAPATH/taxonomy/ncbi_tax
obi grep --require-rank=species --require-rank=genus --require-rank=family --taxonomy $DATAPATH/taxonomy/ncbi_tax $DATAPATH/ rbcl_embl143_ref $DATAPATH/ rbcl_embl143_refDB
obi clean_dms $DATAPATH
obi build_ref_db -t 0.97 --taxonomy $DATAPATH/taxonomy/ncbi_tax $DATAPATH/rbcl_embl143_refDB $DATAPATH/embl143_rbcl
```

The following adjustments have to be made in the script `Obitools_APMG_41_embl143_rbcL.sh`:

```
# give the email address for notifications
#SBATCH --mail-user=
```

```
# change oath to the directory where the output has to be saved
DMS_OUT=/path/to/results/
```

```
# change directory to where the raw sequencing data and the tagfile are saved
DATAPATH=/path/to/raw/data
```

```
# change to path to where rbcL database is saved
cp -r /path/to/embl143_rbcl.obidms /tm
```


## Data analysis with R

### Input
As an input, the diatom amplicon sequencing data after taxonomic assignment through the [ObiTools Pipeline](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12428) are used.
The raw data and the ObiTools script are stored here (add link to data repository).

The folder `input_files` contains the following files: <br>

|file name|description|
|-|-|
|APMG-41_embl143_rbcL.xlsx|tabel of amplicon sequencing data after taxonomic assignment using the ObiTools pipeline|
|AVS_to_Genus_info.xlsx|table assigning every ASV to a diatom genus|
|Age_Period_Group.xlsx|table that assigns sample ages to a period (Holocene, Youger Dryas, Bollong-Allerod, Last Glacial Maximum, Glacial period, Eemian)|
|Envi_for_RDA.xlsx|Input for the RDA, including d18O NGRIP data from [Rasmussen et al. 2006](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2005JD006079) and zooplankton data from [shotgun analysis](https://github.com/StellaZBuchwald/Bering_Sea_Plankton_Sequencing/tree/main/Shotgun)|
|Not_Bacillariophyceae.txt|all taxa not representing Bacillariophyceae in the raw data for diatoms for inital filtering|
|Sample_names.xlsx|table assigning sample names after sequencing to sample age|
|Shape_info.xlsx|table assigning main ASVs to centric or pennate shape|

---

### R script

For running the R script `Bering_Sea2023_Metabarcoding.R` the working directory in the first row must be customitzed to a directory which contains a folder with the `Input_files`.
In this directory, an `Output` folder will be created if it doesn't exist alerady. The `output` is included in `.gitignore` to avoid it being pushed into the repository.

The R script is structured into several sections:
<li> Filtering ObiTools Output
<li> Counts per replicate
<li> NMDS
<li> Sum of replicates
<li> Rarefaction <br>
  &emsp;&emsp;&emsp;&emsp;Rarefaction on ASV level <br>
  &emsp;&emsp;&emsp;&emsp;Mean ASVs per sample <br>
   &emsp;&emsp;&emsp;&emsp;Mean counts oer ASV and sample <br>
<li> Stratigrams
<li> PCA/RDA

---
For each section, output files will be created that might be used as input files in following sections.
As long as all necessary input files are created already, each of the sections can be run independently.

For the `Rarefaction`, the R script by [Stefan Kruse](https://github.com/StefanKruse/R_Rarefaction) was adapted.

For running the script, the following libraries in R must be installed:
<li> tidyverse
<li> readxl
<li> openxlsx
<li> vegan
<li> ggplot2
<li> tidypaleo
<li> ggh4x
<li> olsrr
<li> ggrepel

---
