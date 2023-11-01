# Code for Bering Sea shotgun sequencing

This directory contains the code for filtering, processing and analyzing shotgun sequencing data from a marine sediment core from the Bering Sea for the manuscript:
Buchwald, S. Z., Herzschuh, U., Nürnberg, D., Harms, L., & Stoof-Leichsenring, K. R. Plankton community changes during the last 124 000 years in the subarctic Bering Sea derived from sedimentary ancient DNA. <br>
It contains all necessary `input_files` and the scripts for replicating the results presented in the manuscript.

## Raw sequencing data
The raw sequencing data are stored in the European Nucleotide Archive (Bioproject number PRJE866300).

## Taxonomic classification with Kraken2
### data preparation
Two different sets of samples were sent for shotgun sequencing, resulting in two separate raw data files and two files after the taxonomic assignment with Kraken2. 
The Kraken2 script is stored [here](https://github.com/StellaZBuchwald/Bering_Sea_Plankton_Sequencing/tree/main/Shotgun/Kraken2).
Before taxonomic classification with [Kraken2](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown) [Wood et al. 2019](https://www.biorxiv.org/content/10.1101/762302v1), the raw sequencing data were quality checked using the software FastQC (version 0.11.9), and read duplicates were removed with the program BBmap (version 38.87). For adapter trimming and merging of overlapping reads Fastp (version 0.20.1) was used and the quality of the processed reads was checked again.
These processes are carried out if the script `1_step_shop_Bering_with_clumpify.sl` is run with the following adjustments:

```
# give the email address for notifications
#SBATCH --mail-user=
```

```
# change the directory to the path where the raw data are stored 
INDIR=path/to/raw/shotgun/data
```

### Kraken2 Script
On the quality checked and deduplicated dataset, the Kraken2 script `2_step_shop_Bering.sh` can be run. The following adjustments have to be made to the script:

```
# give the email address for notifications
#SBATCH --mail-user=
```

```
# change directory to the path where the database is stored
DB="/path/to/nt_2021_04_db"
```

The output directory will stay the same to the output directory for the script `1_step_shop_Bering_with_clumpify.sl`.



## Ancient pattern analysis with HOPS
For the analysis of ancient pattern in the DNA sequences, we used the HOPS pipeline ([Hübler et al. 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1903-0)). Scripts and input files can be found in the folder `HOPS`. For creating the database and further information, click [here](https://github.com/rhuebler/HOPS).
We analyzed ancient DNA patterns for different taxa, summarized in the `taxalist.txt`.


We analyszed ancient DNA patterns for six time intervals (ds01-ds06) and samples were grouped accordingly using the `cat` command:

|sample age (ka BP) |time interval||sample age (ka BP)|time interval||sample age (ka BP)|time interval|
|-|-|-|-|-|-|-|-|
|1.82|Holocene (ds01)||16.09|Deglaciation (ds02)||77.48|Glacial II (ds04)|
|3.22|Holocene (ds01)||17.54|Deglaciation (ds02)||86.82|Glacial II (ds04)|
|4.68|Holocene (ds01)||18.75|Deglaciation (ds02)||87.83|Glacial II (ds04)|
|7.21|Holocene (ds01)||20.50|Glacial I (ds03)||94.00|Glacial III (ds05)|
|7.97|Holocene (ds01)||24.55|Glacial I (ds03)||99.59|Glacial III (ds05)|
|9.19|Holocene (ds01)||26.65|Glacial I (ds03)||101.13|Glacial III (ds05)|
|10.15|Holocene (ds01)||32.24|Glacial I (ds03)||103.37|Glacial III (ds05)|
|11.20|Holocene (ds01)||38.11|Glacial I (ds03)||106.79|Glacial III (ds05)|
|11.74|Deglaciation (ds02)||40.03|Glacial I (ds03)||111.63|Glacial III (ds05)|
|12.43|Deglaciation (ds02)||47.09|Glacial I (ds03)||120.38|Eemian (ds06)|
|12.62|Deglaciation (ds02)||52.70|Glacial II (ds04)||121.25|Eemian (ds06)|
|13.09|Deglaciation (ds02)||56.34|Glacial II (ds04)||122.03|Eemian (ds06)|
|14.33|Deglaciation (ds02)||62.68|Glacial II (ds04)||123.34|Eemian (ds06)|
|14.99|Deglaciation (ds02)||68.11|Glacial II (ds04)||123.86|Eemian (ds06)|




Both the `taxalist.txt` and the `configfile.txt` have to be stored in the working directory that has to be indicated in the script `03_hops_nt_KL77_ds01.sh` and the `configfile.txt`.
For running the script, the following modifications have to be made to the `configfile.txt`:

```
# change directory to the path where the database is stored
index=/path/to/database/nt-hops-20-11-03_step8/
resources=/path/to/database/ncbi
```

```
# change directory to where the taxalist is stored
PathToList=/path/to/working/directory/taxalist_2022_10_20.txt
```

Also, the following modifications have to be made to the script `03_hops_nt_KL77_ds01.sh`:

```
# give the email address for notifications
#SBATCH --mail-user=
```

```
# set working directory in the following lines
WORKDIR=/path/to/working/directory/
INPUT=/path/to/working/directory/ds01/ds01.fq.gz
OUTPUT=/path/to/working/directory/output_nt_ds01/
CONFIG=/path/to/working/directory/configfile_nt.txt
TAXFILE=path/to/working/directory/taxalist_2022_09_20.txt
```

```
# set directory to the database in the following lines
INDEXDB=/path/to/database/nt-hops-20-11-03_step8/
NCBIRESC=/path/to/database/ncbi
```

The output file will be stored in the folder `output_nt_ds01` that will be created in the working directory.


## Data analysis with R
### Input
As an input, the shotgun sequencing data after taxonomic assignment through the [Kraken2](https://www.biorxiv.org/content/10.1101/762302v1) pipeline are used.

The folder `input_files` contains the following files: <br>

|file name|description|
|-|-|
|Age_Period_Group.xlsx|table assigning sample ages to a period (Holocene, Youger Dryas, Bolling-Allerod, Last Glacial Maximum, Glacial period, Eemian)|
|Group_infos_phyto.xlsx|table assigning the phytoplankton families to broader taxonomic groups (chlorophyes, phototrophic bacteria, phototrophic protsits, red algae)|
|Group_infos_zoo.xlsx|table assigning the zooplankton families to broader taxonomic groups (crustaceous zooplankton, heterotrophic protsits, gelatinous zooplankton)|
|KL-77-1_nt0.2_lineageDB_age_depth_adj.xlsx|Kraken2 Output of the first batch sent to sequencing|
|KL77-2_nt0.2_lineageDB_age_depth_adj.xlsx|Kraken2 Output of the second batch sent to sequencing|
|NGRIP_data.xlsx|d18O NGRIP data from [Rasmussen et al. 2006](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2005JD006079)|
|Phytoplankton_Family_Genus_Species.xlsx|table assigning each phytoplankton plankton species to a family|
|phyto_families.txt|phytoplankton families that the raw data are filtered for|
|zoo_families.txt|zooplankton families that the raw data are filtered for|
---

### R script

For running the R script `Bering_Sea2023_Shotgun.R` the working directory in the first row must be customized to a directory which contains a folder with the `Input_files`.
In this directory, an `Output` folder will be created if it doesn't exist alerady. The `output` is included in `.gitignore` to avoid it being pushed into the repository.

The R script is structured into several sections:
<li> Creating Subset of data <br>
&emsp;&emsp;&emsp;&emsp;Creating the dataset <br>
&emsp;&emsp;&emsp;&emsp; Group summary <br>
<li>  Rarefaction <br>
&emsp;&emsp;&emsp;&emsp; Rarefaction on family level <br>
&emsp;&emsp;&emsp;&emsp; Mean families per sample <br>
&emsp;&emsp;&emsp;&emsp; Mean counts per family and sample <br>
<li> Stratigrams
<li> Environmental parameters
<li> Zooplankton <br>
&emsp;&emsp;&emsp;&emsp; Creating the dataset <br>
&emsp;&emsp;&emsp;&emsp; Rarefaction <br>
&emsp;&emsp;&emsp;&emsp; Mean total counts per sample <br>
&emsp;&emsp;&emsp;&emsp; Mean families per sample <br>
&emsp;&emsp;&emsp;&emsp; Realtive abundance of main zooplankton families <br>
&emsp;&emsp;&emsp;&emsp; Zooplankton stratigram <br>
<li> PCA/RDA

---
For each section, output files will be created that might be used as input files in following sections.
As long as all necessary input files are created already, each of the sections can be run independently.
For the `Rarefaction`, the R script by [Stefan Kruse](https://github.com/StefanKruse/R_Rarefaction) was adapted.


For running the script, the following libraries must be installed in R:
<li> tidyverse
<li> scales
<li> readxl
<li> openxlsx
<li> vegan
<li> ggplot2
<li> tidypaleo
<li> dplyr
<li> ggh4x
<li> olsrr
<li> ggrepel

---
