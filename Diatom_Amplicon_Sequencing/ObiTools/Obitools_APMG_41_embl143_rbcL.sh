#!/bin/bash
#SBATCH -t48:00:00 ### set the time
#SBATCH --qos=large ## set the range your time demand is in (not mandatory)
#SBATCH -p smp ## set the node
#SBATCH -c 36
#SBATCH --job-name=APMG-41_rbcL
#SBATCH --output=%x_%j.out
#SBATCH --mail-user=
#SBATCH --mail-type=END

# load Obi3
module load bio/OBItools/3.0.0
# new modlue
# module load bio/OBItools/3.0.1b7
######################################################
################ things to adjust ####################
######################################################
#data variables (please adjust)
DMS_OUT=/path/to/results/
DMS=/tmp/APMG-41
# data path and file name
DATAPATH=/path/to/raw/data
FORWARD=210914_NB501850_A_L1-4_APMG-41_R1.fastq
REVERSE=210914_NB501850_A_L1-4_APMG-41_R2.fastq
TAGFILE=$DATAPATH/APMG-41_tagfile_rbcL.txt
ID=APMG-41_
#######################################################
#######################################################
# copy variables to tmp

cp $DATAPATH/$FORWARD /tmp
cp $DATAPATH/$REVERSE /tmp

cp -r /path/to/embl143_rbcl.obidms /tmp

# set tmp variables
DATA1=/tmp/$FORWARD
DATA2=/tmp/$REVERSE
EMBL=/tmp/embl143_rbcl.obidms

# import data
srun obi import --fastq-input $DATA1 $DMS/${ID}reads1
srun obi import --fastq-input $DATA2 $DMS/${ID}reads2
# import tag file
srun obi import --ngsfilter-input $TAGFILE $DMS/${ID}tagfile
# align paired ends
srun obi alignpairedend -R $DMS/${ID}reads1 $DMS/${ID}reads2 $DMS/${ID}aligned_reads
# remove unaligned
srun obi grep -a mode:alignment $DMS/${ID}aligned_reads $DMS/${ID}good_sequences
# assign reads to tag combinations
srun obi ngsfilter -t $DMS/${ID}tagfile -u $DMS/${ID}unidentified_seqs $DMS/${ID}good_sequences $DMS/${ID}identified_sequences
# dereplicate sequences into unique sequences
srun obi uniq --merge sample $DMS/${ID}identified_sequences $DMS/${ID}dereplicated_sequences
# denoise data, only keep COUNT and merged_sample tags
srun obi annotate -k COUNT -k MERGED_sample $DMS/${ID}dereplicated_sequences $DMS/${ID}cleaned_metadata_sequences
#remove sequences with length smaller 10, and count smaller 10
srun obi grep -p "len(sequence)>=10 and sequence['COUNT']>=10" $DMS/${ID}cleaned_metadata_sequences $DMS/${ID}cleaned_10_metadata_sequences
# denoise data, clean from pcr/sequencing errors
srun obi clean -s MERGED_sample -r 0.05 -H $DMS/${ID}cleaned_10_metadata_sequences $DMS/${ID}cleaned_sequences
# taxonomic assignment to embl
srun obi clean_dms $EMBL
srun obi ecotag --taxonomy $EMBL/TAXONOMY/ncbi_tax_12_20 -R $EMBL/VIEWS/rbcl_embl143_db.obiview $DMS/${ID}cleaned_sequences $DMS/${ID}embl143_assigned_sequences
#annonate lineage information
srun obi annotate --with-taxon-at-rank family --with-taxon-at-rank genus --with-taxon-at-rank species --taxonomy $EMBL/TAXONOMY/ncbi_tax_12_20 $DMS/${ID}embl143_assigned_sequences $DMS/${ID}embl143_annotated_assigned_sequences
# export embl assignment as csv file
srun obi export --header --tab-output $DMS/${ID}embl143_annotated_assigned_sequences > $DMS_OUT/${ID}embl143_rbcL.csv
srun cp -r ${DMS}.obidms $DMS_OUT
