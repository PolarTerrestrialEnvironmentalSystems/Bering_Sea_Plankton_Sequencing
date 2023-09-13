#!/bin/bash

#SBATCH --job-name=hopsKL77ds01
#SBATCH -p xfat
#SBATCH --qos=large
#SBATCH -t 96:00:00 
#SBATCH --cpus-per-task=28
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --mem=1400G 
#SBATCH --error="hopsNT-%j.err"
#SBATCH --out="hopsNT-%j.out"



# set variables (requires modification)
#================================================
WORKDIR=/path/to/working/directory/
INPUT=/path/to/working/directory/ds01/ds01.fq.gz
OUTPUT=/path/to/working/directory/output_nt_ds01/
CONFIG=/path/to/working/directory/configfile_nt.txt
INDEXDB=/path/to/database/nt-hops-20-11-03_step8/
TAXFILE=path/to/working/directory/taxalist_2022_09_20.txt
NCBIRESC=/path/to/database/ncbi
MEM=1400

# preparing the working environment
#================================================
cd $WORKDIR
module load bio/hops/0.34

# create Hops config file
#================================================
echo "preProcess=0 
alignment=1 
dupRemOff=0 
filter=def_anc 
index=${INDEXDB}
pathToList=${TAXFILE}
resources=${NCBIRESC}

useSlurm = 0
threadsMalt=${SLURM_CPUS_PER_TASK}
maxMemoryMalt=${MEM}

threadsMaltEx=${SLURM_CPUS_PER_TASK}
maxMemoryMaltEx=${MEM}

threadsPost=${SLURM_CPUS_PER_TASK}
maxMemoryPost=256" > ${CONFIG}

# tasks to be performed
#================================================

srun hops -Xmx${MEM}G -input ${INPUT} -output ${OUTPUT} -m full -c ${CONFIG} 


