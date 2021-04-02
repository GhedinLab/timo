#!/bin/sh

#### USER CHANGES THESE VARIABLES DEPENDING ON RUN ######
### run: bash ./run_IndexVariants.sh
JOB_NAME="variants-run"
RUNDIR="/data/user/timo"
BAM_DIR="/data/user/bamdir" #using bam because of index
STRAIN="cov19" #lowercase strain name
REFSEQ="/data/user/SARS-CoV2.fa"
USER='some.email@email.com' #email for running information
array=0-9 #number of samples! (not number of bams...)


#### DON'T TOUCH BELOW HERE #####
cd ${RUNDIR}
sbatch --mail-type=END --mail-user=$USER --job-name=$JOB_NAME -a ${array} index_variants.sh ${BAM_DIR} ${RUNDIR} ${STRAIN} ${REFSEQ}

