#!/bin/sh

#### USER CHANGES THESE VARIABLES DEPENDING ON RUN ######
### run: bash ./run_IndexVariants.sh
JOB_NAME="OMS524-583-NIHvariants"
RUNDIR="/data/SGSlab/pipelines/timo"
BAM_DIR="/data/SGSlab/NIH-COV2/20210329/out/merged" #using bam because of index
STRAIN="cov19" #lowercase strain name
REFSEQ="/data/SGSlab/SARS-CoV2.fa"
USER='allison.roder@nih.gov' #email for running information
array=0-43 #number of samples! (not number of bams...) originally 0-21


#### DON'T TOUCH BELOW HERE #####
cd ${RUNDIR}
sbatch --mail-type=END --mail-user=$USER --job-name=$JOB_NAME -a ${array} index_variants.sh ${BAM_DIR} ${RUNDIR} ${STRAIN} ${REFSEQ}

