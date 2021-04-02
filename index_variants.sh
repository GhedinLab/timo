#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=3:00:00
#SBATCH --mem=32GB

### DON'T CHANGE ANY CODE BELOW #####

. /data/SGSlab/conda/etc/profile.d/conda.sh
conda activate

### input variables
BAM_DIR=$1
RUNDIR=$2
STRAIN=$3
REFSEQ=$4

F="${BAM_DIR}/*.bam"


### array set up
cd $RUNDIR
bams=($(ls $F))

FS=${bams[$SLURM_ARRAY_TASK_ID]} #read name associated with slurm task id

N=${FS#"${BAM_DIR}/"} #name of sample
NAME=${N%".bam"} #may need to change if not zipped
echo $NAME

module purge
module load samtools/1.11
#module load pysam/intel/0.10.0
module load picard/2.17.11
module load python/3.8

conda 'bioconda::pysam bioconda::pairix'

samtools sort -T ./${NAME}.${STRAIN}.sorted -o ./${NAME}.${STRAIN}.sorted.merged.bam ${BAM_DIR}/${NAME}.bam
samtools index ./${NAME}.${STRAIN}.sorted.merged.bam
java -jar /usr/local/apps/picard/2.17.11/picard.jar MarkDuplicates I=./${NAME}.${STRAIN}.sorted.merged.bam O=./${NAME}.${STRAIN}.rmd.merged.bam M=./${NAME}.${STRAIN}.met.txt REMOVE_DUPLICATES=true
samtools index ./${NAME}.${STRAIN}.rmd.merged.bam
python timo.v1.py --strain ${STRAIN} --infile ./${NAME}.${STRAIN}.rmd.merged.bam --ref ${REFSEQ} --covercutoff 1 --cutoff 0.001
