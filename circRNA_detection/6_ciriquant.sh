#! /bin/bash
#SBATCH -A ghwtcmr
#SBATCH -p bigmem
#SBATCH -c 12
#SBATCH --mem=64G
#SBATCH --job-name=CIRIquant

base_dir=$(pwd)
fastq_dir=$1
bsj_bed=$2
config=$3

# load software
echo "LOADING CONDA ENVIRONMENT to run CIRIquant (ciri) - $(date)"
module load Anaconda3/2021.11
eval "$(conda shell.bash hook)"
conda activate ciri

# get sample ID
echo "Getting sample ID"
cd ${fastq_dir}
samplesheet=samples.txt
sample=`head -n ${SLURM_ARRAY_TASK_ID} $samplesheet | tail -1`
echo "Starting run for ${sample} - $(date)"

# set up
cd ${base_dir}
mkdir -p ciriquant/${sample}

# run CIRIquant
CIRIquant \
    --config ${config} \
    -1 ${fastq_dir}/${sample}_R1.fastq.gz -2 ${fastq_dir}/${sample}_R2.fastq.gz \
    --bed ${bsj_bed} \
    -o ${base_dir}/ciriquant/${sample}/ \
    -t 12 \
    -p ${sample} \
    -l 2 \

echo "Run finished for ${sample} - $(date)"