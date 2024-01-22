#! /bin/bash
#SBATCH -A ghwtcmr
#SBATCH -p bigmem
#SBATCH --mem=32G
#SBATCH -c 8
#SBATCH --job-name=trimGalore

echo "STARTING ADAPTER AND QUALITY TRIMMING: $(date)"
# load software
echo "LOADING CONDA ENVIRONMENT to run Trim Galore! and MultiQC (qc) - $(date)"
module load Anaconda3/2021.11
eval "$(conda shell.bash hook)"
conda activate qc
# get arguments from job submission
base_dir=$(pwd)
fastq_dir=$1
# create samplesheet
echo "Getting sample ID"
cd ${fastq_dir}
samplesheet=samples.txt
sample=`head -n ${SLURM_ARRAY_TASK_ID} $samplesheet | tail -1`
echo "Starting job for ${sample} - $(date)"
# folder set up
mkdir -p ${fastq_dir}/repaired_fastq/singletons
mkdir -p ${fastq_dir}/trimming/
mkdir -p ${fastq_dir}/trimming/trimmed_fastq/
mkdir -p ${fastq_dir}/trimming/qc_preTrim/
mkdir -p ${fastq_dir}/trimming/qc_postTrim/
# repair fastq files
echo "Reparing fastq files - $(date)"
repair.sh \
    in1=${fastq_dir}/${sample}_R1.fastq.gz in2=${fastq_dir}/${sample}_R2.fastq.gz \
    out1=${fastq_dir}/repaired_fastq/${sample}_R1.fastq.gz out2=${fastq_dir}/repaired_fastq/${sample}_R2.fastq.gz \
    outs=${fastq_dir}/repaired_fastq/singletons/${sample}_singletons.fastq.gz \
    repair
# run fastqc on pre-trimmed fastq files
echo "Running fastqc on repaired fastq files - $(date)"
fastqc ${fastq_dir}/repaired_fastq/${sample}* -t 8 --outdir ${fastq_dir}/trimming/qc_preTrim/
# run Trim Galore!
echo "Starting Trim Galore! on ${sample} - $(date)"
trim_galore --cores 8 --paired --quality 15 \
    --fastqc --fastqc_args "--outdir ${fastq_dir}/trimming/qc_postTrim" \
    ${fastq_dir}/repaired_fastq/${sample}_R1.fastq.gz \
    ${fastq_dir}/repaired_fastq/${sample}_R2.fastq.gz \
    -o ${fastq_dir}/trimming/trimmed_fastq
# rename output files
cd ${fastq_dir}/trimming/trimmed_fastq
mv ${sample}_R1_val_1.fq.gz ${sample}_R1.fastq.gz
mv ${sample}_R2_val_2.fq.gz ${sample}_R2.fastq.gz
echo "Job finished - $(date)"