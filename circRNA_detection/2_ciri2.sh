#! /bin/bash
#SBATCH -A ghwtcmr
#SBATCH -p bigmem
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH --job-name=CIRI2

echo "STARTING RUN: $(date)"
# load software
module load Perl/5.34.1-GCCcore-11.3.0
echo "loaded Perl/5.34.1-GCCcore-11.3.0"
bwa_path=/mnt/storage/nobackup/proj/ghwtcmr/ben/software/bwa-0.7.17
ciri_path=/mnt/storage/nobackup/proj/ghwtcmr/ben/software/CIRI_v2.0.6
# get arguments from job submission
base_dir=$(pwd)
fastq_dir=$1
reference_fasta=$2
gtf_anno=$3
# create samplesheet
echo "Getting sample ID"
cd ${fastq_dir}
samplesheet=samples.txt
sample=`head -n ${SLURM_ARRAY_TASK_ID} $samplesheet | tail -1`
echo "Starting run for ${sample} - $(date)"
# set up folders for ciri
cd ${base_dir}
mkdir -p CIRI2/
cd CIRI2/
mkdir -p alignment/
mkdir -p output/
# align with BWA-MEM
echo "Aligning ${sample} with BWA-MEM v0.7.17 - $(date)"
${bwa_path}/bwa mem \
-T 19 \
-t 8 \
${reference_fasta} \
${fastq_dir}/${sample}_R1.fastq.gz \
${fastq_dir}/${sample}_R2.fastq.gz \
> alignment/${sample}.sam
# Detect BSJs with CIRI2
echo "Detecting BSJs in ${sample} with CIRI v2.0.6 - $(date)"
perl ${ciri_path}/CIRI2.pl \
-I alignment/${sample}.sam \
-O output/${sample}.bed \
-F ${reference_fasta} \
-A ${gtf_anno} \
-G ${sample}_ciri.log \
-T 8 \
--max_span 1000000 \
--rel_exp 0.05
echo "Finished run for ${sample} - $(date)"