#! /bin/bash
#SBATCH -A ghwtcmr
#SBATCH -p bigmem
#SBATCH --mem=64G
#SBATCH -c 16
#SBATCH --job-name=pfv2

echo "STARTING RUN: $(date)"
# load dependencies
module load BEDTools/2.27.1-foss-2018b
module load Python/3.8.6-GCCcore-10.2.0
module load SAMtools/1.12-GCC-10.2.0 
export PATH=/nobackup/proj/ghwtcmr/ben/software/STAR-2.7.10a/source:$PATH
export PATH=/nobackup/proj/ghwtcmr/ben/software/bowtie2-2.3.4-linux-x86_64:$PATH
JAVA_HOME=/nobackup/proj/ghwtcmr/ben/software/jdk-16.0.2+7
PATH=$PATH:$HOME/bin:$JAVA_HOME/bin
export JAVA_HOME
export PATH
pfv2_path=/nobackup/proj/ghwtcmr/ben/software/pfv2-ncl-pd-dev
bbmap=/nobackup/proj/ghwtcmr/ben/software/bbmap

# set up env variables
threads=$SLURM_JOB_CPUS_PER_NODE
base_dir=$(pwd)
fastq_dir=$1
reference_dir=$2

# create samplesheet
echo "Getting sample ID "
cd ${fastq_dir}
samplesheet=samples.txt
sample=`head -n ${SLURM_ARRAY_TASK_ID} $samplesheet | tail -1`
echo "Starting run for ${sample} - $(date)"

# folder set up
cd ${base_dir}
mkdir -p pfv2/
cd pfv2/
mkdir -p interleaved/

# interleave PE fastq files into one file
echo "Interleaving ${sample}"
${bbmap}/reformat.sh \
    in1=${fastq_dir}/${sample}_R1.fastq.gz \
    in2=${fastq_dir}/${sample}_R2.fastq.gz \
    out=interleaved/${sample}_interleaved.fastq.gz

# set up PFv2
echo "Running PFv2 set up - $(date)"
sh ${pfv2_path}/setup.sh

# run PTESfinder
echo "running PTESfinder v2 for ${sample} - $(date)"
sh ${pfv2_path}/PFv2.sh \
    -i ${sample} \
    -r interleaved/${sample}_interleaved.fastq.gz \
    -d ./ \
    -S ${reference_dir}/star_index_length_100 \
    -t ${reference_dir}/transcriptome \
    -g ${reference_dir}/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    -b ${reference_dir}/genome \
    -l 100 \
    -c ${pfv2_path}
echo "PFv2 completed for ${sample} - $(date)"