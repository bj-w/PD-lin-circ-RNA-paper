#! /bin/bash
#SBATCH -A ghwtcmr
#SBATCH --mem=16G
#SBATCH --exclude sb107,sb021,sb012
#SBATCH -c 8

##Fastq files QCs are running in parallel, submitted to the array of jobs
##Fastq file QC includes fastQC screen for the sequencing quality and fastq_screen to check the contamination from other organisms


Out=$1
Job=$2


module load FastQC/0.11.7-Java-1.8.0_144 
module load Bowtie2/2.3.4.2-foss-2018b
module load SAMtools/1.9-foss-2018b


fastq_path="/nobackup/proj/ghwtcmr/Programs_DD/FastQ-Screen-0.14.1"
fastq_folder=${Out}/${Job}_results/${Job}_fastqc


cd ${Out}/${Job}_rawfiles


##Making the list of fastq files
(ls -d $PWD/*.fastq*) > ${Out}/sample_list.txt
sample_list=${Out}/sample_list.txt



##Assign the sample to the job in the array

sample=`sed -n "$SLURM_ARRAY_TASK_ID"p $sample_list |  awk '{print $1}'` 

fastqc ${sample} --outdir=${fastq_folder} 2>> ${Out}/${Job}_logs/${Job}_fastqc.log

echo "Finished Fastqc analysis"

${fastq_path}/fastq_screen ${sample} --aligner bowtie2 --force \
--threads 8 --outdir ${fastq_folder} 2>> ${Out}/${Job}_logs/${Job}_fastq_screen.log


echo "Finished fastq screen analysis"

