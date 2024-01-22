#! /bin/bash
#SBATCH -A ghwtcmr
#SBATCH --mem=10G


Out=$1
Job=$2


export LC_ALL=en_GB.utf-8
export LANG=en_GB.utf-8
module load MultiQC/1.7-foss-2018b-Python-3.6.6


fastq_folder=${Out}/${Job}_results/${Job}_fastqc
hisat_folder=${Out}/${Job}_results/${Job}_hisat2_bam
salmon_folder=${Out}/${Job}_results/${Job}_salmon
qualimap_folder=${Out}/${Job}_results/${Job}_qualimap
featurecount_folder=${Out}/${Job}_results/${Job}_featurecount


multiqc ${fastq_folder} -o ${Out}/${Job}_report -n ${Job}_fastq_report
multiqc ${hisat_folder} ${salmon_folder} -o ${Out}/${Job}_report -n ${Job}_alignment_report

multiqc ${qualimap_folder}/*/ -o ${Out}/${Job}_report -n ${Job}_qualimap_report
multiqc ${featurecount_folder} -o ${Out}/${Job}_report -n ${Job}_featurecount_report


