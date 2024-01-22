#! /bin/bash
#SBATCH -A ghwtcmr
#SBATCH -p bigmem
#SBATCH --mem=64G
#SBATCH -c 12


##This script assumes that strand-specificity is based on dUTP; 

##HISAT2 aligns reads to the known transcriptome, not whole genome - the novel transcripts will not be detected.  



##Input arguments
Out=$1
Job=$2
Species=$3
paired=$4

##Loading required modules
module load SAMtools/1.9-foss-2018b
module load HISAT2/2.1.0-foss-2017b

##Defining variables
threads=$SLURM_JOB_CPUS_PER_NODE
hisat_folder=${Out}/${Job}_results/${Job}_hisat2_bam
samplesheet=${Out}/samplesheet.txt


cd ${Out}/${Job}_rawfiles
name=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet |  awk '{print $1}'` 

if [[ "${Species}" == 'human' ]]; then
    echo "Samples are from humans"
    index="/nobackup/proj/ghwtcmr/Genomes_DD/HISAT2/GRCh38/grch38_tran/genome_tran"


elif [[ "${Species}" == 'mouse' ]]; then
    echo "Samples are from mouse"
    index="/nobackup/proj/ghwtcmr/Genomes_DD/HISAT2/grcm38_tran/genome_tran"
    
fi




if [[ "${paired}" == '2' ]]; then
    echo "Samples are paired-end"


#assigning the sample name and forward and reverse reads
    r1=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet |  awk '{print $2}'` 
    r2=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet |  awk '{print $3}'`

#generate sorted and indexed bam outputs, produces statistics 
    hisat2 -p $threads \
    --dta --rna-strandness RF -x $index \
    -1 $r1 -2 $r2 2>> ${hisat_folder}/${name}_stats_hisat2.txt | \
    tee >(samtools flagstat - > ${hisat_folder}/${name}_hisat2_output.flagstat) | \
    samtools sort -O bam | \
    tee ${hisat_folder}/${name}_hisat2_output.bam | \
    samtools index - ${hisat_folder}/${name}_hisat2_output.bam.bai >>2 ${Out}/${Job}_logs/${Job}_samtools.log


elif [[ "${paired}" == '1' ]]; then
    echo "Samples are single-ended" 
    
#assigning the sample name and forward read

    r1=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet |  awk '{print $2}'`

#generate sorted and indexed bam outputs, produces statistics
    hisat2 -p $threads \
    --dta  -x $index \
    -U $r1  2>> ${hisat_folder}/${name}_stats_hisat2.txt | \
    tee >(samtools flagstat - > ${hisat_folder}/${name}_hisat2_output.flagstat) | \
    samtools sort -O bam | \
    tee ${hisat_folder}/${name}_hisat2_output.bam | \
    samtools index - ${hisat_folder}/${name}_hisat2_output.bam.bai >>2 ${Out}/${Job}_logs/${Job}_samtools.log

fi