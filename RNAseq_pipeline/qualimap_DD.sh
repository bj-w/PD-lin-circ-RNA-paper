#! /bin/bash
#SBATCH -A ghwtcmr
#SBATCH --mem=32G
#SBATCH --exclude sb107,sb021,sb012
#SBATCH -c 8

Out=$1
Job=$2
Species=$3
paired=$4

module load Java/11.0.2

qualimap_path="/nobackup/proj/ghwtcmr/Programs_DD/qualimap_v2.2.1"
featurecount_path="/nobackup/proj/ghwtcmr/Programs_DD/subread-2.0.1-Linux-x86_64/bin"



samplesheet=${Out}/samplesheet.txt

name=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet |  awk '{print $1}'` 

hisat_folder=${Out}/${Job}_results/${Job}_hisat2_bam

qualimap_folder=${Out}/${Job}_results/${Job}_qualimap
featurecount_folder=${Out}/${Job}_results/${Job}_featurecount


mkdir ${qualimap_folder}/${name}


#Defining the reference for mouse or human

if [[ "${Species}" == 'human' ]]; then
    echo "Samples are from humans"
    gtf_encode="/nobackup/proj/ghwtcmr/Genomes_DD/qualimap/Homo_sapiens.GRCh38.101.gtf"

elif [[ "${Species}" == 'mouse' ]]; then
    echo "Samples are from mice"
    gtf_encode="/nobackup/proj/ghwtcmr/Genomes_DD/qualimap/Mus_musculus.GRCm39.103.gtf"
    
fi

#Running qualimap and featureCount for paired and single ended files

if [[ "${paired}" == '2' ]]; then
    echo "Samples are paired-end"
    ${qualimap_path}/qualimap rnaseq -a proportional -bam ${hisat_folder}/${name}_hisat2_output.bam -gtf ${gtf_encode} -oc ${qualimap_folder}/${name}_count.txt \
    -outdir ${qualimap_folder}/${name} -outformat HTML  -pe -s -p strand-specific-reverse --java-mem-size=50G 2>> ${Out}/${Job}_logs/${Job}_qualimap.log

    ${featurecount_path}/featureCounts -a ${gtf_encode} -o ${featurecount_folder}/${name} \
    --minOverlap 1 -M --fraction -s 0 -p -T 8 -C \
    ${hisat_folder}/${name}_hisat2_output.bam


elif [[ "${paired}" == '1' ]]; then
    ${qualimap_path}/qualimap rnaseq -a proportional -bam ${hisat_folder}/${name}_hisat2_output.bam -gtf ${gtf_encode} -oc ${qualimap_folder}/${name}_count.txt -outdir ${qualimap_folder}/${name} -outformat HTML  -s  --java-mem-size=50G 2>> ${Out}/${Job}_logs/${Job}_qualimap.log
    
    ${featurecount_path}/featureCounts -a ${gtf_encode} -o ${featurecount_folder}/${name} \
    --minOverlap 1 -M --fraction -s 0 -T 8 \
    ${hisat_folder}/${name}_hisat2_output.bam
fi