#! /bin/bash
#SBATCH -A ghwtcmr
#SBATCH -p bigmem
#SBATCH --mem=40G
#SBATCH -c 8


Out=$1
Job=$2
Species=$3
paired=$4



threads=$SLURM_JOB_CPUS_PER_NODE
salmon_folder=${Out}/${Job}_results/${Job}_salmon
salmon_path=/mnt/storage/nobackup/proj/ghwtcmr/Programs_DD/salmon-latest_linux_x86_64/bin

cd ${Out}/${Job}_rawfiles

samplesheet=${Out}/samplesheet.txt
name=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet |  awk '{print $1}'` 



if [[ "${Species}" == 'human' ]]; then
    echo "Samples are from humans"
    index=/nobackup/proj/ghwtcmr/Genomes_DD/salmon/salmon_index/salmon_index_human
    

elif [[ "${Species}" == 'mouse' ]]; then
    echo "Samples are from mouse"
    index=/nobackup/proj/ghwtcmr/Genomes_DD/salmon/salmon_index/salmon_index_mouse
    
fi




if [[ "${paired}" == '2' ]]; then
    echo "Samples are paired-end"
    r1=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet |  awk '{print $2}'` 
    r2=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet |  awk '{print $3}'`

    ${salmon_path}/salmon quant -i ${index} \
    -p $threads \
    --libType A \
    -1 ${r1} -2 ${r2} \
    --seqBias --gcBias --posBias --validateMappings \
    -o ${salmon_folder}/${name} 2>> ${Out}/${Job}_logs/${Job}_salmon.log

elif [[ "${paired}" == '1' ]]; then
    echo "Samples are single-ended"
    r1=`sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet |  awk '{print $2}'` 
    ${salmon_path}/salmon quant -i ${index} \
    -p $threads \
    --libType A \
    -r ${r1}  \
    --seqBias --gcBias --posBias --validateMappings \
    -o ${salmon_folder}/${name} 2>> ${Out}/${Job}_logs/${Job}_salmon.log

fi
