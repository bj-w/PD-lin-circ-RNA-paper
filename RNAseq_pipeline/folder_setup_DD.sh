#! /bin/bash
#SBATCH -A ghwtcmr
#SBATCH --mem=12G


In=$1
Out=$2
Job=$3
paired=$4

cd $Out
mkdir ${Job}_logs ${Job}_report ${Job}_results ${Job}_rawfiles

#Making separate folders for the results of the analysis
fastq_folder=${Out}/${Job}_results/${Job}_fastqc
mkdir ${fastq_folder}

hisat_folder=${Out}/${Job}_results/${Job}_hisat2_bam
mkdir ${hisat_folder}

salmon_folder=${Out}/${Job}_results/${Job}_salmon
mkdir ${salmon_folder}

qualimap_folder=${Out}/${Job}_results/${Job}_qualimap
mkdir ${qualimap_folder}

featurecount_folder=${Out}/${Job}_results/${Job}_featurecount
mkdir ${featurecount_folder}

 

cd $In


##Sample concatenation; assumes that if the files are ending in _L*_R[1-2]_001.fastq.gz, then concatenation is required, if the ending is "*_R[1-2].fastq.gz", then concatenation is not reqiored

if ls $In/*R1.fastq.gz  
then
    echo 'No concatentation needed, files are ending in R*.fastq.gz'
    cp *.fastq.gz ${Out}/${Job}_rawfiles
    echo 'Fastq files copied'
else 
    echo "Concatenation required"
    ls -1 *R*_001.fastq.gz | sed 's/_L[0-9]*_R[1-2]_001.fastq.gz//' | uniq > ${Out}/${Job}_rawfiles/fastq_combine.txt

    samplesheet=${Out}/${Job}_rawfiles/fastq_combine.txt
    END=$(wc -l < ${samplesheet})

    if [[ "${paired}" == '2' ]]; then
        for i in $(seq 1 $END)
        do
        sample=$(sed -n "$i"p $samplesheet |  awk '{print $1}') ;
        cat $(ls ${sample}_L*_R1_001.fastq.gz) >  ${Out}/${Job}_rawfiles/${sample}_R1.fastq.gz
        cat $(ls ${sample}_L*_R2_001.fastq.gz) >  ${Out}/${Job}_rawfiles/${sample}_R2.fastq.gz
        done

    elif [[ "${paired}" == '1' ]]; then

        for i in $(seq 1 $END)
        do
        sample=$(sed -n "$i"p $samplesheet |  awk '{print $1}') ;
        cat $(ls ${sample}_L*_R1_001.fastq.gz) >  ${Out}/${Job}_rawfiles/${sample}_R1.fastq.gz
        done

    fi

    rm ${Out}/${Job}_rawfiles/fastq_combine.txt

    echo "The fastq files are concatenated and copied to the working directory"

fi





#make a samplesheet matrix; will generate name name_R1.fastq.gz name_R2.fastq.gz; number of lines will be equal to the number of files. 
cd ${Out}/${Job}_rawfiles


for f in `ls -1 *_R1.fastq.gz | sed 's/_R1.fastq.gz//' `
do
echo ${f} ${f}_R1.fastq.gz ${f}_R2.fastq.gz >> ${Out}/samples.txt
done

#To prevent the duplicates when the pipeline is rerun
sort ${Out}/samples.txt | uniq > ${Out}/samplesheet.txt

rm ${Out}/samples.txt

echo "The set up of the folders and generation of the samplesheet are done"
