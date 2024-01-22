#! /bin/bash
#SBATCH -A ghwtcmr
#SBATCH -p bigmem
#SBATCH --mem=10G


###This is the master script for RNASeq analysis of paired end or single-end human or mouse RNASeq.####
### It requires 5 arguments: 
##In - path to fastq files; if they are not concatenated, the endings should be  *_R[1-2]_001.fastq.gz, alternatively for the alreaday concatenated samples the ending should be *_R[1-2].fastq.gz
##Out - path to the output, the folder must exist
##Job - the name of the batch analysis, the folder names and files will start with this job name; e.g. if the job name is march10_2021, then the folder will be named march10_2021_report
##Species - available options are mouse and human 
##Paired or single end: for paired the argument is 2, for single-end is 1

In=$1
Out=$2
Job=$3
Species=$4
paired=$5


script_path="/mnt/nfs/home/b8057117/scripts/rnaseq_dd_bw"


JobID1="Set up folder and concatenating the files"
JobID2="Running fastq QC"
JobID3="HISAT2 alignment"
JobID4="Salmon mapping"
JobID5="QC controls"
JobID6="Make a report"


srun --job-name="${JobID1}" ${script_path}/folder_setup_DD.sh "${In}" "${Out}" "${Job}" "${paired}" 

echo "setting up folder system"
sleep 2

wait

##Calculates 'number of samples' variable for the array job submission
no_of_lines=$(wc -l < $2/samplesheet.txt)
echo "${no_of_lines}"

job2=$(sbatch --parsable --job-name="${JobID2}" --array=1-"$((${no_of_lines}*2))"%24 ${script_path}/fastqc_rna_DD.sh "${Out}" "${Job}")

sleep 2
echo ${job2}


job3=$(sbatch --parsable --job-name="${JobID3}" --array=1-"${no_of_lines}"%24 ${script_path}/hisat_DD.sh "${Out}" "${Job}" "${Species}" "${paired}")

sleep 2
echo ${job3}
j3=${job3} 

job4=$(sbatch --parsable --job-name="${JobID4}" --array=1-"${no_of_lines}"%24 ${script_path}/salmon_DD.sh "${Out}" "${Job}" "${Species}" "${paired}")

sleep 2
echo ${job4}


job5=$(sbatch --parsable --dependency=aftercorr:$j3 --job-name="${JobID5}" --array=1-"${no_of_lines}"%24 ${script_path}/qualimap_DD.sh "${Out}" "${Job}" "${Species}" "${paired}")

sleep 2
echo ${job5}
j5=${job5}

job6=$(sbatch --parsable --dependency=aftercorr:$j5 --job-name="${JobID6}" ${script_path}/report_DD.sh "${Out}" "${Job}" )

sleep 2
echo ${job6} 