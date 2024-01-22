#! /bin/bash
#SBATCH -A ghwtcmr
#SBATCH -p short
#SBATCH --time=0-00:10:00
#SBATCH --mem=8G
#SBATCH --job-name=CIRCexplorer2

base_dir=$(pwd)
fastq_dir=$1
pf_dir=$2
references=$3

# load software
echo "LOADING CIRCexplorer2 CONDA ENVIRONMENT (CIRCexplorer) - $(date)"
module load Anaconda3/2021.11
eval "$(conda shell.bash hook)"
conda activate CIRCexplorer

# create samplesheet
echo "Getting sample ID"
cd ${fastq_dir}
samplesheet=samples.txt
sample=`head -n ${SLURM_ARRAY_TASK_ID} $samplesheet | tail -1`
echo "Starting run for ${sample}"

# set up
cd ${base_dir}
mkdir -p CIRCexplorer2
cd CIRCexplorer2
mkdir -p parse

# parse junction
echo "PARSING STAR CHIMERIC JUNCTIONS FOR ${sample} - $(date)"
CIRCexplorer2 parse \
    -t STAR ${pf_dir}/${sample}/star_Chimeric.out.junction \
    -b parse/${sample}_parsed.bed

# annotate detected circRNAs
echo "ANNOTATING CIRCS IN ${sample} - $(date)"
CIRCexplorer2 annotate \
    -r ${references}/*.genePred \
    -g ${references}/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    -b parse/${sample}_parsed.bed \
    -o ${sample}.bed

echo "Run finished for ${sample} - $(date)"