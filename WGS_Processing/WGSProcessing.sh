#!/bin/bash

#SBATCH --output=/mnt/shared/projects/sruc/dairy_genetics/bbsrc_feed_efficiency/WGS/_slurmfiles/WGSProcessing_%j.out
#SBATCH --array=1-121
#SBATCH --cpus-per-task=8
#SBATCH --mem=58G
#SBATCH --partition=long
#SBATCH --chdir /mnt/shared/projects/sruc/dairy_genetics/bbsrc_feed_efficiency/WGS/
#SBATCH --time 13-00:00:00

##### Assigning the number of requested cores to the variable CORES #####
CORES=$SLURM_CPUS_PER_TASK

##### Activating the conda environment #####
source activate WGS

##### Defining the reference files as variables #####
REFERENCE_FA=Reference/Bos_taurus.ARS-UCD1.3.dna_rm.toplevel.fa
REFERENCE_VCF=Reference/bos_taurus.vcf.gz

##### Using the array number to extract the sample ID #####
i=$((SLURM_ARRAY_TASK_ID)) #Assigning the array number to the variable i
echo $i #Making sure that the variable $i has been set properly
sample=$(sed "${i}q;d" Array.txt ) #Array.txt is a .txt file with each sample ID on a new line, this takes the ith line of Array.txt and assigns the output to the variable sample

##### Creating the directory and file varaibles #####
DIR=test/clean/${sample} #This is the path to where the raw files for this sample are stored

if [[ -f ${DIR}/Merged*fq.gz ]]
then
files=(${DIR}/Merged_1.fq.gz) #Some samples came with more than 2 raw files, for those samples before running this script I merged the forwards and back reads into one file each which started with "Merged"
else
files=(${DIR}/*_1.fq.gz) #If I didn't need to merge them then we're happy just to take the raw file name
fi

file=${files[0]}
FILE=${file::-8} #Removes _1.fq.gz from the file variable

echo $FILE #Checking to make sure the file prefix has been assigned to the FILE variable correctly

tempDIR=/mnt/shared/scratch/cjames/private/WGS/${sample} #So that the intermediate files don't clog up /projects/, I've made a tmp folder in /scratch/
mkdir -p ${tempDIR} #Making the tmp folder if it doesn't already exist
tempID=${tempDIR}/${sample} #This means I don't have to type out ${tempDIR}/${sample} at each step

##### STARTING THE WGS PROCESSING! #####

#### Fastp to check quality of reads and potentially do any trimming if needed ####
echo "FASTP RUNNING
" #To make this section of the output file easier to find
fastp -i ${FILE}_1.fq.gz -I ${FILE}_2.fq.gz -o ${tempID}_1.QC.fastq.gz -O ${tempID}_2.QC.fastq.gz -j ${tempID}.fastp.json -h ${tempID}.fastp.html 

#### SAM allignment and sorting ####
echo "ALLIGNING SAM
"
bwa mem -M -t ${CORES} -R "@RG\tID:${tempID}\tSM:${tempID}\tPL:ILLUMINA" ${REFERENCE_FA} ${tempID}_1.QC.fastq.gz ${tempID}_2.QC.fastq.gz > ${tempID}_aligned.sam
samtools view -b -@ ${CORES}  ${tempID}_aligned.sam -o  ${tempID}_aligned.bam
samtools sort -@ ${CORES} ${tempID}_aligned.bam -o ${tempID}_aligned_sorted.bam 
picard CollectAlignmentSummaryMetrics I=${tempID}_aligned_sorted.bam O=${tempID}_metrics.txt R=${REFERENCE_FA} # Creates a metric file to check how the allignment went

#### Markign duplicates ####
echo "
MarkDuplicates
"
picard -Xmx50g   MarkDuplicates I=${tempID}_aligned_sorted.bam O=${tempID}_aligned_sorted_dedup.bam M=${tempID}_aligned_sorted_dedup.metrics.txt

#### Samtolls sorting ####
echo "
Sort it again
"
samtools sort ${tempID}_aligned_sorted_dedup.bam -o ${tempID}_aligned_sorted_dedup_sorted.bam -@ ${CORES}

#### Base quality score recalibration ####
gatk BaseRecalibrator -R ${REFERENCE_FA} -I ${tempID}_aligned_sorted_dedup_sorted.bam --known-sites ${REFERENCE_VCF} -O ${tempID}_recal_data.table #Creates a data table for the next steps
gatk ApplyBQSR -R ${REFERENCE_FA} -I ${tempID}_aligned_sorted_dedup_sorted.bam -bqsr ${tempID}_recal_data.table -O ${tempID}_aligned_sorted_dedup_sorted_bqsr.bam #Uses the data table for recallibration

#### Identify SNPs and indels ####
gatk HaplotypeCaller -R ${REFERENCE_FA} -I ${tempID}_aligned_sorted_dedup_sorted_bqsr.bam -dbsnp ${REFERENCE_VCF} -stand-call-conf 30 -O ${tempID}_genomic_variants.vcf.gz -ERC GVCF

##### Copying processed WGS to /projects/ #####
cp ${tempID}_genomic_variants.vcf.gz CleanData/${sample}_genomic_variants.vcf.gz
cp ${tempID}_genomic_variants.vcf.gz.tbi CleanData/${sample}_genomic_variants.vcf.gz.tbi
