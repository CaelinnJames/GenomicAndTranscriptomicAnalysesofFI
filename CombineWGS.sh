#!/bin/bash

#SBATCH --output=/mnt/shared/projects/sruc/dairy_genetics/bbsrc_feed_efficiency/WGS/_slurmfiles/%j.out
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=1-29
#SBATCH --partition=long
#SBATCH --chdir /mnt/shared/projects/sruc/dairy_genetics/bbsrc_feed_efficiency/WGS/CleanData
#SBATCH --time 14-00:00:00

##### Assigning the number of requested cores to the variable CORES #####
CORES=$SLURM_CPUS_PER_TASK

##### Activating the conda environment #####
source activate WGS

##### Using the array number to define the focal chromosome ##### 
i=$((SLURM_ARRAY_TASK_ID))
echo $i

##### Defining the reference files as variables #####
REFERENCE_FA=/mnt/shared/projects/sruc/dairy_genetics/bbsrc_feed_efficiency/WGS/Reference/Bos_taurus.ARS-UCD1.3.dna_rm.toplevel.fa
REFERENCE_VCF=/mnt/shared/projects/sruc/dairy_genetics/bbsrc_feed_efficiency/WGS/Reference/bos_taurus.vcf.gz
REFERENCE_VCF_CHROM=/mnt/shared/projects/sruc/dairy_genetics/bbsrc_feed_efficiency/WGS/Reference/Chr${i}-Run8-TAUIND-public.vcf.gz

##### Creating the directory and file varaibles #####
tempDIR=/mnt/shared/scratch/cjames/private/WGS/WGSCombined
mkdir -p ${tempDIR}

##### Creating the sample map for GenomicsDBImport #####
echo "CREATING SAMPLE MAP
"
awk -F' ' '{print $0 " "$1}' ../Finished.txt > ${tempDIR}/tmp_$i.txt #Takes Finished.txt (list of sample IDs) and appends each sample ID as a second column. GenomicsDBImport requires a two‑column sample map: <sampleID> <path_to_gvcf>
tr ' ' '\t' <  ${tempDIR}/tmp_$i.txt > ${tempDIR}/SampleMap_$i.txt #Converts spaces to tabs
sed -e 's/$/_genomic_variants.vcf.gz/' -i ${tempDIR}/SampleMap_$i.txt #Appends the suffix of each sample’s GVCF file.

nom=$(sed -n '$=' ${tempDIR}/SampleMap_$i.txt )
echo  " $nom SAMPLES IN SAMPLE MAP FOR CHROMOSOME $i 
" #Checking that all of the samples are there

##### Running GATK Joint Genotyping #####
echo "GATK
"
SampleMap=${tempDIR}/SampleMap_$i.txt 



gatk GenomicsDBImport \  #Creates genomics database for future steps
--genomicsdb-workspace-path ${tempDIR}/database_${i} \   #Imports all GVCFs for chromosome i into a GenomicsDB workspace
--batch-size 30 \   #controls how many samples load at once
--sample-name-map $SampleMap \  #Sample map made earlier
--reader-threads ${CORES} \   #Assigns the number of threads/cores
-L ${i} \    #restricts to the chromosome i
--overwrite-existing-genomicsdb-workspace true   

gatk GenotypeGVCFs \
-R ${REFERENCE_FA} \  #Reference FA file
-V gendb://${tempDIR}/database_${i} \ #Database made in the previous step
-O ${tempDIR}/joint_${i}.vcf.gz #Output

##### SNPs #####

gatk VariantRecalibrator \  #Builds a Gaussian mixture model to score SNPs.
-mode SNP \  #SNPs and indels have to be run seperately
-V ${tempDIR}/joint_${i}.vcf.gz \ #Input VCF
-R ${REFERENCE_FA} \ #Refernce FA
-O ${tempDIR}/joint_${i}_SNPs_recal \ Output
--resource:1000bulls,known=true,training=true,truth=true,prior=10.0 ${REFERENCE_VCF_CHROM} \  #Adds a training resource
#known=true → variants are known but not used for training
#training=true → used to train the Gaussian mixture model
#truth=true → high‑confidence variants
#prior=10.0 → confidence weight (higher = more trusted)
--resource:ensembl,known=true,training=true,truth=true,prior=10.0 ${REFERENCE_VCF} \
-an QD \ #Annotation: Quality by Depth
-an MQ \ #Mapping Quality
-an MQRankSum \ #Mapping quality rank sum test
-an ReadPosRankSum \ #Read position bias
-an FS \ #Fisher strand bias
-an SOR \ #Strand odds ratio
--tranches-file ${tempDIR}/joint_${i}_SNP_tranches #Output file describing sensitivity/specificity tradeoffs

gatk ApplyVQSR \ #Applies the recalibration model to filter SNPs.
-mode SNP \
-V ${tempDIR}/joint_${i}.vcf.gz \
-O ${tempDIR}/joint_${i}_SNP.vcf.gz \
--recal-file ${tempDIR}/joint_${i}_SNPs_recal \ #The recalibration model produced by VariantRecalibrator
--tranches-file ${tempDIR}/joint_${i}_SNP_tranches \ #Tranches file describing cutoffs.
--truth-sensitivity-filter-level 99.0 \ #Keeps variants up to the 99% truth sensitivity threshold. Higher = more variants retained, lower = stricter filtering.
-create-output-variant-index true #Automatically generates a .tbi index.

##### Indels #####

gatk VariantRecalibrator \
-mode INDEL \
-V ${tempDIR}/joint_${i}.vcf.gz \
-R ${REFERENCE_FA} \
-O ${tempDIR}/joint_${i}_INDEL_recal \
--resource:1000bulls,known=true,training=true,truth=true,prior=10.0 ${REFERENCE_VCF_CHROM} \
--resource:ensembl,known=true,training=true,truth=true,prior=10.0 ${REFERENCE_VCF} \
-an QD \
-an MQ \
-an MQRankSum \
-an ReadPosRankSum \
-an FS \
-an SOR \
--tranches-file ${tempDIR}/joint_${i}_INDEL_tranches

gatk ApplyVQSR \
-mode INDEL \
-V ${tempDIR}/joint_${i}.vcf.gz \
-O ${tempDIR}/joint_${i}_INDEL.vcf.gz \
--recal-file ${tempDIR}/joint_${i}_INDEL_recal \
--tranches-file ${tempDIR}/joint_${i}_INDEL_tranches \
--truth-sensitivity-filter-level 99.0 \
-create-output-variant-index true