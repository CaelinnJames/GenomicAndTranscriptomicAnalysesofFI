#!/bin/bash

#SBATCH --output=/mnt/shared/projects/sruc/dairy_genetics/bbsrc_feed_efficiency/WGS/_slurmfiles/%j.out
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --chdir /mnt/shared/projects/sruc/dairy_genetics/bbsrc_feed_efficiency/WGS/CleanData

##### Assigning the number of requested cores to the variable CORES #####
CORES=$SLURM_CPUS_PER_TASK

##### Activating the conda environment #####
source activate WGS

##### Creating the directory and file varaibles #####
tempDIR=/mnt/shared/scratch/cjames/private/WGS/WGSCombined
mkdir -p ${tempDIR}

##### Combining the full VCF files #####
ls ${tempDIR}/joint_{[0-9],[0-9][0-9]}.vcf.gz > variant_files.list #Lists all per‑chromosome joint VCFs (joint_1.vcf.gz … joint_29.vcf.gz)
picard MergeVcfs I=variant_files.list O=${tempDIR}/joint.vcf.gz #Merges all chromosome‑specific VCFs into a single whole‑genome VCF. MergeVcfs reads the list file and concatenates VCFs in order.

##### Combining the SNP VCF files #####
ls ${tempDIR}/joint_{[0-9],[0-9][0-9]}_SNP.vcf.gz > SNP_variant_files.list 
picard MergeVcfs I=SNP_variant_files.list O=${tempDIR}/joint_SNP.vcf.gz 

##### Combining the Indel VCF files #####
ls ${tempDIR}/joint_{[0-9],[0-9][0-9]}_INDEL.vcf.gz > INDEL_variant_files.list 
picard MergeVcfs I=INDEL_variant_files.list O=${tempDIR}/joint_INDEL.vcf.gz 

##### Move final output files to /projects/ #####
cp ${tempDIR}/joint.vcf.gz ./joint.vcf.gz
cp ${tempDIR}/joint_SNP.vcf.gz ./joint_SNP.vcf.gz
cp ${tempDIR}/joint_INDEL.vcf.gz ./joint_INDEL.vcf.gz


##### Change to conda env with plink #####
conda deactivate WGS
source activate GWAS

##### Use plink to convert to bfile format #####
plink --vcf joint.vcf.gz --recode --make-bed --cow --set-missing-id-vars --out joint_06_24
plink --vcf joint_SNP.vcf.gz --recode --make-bed --cow --set-missing-id-vars --out joint_06_24_SNPs
plink --vcf joint_INDEL.vcf.gz --recode --make-bed --cow --set-missing-id-vars --out joint_06_24_INDELs
