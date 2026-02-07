Scripts for processing whole genome sequence data for holstein dairy cattle 

1) WGSProcessing.sh -> Processes raw WGS files individually
2) CombineWGS.sh -> Combines processed WGS files for all individuals and extracts SNPs and Indels using Ensembl and 1000 bulls as reference genomes. This is done as an array job with each array being 1 chromosome
3) CombineWGS2.sh -> Combines the SNP and Indel files for each chromosome and converts them to Plink Bfile format
