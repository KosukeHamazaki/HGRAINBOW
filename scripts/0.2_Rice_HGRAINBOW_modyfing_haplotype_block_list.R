##############################################################
###### 0.2_Rice_HGRAINBOW_modyfing_haplotype_block_list ######
##############################################################

###### 1. Settings ######
##### 1.0. Reset workspace ######
#rm(list=ls())



##### 1.1. Setting working directory to the "HGRAINBOW" directory #####
project <- "HGRAINBOW"
setwd(paste0("/media/hamazaki/d4000953-5e56-40ce-97ea-4cfee57fc91a/research/rice/Project/", project))



##### 1.2. Import packages #####
require(data.table)



##### 1.3. Setting some parameters #####





###### 2. Modification of raw data ######
##### 2.1. Read original files into R #####
geno <- data.frame(fread("data/genotype/L3024_core_extract_L414_ind1A_ind1B_MAF_0.025_geno.tsv"))
block_raw_data <- data.frame(fread("raw_data/genotype/ind1A_ind1B_haplotype_block_min_maf_0.025_list.blocks.det"))



##### 2.2. Extract haplotype block information #####
block.name <- paste0("haploblock_", rep(1:length(block_raw_data[, 5]), block_raw_data[, 5]))

chr <- geno[, 1]
pos <- geno[, 2]
mark <- rownames(geno)
chr.pos <- paste(chr, pos, sep = "_")

n.block <- nrow(block_raw_data)
block.SNPs <- NULL

chrom.old <- 1
mark.correct <- 0
for(i in 1:n.block){
  chrom.now <- block_raw_data[i, 1]
  
  if(chrom.now != chrom.old){
    marks.now.0.1st <- as.numeric(strsplit(block_raw_data[i, 6], "\\|")[[1]])[1]
    mark.real.1st <- block_raw_data[i, 2]
    mark.correct <- marks.now.0.1st - mark.real.1st
  }

  marks.now.0 <- as.numeric(strsplit(block_raw_data[i, 6], "\\|")[[1]])
  marks.now <- marks.now.0 - mark.correct
  mark.block.now <- mark[match(paste(chrom.now, marks.now, sep = "_"), chr.pos)]
  
  block.SNPs <- c(block.SNPs, mark.block.now)
  chrom.old <- chrom.now
}



##### 2.3. Make haplotype block list and save it to .csv #####
block.list <- data.frame(block = block.name, marker = block.SNPs)

file.block <- paste0("data/extra/haplotype_block_list.csv")
fwrite(block.list, file = file.block)
