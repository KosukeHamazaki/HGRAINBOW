############################################################
###### 3.2_Rice_HGPAG_Manhattan_plot_from_the_results ######
############################################################

###### 1. Settings ######
##### 1.0. Reset workspace ######
#rm(list=ls())



##### 1.1. Setting working directory to the "HGPAG" directory #####
project <- "HGPAG"
setwd(paste0("/media/hamazaki/d4000953-5e56-40ce-97ea-4cfee57fc91a/research/rice/Project/", project))



##### 1.2. Import packages #####
require(data.table)
require(MASS)
require(rrBLUP)
require(SKAT)
require(RAINBOW)

dir.source <- paste0("/home/hamazaki/research_secret/rice/Project/", project)
source(paste0(dir.source, "/1.1_Rice_HGPAG_score_SKAT_geneset.R"))
source(paste0(dir.source, "/1.2_Rice_HGPAG_haplotype_group_fixed_GWAS.R"))
source(paste0(dir.source, "/1.3_Rice_HGPAG_SS_gwas.R"))



##### 1.3. Setting some parameters #####
trial.no <- 2
n.causal.site <- 2
direction.eff <- "plus"

n.haplo.mark.min <- 5
n.gene <- 2
h2 <- 0.6
weight.haplo <- 2
prop <- 1
n.iter <- 100


dir.base <- paste0("midstream/number_of_causals=", n.causal.site,
                   "_direction_of_effect=", direction.eff,
                   "_trial_no=", trial.no)




###### 2. Modification of raw data ######
##### 2.1. Read original files into R #####
geno <- data.frame(fread("data/genotype/L3024_core_extract_L209_ind1A_MAF_0.05_geno.tsv"))
block.list <- data.frame(fread("data/extra/haplotype_block_list.csv"))

file.pheno <- paste0("data/phenotype/number_of_causals=", n.causal.site,
                     "_direction_of_effect=", direction.eff,
                     "_trial_no=", trial.no, ".csv")
pheno <- as.matrix(read.csv(file = file.pheno))

##### 2.2. Some modification #####
chr <- geno[, 1]
pos <- geno[, 2]
marker <- rownames(geno) 
map <- data.frame(marker, chr, pos)

chr.tab <- table(chr)
chr.max <- max(chr)
chr.cum <- cumsum(chr.tab)
cum.pos <- pos
if (length(chr.tab) != 1) {
  for (i in 1:(chr.max - 1)) {
    cum.pos[(chr.cum[i] + 1):(chr.cum[i + 1])] <- pos[(chr.cum[i] + 
                                                         1):(chr.cum[i + 1])] + cum.pos[chr.cum[i]]
  }
}


gene.set <- block.list
gene.names <- as.character(gene.set[, 1])
mark.id <- as.character(gene.set[, 2])
gene.name <- as.character(unique(gene.names))
n.scores <- length(unique(gene.set[, 1]))
chr.set.mean <- pos.set.mean <- cum.pos.set.mean <- rep(NA, 
                                                        n.scores)
ids <- as.data.frame(matrix(rep(NA, n.scores * 2), 
                            ncol = 2))
marker.now <- gene.name
for (k in 1:n.scores) {
  id <- mark.id[gene.names == gene.name[k]]
  ids[k, ] <- c(as.character(id[1]), as.character(id[length(id)]))
  num.sel <- match(id, map[, 1])
  chr.sel <- map[num.sel, 2]
  pos.sel <- map[num.sel, 3]
  cum.pos.sel <- cum.pos[num.sel]
  chr.set.mean[k] <- chr.sel[1]
  pos.set.mean[k] <- mean(pos.sel)
  cum.pos.set.mean[k] <- mean(cum.pos.sel)
}
map2 <- data.frame(marker = marker.now, chr = chr.set.mean, 
                   pos = pos.set.mean)


x <- t(geno[, -c(1:2)])
n.line <- nrow(x)
n.marker <- ncol(x)
line.names <- rownames(x)

n.haplo <- length(unique(block.list[, 1]))
n.marker.haplo <- nrow(block.list)




###### 3. Perform GWAS by 4 methods ######
##### 3.0. Some preparations #####
#### 3.0.0. Define directory name to save the results ####
dir.iter.0 <- paste0(dir.base, "/iteration_")



#### 3.0.1. Estimate additive genetic relationship matrix ####
file.amat <- paste0("data/extra/additive_relationship_matrix.RData")
load.amat <- try(load(file.amat), silent = TRUE)
if(class(load.amat) == "try-error"){
  K.A <- A.mat(x)
}
#save(K.A, file = file.amat)




#### 3.0.2. Modify genotype and phenotype data to GWAS format ####
### 3.0.2.1. genotype data ###
geno.GWAS <- data.frame(map, t(x))
rownames(geno.GWAS) <- 1:nrow(geno)
colnames(geno.GWAS) <- c("marker", "chrom", "pos", rownames(x))


### 3.0.2.2. additive relationship matrix ###
colnames(K.A) <- rownames(K.A) <- rownames(x)
Z <- diag(nrow(x))
ZETA <- list(A = list(Z = Z, K = K.A))  ### Make a list for input of the function ; in this case, additive only.


for(i in 1:n.iter){
  print(paste0("iteration_", i))
  dir.iter <- paste0(dir.iter.0, i)
  
  ### 3.0.2.3. phenotype data ###
  y <- pheno[i, ]
  trait.name <- paste0("simulation_", i)
  
  pheno.GWAS <- data.frame(names(y), y)
  colnames(pheno.GWAS) <- c("Line_names", trait.name)
  
  file.gene <- paste0(dir.iter, "/causal_information.csv")
  causal.info <- data.frame(fread(file.gene))
  causals <- causal.info[, 1]
  
  gene.name <- unique(block.list[, 1])
  gene.names <- block.list[, 1]
  mark.id <- block.list[, 2]
  
  haplo.name <- block.list[match(marker[causals], mark.id), 1]
  
  for(j in 1:length(causals)){
    if(is.na(haplo.name)[j]){
      cum.pos.geneset <- cum.pos[match(mark.id, marker)]
      nearest <- which.min(abs(cum.pos[causals[j]] - cum.pos.geneset))
      
      haplo.name[j] <- block.list[nearest, 1]
    }
  }
  
  
  haplo.no <- match(haplo.name, gene.name)
  gene.candidate <- unique(haplo.no)
  
  ##### 3.1. Perform Normal GWAS by RAINBOW #####
  save.normal <- paste0(dir.iter, "/normal_")
  file.normal <- paste0(save.normal, "simple_result.RData")
  load(file = file.normal)
  
  png(paste0(save.normal, trait.name, "_manhattan_2.png"), width = 800)
  manhattan(input = normal.res)
  points(cum.pos[causals], normal.res[causals, 4], pch = 16, col = c(2, 2, 3), cex = 2)
  abline(v = cum.pos[causals], lwd = 2, col = c(2, 2, 3), lty = 2)
  title(main = paste0(trait.name, " : normal"))
  dev.off()
  
  
  
  
  
  ##### 3.2. Perform haplo-group-fixed (HGF) GWAS #####
  num.haps <- c(2, 3)
  
  for(num.hap in num.haps){
    save.HGF <- paste0(dir.iter, "/HGF", num.hap, "_")
    file.HGF <- paste0(save.HGF, "simple_result.RData")
    load(file = file.HGF)
    
    png(paste0(save.HGF, trait.name, "_manhattan_2.png"), width = 800)
    manhattan(input = HGF.res)
    points(cum.pos[causals], HGF.res[haplo.no, 4], pch = 16, col = c(2, 2, 3), cex = 2)
    abline(v = cum.pos[causals], lwd = 2, col = c(2, 2, 3), lty = 2)
    title(main = paste0(trait.name, " : HGF", num.hap))
    dev.off()
    
  }
  
  
  
  
  ##### 3.3. Perform SKAT #####
  save.SKAT <- paste0(dir.iter, "/SKAT_")
  file.SKAT <- paste0(save.SKAT, "simple_result.RData")
  load(file = file.SKAT)


  png(paste0(save.SKAT, trait.name, "_manhattan_2.png"), width = 800)
  manhattan(input = SKAT.res)
  points(cum.pos[causals], SKAT.res[haplo.no, 4], pch = 16, col = c(2, 2, 3), cex = 2)
  abline(v = cum.pos[causals], lwd = 2, col = c(2, 2, 3), lty = 2)
  title(main = paste0(trait.name, " : SKAT"))
  dev.off()
  
  
  
  
  
  ##### 3.4. Perform Haplotype-based multi SNP GWAS by RAINBOW #####
  save.multisnp <- paste0(dir.iter, "/multisnp_")
  file.multisnp <- paste0(save.multisnp, "simple_result.RData")
  load(file = file.multisnp)

  
  png(paste0(save.multisnp, trait.name, "_manhattan_2.png"), width = 800)
  manhattan(input = multisnp.res)
  points(cum.pos[causals], multisnp.res[haplo.no, 4], pch = 16, col = c(2, 2, 3), cex = 2)
  abline(v = cum.pos[causals], lwd = 2, col = c(2, 2, 3), lty = 2)
  title(main = paste0(trait.name, " : multisnp"))
  dev.off()
}
