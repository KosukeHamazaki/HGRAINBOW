###############################################################################
###### 2.1_Rice_HGRAINBOW_Haplotype_based_GWAS_for_RAINBOW_package_paper ######
###############################################################################

###### 1. Settings ######
##### 1.0. Reset workspace ######
#rm(list=ls())



##### 1.1. Setting working directory to the "HGRAINBOW" directory #####
project <- "HGRAINBOW"
setwd(paste0("/media/hamazaki/d4000953-5e56-40ce-97ea-4cfee57fc91a/research/rice/Project/", project))



##### 1.2. Import packages #####
require(data.table)
require(MASS)
require(rrBLUP)
require(SKAT)
require(RAINBOW)

dir.source <- paste0("/home/hamazaki/research_secret/rice/Project/", project)
source(paste0(dir.source, "/1.1_Rice_HGRAINBOW_score_SKAT_geneset.R"))
source(paste0(dir.source, "/1.2_Rice_HGRAINBOW_haplotype_group_fixed_GWAS.R"))
source(paste0(dir.source, "/1.3_Rice_HGRAINBOW_SS_gwas.R"))



##### 1.3. Setting some parameters #####
trial.no <- 1
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
geno <- data.frame(fread("data/genotype/L3024_core_extract_L414_ind1A_ind1B_MAF_0.025_geno.tsv"))
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
map2 <- try(read.csv("data/extra/map2.csv"), silent = TRUE)

if(class(map2) == "error"){
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
  write.csv(map2, "data/extra/map2.csv", row.names = FALSE)
}


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
  save(K.A, file = file.amat)
}




#### 3.0.2. Modify genotype and phenotype data to GWAS format ####
### 3.0.2.1. genotype data ###
geno.GWAS <- data.frame(map, t(x))
rownames(geno.GWAS) <- 1:nrow(geno)
colnames(geno.GWAS) <- c("marker", "chrom", "pos", rownames(x))


### 3.0.2.2. additive relationship matrix ###
colnames(K.A) <- rownames(K.A) <- rownames(x)
Z <- diag(nrow(x))
rownames(Z) <- colnames(Z) <- rownames(x)
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
  print(paste0("iteration_", i, " : normal"))
  save.normal <- paste0(dir.iter, "/normal_")
  normal.res <- RGWAS.normal(pheno = pheno.GWAS, geno = geno.GWAS, ZETA = ZETA, covariate =  NULL, covariate.factor = NULL,
                             structure.matrix = NULL, n.PC = 2, min.MAF = 0.02, P3D = TRUE, n.core = 1,
                             sig.level = 0.05, plot.qq = TRUE, plot.Manhattan = TRUE, plot.method = 1,
                             plot.col1 = c("dark blue", "cornflowerblue"), plot.col2 = 1,
                             plot.type = "p", plot.pch = 16, saveName = save.normal,
                             main.qq = paste0(trait.name, " : normal"), main.man = paste0(trait.name, " : normal"),
                             plot.add.last = TRUE, return.EMM.res = FALSE,
                             thres = FALSE, verbose = FALSE, count = TRUE, time = TRUE)
  points(cum.pos[causals], normal.res[causals, 4], pch = 16, col = c(rep(2, n.causal.site), rep(3, n.gene - 1)), cex = 3)
  abline(v = cum.pos[causals], lwd = 2, col = c(rep(2, n.causal.site), rep(3, n.gene - 1)), lty = 2)
  dev.off()
  
  file.normal <- paste0(save.normal, "simple_result.RData")
  save(normal.res, file = file.normal)
  
  
  
  
  ##### 3.2. Perform haplo-group-fixed (HGF) GWAS #####
  num.haps <- 2:4
  grouping.methods <- c("kmed", "phylo")
  
  for(num.hap in num.haps){
    for(grouping.method in grouping.methods){
      print(paste0("iteration_", i, " : HGF", num.hap, "_", grouping.method))
      save.HGF <- paste0(dir.iter, "/HGF", num.hap, "_", grouping.method, "_")
      HGF.score <- score_HGF_geneset(M = x, y = y, gene.set = block.list, map = map,
                                     ZETA = ZETA, num.hap = num.hap, n.PC = 2, grouping.method = grouping.method)
      HGF.res <- cbind(map2, HGF.score)
      colnames(HGF.res)[4] <- trait.name
      
      png(paste0(save.HGF, trait.name, "_qq.png"))
      qq(scores = HGF.score)
      title(main = paste0(trait.name, " : HGF", num.hap, "_", grouping.method))
      dev.off()
      
      
      png(paste0(save.HGF, trait.name, "_manhattan.png"), width = 800)
      manhattan(input = HGF.res)
      points(cum.pos[causals], HGF.res[haplo.no, 4], pch = 16, col = c(rep(2, n.causal.site), rep(3, n.gene - 1)), cex = 3)
      abline(v = cum.pos[causals], lwd = 2, col = c(rep(2, n.causal.site), rep(3, n.gene - 1)), lty = 2)
      title(main = paste0(trait.name, " : HGF", num.hap, "_", grouping.method))
      dev.off()
      
      file.HGF <- paste0(save.HGF, "simple_result.RData")
      save(HGF.res, file = file.HGF)
      cat("\n")
    }
  }
  
  
  
  
  ##### 3.3. Perform SKAT #####
  print(paste0("iteration_", i, " : SKAT"))
  save.SKAT <- paste0(dir.iter, "/SKAT_")
  SKAT.score <- score_SKAT_geneset(M = x, y = y, gene.set = block.list, map = map)
  cat("\n")
  SKAT.res <- cbind(map2, SKAT.score)
  colnames(SKAT.res)[4] <- trait.name
  
  png(paste0(save.SKAT, trait.name, "_qq.png"))
  qq(scores = SKAT.score)
  title(main = paste0(trait.name, " : SKAT"))
  dev.off()
  
  
  png(paste0(save.SKAT, trait.name, "_manhattan.png"), width = 800)
  manhattan(input = SKAT.res)
  points(cum.pos[causals], SKAT.res[haplo.no, 4], pch = 16, col = c(rep(2, n.causal.site), rep(3, n.gene - 1)), cex = 3)
  abline(v = cum.pos[causals], lwd = 2, col = c(rep(2, n.causal.site), rep(3, n.gene - 1)), lty = 2)
  title(main = paste0(trait.name, " : SKAT"))
  dev.off()
  
  file.SKAT <- paste0(save.SKAT, "simple_result.RData")
  save(SKAT.res, file = file.SKAT)
  
  
  
  
  ##### 3.4. Perform Haplotype-based multi SNP GWAS by RAINBOW #####
  print(paste0("iteration_", i, " : multisnp"))
  save.multisnp <- paste0(dir.iter, "/multisnp_")
  multisnp.res <- RGWAS.multisnp(pheno.GWAS, geno.GWAS, ZETA = ZETA, gene.set = block.list, covariate = NULL, covariate.factor = NULL,
                                 structure.matrix = NULL, n.PC = 2, min.MAF = 0.02, test.method = "LR", n.core = 1,
                                 kernel.method = "linear", kernel.h = "tuned", haplotype = TRUE, num.hap = NULL,
                                 test.effect = "additive", window.size.half = 5, window.slide = 11,
                                 chi0.mixture = 0.5, weighting.center = FALSE,  weighting.other = NULL,
                                 sig.level = 0.05, plot.qq = TRUE, plot.Manhattan = TRUE, plot.method = 1,
                                 plot.col1 = c("dark blue", "cornflowerblue"), plot.col2 = 1,  
                                 plot.type = "p", plot.pch = 16, saveName = save.multisnp,
                                 main.qq = paste0(trait.name, " : multisnp"), main.man = paste0(trait.name, " : multisnp"),
                                 plot.add.last = TRUE, return.EMM.res = FALSE, thres = FALSE,
                                 verbose = FALSE, count = TRUE, time = TRUE)
  points(cum.pos[causals], multisnp.res[haplo.no, 4], pch = 16, col = c(rep(2, n.causal.site), rep(3, n.gene - 1)), cex = 3)
  abline(v = cum.pos[causals], lwd = 2, col = c(rep(2, n.causal.site), rep(3, n.gene - 1)), lty = 2)
  dev.off()
  
  file.multisnp <- paste0(save.multisnp, "simple_result.RData")
  save(multisnp.res, file = file.multisnp)
}
