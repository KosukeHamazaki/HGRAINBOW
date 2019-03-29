#########################################################################
###### 3.4_Rice_HG RAINBOW_Manhattan_plot_for_overwhelming_results ######
#########################################################################

###### 1. Settings ######
##### 1.0. Reset workspace ######
#rm(list=ls())



##### 1.1. Setting working directory to the "HGPAG" directory #####
project <- "HGRAINBOW"
setwd(paste0("/media/hamazaki/d4000953-5e56-40ce-97ea-4cfee57fc91a/research/rice/Project/", project))



##### 1.2. Import packages #####
require(data.table)
require(MASS)
require(rrBLUP)
require(SKAT)
require(RAINBOW)

dir.source <- paste0("~/research_secret/rice/Project/", project)
source(paste0(dir.source, "/1.1_Rice_HGRAINBOW_score_SKAT_geneset.R"))
source(paste0(dir.source, "/1.2_Rice_HGRAINBOW_haplotype_group_fixed_GWAS.R"))
source(paste0(dir.source, "/1.3_Rice_HGRAINBOW_SS_gwas.R"))



##### 1.3. Setting some parameters #####
trial.no <- 1
n.causal.site <- 2
direction.eff <- "minus"

n.haplo.mark.min <- 5
n.gene <- 2
h2 <- 0.6
weight.haplo <- 2
prop <- 1
n.iter <- 100

n.top.false <- 10
n.inflator <- length((n.top.false))

method.names <- c("Single-SNP", "HGF2_kmed", "HGF2_phylo", "HGF3_kmed",
                  "HGF3_phylo", "HGF4_kmed", "HGF4_phylo", "SKAT", "RAINBOW")
method.levels <- c("RAINBOW", "Single-SNP", "HGF2_kmed", "HGF2_phylo", "HGF3_kmed",
                   "HGF3_phylo", "HGF4_kmed", "HGF4_phylo", "SKAT")
method.fac <- factor(rep(method.names, each = n.iter),
                     levels = method.levels)
n.method <- length(method.names)

method.names.abb <- c("SS", "H2k", "H2p", "H3k", "H3p", "H4k", "H4p", "SK", "R")
method.levels.abb <- c("R", "SS", "H2k", "H2p", "H3k", "H3p", "H4k", "H4p", "SK")
method.fac.abb <- factor(rep(method.names.abb, each = n.iter),
                         levels = method.levels.abb)


dir.base <- paste0("midstream/number_of_causals=", n.causal.site,
                   "_direction_of_effect=", direction.eff,
                   "_trial_no=", trial.no)


dir.res.base <- paste0("results/number_of_causals=", n.causal.site,
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
no.now <- 9
thres.now <- "BH_0.01"
dir.overwhelm <- paste0(dir.res.base, "/overwhelming_results_", method.names.abb[no.now], "_", thres.now,
                        "_#_of_top_false_", n.top.false)
dir.overwhelm.manhattan2 <- paste0(dir.res.base, "/overwhelming_results_", method.names.abb[no.now],  "_", thres.now,
                                   "_manhattan2")
dir.create(dir.overwhelm.manhattan2)
check.no <- as.numeric(substr(unique(list.files(dir.overwhelm)), 11, 12))

for(i in check.no){
  print(paste0("iteration_", i))
  dir.iter <- paste0(dir.iter.0, i)
  dir.man.iter <- paste0(dir.overwhelm.manhattan2, "/iteration_", i)
  dir.create(dir.man.iter)
  
  ### 3.0.2.3. phenotype data ###
  y <- pheno[i, ]
  trait.name <- paste0("Iteration_", i)
  
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
  input.normal <- paste0(dir.iter, "/normal_")
  save.normal <- paste0(dir.man.iter, "/SS_")
  file.normal <- paste0(input.normal, "simple_result.RData")
  load(file = file.normal)
  
  
  ##### 3.2. Perform haplo-group-fixed (HGF) GWAS #####
  input.HGF <- paste0(dir.iter, "/HGF2_phylo_")
  save.HGF <- paste0(dir.man.iter, "/HGF2p_")
  file.HGF <- paste0(input.HGF, "simple_result.RData")
  load(file = file.HGF)
  
  
  ##### 3.3. Perform SKAT #####
  input.SKAT <- paste0(dir.iter, "/SKAT_")
  save.SKAT <- paste0(dir.man.iter, "/SK_")
  file.SKAT <- paste0(input.SKAT, "simple_result.RData")
  load(file = file.SKAT)
  
  
  ##### 3.4. Perform Haplotype-based multi SNP GWAS by RAINBOW #####
  input.multisnp <- paste0(dir.iter, "/multisnp_")
  save.multisnp <- paste0(dir.man.iter, "/R_")
  file.multisnp <- paste0(input.multisnp, "simple_result.RData")
  load(file = file.multisnp)
  
  
  y.max <- max(normal.res[, 4], HGF.res[, 4],
               SKAT.res[, 4], multisnp.res[, 4]) + 1
  
  
  
  png(paste0(save.multisnp, trait.name, "_manhattan_2.png"), width = 850)
  op <- par(mar = c(7, 8, 4, 1), mgp = c(4.5, 1.8, 0))
  manhattan(input = multisnp.res, cex.lab = 2.8, cex.axis.x = 2,
            cex.axis.y = 3, y.max = y.max, lwd.thres = 3, sig.level = 0.01)
  points(cum.pos[causals], multisnp.res[haplo.no, 4], pch = 16, col = c("red", "red", "purple"), cex = 4)
  abline(v = cum.pos[causals], lwd = 4, col = c("red", "red", "purple"), lty = 2)
  title(main = paste0(trait.name, " : RAINBOW"), cex.main = 3.5)
  dev.off()
  par(op)
  
  
  png(paste0(save.normal, trait.name, "_manhattan_2.png"), width = 850)
  op <- par(mar = c(7, 8, 4, 1), mgp = c(4.5, 1.8, 0))
  manhattan(input = normal.res, cex.lab = 2.8, cex.axis.x = 2,
            cex.axis.y = 3, y.max = y.max, lwd.thres = 3, sig.level = 0.01)
  points(cum.pos[causals], normal.res[causals, 4], pch = 16, col = c("red", "red", "purple"), cex = 4)
  abline(v = cum.pos[causals], lwd = 4, col = c("red", "red", "purple"), lty = 2)
  title(main = paste0(trait.name, " : Single-SNP"), cex.main = 3.5)
  dev.off()
  par(op)
  
  
  
  
  
  png(paste0(save.HGF, trait.name, "_manhattan_2.png"), width = 850)
  op <- par(mar = c(7, 8, 4, 1), mgp = c(4.5, 1.8, 0))
  manhattan(input = HGF.res, cex.lab = 2.8, cex.axis.x = 2,
            cex.axis.y = 3, y.max = y.max, lwd.thres = 3, sig.level = 0.01)
  points(cum.pos[causals], HGF.res[haplo.no, 4], pch = 16, col = c("red", "red", "purple"), cex = 4)
  abline(v = cum.pos[causals], lwd = 4, col = c("red", "red", "purple"), lty = 2)
  title(main = paste0(trait.name, " : HGF2p"), cex.main = 3.5)
  dev.off()
  par(op)
  
  
  
  
  
  png(paste0(save.SKAT, trait.name, "_manhattan_2.png"), width = 850)
  op <- par(mar = c(7, 8, 4, 1), mgp = c(4.5, 1.8, 0))
  manhattan(input = SKAT.res, cex.lab = 2.8, cex.axis.x = 2,
            cex.axis.y = 3, y.max = y.max, lwd.thres = 3, sig.level = 0.01)
  points(cum.pos[causals], SKAT.res[haplo.no, 4], pch = 16, col = c("red", "red", "purple"), cex = 4)
  abline(v = cum.pos[causals], lwd = 4, col = c("red", "red", "purple"), lty = 2)
  title(main = paste0(trait.name, " : SKAT"), cex.main = 3.5)
  dev.off()
  par(op)

  
  
  
  
  
  
  
  save.all <- paste0(dir.overwhelm.manhattan2, "/ALL_")
  png(paste0(save.all, trait.name, "_manhattan_2.png"), height = 950, width = 1600)
  op <- par(mfrow = c(2, 2), mar = c(7, 8, 5, 1), mgp = c(4.5, 1.8, 0))
  manhattan(input = multisnp.res, cex.lab = 2.8, cex.axis.x = 2.4,
            cex.axis.y = 3, y.max = y.max, lwd.thres = 3, sig.level = 0.01)
  points(cum.pos[causals], multisnp.res[haplo.no, 4], pch = 16, col = c("red", "red", "purple"), cex = 5)
  abline(v = cum.pos[causals], lwd = 4, col = c("red", "red", "purple"), lty = 2)
  title(main = paste0("RAINBOW"), cex.main = 4.5)

  manhattan(input = normal.res, cex.lab = 2.8, cex.axis.x = 2.4,
            cex.axis.y = 3, y.max = y.max, lwd.thres = 3, sig.level = 0.01)
  points(cum.pos[causals], normal.res[causals, 4], pch = 16, col = c("red", "red", "purple"), cex = 5)
  abline(v = cum.pos[causals], lwd = 4, col = c("red", "red", "purple"), lty = 2)
  title(main = paste0("Single-SNP"), cex.main = 4.5)

  manhattan(input = HGF.res, cex.lab = 2.8, cex.axis.x = 2.4,
            cex.axis.y = 3, y.max = y.max, lwd.thres = 3, sig.level = 0.01)
  points(cum.pos[causals], HGF.res[haplo.no, 4], pch = 16, col = c("red", "red", "purple"), cex = 5)
  abline(v = cum.pos[causals], lwd = 4, col = c("red", "red", "purple"), lty = 2)
  title(main = paste0("HGF2p (k=2, UPGMA)"), cex.main = 4.5)

  
  manhattan(input = SKAT.res, cex.lab = 2.8, cex.axis.x = 2.4,
            cex.axis.y = 3, y.max = y.max, lwd.thres = 3, sig.level = 0.01)
  points(cum.pos[causals], SKAT.res[haplo.no, 4], pch = 16, col = c("red", "red", "purple"), cex = 5)
  abline(v = cum.pos[causals], lwd = 4, col = c("red", "red", "purple"), lty = 2)
  title(main = paste0("SKAT"), cex.main = 4.5)
  dev.off()
  par(op)  
}
