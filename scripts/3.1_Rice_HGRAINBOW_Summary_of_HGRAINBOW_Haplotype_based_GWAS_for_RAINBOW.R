######################################################################################
###### 3.1_Rice_HGRAINBOW_Summary_of_HGRAINBOW_Haplotype_based_GWAS_for_RAINBOW ######
######################################################################################

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
require(ggplot2)
require(gridExtra)

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

n.top.false <- c(1, 10, 20, 50)
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
dir.create(dir.res.base)

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




###### 3. Summary results of GWAS by 4 methods ######
##### 3.0. Some preparations #####
#### 3.0.0. Define directory name to save the results ####
dir.iter.0 <- paste0(dir.base, "/iteration_")
method.now <- c("BH", "Bonf")
sig.levels <- c(0.05, 0.01)
n.thres <- length(method.now) * length(sig.levels)
thres.names <- paste(rep(method.now, each = length(sig.levels)),
                     rep(sig.levels, length(method.now)), sep = "_")

MAFs <- causals.mat <- matrix(nrow = n.iter, ncol = n.causal.site + (n.gene - 1))
FDRs <- FPRs <- TPRs <- FNRs <- Precisions <- Recalls <- Accuracies <- Hms <- 
  FDRs.QTN12 <- FPRs.QTN12 <- TPRs.QTN12 <- Hms.QTN12 <-  
  FNRs.QTN12 <- Recalls.QTN12 <- Accuracies.QTN12 <- Precisions.QTN12 <- 
  FDRs.QTN3 <- FPRs.QTN3 <- TPRs.QTN3 <- Hms.QTN3 <-  
  FNRs.QTN3 <- Recalls.QTN3 <- Accuracies.QTN3 <- Precisions.QTN3 <- 
  array(NA, dim = c(n.iter, n.method, n.thres),
        dimnames = list(iter = 1:n.iter, method = method.names, thres = thres.names))
FDRs.other <- FPRs.other <- TPRs.other <- Hms.other <-  
  FNRs.other <- Recalls.other <- Accuracies.other <- Precisions.other <- 
  inflation.levels <- array(NA, dim = c(n.iter, n.method, n.inflator),
                            dimnames = list(iter = 1:n.iter, method = method.names, inflation = n.top.false))
AUCs <- AUC.relaxes <- AUC.causals <- log10p1s <- log10p2s <-
  log10p1.relaxes <- log10p2.relaxes <- matrix(rep(NA, n.iter * n.method), nrow = n.iter,
                                               dimnames = list(iter = 1:n.iter, method = method.names))
cor.causals <- rep(NA, n.iter)
for(i in 1:n.iter){
  print(paste0("iteration_", i))
  dir.iter <- paste0(dir.iter.0, i)
  
  ### 3.0.2.3. phenotype data ###
  file.gene <- paste0(dir.iter, "/causal_information.csv")
  causal.info <- data.frame(fread(file.gene))
  causals <- causal.info[, 1]
  causals.mat[i, ] <- causals
  freq.now <- apply(x[, causals], 2,  function(x) mean(x + 1) / 2)
  MAF <- pmin(freq.now, 1 - freq.now)
  MAFs[i, ] <- MAF
  
  cor.causals[i] <- cor(x[, causals[1]], x[, causals[2]])
  
  
  
  ##### 3.1. Summarize Normal GWAS by RAINBOW #####
  save.normal <- paste0(dir.iter, "/normal_")
  file.normal <- paste0(save.normal, "simple_result.RData")
  load(file = file.normal)
  
  SS.normal <- SS_gwas(res = normal.res, x = x, qtn.candidate = causals, method.thres = method.now,
                       map.x = map, gene.set = NULL, saveName = save.normal, n.top.false.block = n.top.false)
  AUC.causals[i, 1] <- AUC_causal(res = normal.res, x = x, map.x = map, qtn.candidate = causals,
                                  LD_length = 15000, cor.thres = 0.35, window.size = 0, gene.set = NULL)
  
  
  
  ##### 3.2. Summarize haplo-group-fixed (HGF) GWAS #####
  num.haps <- c(2, 3, 4)
  HGF.methods <- c("kmed", "phylo")
  
  
  SS.HGFs <- NULL
  AUC.causal.nows <- NULL
  for(num.hap in num.haps){
    for(HGF.method in HGF.methods){
      save.HGF <- paste0(dir.iter, "/HGF", num.hap, "_", HGF.method, "_")
      file.HGF <- paste0(save.HGF, "simple_result.RData")
      load(file = file.HGF)
      
      SS.HGF <- SS_gwas(res = HGF.res, x = x, qtn.candidate = causals, method.thres = method.now,
                        map.x = map, gene.set = block.list, saveName = save.HGF, n.top.false.block = n.top.false)
      SS.HGFs <- c(SS.HGFs, list(SS.HGF))
      
      AUC.causal.now <- AUC_causal(res = HGF.res, x = x, map.x = map, qtn.candidate = causals,
                                   LD_length = 15000, cor.thres = 0.35, window.size = 0, gene.set = gene.set)
      AUC.causal.nows <- c(AUC.causal.nows, AUC.causal.now)
    }
  }
  
  SS.HGF2.kmed <- SS.HGFs[[1]]
  SS.HGF2.phylo <- SS.HGFs[[2]]
  SS.HGF3.kmed <- SS.HGFs[[3]]
  SS.HGF3.phylo <- SS.HGFs[[4]]
  SS.HGF4.kmed <- SS.HGFs[[5]]
  SS.HGF4.phylo <- SS.HGFs[[6]]
  
  
  AUC.causals[i, 2:7] <- AUC.causal.nows
  
  
  ##### 3.3. Summarize SKAT #####
  save.SKAT <- paste0(dir.iter, "/SKAT_")
  file.SKAT <- paste0(save.SKAT, "simple_result.RData")
  load(file = file.SKAT)
  
  SS.SKAT <- SS_gwas(res = SKAT.res, x = x, qtn.candidate = causals, method.thres = method.now,
                     map.x = map, gene.set = block.list, saveName = save.SKAT, n.top.false.block = n.top.false)
  
  AUC.causals[i, 8] <- AUC_causal(res = SKAT.res, x = x, map.x = map, qtn.candidate = causals,
                                  LD_length = 15000, cor.thres = 0.35, window.size = 0, gene.set = gene.set)
  
  
  
  ##### 3.4. Summarize Haplotype-based multi-SNP GWAS by RAINBOW #####
  save.multisnp <- paste0(dir.iter, "/multisnp_")
  file.multisnp <- paste0(save.multisnp, "simple_result.RData")
  load(file = file.multisnp)
  
  SS.multisnp <- SS_gwas(res = multisnp.res, x = x, qtn.candidate = causals, method.thres = method.now,
                         map.x = map, gene.set = block.list, saveName = save.multisnp, n.top.false.block = n.top.false)
  
  AUC.causals[i, 9] <- AUC_causal(res = multisnp.res, x = x, map.x = map, qtn.candidate = causals,
                                  LD_length = 15000, cor.thres = 0.35, window.size = 0, gene.set = gene.set)
  
  SS.all <- c(list(normal = SS.normal), list(HGF2.kmed = SS.HGF2.kmed), list(HGF2.phylo = SS.HGF2.phylo),
              list(HGF3.kmed = SS.HGF3.kmed), list(HGF3.phylo = SS.HGF3.phylo),
              list(HGF4.kmed = SS.HGF4.kmed), list(HGF4.phylo = SS.HGF4.phylo),
              list(SKAT = SS.SKAT), list(multisnp = SS.multisnp))
  for(j in 1:n.method){
    SS.now <- SS.all[[j]]
    
    FDR <- SS.now$FDR
    FPR <- SS.now$FPR
    FNR <- SS.now$FNR
    Recall <- SS.now$Recall
    Precision <- SS.now$Precision
    Accuracy <- SS.now$Accuracy
    Hm <- SS.now$Hm
    if(nrow(SS.now$FDR.each) == 2){
    FDR.QTN12 <- SS.now$FDR.each[1, ]
    FPR.QTN12 <- SS.now$FPR.each[1, ]
    FNR.QTN12 <- SS.now$FNR.each[1, ]
    Recall.QTN12 <- SS.now$Recall.each[1, ]
    Precision.QTN12 <- SS.now$Precision.each[1, ]
    Accuracy.QTN12 <- SS.now$Accuracy.each[1, ]
    Hm.QTN12 <- SS.now$Hm.each[1, ]
    FDR.QTN3 <- SS.now$FDR.each[2, ]
    FPR.QTN3 <- SS.now$FPR.each[2, ]
    FNR.QTN3 <- SS.now$FNR.each[2, ]
    Recall.QTN3 <- SS.now$Recall.each[2, ]
    Precision.QTN3 <- SS.now$Precision.each[2, ]
    Accuracy.QTN3 <- SS.now$Accuracy.each[2, ]
    Hm.QTN3 <- SS.now$Hm.each[2, ]
    }else{
      FDR.QTN12 <- apply(SS.now$FDR.each[1:2, ], 2, mean, na.rm = TRUE)
      FPR.QTN12 <- apply(SS.now$FPR.each[1:2, ], 2, mean, na.rm = TRUE)
      FNR.QTN12 <- apply(SS.now$FNR.each[1:2, ], 2, mean, na.rm = TRUE)
      Recall.QTN12 <- apply(SS.now$Recall.each[1:2, ], 2, mean, na.rm = TRUE)
      Precision.QTN12 <- apply(SS.now$Precision.each[1:2, ], 2, mean, na.rm = TRUE)
      Accuracy.QTN12 <- apply(SS.now$Accuracy.each[1:2, ], 2, mean, na.rm = TRUE)
      Hm.QTN12 <- apply(SS.now$Hm.each[1:2, ], 2, mean, na.rm = TRUE)
      FDR.QTN3 <- SS.now$FDR.each[3, ]
      FPR.QTN3 <- SS.now$FPR.each[3, ]
      FNR.QTN3 <- SS.now$FNR.each[3, ]
      Recall.QTN3 <- SS.now$Recall.each[3, ]
      Precision.QTN3 <- SS.now$Precision.each[3, ]
      Accuracy.QTN3 <- SS.now$Accuracy.each[3, ]
      Hm.QTN3 <- SS.now$Hm.each[3, ]
    }
    FDR.other <- SS.now$FDR.other
    FPR.other <- SS.now$FPR.other
    FNR.other <- SS.now$FNR.other
    Recall.other <- SS.now$Recall.other
    Precision.other <- SS.now$Precision.other
    Accuracy.other <- SS.now$Accuracy.other
    Hm.other <- SS.now$Hm.other
    AUC <- SS.now$AUC
    AUC.relax <- SS.now$AUC.relax
    log10p1 <- SS.now$log.p[1]
    if(j == 1){
      log10p1 <- mean(SS.now$log.p[1:n.causal.site], na.rm = TRUE)
      log10p2 <- SS.now$log.p[n.causal.site + 1]
      log10p1.relax <- mean(SS.now$max.trues[1:n.causal.site], na.rm = TRUE)
      log10p2.relax <- SS.now$max.trues[n.causal.site + 1]
    }else{
      log10p1 <- SS.now$log.p[1]
      log10p2 <- SS.now$log.p[2]
      log10p1.relax <- SS.now$max.trues[1]
      log10p2.relax <- SS.now$max.trues[2]
    }
    inflation.level <- SS.now$mean.false
    FDRs[i, j, ] <- FDR
    FPRs[i, j, ] <- FPR
    FNRs[i, j, ] <- FNR
    Recalls[i, j, ] <- Recall
    Precisions[i, j, ] <- Precision
    Accuracies[i, j, ] <- Accuracy
    Hms[i, j, ] <- Hm
    FDRs.QTN12[i, j, ] <- FDR.QTN12
    FPRs.QTN12[i, j, ] <- FPR.QTN12
    FNRs.QTN12[i, j, ] <- FNR.QTN12
    Recalls.QTN12[i, j, ] <- Recall.QTN12
    Precisions.QTN12[i, j, ] <- Precision.QTN12
    Accuracies.QTN12[i, j, ] <- Accuracy.QTN12
    Hms.QTN12[i, j, ] <- Hm.QTN12
    FDRs.QTN3[i, j, ] <- FDR.QTN3
    FPRs.QTN3[i, j, ] <- FPR.QTN3
    FNRs.QTN3[i, j, ] <- FNR.QTN3
    Recalls.QTN3[i, j, ] <- Recall.QTN3
    Precisions.QTN3[i, j, ] <- Precision.QTN3
    Accuracies.QTN3[i, j, ] <- Accuracy.QTN3
    Hms.QTN3[i, j, ] <- Hm.QTN3
    FDRs.other[i, j, ] <- FDR.other
    FPRs.other[i, j, ] <- FPR.other
    FNRs.other[i, j, ] <- FNR.other
    Recalls.other[i, j, ] <- Recall.other
    Precisions.other[i, j, ] <- Precision.other
    Accuracies.other[i, j, ] <- Accuracy.other
    Hms.other[i, j, ] <- Hm.other
    AUCs[i, j] <- AUC
    AUC.relaxes[i, j] <- AUC.relax
    log10p1s[i, j] <- log10p1
    log10p2s[i, j] <- log10p2
    inflation.levels[i, j, ] <- inflation.level
    log10p1.relaxes[i, j] <- log10p1.relax
    log10p2.relaxes[i, j] <- log10p2.relax
  }
}



##### 3.5. See results quickly #####
apply(FDRs, c(2, 3), mean, na.rm = TRUE)
apply(FPRs, c(2, 3), mean, na.rm = TRUE)
apply(FNRs, c(2, 3), mean, na.rm = TRUE)
apply(Recalls, c(2, 3), mean, na.rm = TRUE)
apply(Precisions, c(2, 3), mean, na.rm = TRUE)
apply(Accuracies, c(2, 3), mean, na.rm = TRUE)
apply(Hms, c(2, 3), mean, na.rm = TRUE)
apply(FDRs.other, c(2, 3), mean, na.rm = TRUE)
apply(FPRs.other, c(2, 3), mean, na.rm = TRUE)
apply(FNRs.other, c(2, 3), mean, na.rm = TRUE)
apply(Recalls.other, c(2, 3), mean, na.rm = TRUE)
apply(Precisions.other, c(2, 3), mean, na.rm = TRUE)
apply(Accuracies.other, c(2, 3), mean, na.rm = TRUE)
apply(Hms.other, c(2, 3), mean, na.rm = TRUE)
apply(AUCs, 2, mean, na.rm = TRUE)
apply(AUC.relaxes, 2, mean, na.rm = TRUE)
apply(log10p1s, 2, mean, na.rm = TRUE)
apply(log10p2s, 2, mean, na.rm = TRUE)
apply(inflation.levels, c(2, 3), mean, na.rm = TRUE)
apply(log10p1.relaxes, 2, mean, na.rm = TRUE)
apply(log10p2.relaxes, 2, mean, na.rm = TRUE)

adjusted.log10p1s <- array(rep(log10p1s, 4), dim = c(n.iter, n.method, n.inflator)) - inflation.levels
adjusted.log10p2s <- array(rep(log10p2s, 4), dim = c(n.iter, n.method, n.inflator)) - inflation.levels

apply(adjusted.log10p1s, c(2, 3), mean, na.rm = TRUE)
apply(adjusted.log10p2s, c(2, 3), mean, na.rm = TRUE)

adjusted.log10p1.relaxes <- array(rep(log10p1.relaxes, 4), dim = c(n.iter, n.method, n.inflator)) - inflation.levels
adjusted.log10p2.relaxes <- array(rep(log10p2.relaxes, 4), dim = c(n.iter, n.method, n.inflator)) - inflation.levels
apply(adjusted.log10p1.relaxes, c(2, 3), mean, na.rm = TRUE)
apply(adjusted.log10p2.relaxes, c(2, 3), mean, na.rm = TRUE)

apply(AUC.causals, 2, mean, na.rm = TRUE)
apply(AUC.causals[apply(adjusted.log10p1.relaxes >= 2, 1, any), ], 2, mean, na.rm = T)




###### 4. Summarize results and draw boxplots ###### 
##### 4.1. Summarize results #####
for(thres.no in 1:n.thres){
  for(inflator.no in 1:n.inflator){
    thres.now <- thres.names[thres.no]
    n.top.false.now <- n.top.false[inflator.no]
    
    df.summary <- data.frame(method = method.fac.abb,
                             FDR = c(FDRs[, , thres.no]), FPR = c(FPRs[, , thres.no]),
                             FNR = c(FNRs[, , thres.no]), Recall = c(Recalls[, , thres.no]),
                             Precision = c(Precisions[, , thres.no]),
                             Accuracy = c(Accuracies[, , thres.no]), F_value = c(Hms[, , thres.no]),
                             FDR.QTN12 = c(FDRs.QTN12[, , thres.no]), FPR.QTN12 = c(FPRs.QTN12[, , thres.no]),
                             FNR.QTN12 = c(FNRs.QTN12[, , thres.no]), Recall.QTN12 = c(Recalls.QTN12[, , thres.no]),
                             Precision.QTN12 = c(Precisions.QTN12[, , thres.no]),
                             Accuracy.QTN12 = c(Accuracies.QTN12[, , thres.no]),
                             F_value.QTN12 = c(Hms.QTN12[, , thres.no]),
                             FDR.QTN3 = c(FDRs.QTN3[, , thres.no]), FPR.QTN3 = c(FPRs.QTN3[, , thres.no]),
                             FNR.QTN3 = c(FNRs.QTN3[, , thres.no]), Recall.QTN3 = c(Recalls.QTN3[, , thres.no]),
                             Precision.QTN3 = c(Precisions.QTN3[, , thres.no]),
                             Accuracy.QTN3 = c(Accuracies.QTN3[, , thres.no]),
                             F_value.QTN3 = c(Hms.QTN3[, , thres.no]),
                             FDR.other = c(FDRs.other[, , inflator.no]), FPR.other = c(FPRs.other[, , inflator.no]),
                             FNR.other = c(FNRs.other[, , inflator.no]), Recall.other = c(Recalls.other[, , inflator.no]),
                             Precision.other = c(Precisions.other[, , inflator.no]),
                             Accuracy.other = c(Accuracies.other[, , inflator.no]),
                             F_value.other = c(Hms.other[, , inflator.no]),
                             AUC = c(AUCs), AUC.causal = c(AUC.causals), AUC.relax = c(AUC.relaxes),
                             log10p.1 = c(log10p1s), log10p.2 = c(log10p2s),
                             log10p.1.relax = c(log10p1.relaxes),
                             log10p.2.relax = c(log10p2.relaxes),
                             Inflator = c(inflation.levels[, , inflator.no]),
                             Adjusted.log10p.1 = c(adjusted.log10p1s[, , inflator.no]),
                             Adjusted.log10p.2 = c(adjusted.log10p2s[, , inflator.no]),
                             Adjusted.log10p.1.relax = c(adjusted.log10p1.relaxes[, , inflator.no]),
                             Adjusted.log10p.2.relax = c(adjusted.log10p2.relaxes[, , inflator.no]))
    
    file.df.summary <- paste0(dir.base, "/Summary_statistics_of_", n.iter, "simulations_", thres.now,
                              "_#_of_top_false_", n.top.false.now, ".csv")
    write.csv(df.summary, file = file.df.summary)
    
    
    
    
    
  #   ##### 4.2. Draw boxplots #####
  #   y.range <- range(log10p1s, log10p1.relaxes)
  #   
  #   p.1 <- ggplot(data = df.summary, aes(x = method, y = log10p.1)) + 
  #     ylab(expression(-log[10](italic(p)))) +
  #     ylim(y.range[1], y.range[2]) +
  #     geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
  #     ggtitle(expression(paste(bold(-log[bold("10")](bolditalic(p))), bold(" in terms of haplotype block")))) +
  #     theme(plot.title = element_text(size = 22, hjust = 0.5, family = "Verdana"),
  #           axis.title = element_text(size = 18),
  #           axis.text.x = element_text(size = 20, face = "bold"),
  #           axis.text.y = element_text(size = 16))
  #   
  #   p.2 <- ggplot(data = df.summary, aes(x = method, y = log10p.1.relax)) + 
  #     ylab(expression(-log[10](italic(p)))) +
  #     ylim(y.range[1], y.range[2]) +
  #     geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
  #     ggtitle(expression(paste(bold(-log[bold("10")](bolditalic(p))), bold(" in terms of LD block")))) +
  #     theme(plot.title = element_text(size = 22, hjust = 0.5, family = "Verdana"),
  #           axis.title = element_text(size = 18),
  #           axis.text.x = element_text(size = 20, face = "bold"),
  #           axis.text.y = element_text(size = 16))
  #   
  #   y.range.adjusted <- range(adjusted.log10p1s, adjusted.log10p1.relaxes)
  #   
  #   p.3 <- ggplot(data = df.summary, aes(x = method, y = Adjusted.log10p.1)) + 
  #     ylab(expression(-log[10](italic(p)))) +
  #     ylim(y.range.adjusted[1], y.range.adjusted[2]) +
  #     geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
  #     ggtitle(expression(paste(bold("Adjusted "), bold(-log[bold("10")](bolditalic(p))), bold(" in terms of haplotype block")))) +
  #     theme(plot.title = element_text(size = 22, hjust = 0.5, family = "Verdana"),
  #           axis.title = element_text(size = 18),
  #           axis.text.x = element_text(size = 20, face = "bold"),
  #           axis.text.y = element_text(size = 16))
  #   
  #   
  #   p.4 <- ggplot(data = df.summary, aes(x = method, y = Adjusted.log10p.1.relax)) + 
  #     ylab(expression(-log[10](italic(p)))) +
  #     ylim(y.range.adjusted[1], y.range.adjusted[2]) +
  #     geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
  #     ggtitle(expression(paste(bold("Adjusted "), bold(-log[bold("10")](bolditalic(p))), bold(" in terms of LD block")))) +
  #     theme(plot.title = element_text(size = 22, hjust = 0.5, family = "Verdana"),
  #           axis.title = element_text(size = 18),
  #           axis.text.x = element_text(size = 20, face = "bold"),
  #           axis.text.y = element_text(size = 16))
  #   
  #   file.boxplot <- paste0(dir.res.base, "/Results_of_-log10(p)_boxplot_#_of_top_false_", n.top.false.now, ".png")
  #   png(file.boxplot, height = 1000, width = 1400)
  #   grid.arrange(p.1, p.2, p.3, p.4, ncol = 2)
  #   dev.off()
  #   
  #   
  #   
  #   
  #   
  #   ###### 5. Draw barplot for Recall, Precision and Hm ###### 
  #   method.mean.fac <- factor(rep(method.names.abb, 3), 
  #                             levels = method.levels.abb)
  #   SS.mean.fac <- factor(rep(c("Recall", "Precision", "F-measure"), each = n.method),
  #                         levels = c("Recall", "Precision", "F-measure"))
  #   SS.mean.values <- c(apply(Recalls[, , thres.no], 2, mean, na.rm = TRUE),
  #                       apply(Precisions[, , thres.no], 2, mean, na.rm = TRUE),
  #                       apply(Hms[, , thres.no], 2, mean, na.rm = TRUE))
  #   SS.sd.values <- c(apply(Recalls[, , thres.no], 2, sd, na.rm = TRUE),
  #                     apply(Precisions[, , thres.no], 2, sd, na.rm = TRUE),
  #                     apply(Hms[, , thres.no], 2, sd, na.rm = TRUE))
  #   df.summary.mean <- data.frame(method = method.mean.fac,
  #                                 summary_statistics = SS.mean.fac,
  #                                 values = SS.mean.values,
  #                                 sd = SS.sd.values)
  #   
  #   p <- ggplot(data = df.summary.mean, aes(x = method, y = values, fill = summary_statistics)) +
  #     geom_bar(stat = "identity", width = 0.6, position = position_dodge(0.65)) +
  #     scale_color_manual(values = c("red", "green", "blue")) + 
  #     ggtitle(paste0("Summary statistics of GWAS results")) +
  #     theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Verdana", face = "bold"),
  #           axis.title = element_text(size = 22),
  #           axis.text.x = element_text(size = 25, face = "bold"),
  #           axis.text.y = element_text(size = 19),
  #           legend.text = element_text(size = 23, face = "bold"),
  #           legend.title = element_text(size = 18),
  #           legend.key.height = unit(31, units = "pt"))
  #   # p <- p + geom_errorbar(aes(ymin = values - sd, ymax = values + sd, width = 0.25),
  #   #                          position = position_dodge(0.65)) 
  #   file.barplot <- paste0(dir.res.base, "/Results_of_summary_statistics_barplot_", thres.now, ".png")
  #   png(file.barplot, height = 700, width = 1100)
  #   print(p)
  #   dev.off()
  #   
  #   
  #   
  #   SS.mean.fac.other <- factor(rep(c("Recall", "Precision", "F-measure"), each = n.method),
  #                               levels = c("Recall", "Precision", "F-measure"))
  #   SS.mean.values.other <- c(apply(Recalls.other[, , inflator.no], 2, mean, na.rm = TRUE),
  #                             apply(Precisions.other[, , inflator.no], 2, mean, na.rm = TRUE),
  #                             apply(Hms.other[, , inflator.no], 2, mean, na.rm = TRUE))
  #   SS.sd.values.other <- c(apply(Recalls.other[, , inflator.no], 2, sd, na.rm = TRUE),
  #                           apply(Precisions.other[, , inflator.no], 2, sd, na.rm = TRUE),
  #                           apply(Hms.other[, , inflator.no], 2, sd, na.rm = TRUE))
  #   df.summary.mean.other <- data.frame(method = method.mean.fac,
  #                                       summary_statistics = SS.mean.fac.other,
  #                                       values = SS.mean.values.other,
  #                                       sd = SS.sd.values.other)
  #   
  #   p.other <- ggplot(data = df.summary.mean.other, aes(x = method, y = values, fill = summary_statistics)) +
  #     geom_bar(stat = "identity", width = 0.6, position = position_dodge(0.65)) +
  #     scale_color_manual(values = c("red", "green", "blue")) + 
  #     ggtitle(paste0("Summary statistics of GWAS results")) +
  #     theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Verdana", face = "bold"),
  #           axis.title = element_text(size = 22),
  #           axis.text.x = element_text(size = 25, face = "bold"),
  #           axis.text.y = element_text(size = 19),
  #           legend.text = element_text(size = 23, face = "bold"),
  #           legend.title = element_text(size = 18),
  #           legend.key.height = unit(31, units = "pt"))
  #   # p <- p + geom_errorbar(aes(ymin = values - sd, ymax = values + sd, width = 0.25),
  #   #                          position = position_dodge(0.65)) 
  #   file.barplot.other <- paste0(dir.res.base, "/Results_of_summary_statistics_barplot_inflator_threshold_",
  #                                "#_of_top_false_", n.top.false.now, ".png")
  #   png(file.barplot.other, height = 700, width = 1100)
  #   print(p.other)
  #   dev.off()
  #   
  #   
  #   
  #   
  #   
  #   
  #   ###### 6. Overwhelming results ###### 
  #   greatest.method.no.0 <- apply(adjusted.log10p1s[, , inflator.no], 1, which.max)
  #   greatest.method.no <- apply(adjusted.log10p1.relaxes[, , inflator.no], 1, which.max)
  #   greatest.method.count.0 <- table(greatest.method.no.0)
  #   greatest.method.count <- table(greatest.method.no)
  #   
  #   if(length(greatest.method.count.0) <= (n.method - 1)){
  #     missing.no <- which(!((1:n.method) %in% names(greatest.method.count.0)))
  #     greatest.method.count.0 <- c(greatest.method.count.0, rep(0, length(missing.no)))
  #     names(greatest.method.count.0) <-
  #       c(names(greatest.method.count.0)[1:(n.method - length(missing.no))], missing.no)
  #     greatest.method.count.0 <- greatest.method.count.0[order(names(greatest.method.count.0))]
  #   }
  #   
  #   if(length(greatest.method.count) <= (n.method - 1)){
  #     missing.no <- which(!((1:n.method) %in% names(greatest.method.count)))
  #     greatest.method.count <- c(greatest.method.count, rep(0, length(missing.no)))
  #     names(greatest.method.count) <-
  #       c(names(greatest.method.count)[1:(n.method - length(missing.no))], missing.no)
  #     greatest.method.count <- greatest.method.count[order(names(greatest.method.count))]
  #   }
  #   
  #   
  #   greatest.method.count.all <- rbind(greatest.method.count.0,
  #                                      greatest.method.count)
  #   rownames(greatest.method.count.all) <- c("Restrict", "Relax")
  #   
  #   file.greate.method <- paste0(dir.res.base, "/greatest_method_count_", 
  #                                "#_of_top_false_", n.top.false.now, ".csv")
  #   write.csv(greatest.method.count.all, file = file.greate.method)
  #   
  #   
  #   for(no.now in 1:n.method){
  #     adjusted.log10p1.relaxes[(greatest.method.no == no.now) & 
  #                                (adjusted.log10p1.relaxes[, no.now, inflator.no] >= 2) &
  #                                apply(adjusted.log10p1.relaxes[, - no.now, inflator.no] < 1, 1, all), , inflator.no, drop = FALSE]
  #     
  #     overwhelm.no <- which((greatest.method.no == no.now) & 
  #                             (adjusted.log10p1.relaxes[, no.now, inflator.no] >= 2) &
  #                             apply(adjusted.log10p1.relaxes[, - no.now, inflator.no] < 1, 1, all))
  #     
  #     if(length(overwhelm.no) >= 1){
  #       dir.overwhelm <- paste0(dir.res.base, "/overwhelming_results_", method.names.abb[no.now], 
  #                               "_#_of_top_false_", n.top.false.now)
  #       dir.create(dir.overwhelm)
  #       
  #       for(iter.no in overwhelm.no){
  #         file.overwhelms.0 <- list.files(paste0(dir.iter.0, iter.no, "/"))
  #         file.overwhelms <- paste0(dir.iter.0, iter.no, "/", file.overwhelms.0)
  #         
  #         dir.overwhelm.now <- paste0(dir.overwhelm, "/iteration_", iter.no)
  #         dir.create(dir.overwhelm.now)
  #         file.overwhelm.news <- paste0(dir.overwhelm.now, "/", file.overwhelms.0)
  #         
  #         file.copy(file.overwhelms, file.overwhelm.news, overwrite = TRUE)
  #       }
  #     }
  #   }
  #   
  #   
  #   ######  7. Whether each method can detect causal haplotype block itself ###### 
  #   method.fac.causal <- factor(rep(method.names.abb, 2), 
  #                               levels = method.levels.abb)
  #   AUC.causal.mean.fac <- factor(rep(c("all", "detected"), each = n.method),
  #                                 levels = c("all", "detected"))
  #   AUC.causal.mean.values <- c(apply(AUC.causals, 2, mean, na.rm = TRUE),
  #                               apply(AUC.causals[apply(adjusted.log10p1.relaxes[, , inflator.no] >= 2, 1, any), ], 2, mean, na.rm = T))
  #   AUC.causal.sd.values <- c(apply(AUC.causals, 2, sd, na.rm = TRUE),
  #                             apply(AUC.causals[apply(adjusted.log10p1.relaxes[, , inflator.no] >= 2, 1, any), ], 2, sd, na.rm = T))
  #   df.causal.mean <- data.frame(method = method.fac.causal,
  #                                AUCs = AUC.causal.mean.fac,
  #                                values = AUC.causal.mean.values,
  #                                sd = AUC.causal.sd.values)
  #   
  #   p.causal <- ggplot(data = df.causal.mean, aes(x = method, y = values, fill = AUCs)) +
  #     geom_bar(stat = "identity", width = 0.6, position = position_dodge(0.65)) +
  #     scale_color_manual(values = c("red", "green")) + 
  #     ggtitle(paste0("AUC around causal haplotype block")) +
  #     theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Verdana", face = "bold"),
  #           axis.title = element_text(size = 22),
  #           axis.text.x = element_text(size = 25, face = "bold"),
  #           axis.text.y = element_text(size = 19),
  #           legend.text = element_text(size = 23, face = "bold"),
  #           legend.title = element_text(size = 18),
  #           legend.key.height = unit(31, units = "pt"))
  #   # p <- p + geom_errorbar(aes(ymin = values - sd, ymax = values + sd, width = 0.25),
  #   #                          position = position_dodge(0.65)) 
  #   file.barplot <- paste0(dir.res.base, "/Results_of_AUC_causals_barplot_",
  #                          "#_of_top_false_", n.top.false.now, ".png")
  #   png(file.barplot, height = 700, width = 1100)
  #   print(p.causal)
  #   dev.off()
  #   
  }
}













###### 8. Relationship with MAF ###### 
# order.min.MAF <- order(apply(MAFs, 1, function(x) mean(x[1:2])))
# 
# FDRs.res.ord <- cbind(order.min.MAF, cor.causals[order.min.MAF],
#                           MAFs[order.min.MAF, 1:2],
#                           round(FDRs[order.min.MAF, ], 3))
# rownames(FDRs.res.ord) <- paste0("iteration_", order.min.MAF)
# colnames(FDRs.res.ord) <- c("iter.no", "cor.causals", "MAF1", "MAF2",
#                                 "Single-SNP", "HGF2", "HGF3", "SKAT", "multisnp")
# 
# 
# Recalls.res.ord <- cbind(order.min.MAF, cor.causals[order.min.MAF],
#                       MAFs[order.min.MAF, 1:2],
#                       round(Recalls[order.min.MAF, ], 3))
# rownames(Recalls.res.ord) <- paste0("iteration_", order.min.MAF)
# colnames(Recalls.res.ord) <- c("iter.no", "cor.causals", "MAF1", "MAF2",
#                             "Single-SNP", "HGF2", "HGF3", "SKAT", "multisnp")






# log10p1s.res.ord <- cbind(order.min.MAF, cor.causals[order.min.MAF],
#                           MAFs[order.min.MAF, 1:2],
#                           round(log10p1s[order.min.MAF, ], 3))
# rownames(log10p1s.res.ord) <- paste0("iteration_", order.min.MAF)
# colnames(log10p1s.res.ord) <- c("iter.no", "cor.causals", "MAF1", "MAF2",
#                                 "Single-SNP", "HGF2", "HGF3", "SKAT", "multisnp")
# 
# adjusted.log10p1s.res.ord <- cbind(order.min.MAF, cor.causals[order.min.MAF],
#                           MAFs[order.min.MAF, 1:2],
#                           round(adjusted.log10p1s[order.min.MAF, ], 3))
# rownames(adjusted.log10p1s.res.ord) <- paste0("iteration_", order.min.MAF)
# colnames(adjusted.log10p1s.res.ord) <- c("iter.no", "cor.causals", "MAF1", "MAF2",
#                                 "Single-SNP", "HGF2", "HGF3", "SKAT", "multisnp")
# adjusted.log10p1s.res.ord
# 
# apply(adjusted.log10p1s.res.ord[, 5:9] >= 3, 2, sum)
# 
# 
# #####
# log10p1.relaxes.res.ord <- cbind(order.min.MAF, cor.causals[order.min.MAF],
#                           MAFs[order.min.MAF, 1:2],
#                           round(log10p1.relaxes[order.min.MAF, ], 3))
# rownames(log10p1.relaxes.res.ord) <- paste0("iteration_", order.min.MAF)
# colnames(log10p1.relaxes.res.ord) <- c("iter.no", "cor.causals", "MAF1", "MAF2",
#                                 "Single-SNP", "HGF2", "HGF3", "SKAT", "multisnp")
# 
# adjusted.log10p1.relaxes.res.ord <- cbind(order.min.MAF, cor.causals[order.min.MAF],
#                                    MAFs[order.min.MAF, 1:2],
#                                    round(adjusted.log10p1.relaxes[order.min.MAF, ], 3))
# rownames(adjusted.log10p1.relaxes.res.ord) <- paste0("iteration_", order.min.MAF)
# colnames(adjusted.log10p1.relaxes.res.ord) <- c("iter.no", "cor.causals", "MAF1", "MAF2",
#                                          "Single-SNP", "HGF2", "HGF3", "SKAT", "multisnp")
# 
# apply(adjusted.log10p1.relaxes.res.ord[, 9] - adjusted.log10p1.relaxes.res.ord[, 5:8], 1, sum)








save(inflation.levels, file = paste0(dir.base, "/inflators.RData"))





