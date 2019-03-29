###########################################################################################
###### 3.3_Rice_HGRAINBOW_Summary_plot_of_HGRAINBOW_Haplotype_based_GWAS_for_RAINBOW ######
###########################################################################################

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

manhattan.no <- c(9, 34, 39, 44)

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
}




###### 4. Summarize results and draw boxplots ###### 
##### 4.1. Summarize results #####
direction.effs <- c("plus", "minus")
ranges <- ranges.adjusted <- ranges.2 <- ranges.adjusted.2 <- NULL
for(eff.now in direction.effs){
  for(thres.no in 1:n.thres){
    for(inflator.no in 1:n.inflator){
      thres.now <- thres.names[thres.no]
      n.top.false.now <- n.top.false[inflator.no]
      dir.base.now <- paste0("midstream/number_of_causals=", n.causal.site,
                             "_direction_of_effect=", eff.now,
                             "_trial_no=", trial.no)
      file.df.summary <- paste0(dir.base.now, "/Summary_statistics_of_", n.iter, "simulations_", thres.now,
                                "_#_of_top_false_", n.top.false.now, ".csv")
      df.summary <- read.csv(file = file.df.summary, row.names = 1)
      range.now <- range(df.summary$log10p.1, df.summary$log10p.1.relax)
      range.adjusted.now <- range(df.summary$Adjusted.log10p.1, df.summary$Adjusted.log10p.1.relax)
      
      range.now.2 <- range(df.summary$log10p.2, df.summary$log10p.2.relax)
      range.adjusted.now.2 <- range(df.summary$Adjusted.log10p.2, df.summary$Adjusted.log10p.2.relax)
      
      ranges <- c(ranges, range.now)
      ranges.adjusted <- c(ranges.adjusted, range.adjusted.now)
      
      ranges.2 <- c(ranges.2, range.now.2)
      ranges.adjusted.2 <- c(ranges.adjusted.2, range.adjusted.now.2)
    }
  }
}

y.range.plus <- range(ranges)
y.range.adjusted.plus <- range(ranges.adjusted)
y.range.plus.2 <- range(ranges.2)
y.range.adjusted.plus.2 <- range(ranges.adjusted.2)

inflator.thres.base <- 1.5
inflator.thres.diff <- 0.3
load(paste0(dir.base, "/inflators.RData"))
inflator.thress <- t(apply(inflation.levels, c(3, 2), mean)) + 
  inflator.thres.base - apply(inflation.levels, c(3, 2), mean)[2, ]

for(thres.no in 1:n.thres){
  for(inflator.no in 1:n.inflator){
    thres.now <- thres.names[thres.no]
    n.top.false.now <- n.top.false[inflator.no]
    
    file.df.summary <- paste0(dir.base, "/Summary_statistics_of_", n.iter, "simulations_", thres.now,
                              "_#_of_top_false_", n.top.false.now, ".csv")
    df.summary <- read.csv(file = file.df.summary, row.names = 1)
    df.summary[, 1] <- method.fac.abb
    
    log10p1s.now <- matrix(df.summary$log10p.1, nrow = n.iter)
    log10p1.relaxes.now <- matrix(df.summary$log10p.1.relax, nrow = n.iter)
    adjusted.log10p1s.now <- matrix(df.summary$Adjusted.log10p.1, nrow = n.iter)
    adjusted.log10p1.relaxes.now <- matrix(df.summary$Adjusted.log10p.1.relax, nrow = n.iter)
    
    log10p2s.now <- matrix(df.summary$log10p.2, nrow = n.iter)
    log10p2.relaxes.now <- matrix(df.summary$log10p.2.relax, nrow = n.iter)
    adjusted.log10p2s.now <- matrix(df.summary$Adjusted.log10p.2, nrow = n.iter)
    adjusted.log10p2.relaxes.now <- matrix(df.summary$Adjusted.log10p.2.relax, nrow = n.iter)
    
    Recalls.now <- matrix(df.summary$Recall, nrow = n.iter)
    Precisions.now <- matrix(df.summary$Precision, nrow = n.iter)
    Hms.now <- matrix(df.summary$F_value, nrow = n.iter)
    Recalls.QTN12.now <- matrix(df.summary$Recall.QTN12, nrow = n.iter)
    Precisions.QTN12.now <- matrix(df.summary$Precision.QTN12, nrow = n.iter)
    Hms.QTN12.now <- matrix(df.summary$F_value.QTN12, nrow = n.iter)
    Recalls.QTN3.now <- matrix(df.summary$Recall.QTN3, nrow = n.iter)
    Precisions.QTN3.now <- matrix(df.summary$Precision.QTN3, nrow = n.iter)
    Hms.QTN3.now <- matrix(df.summary$F_value.QTN3, nrow = n.iter)
    Recalls.other.now <- matrix(df.summary$Recall.other, nrow = n.iter)
    Precisions.other.now <- matrix(df.summary$Precision.other, nrow = n.iter)
    Hms.other.now <- matrix(df.summary$F_value.other, nrow = n.iter)
    AUC.causals <- matrix(df.summary$AUC.causal, nrow = n.iter)
    
    ##### 4.2.1. Draw boxplots for QTN1 with labels #####
    p.1 <- ggplot(data = df.summary, aes(x = method, y = log10p.1)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.plus[1], y.range.plus[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold(-log[bold("10")](bolditalic(p))), bold(" in terms of haplotype block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 30),
            axis.text.x = element_text(size = 22, face = "bold"),
            axis.text.y = element_text(size = 29))
    
    p.2 <- ggplot(data = df.summary, aes(x = method, y = log10p.1.relax)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.plus[1], y.range.plus[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold(-log[bold("10")](bolditalic(p))), bold(" in terms of LD block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 30),
            axis.text.x = element_text(size = 22, face = "bold"),
            axis.text.y = element_text(size = 29))
    
    
    p.3 <- ggplot(data = df.summary, aes(x = method, y = Adjusted.log10p.1)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.adjusted.plus[1], y.range.adjusted.plus[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold("Adjusted "), bold(-log[bold("10")](bolditalic(p))), bold(" in terms of haplotype block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 30),
            axis.text.x = element_text(size = 22, face = "bold"),
            axis.text.y = element_text(size = 29))
    
    
    p.4 <- ggplot(data = df.summary, aes(x = method, y = Adjusted.log10p.1.relax)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.adjusted.plus[1], y.range.adjusted.plus[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold("Adjusted "), bold(-log[bold("10")](bolditalic(p))), bold(" in terms of LD block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 30),
            axis.text.x = element_text(size = 22, face = "bold"),
            axis.text.y = element_text(size = 29))
    
    file.boxplot <- paste0(dir.res.base, "/Results_of_-log10(p)_QTN12_boxplot_#_of_top_false_", n.top.false.now, "_with_label.png")
    png(file.boxplot, height = 1000, width = 1300)
    grid.arrange(p.1, p.2, p.3, p.4, ncol = 2)
    dev.off()
    
    
    
    ##### 4.2.2 Draw boxplots of QTN3 with label #####
    p.1 <- ggplot(data = df.summary, aes(x = method, y = log10p.2)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.plus.2[1], y.range.plus.2[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold(-log[bold("10")](bolditalic(p))), bold(" in terms of haplotype block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 30),
            axis.text.x = element_text(size = 22, face = "bold"),
            axis.text.y = element_text(size = 29))
    
    p.2 <- ggplot(data = df.summary, aes(x = method, y = log10p.2.relax)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.plus.2[1], y.range.plus.2[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold(-log[bold("10")](bolditalic(p))), bold(" in terms of LD block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 30),
            axis.text.x = element_text(size = 22, face = "bold"),
            axis.text.y = element_text(size = 29))
    
    
    p.3 <- ggplot(data = df.summary, aes(x = method, y = Adjusted.log10p.2)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.adjusted.plus.2[1], y.range.adjusted.plus.2[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold("Adjusted "), bold(-log[bold("10")](bolditalic(p))), bold(" in terms of haplotype block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 30),
            axis.text.x = element_text(size = 22, face = "bold"),
            axis.text.y = element_text(size = 29))
    
    
    p.4 <- ggplot(data = df.summary, aes(x = method, y = Adjusted.log10p.2.relax)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.adjusted.plus.2[1], y.range.adjusted.plus.2[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold("Adjusted "), bold(-log[bold("10")](bolditalic(p))), bold(" in terms of LD block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 30),
            axis.text.x = element_text(size = 22, face = "bold"),
            axis.text.y = element_text(size = 29))
    
    file.boxplot <- paste0(dir.res.base, "/Results_of_-log10(p)_QTN3_boxplot_#_of_top_false_", n.top.false.now, "_with_label.png")
    png(file.boxplot, height = 1000, width = 1300)
    grid.arrange(p.1, p.2, p.3, p.4, ncol = 2)
    dev.off()
    
    
    
    
    ##### 4.2.3. Draw boxplots of QTN1 with no label #####
    p.1 <- ggplot(data = df.summary, aes(x = method, y = log10p.1)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.plus[1], y.range.plus[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold(-log[bold("10")](bolditalic(p))), bold(" in terms of haplotype block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 23, face = "bold"),
            axis.text.y = element_blank())
    
    p.2 <- ggplot(data = df.summary, aes(x = method, y = log10p.1.relax)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.plus[1], y.range.plus[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold(-log[bold("10")](bolditalic(p))), bold(" in terms of LD block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 23, face = "bold"),
            axis.text.y = element_blank())
    
    
    p.3 <- ggplot(data = df.summary, aes(x = method, y = Adjusted.log10p.1)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.adjusted.plus[1], y.range.adjusted.plus[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold("Adjusted "), bold(-log[bold("10")](bolditalic(p))), bold(" in terms of haplotype block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 23, face = "bold"),
            axis.text.y = element_blank())
    
    
    p.4 <- ggplot(data = df.summary, aes(x = method, y = Adjusted.log10p.1.relax)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.adjusted.plus[1], y.range.adjusted.plus[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold("Adjusted "), bold(-log[bold("10")](bolditalic(p))), bold(" in terms of LD block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 23, face = "bold"),
            axis.text.y = element_blank())
    
    file.boxplot <- paste0(dir.res.base, "/Results_of_-log10(p)_QTN12_boxplot_#_of_top_false_", n.top.false.now, "_with_no_label.png")
    png(file.boxplot, height = 1000, width = 1100)
    grid.arrange(p.1, p.2, p.3, p.4, ncol = 2)
    dev.off()
    
    
    
    ##### 4.2.4 Draw boxplots of QTN3 with no label #####
    p.1 <- ggplot(data = df.summary, aes(x = method, y = log10p.2)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.plus.2[1], y.range.plus.2[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold(-log[bold("10")](bolditalic(p))), bold(" in terms of haplotype block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 23, face = "bold"),
            axis.text.y = element_blank())
    
    p.2 <- ggplot(data = df.summary, aes(x = method, y = log10p.2.relax)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.plus.2[1], y.range.plus.2[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold(-log[bold("10")](bolditalic(p))), bold(" in terms of LD block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 23, face = "bold"),
            axis.text.y = element_blank())
    
    
    p.3 <- ggplot(data = df.summary, aes(x = method, y = Adjusted.log10p.2)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.adjusted.plus.2[1], y.range.adjusted.plus.2[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold("Adjusted "), bold(-log[bold("10")](bolditalic(p))), bold(" in terms of haplotype block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 23, face = "bold"),
            axis.text.y = element_blank())
    
    
    p.4 <- ggplot(data = df.summary, aes(x = method, y = Adjusted.log10p.2.relax)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.adjusted.plus.2[1], y.range.adjusted.plus.2[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold("Adjusted "), bold(-log[bold("10")](bolditalic(p))), bold(" in terms of LD block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 23, face = "bold"),
            axis.text.y = element_blank())
    
    file.boxplot <- paste0(dir.res.base, "/Results_of_-log10(p)_QTN3_boxplot_#_of_top_false_", n.top.false.now, "_with_no_label.png")
    png(file.boxplot, height = 1000, width = 1100)
    grid.arrange(p.1, p.2, p.3, p.4, ncol = 2)
    dev.off()
    
    
    
    ###### 5. Draw barplot for Recall, Precision and Hm ###### 
    method.mean.fac <- factor(rep(method.names.abb, 3), 
                              levels = method.levels.abb)
    SS.mean.fac <- factor(rep(c("Recall", "Precision", "F-measure"), each = n.method),
                          levels = c("Recall", "Precision", "F-measure"))
    SS.mean.values <- c(apply(Recalls.now, 2, mean, na.rm = TRUE),
                        apply(Precisions.now, 2, mean, na.rm = TRUE),
                        apply(Hms.now, 2, mean, na.rm = TRUE))
    SS.sd.values <- c(apply(Recalls.now, 2, sd, na.rm = TRUE),
                      apply(Precisions.now, 2, sd, na.rm = TRUE),
                      apply(Hms.now, 2, sd, na.rm = TRUE))
    df.summary.mean <- data.frame(method = method.mean.fac,
                                  summary_statistics = SS.mean.fac,
                                  values = SS.mean.values,
                                  sd = SS.sd.values)
    
    p <- ggplot(data = df.summary.mean, aes(x = method, y = values, fill = summary_statistics)) +
      geom_bar(stat = "identity", width = 0.6, position = position_dodge(0.65)) +
      ylim(0, 1) +
      scale_color_manual(values = c("red", "green", "blue")) + 
      ggtitle(paste0("Summary statistics of GWAS results")) +
      theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Verdana", face = "bold"),
            axis.title = element_text(size = 22),
            axis.text.x = element_text(size = 25, face = "bold"),
            axis.text.y = element_text(size = 26),
            legend.text = element_text(size = 23, face = "bold"),
            legend.title = element_text(size = 18),
            legend.key.height = unit(31, units = "pt"))
    # p <- p + geom_errorbar(aes(ymin = values - sd, ymax = values + sd, width = 0.25),
    #                          position = position_dodge(0.65)) 
    file.barplot <- paste0(dir.res.base, "/Results_of_summary_statistics_barplot_", thres.now, ".png")
    png(file.barplot, height = 700, width = 800)
    print(p)
    dev.off()
    
    
    
    SS.mean.fac.QTN12 <- factor(rep(c("Recall", "Precision", "F-measure"), each = n.method),
                                levels = c("Recall", "Precision", "F-measure"))
    SS.mean.values.QTN12 <- c(apply(Recalls.QTN12.now, 2, mean, na.rm = TRUE),
                              apply(Precisions.QTN12.now, 2, mean, na.rm = TRUE),
                              apply(Hms.QTN12.now, 2, mean, na.rm = TRUE))
    SS.sd.values.QTN12 <- c(apply(Recalls.QTN12.now, 2, sd, na.rm = TRUE),
                            apply(Precisions.QTN12.now, 2, sd, na.rm = TRUE),
                            apply(Hms.QTN12.now, 2, sd, na.rm = TRUE))
    df.summary.mean.QTN12 <- data.frame(method = method.mean.fac,
                                        summary_statistics = SS.mean.fac.QTN12,
                                        values = SS.mean.values.QTN12,
                                        sd = SS.sd.values.QTN12)
    
    p.QTN12 <- ggplot(data = df.summary.mean.QTN12, aes(x = method, y = values, fill = summary_statistics)) +
      geom_bar(stat = "identity", width = 0.6, position = position_dodge(0.65)) +
      ylim(0, 1) +
      scale_color_manual(values = c("red", "green", "blue")) + 
      ggtitle(paste0("Summary statistics of GWAS results")) +
      theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Verdana", face = "bold"),
            axis.title = element_text(size = 22),
            axis.text.x = element_text(size = 25, face = "bold"),
            axis.text.y = element_text(size = 26),
            legend.text = element_text(size = 23, face = "bold"),
            legend.title = element_text(size = 18),
            legend.key.height = unit(31, units = "pt"))
    # p.QTN12 <- p.QTN12 + geom_errorbar(aes(ymin = values - sd, ymax = values + sd, width = 0.25),
    #                          position = position_dodge(0.65)) 
    file.barplot.QTN12 <- paste0(dir.res.base, "/Results_of_summary_statistics_barplot_", thres.now, "_QTN12.png")
    png(file.barplot.QTN12, height = 700, width = 800)
    print(p.QTN12)
    dev.off()
    
    
    
    SS.mean.fac.QTN3 <- factor(rep(c("Recall", "Precision", "F-measure"), each = n.method),
                               levels = c("Recall", "Precision", "F-measure"))
    SS.mean.values.QTN3 <- c(apply(Recalls.QTN3.now, 2, mean, na.rm = TRUE),
                             apply(Precisions.QTN3.now, 2, mean, na.rm = TRUE),
                             apply(Hms.QTN3.now, 2, mean, na.rm = TRUE))
    SS.sd.values.QTN3 <- c(apply(Recalls.QTN3.now, 2, sd, na.rm = TRUE),
                           apply(Precisions.QTN3.now, 2, sd, na.rm = TRUE),
                           apply(Hms.QTN3.now, 2, sd, na.rm = TRUE))
    df.summary.mean.QTN3 <- data.frame(method = method.mean.fac,
                                       summary_statistics = SS.mean.fac.QTN3,
                                       values = SS.mean.values.QTN3,
                                       sd = SS.sd.values.QTN3)
    
    p.QTN3 <- ggplot(data = df.summary.mean.QTN3, aes(x = method, y = values, fill = summary_statistics)) +
      geom_bar(stat = "identity", width = 0.6, position = position_dodge(0.65)) +
      ylim(0, 1) +
      scale_color_manual(values = c("red", "green", "blue")) + 
      ggtitle(paste0("Summary statistics of GWAS results")) +
      theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Verdana", face = "bold"),
            axis.title = element_text(size = 22),
            axis.text.x = element_text(size = 25, face = "bold"),
            axis.text.y = element_text(size = 26),
            legend.text = element_text(size = 23, face = "bold"),
            legend.title = element_text(size = 18),
            legend.key.height = unit(31, units = "pt"))
    # p.QTN3 <- p.QTN3 + geom_errorbar(aes(ymin = values - sd, ymax = values + sd, width = 0.25),
    #                          position = position_dodge(0.65)) 
    file.barplot.QTN3 <- paste0(dir.res.base, "/Results_of_summary_statistics_barplot_", thres.now, "_QTN3.png")
    png(file.barplot.QTN3, height = 700, width = 800)
    print(p.QTN3)
    dev.off()
    
    
    
    
    SS.mean.fac.other <- factor(rep(c("Recall", "Precision", "F-measure"), each = n.method),
                                levels = c("Recall", "Precision", "F-measure"))
    SS.mean.values.other <- c(apply(Recalls.other.now, 2, mean, na.rm = TRUE),
                              apply(Precisions.other.now, 2, mean, na.rm = TRUE),
                              apply(Hms.other.now, 2, mean, na.rm = TRUE))
    SS.sd.values.other <- c(apply(Recalls.other.now, 2, sd, na.rm = TRUE),
                            apply(Precisions.other.now, 2, sd, na.rm = TRUE),
                            apply(Hms.other.now, 2, sd, na.rm = TRUE))
    df.summary.mean.other <- data.frame(method = method.mean.fac,
                                        summary_statistics = SS.mean.fac.other,
                                        values = SS.mean.values.other,
                                        sd = SS.sd.values.other)
    
    p.other <- ggplot(data = df.summary.mean.other, aes(x = method, y = values, fill = summary_statistics)) +
      geom_bar(stat = "identity", width = 0.6, position = position_dodge(0.65)) +
      ylim(0, 1) +
      scale_color_manual(values = c("red", "green", "blue")) + 
      ggtitle(paste0("Summary statistics of GWAS results")) +
      theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Verdana", face = "bold"),
            axis.title = element_text(size = 22),
            axis.text.x = element_text(size = 25, face = "bold"),
            axis.text.y = element_text(size = 26),
            legend.text = element_text(size = 23, face = "bold"),
            legend.title = element_text(size = 18),
            legend.key.height = unit(31, units = "pt"))
    # p.other <- p.other + geom_errorbar(aes(ymin = values - sd, ymax = values + sd, width = 0.25),
    #                          position = position_dodge(0.65)) 
    file.barplot.other <- paste0(dir.res.base, "/Results_of_summary_statistics_barplot_inflator_threshold_",
                                 "#_of_top_false_", n.top.false.now, ".png")
    png(file.barplot.other, height = 700, width = 800)
    print(p.other)
    dev.off()
    
    
    
    
    
    
    ###### 6. Overwhelming results ###### 
    greatest.method.no.0 <- apply(adjusted.log10p1s.now, 1, which.max)
    greatest.method.no <- apply(adjusted.log10p1.relaxes.now, 1, which.max)
    greatest.method.count.0 <- table(greatest.method.no.0)
    greatest.method.count <- table(greatest.method.no)
    
    if(length(greatest.method.count.0) <= (n.method - 1)){
      missing.no <- which(!((1:n.method) %in% names(greatest.method.count.0)))
      greatest.method.count.0 <- c(greatest.method.count.0, rep(0, length(missing.no)))
      names(greatest.method.count.0) <-
        c(names(greatest.method.count.0)[1:(n.method - length(missing.no))], missing.no)
      greatest.method.count.0 <- greatest.method.count.0[order(names(greatest.method.count.0))]
    }
    
    if(length(greatest.method.count) <= (n.method - 1)){
      missing.no <- which(!((1:n.method) %in% names(greatest.method.count)))
      greatest.method.count <- c(greatest.method.count, rep(0, length(missing.no)))
      names(greatest.method.count) <-
        c(names(greatest.method.count)[1:(n.method - length(missing.no))], missing.no)
      greatest.method.count <- greatest.method.count[order(names(greatest.method.count))]
    }
    
    
    greatest.method.count.all <- rbind(greatest.method.count.0,
                                       greatest.method.count)
    rownames(greatest.method.count.all) <- c("Restrict", "Relax")
    
    file.greate.method <- paste0(dir.res.base, "/greatest_method_count_", 
                                 "#_of_top_false_", n.top.false.now, ".csv")
    write.csv(greatest.method.count.all, file = file.greate.method)
    
    
    for(no.now in 1:n.method){
      adjusted.log10p1.relaxes.now[(greatest.method.no == no.now) & (Recalls.QTN12.now[, no.now] > 0) &
                                     (adjusted.log10p1.relaxes.now[, no.now] >= inflator.thress[no.now, inflator.no]) &
                                     apply(adjusted.log10p1.relaxes.now[, - no.now] < (inflator.thress[no.now, inflator.no] - inflator.thres.diff), 1, all), , drop = FALSE]
      
      overwhelm.no <- which((greatest.method.no == no.now) & (Recalls.QTN12.now[, no.now] > 0) &
                              (adjusted.log10p1.relaxes.now[, no.now] >= inflator.thress[no.now, inflator.no]) &
                              apply(adjusted.log10p1.relaxes.now[, - no.now] < (inflator.thress[no.now, inflator.no] - inflator.thres.diff), 1, all))
      
      if(length(overwhelm.no) >= 1){
        dir.overwhelm <- paste0(dir.res.base, "/overwhelming_results_", method.names.abb[no.now], "_", thres.now,
                                "_#_of_top_false_", n.top.false.now)
        dir.create(dir.overwhelm)
        
        dir.overwhelm.manhattan <- paste0(dir.res.base, "/overwhelming_results_", method.names.abb[no.now], "_", thres.now,
                                          "_#_of_top_false_", n.top.false.now, "_manhattan")
        dir.create(dir.overwhelm.manhattan)
        
        for(iter.no in overwhelm.no){
          file.overwhelms.0 <- list.files(paste0(dir.iter.0, iter.no, "/"))
          file.overwhelms.manhattan.0 <- file.overwhelms.0[manhattan.no]
          file.overwhelms <- paste0(dir.iter.0, iter.no, "/", file.overwhelms.0)
          file.overwhelms.manhattan <- file.overwhelms[manhattan.no]
          
          dir.overwhelm.now <- paste0(dir.overwhelm, "/iteration_", iter.no)
          dir.create(dir.overwhelm.now)
          file.overwhelm.news <- paste0(dir.overwhelm.now, "/", file.overwhelms.0)
          
          dir.overwhelm.manhattan.now <- paste0(dir.overwhelm.manhattan, "/iteration_", iter.no)
          dir.create(dir.overwhelm.manhattan.now)
          file.overwhelm.manhattan.news <- paste0(dir.overwhelm.manhattan.now, "/", file.overwhelms.manhattan.0)
          
          file.copy(file.overwhelms, file.overwhelm.news, overwrite = TRUE)
          file.copy(file.overwhelms.manhattan, file.overwhelm.manhattan.news, overwrite = TRUE)
        }
      }
    }
    
    ######  7. Whether each method can detect causal haplotype block itself ###### 
    method.fac.causal <- factor(rep(method.names.abb, 2), 
                                levels = method.levels.abb)
    AUC.causal.mean.fac <- factor(rep(c("all", "detected"), each = n.method),
                                  levels = c("all", "detected"))
    AUC.detected.mean <- AUC.detected.sd <- rep(NA, n.method)
    
    for(method.no in 1:n.method){
      AUC.causals.detected <- AUC.causals[adjusted.log10p1.relaxes.now[, method.no] >= inflator.thress[method.no, inflator.no], method.no]
      AUC.detected.mean[method.no] <- mean(AUC.causals.detected, na.rm = TRUE)
      AUC.detected.sd[method.no] <- sd(AUC.causals.detected, na.rm = TRUE)
    }
    
    AUC.causal.mean.values <- c(apply(AUC.causals, 2, mean, na.rm = TRUE), AUC.detected.mean)
    AUC.causal.sd.values <- c(apply(AUC.causals, 2, sd, na.rm = TRUE), AUC.detected.sd)
    # AUC.causal.mean.values <- c(apply(AUC.causals, 2, mean, na.rm = TRUE),
    #                             apply(AUC.causals[apply(adjusted.log10p1.relaxes.now >= inflator.thress[thres.no], 1, any), ], 2, mean, na.rm = T))
    # AUC.causal.sd.values <- c(apply(AUC.causals, 2, sd, na.rm = TRUE),
    #                           apply(AUC.causals[apply(adjusted.log10p1.relaxes.now >= inflator.thress[thres.no], 1, any), ], 2, sd, na.rm = T))
    df.causal.mean <- data.frame(method = method.fac.causal,
                                 AUCs = AUC.causal.mean.fac,
                                 values = AUC.causal.mean.values,
                                 sd = AUC.causal.sd.values)
    
    p.causal <- ggplot(data = df.causal.mean, aes(x = method, y = values, fill = AUCs)) +
      geom_bar(stat = "identity", width = 0.6, position = position_dodge(0.65)) +
      ylim(0, 1) +
      scale_color_manual(values = c("red", "green")) + 
      ggtitle(paste0("AUC around causal haplotype block")) +
      theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Verdana", face = "bold"),
            axis.title = element_text(size = 22),
            axis.text.x = element_text(size = 25, face = "bold"),
            axis.text.y = element_text(size = 26),
            legend.text = element_text(size = 23, face = "bold"),
            legend.title = element_text(size = 18),
            legend.key.height = unit(31, units = "pt"))
    # p.causal <- p.causal + geom_errorbar(aes(ymin = values - sd, ymax = values + sd, width = 0.25),
    #                          position = position_dodge(0.65)) 
    file.barplot <- paste0(dir.res.base, "/Results_of_AUC_causals_barplot_",
                           "#_of_top_false_", n.top.false.now, ".png")
    png(file.barplot, height = 700, width = 800)
    print(p.causal)
    dev.off()
    
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











direction.eff <- "minus"
dir.base <- paste0("midstream/number_of_causals=", n.causal.site,
                   "_direction_of_effect=", direction.eff,
                   "_trial_no=", trial.no)


dir.res.base <- paste0("results/number_of_causals=", n.causal.site,
                       "_direction_of_effect=", direction.eff,
                       "_trial_no=", trial.no)
dir.create(dir.res.base)

file.pheno <- paste0("data/phenotype/number_of_causals=", n.causal.site,
                     "_direction_of_effect=", direction.eff,
                     "_trial_no=", trial.no, ".csv")
pheno <- as.matrix(read.csv(file = file.pheno))

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
}




###### 4. Summarize results and draw boxplots ###### 
##### 4.1. Summarize results #####
load(paste0(dir.base, "/inflators.RData"))
inflator.thress <- t(apply(inflation.levels, c(3, 2), mean)) + 
  inflator.thres.base - apply(inflation.levels, c(3, 2), mean)[2, ]

for(thres.no in 1:n.thres){
  for(inflator.no in 1:n.inflator){
    thres.now <- thres.names[thres.no]
    n.top.false.now <- n.top.false[inflator.no]
    
    file.df.summary <- paste0(dir.base, "/Summary_statistics_of_", n.iter, "simulations_", thres.now,
                              "_#_of_top_false_", n.top.false.now, ".csv")
    df.summary <- read.csv(file = file.df.summary, row.names = 1)
    df.summary[, 1] <- method.fac.abb
    
    log10p1s.now <- matrix(df.summary$log10p.1, nrow = n.iter)
    log10p1.relaxes.now <- matrix(df.summary$log10p.1.relax, nrow = n.iter)
    adjusted.log10p1s.now <- matrix(df.summary$Adjusted.log10p.1, nrow = n.iter)
    adjusted.log10p1.relaxes.now <- matrix(df.summary$Adjusted.log10p.1.relax, nrow = n.iter)
    
    log10p2s.now <- matrix(df.summary$log10p.2, nrow = n.iter)
    log10p2.relaxes.now <- matrix(df.summary$log10p.2.relax, nrow = n.iter)
    adjusted.log10p2s.now <- matrix(df.summary$Adjusted.log10p.2, nrow = n.iter)
    adjusted.log10p2.relaxes.now <- matrix(df.summary$Adjusted.log10p.2.relax, nrow = n.iter)
    
    Recalls.now <- matrix(df.summary$Recall, nrow = n.iter)
    Precisions.now <- matrix(df.summary$Precision, nrow = n.iter)
    Hms.now <- matrix(df.summary$F_value, nrow = n.iter)
    Recalls.QTN12.now <- matrix(df.summary$Recall.QTN12, nrow = n.iter)
    Precisions.QTN12.now <- matrix(df.summary$Precision.QTN12, nrow = n.iter)
    Hms.QTN12.now <- matrix(df.summary$F_value.QTN12, nrow = n.iter)
    Recalls.QTN3.now <- matrix(df.summary$Recall.QTN3, nrow = n.iter)
    Precisions.QTN3.now <- matrix(df.summary$Precision.QTN3, nrow = n.iter)
    Hms.QTN3.now <- matrix(df.summary$F_value.QTN3, nrow = n.iter)
    Recalls.other.now <- matrix(df.summary$Recall.other, nrow = n.iter)
    Precisions.other.now <- matrix(df.summary$Precision.other, nrow = n.iter)
    Hms.other.now <- matrix(df.summary$F_value.other, nrow = n.iter)
    AUC.causals <- matrix(df.summary$AUC.causal, nrow = n.iter)
    
    ##### 4.2.1. Draw boxplots for QTN1 with labels #####
    p.1 <- ggplot(data = df.summary, aes(x = method, y = log10p.1)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.plus[1], y.range.plus[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold(-log[bold("10")](bolditalic(p))), bold(" in terms of haplotype block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 30),
            axis.text.x = element_text(size = 22, face = "bold"),
            axis.text.y = element_text(size = 29))
    
    p.2 <- ggplot(data = df.summary, aes(x = method, y = log10p.1.relax)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.plus[1], y.range.plus[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold(-log[bold("10")](bolditalic(p))), bold(" in terms of LD block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 30),
            axis.text.x = element_text(size = 22, face = "bold"),
            axis.text.y = element_text(size = 29))
    
    
    p.3 <- ggplot(data = df.summary, aes(x = method, y = Adjusted.log10p.1)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.adjusted.plus[1], y.range.adjusted.plus[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold("Adjusted "), bold(-log[bold("10")](bolditalic(p))), bold(" in terms of haplotype block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 30),
            axis.text.x = element_text(size = 22, face = "bold"),
            axis.text.y = element_text(size = 29))
    
    
    p.4 <- ggplot(data = df.summary, aes(x = method, y = Adjusted.log10p.1.relax)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.adjusted.plus[1], y.range.adjusted.plus[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold("Adjusted "), bold(-log[bold("10")](bolditalic(p))), bold(" in terms of LD block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 30),
            axis.text.x = element_text(size = 22, face = "bold"),
            axis.text.y = element_text(size = 29))
    
    file.boxplot <- paste0(dir.res.base, "/Results_of_-log10(p)_QTN12_boxplot_#_of_top_false_", n.top.false.now, "_with_label.png")
    png(file.boxplot, height = 1000, width = 1300)
    grid.arrange(p.1, p.2, p.3, p.4, ncol = 2)
    dev.off()
    
    
    
    ##### 4.2.2 Draw boxplots of QTN3 with label #####
    p.1 <- ggplot(data = df.summary, aes(x = method, y = log10p.2)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.plus.2[1], y.range.plus.2[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold(-log[bold("10")](bolditalic(p))), bold(" in terms of haplotype block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 30),
            axis.text.x = element_text(size = 22, face = "bold"),
            axis.text.y = element_text(size = 29))
    
    p.2 <- ggplot(data = df.summary, aes(x = method, y = log10p.2.relax)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.plus.2[1], y.range.plus.2[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold(-log[bold("10")](bolditalic(p))), bold(" in terms of LD block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 30),
            axis.text.x = element_text(size = 22, face = "bold"),
            axis.text.y = element_text(size = 29))
    
    
    p.3 <- ggplot(data = df.summary, aes(x = method, y = Adjusted.log10p.2)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.adjusted.plus.2[1], y.range.adjusted.plus.2[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold("Adjusted "), bold(-log[bold("10")](bolditalic(p))), bold(" in terms of haplotype block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 30),
            axis.text.x = element_text(size = 22, face = "bold"),
            axis.text.y = element_text(size = 29))
    
    
    p.4 <- ggplot(data = df.summary, aes(x = method, y = Adjusted.log10p.2.relax)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.adjusted.plus.2[1], y.range.adjusted.plus.2[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold("Adjusted "), bold(-log[bold("10")](bolditalic(p))), bold(" in terms of LD block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 30),
            axis.text.x = element_text(size = 22, face = "bold"),
            axis.text.y = element_text(size = 29))
    
    file.boxplot <- paste0(dir.res.base, "/Results_of_-log10(p)_QTN3_boxplot_#_of_top_false_", n.top.false.now, "_with_label.png")
    png(file.boxplot, height = 1000, width = 1300)
    grid.arrange(p.1, p.2, p.3, p.4, ncol = 2)
    dev.off()
    
    
    
    
    ##### 4.2.3. Draw boxplots of QTN1 with no label #####
    p.1 <- ggplot(data = df.summary, aes(x = method, y = log10p.1)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.plus[1], y.range.plus[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold(-log[bold("10")](bolditalic(p))), bold(" in terms of haplotype block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 23, face = "bold"),
            axis.text.y = element_blank())
    
    p.2 <- ggplot(data = df.summary, aes(x = method, y = log10p.1.relax)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.plus[1], y.range.plus[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold(-log[bold("10")](bolditalic(p))), bold(" in terms of LD block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 23, face = "bold"),
            axis.text.y = element_blank())
    
    
    p.3 <- ggplot(data = df.summary, aes(x = method, y = Adjusted.log10p.1)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.adjusted.plus[1], y.range.adjusted.plus[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold("Adjusted "), bold(-log[bold("10")](bolditalic(p))), bold(" in terms of haplotype block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 23, face = "bold"),
            axis.text.y = element_blank())
    
    
    p.4 <- ggplot(data = df.summary, aes(x = method, y = Adjusted.log10p.1.relax)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.adjusted.plus[1], y.range.adjusted.plus[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold("Adjusted "), bold(-log[bold("10")](bolditalic(p))), bold(" in terms of LD block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 23, face = "bold"),
            axis.text.y = element_blank())
    
    file.boxplot <- paste0(dir.res.base, "/Results_of_-log10(p)_QTN12_boxplot_#_of_top_false_", n.top.false.now, "_with_no_label.png")
    png(file.boxplot, height = 1000, width = 1100)
    grid.arrange(p.1, p.2, p.3, p.4, ncol = 2)
    dev.off()
    
    
    
    ##### 4.2.4 Draw boxplots of QTN3 with no label #####
    p.1 <- ggplot(data = df.summary, aes(x = method, y = log10p.2)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.plus.2[1], y.range.plus.2[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold(-log[bold("10")](bolditalic(p))), bold(" in terms of haplotype block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 23, face = "bold"),
            axis.text.y = element_blank())
    
    p.2 <- ggplot(data = df.summary, aes(x = method, y = log10p.2.relax)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.plus.2[1], y.range.plus.2[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold(-log[bold("10")](bolditalic(p))), bold(" in terms of LD block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 23, face = "bold"),
            axis.text.y = element_blank())
    
    
    p.3 <- ggplot(data = df.summary, aes(x = method, y = Adjusted.log10p.2)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.adjusted.plus.2[1], y.range.adjusted.plus.2[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold("Adjusted "), bold(-log[bold("10")](bolditalic(p))), bold(" in terms of haplotype block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 23, face = "bold"),
            axis.text.y = element_blank())
    
    
    p.4 <- ggplot(data = df.summary, aes(x = method, y = Adjusted.log10p.2.relax)) + 
      ylab(expression(-log[10](italic(p)))) +
      ylim(y.range.adjusted.plus.2[1], y.range.adjusted.plus.2[2]) +
      geom_boxplot(fill = c("pink", "skyblue", "green", "green", "lightgreen", "lightgreen", "yellow", "yellow", "gray")) +
      ggtitle(expression(paste(bold("Adjusted "), bold(-log[bold("10")](bolditalic(p))), bold(" in terms of LD block")))) +
      theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Verdana"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 23, face = "bold"),
            axis.text.y = element_blank())
    
    file.boxplot <- paste0(dir.res.base, "/Results_of_-log10(p)_QTN3_boxplot_#_of_top_false_", n.top.false.now, "_with_no_label.png")
    png(file.boxplot, height = 1000, width = 1100)
    grid.arrange(p.1, p.2, p.3, p.4, ncol = 2)
    dev.off()
    
    
    
    ###### 5. Draw barplot for Recall, Precision and Hm ###### 
    method.mean.fac <- factor(rep(method.names.abb, 3), 
                              levels = method.levels.abb)
    SS.mean.fac <- factor(rep(c("Recall", "Precision", "F-measure"), each = n.method),
                          levels = c("Recall", "Precision", "F-measure"))
    SS.mean.values <- c(apply(Recalls.now, 2, mean, na.rm = TRUE),
                        apply(Precisions.now, 2, mean, na.rm = TRUE),
                        apply(Hms.now, 2, mean, na.rm = TRUE))
    SS.sd.values <- c(apply(Recalls.now, 2, sd, na.rm = TRUE),
                      apply(Precisions.now, 2, sd, na.rm = TRUE),
                      apply(Hms.now, 2, sd, na.rm = TRUE))
    df.summary.mean <- data.frame(method = method.mean.fac,
                                  summary_statistics = SS.mean.fac,
                                  values = SS.mean.values,
                                  sd = SS.sd.values)
    
    p <- ggplot(data = df.summary.mean, aes(x = method, y = values, fill = summary_statistics)) +
      geom_bar(stat = "identity", width = 0.6, position = position_dodge(0.65)) +
      ylim(0, 1) +
      scale_color_manual(values = c("red", "green", "blue")) + 
      ggtitle(paste0("Summary statistics of GWAS results")) +
      theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Verdana", face = "bold"),
            axis.title = element_text(size = 22),
            axis.text.x = element_text(size = 25, face = "bold"),
            axis.text.y = element_text(size = 26),
            legend.text = element_text(size = 23, face = "bold"),
            legend.title = element_text(size = 18),
            legend.key.height = unit(31, units = "pt"))
    # p <- p + geom_errorbar(aes(ymin = values - sd, ymax = values + sd, width = 0.25),
    #                          position = position_dodge(0.65)) 
    file.barplot <- paste0(dir.res.base, "/Results_of_summary_statistics_barplot_", thres.now, ".png")
    png(file.barplot, height = 700, width = 800)
    print(p)
    dev.off()
    
    
    
    SS.mean.fac.QTN12 <- factor(rep(c("Recall", "Precision", "F-measure"), each = n.method),
                                levels = c("Recall", "Precision", "F-measure"))
    SS.mean.values.QTN12 <- c(apply(Recalls.QTN12.now, 2, mean, na.rm = TRUE),
                              apply(Precisions.QTN12.now, 2, mean, na.rm = TRUE),
                              apply(Hms.QTN12.now, 2, mean, na.rm = TRUE))
    SS.sd.values.QTN12 <- c(apply(Recalls.QTN12.now, 2, sd, na.rm = TRUE),
                            apply(Precisions.QTN12.now, 2, sd, na.rm = TRUE),
                            apply(Hms.QTN12.now, 2, sd, na.rm = TRUE))
    df.summary.mean.QTN12 <- data.frame(method = method.mean.fac,
                                        summary_statistics = SS.mean.fac.QTN12,
                                        values = SS.mean.values.QTN12,
                                        sd = SS.sd.values.QTN12)
    
    p.QTN12 <- ggplot(data = df.summary.mean.QTN12, aes(x = method, y = values, fill = summary_statistics)) +
      geom_bar(stat = "identity", width = 0.6, position = position_dodge(0.65)) +
      ylim(0, 1) +
      scale_color_manual(values = c("red", "green", "blue")) + 
      ggtitle(paste0("Summary statistics of GWAS results")) +
      theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Verdana", face = "bold"),
            axis.title = element_text(size = 22),
            axis.text.x = element_text(size = 25, face = "bold"),
            axis.text.y = element_text(size = 26),
            legend.text = element_text(size = 23, face = "bold"),
            legend.title = element_text(size = 18),
            legend.key.height = unit(31, units = "pt"))
    # p.QTN12 <- p.QTN12 + geom_errorbar(aes(ymin = values - sd, ymax = values + sd, width = 0.25),
    #                          position = position_dodge(0.65)) 
    file.barplot.QTN12 <- paste0(dir.res.base, "/Results_of_summary_statistics_barplot_", thres.now, "_QTN12.png")
    png(file.barplot.QTN12, height = 700, width = 800)
    print(p.QTN12)
    dev.off()
    
    
    
    SS.mean.fac.QTN3 <- factor(rep(c("Recall", "Precision", "F-measure"), each = n.method),
                               levels = c("Recall", "Precision", "F-measure"))
    SS.mean.values.QTN3 <- c(apply(Recalls.QTN3.now, 2, mean, na.rm = TRUE),
                             apply(Precisions.QTN3.now, 2, mean, na.rm = TRUE),
                             apply(Hms.QTN3.now, 2, mean, na.rm = TRUE))
    SS.sd.values.QTN3 <- c(apply(Recalls.QTN3.now, 2, sd, na.rm = TRUE),
                           apply(Precisions.QTN3.now, 2, sd, na.rm = TRUE),
                           apply(Hms.QTN3.now, 2, sd, na.rm = TRUE))
    df.summary.mean.QTN3 <- data.frame(method = method.mean.fac,
                                       summary_statistics = SS.mean.fac.QTN3,
                                       values = SS.mean.values.QTN3,
                                       sd = SS.sd.values.QTN3)
    
    p.QTN3 <- ggplot(data = df.summary.mean.QTN3, aes(x = method, y = values, fill = summary_statistics)) +
      geom_bar(stat = "identity", width = 0.6, position = position_dodge(0.65)) +
      ylim(0, 1) +
      scale_color_manual(values = c("red", "green", "blue")) + 
      ggtitle(paste0("Summary statistics of GWAS results")) +
      theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Verdana", face = "bold"),
            axis.title = element_text(size = 22),
            axis.text.x = element_text(size = 25, face = "bold"),
            axis.text.y = element_text(size = 26),
            legend.text = element_text(size = 23, face = "bold"),
            legend.title = element_text(size = 18),
            legend.key.height = unit(31, units = "pt"))
    # p.QTN3 <- p.QTN3 + geom_errorbar(aes(ymin = values - sd, ymax = values + sd, width = 0.25),
    #                          position = position_dodge(0.65)) 
    file.barplot.QTN3 <- paste0(dir.res.base, "/Results_of_summary_statistics_barplot_", thres.now, "_QTN3.png")
    png(file.barplot.QTN3, height = 700, width = 800)
    print(p.QTN3)
    dev.off()
    
    
    
    
    SS.mean.fac.other <- factor(rep(c("Recall", "Precision", "F-measure"), each = n.method),
                                levels = c("Recall", "Precision", "F-measure"))
    SS.mean.values.other <- c(apply(Recalls.other.now, 2, mean, na.rm = TRUE),
                              apply(Precisions.other.now, 2, mean, na.rm = TRUE),
                              apply(Hms.other.now, 2, mean, na.rm = TRUE))
    SS.sd.values.other <- c(apply(Recalls.other.now, 2, sd, na.rm = TRUE),
                            apply(Precisions.other.now, 2, sd, na.rm = TRUE),
                            apply(Hms.other.now, 2, sd, na.rm = TRUE))
    df.summary.mean.other <- data.frame(method = method.mean.fac,
                                        summary_statistics = SS.mean.fac.other,
                                        values = SS.mean.values.other,
                                        sd = SS.sd.values.other)
    
    p.other <- ggplot(data = df.summary.mean.other, aes(x = method, y = values, fill = summary_statistics)) +
      geom_bar(stat = "identity", width = 0.6, position = position_dodge(0.65)) +
      ylim(0, 1) +
      scale_color_manual(values = c("red", "green", "blue")) + 
      ggtitle(paste0("Summary statistics of GWAS results")) +
      theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Verdana", face = "bold"),
            axis.title = element_text(size = 22),
            axis.text.x = element_text(size = 25, face = "bold"),
            axis.text.y = element_text(size = 26),
            legend.text = element_text(size = 23, face = "bold"),
            legend.title = element_text(size = 18),
            legend.key.height = unit(31, units = "pt"))
    # p.other <- p.other + geom_errorbar(aes(ymin = values - sd, ymax = values + sd, width = 0.25),
    #                          position = position_dodge(0.65)) 
    file.barplot.other <- paste0(dir.res.base, "/Results_of_summary_statistics_barplot_inflator_threshold_",
                                 "#_of_top_false_", n.top.false.now, ".png")
    png(file.barplot.other, height = 700, width = 800)
    print(p.other)
    dev.off()
    
    
    
    
    
    
    ###### 6. Overwhelming results ###### 
    greatest.method.no.0 <- apply(adjusted.log10p1s.now, 1, which.max)
    greatest.method.no <- apply(adjusted.log10p1.relaxes.now, 1, which.max)
    greatest.method.count.0 <- table(greatest.method.no.0)
    greatest.method.count <- table(greatest.method.no)
    
    if(length(greatest.method.count.0) <= (n.method - 1)){
      missing.no <- which(!((1:n.method) %in% names(greatest.method.count.0)))
      greatest.method.count.0 <- c(greatest.method.count.0, rep(0, length(missing.no)))
      names(greatest.method.count.0) <-
        c(names(greatest.method.count.0)[1:(n.method - length(missing.no))], missing.no)
      greatest.method.count.0 <- greatest.method.count.0[order(names(greatest.method.count.0))]
    }
    
    if(length(greatest.method.count) <= (n.method - 1)){
      missing.no <- which(!((1:n.method) %in% names(greatest.method.count)))
      greatest.method.count <- c(greatest.method.count, rep(0, length(missing.no)))
      names(greatest.method.count) <-
        c(names(greatest.method.count)[1:(n.method - length(missing.no))], missing.no)
      greatest.method.count <- greatest.method.count[order(names(greatest.method.count))]
    }
    
    
    greatest.method.count.all <- rbind(greatest.method.count.0,
                                       greatest.method.count)
    rownames(greatest.method.count.all) <- c("Restrict", "Relax")
    
    file.greate.method <- paste0(dir.res.base, "/greatest_method_count_", 
                                 "#_of_top_false_", n.top.false.now, ".csv")
    write.csv(greatest.method.count.all, file = file.greate.method)
    
    
    for(no.now in 1:n.method){
      adjusted.log10p1.relaxes.now[(greatest.method.no == no.now) & (Recalls.QTN12.now[, no.now] > 0) &
                                     (adjusted.log10p1.relaxes.now[, no.now] >= inflator.thress[no.now, inflator.no]) &
                                     apply(adjusted.log10p1.relaxes.now[, - no.now] < (inflator.thress[no.now, inflator.no] - inflator.thres.diff), 1, all), , drop = FALSE]
      
      overwhelm.no <- which((greatest.method.no == no.now) & (Recalls.QTN12.now[, no.now] > 0) &
                              (adjusted.log10p1.relaxes.now[, no.now] >= inflator.thress[no.now, inflator.no]) &
                              apply(adjusted.log10p1.relaxes.now[, - no.now] < (inflator.thress[no.now, inflator.no] - inflator.thres.diff), 1, all))
      
      if(length(overwhelm.no) >= 1){
        dir.overwhelm <- paste0(dir.res.base, "/overwhelming_results_", method.names.abb[no.now], "_", thres.now,
                                "_#_of_top_false_", n.top.false.now)
        dir.create(dir.overwhelm)
        
        dir.overwhelm.manhattan <- paste0(dir.res.base, "/overwhelming_results_", method.names.abb[no.now], "_", thres.now,
                                          "_#_of_top_false_", n.top.false.now, "_manhattan")
        dir.create(dir.overwhelm.manhattan)
        
        for(iter.no in overwhelm.no){
          file.overwhelms.0 <- list.files(paste0(dir.iter.0, iter.no, "/"))
          file.overwhelms.manhattan.0 <- file.overwhelms.0[manhattan.no]
          file.overwhelms <- paste0(dir.iter.0, iter.no, "/", file.overwhelms.0)
          file.overwhelms.manhattan <- file.overwhelms[manhattan.no]
          
          dir.overwhelm.now <- paste0(dir.overwhelm, "/iteration_", iter.no)
          dir.create(dir.overwhelm.now)
          file.overwhelm.news <- paste0(dir.overwhelm.now, "/", file.overwhelms.0)
          
          dir.overwhelm.manhattan.now <- paste0(dir.overwhelm.manhattan, "/iteration_", iter.no)
          dir.create(dir.overwhelm.manhattan.now)
          file.overwhelm.manhattan.news <- paste0(dir.overwhelm.manhattan.now, "/", file.overwhelms.manhattan.0)
          
          file.copy(file.overwhelms, file.overwhelm.news, overwrite = TRUE)
          file.copy(file.overwhelms.manhattan, file.overwhelm.manhattan.news, overwrite = TRUE)
        }
      }
    }
    
    ######  7. Whether each method can detect causal haplotype block itself ###### 
    method.fac.causal <- factor(rep(method.names.abb, 2), 
                                levels = method.levels.abb)
    AUC.causal.mean.fac <- factor(rep(c("all", "detected"), each = n.method),
                                  levels = c("all", "detected"))
    AUC.detected.mean <- AUC.detected.sd <- rep(NA, n.method)
    
    for(method.no in 1:n.method){
      AUC.causals.detected <- AUC.causals[adjusted.log10p1.relaxes.now[, method.no] >= inflator.thress[method.no, inflator.no], method.no]
      AUC.detected.mean[method.no] <- mean(AUC.causals.detected, na.rm = TRUE)
      AUC.detected.sd[method.no] <- sd(AUC.causals.detected, na.rm = TRUE)
    }
    
    AUC.causal.mean.values <- c(apply(AUC.causals, 2, mean, na.rm = TRUE), AUC.detected.mean)
    AUC.causal.sd.values <- c(apply(AUC.causals, 2, sd, na.rm = TRUE), AUC.detected.sd)
    # AUC.causal.mean.values <- c(apply(AUC.causals, 2, mean, na.rm = TRUE),
    #                             apply(AUC.causals[apply(adjusted.log10p1.relaxes.now >= inflator.thress[thres.no], 1, any), ], 2, mean, na.rm = T))
    # AUC.causal.sd.values <- c(apply(AUC.causals, 2, sd, na.rm = TRUE),
    #                           apply(AUC.causals[apply(adjusted.log10p1.relaxes.now >= inflator.thress[thres.no], 1, any), ], 2, sd, na.rm = T))
    df.causal.mean <- data.frame(method = method.fac.causal,
                                 AUCs = AUC.causal.mean.fac,
                                 values = AUC.causal.mean.values,
                                 sd = AUC.causal.sd.values)
    
    p.causal <- ggplot(data = df.causal.mean, aes(x = method, y = values, fill = AUCs)) +
      geom_bar(stat = "identity", width = 0.6, position = position_dodge(0.65)) +
      ylim(0, 1) +
      scale_color_manual(values = c("red", "green")) + 
      ggtitle(paste0("AUC around causal haplotype block")) +
      theme(plot.title = element_text(size = 30, hjust = 0.5, family = "Verdana", face = "bold"),
            axis.title = element_text(size = 22),
            axis.text.x = element_text(size = 25, face = "bold"),
            axis.text.y = element_text(size = 26),
            legend.text = element_text(size = 23, face = "bold"),
            legend.title = element_text(size = 18),
            legend.key.height = unit(31, units = "pt"))
    # p.causal <- p.causal + geom_errorbar(aes(ymin = values - sd, ymax = values + sd, width = 0.25),
    #                          position = position_dodge(0.65)) 
    file.barplot <- paste0(dir.res.base, "/Results_of_AUC_causals_barplot_",
                           "#_of_top_false_", n.top.false.now, ".png")
    png(file.barplot, height = 700, width = 800)
    print(p.causal)
    dev.off()
    
  }
}