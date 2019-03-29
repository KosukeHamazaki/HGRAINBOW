#####################################################################################
###### 0.3_Rice_HGRAINBOW_Simulation_of_phenotypic_values_and_some_preparation ######
#####################################################################################

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
dir.create(dir.base)

seeds <- matrix(sample(1:1000000, n.iter * 4), nrow = n.iter)
file.seeds <- paste0(dir.base, "/seeds.csv")
write.csv(seeds, file = file.seeds)




###### 2. Modification of raw data ######
##### 2.1. Read original files into R #####
geno <- data.frame(fread("data/genotype/L3024_core_extract_L414_ind1A_ind1B_MAF_0.025_geno.tsv"))
block.list <- data.frame(fread("data/extra/haplotype_block_list.csv"))


##### 2.2. Some modification #####
chr <- geno[, 1]
pos <- geno[, 2]
marker <- rownames(geno) 
map <- data.frame(marker, chr, pos)

x <- t(geno[, -c(1:2)])
n.line <- nrow(x)
n.marker <- ncol(x)
line.names <- rownames(x)

n.haplo <- length(unique(block.list[, 1]))
n.marker.haplo <- nrow(block.list)




###### 3. Simulation of phenotypic values ######
##### 3.1. Make directory to "midstream" #####
dir.iter.0 <- paste0(dir.base, "/iteration_")
phenotypes <- matrix(NA, nrow = n.iter, ncol = n.line)
K.A <- NULL

for(i in 1:n.iter){
  seed <- seeds[i, ]
  dir.iter <- paste0(dir.iter.0, i)
  dir.create(dir.iter)
  
  
  ##### 3.2. Fixed effects #####
  #### 3.2.1. Effect from haplotype ####
  block.list.tab <- table(block.list[, 1])
  block.list.tab.sort <- block.list.tab[match(unique(block.list[, 1]), names(block.list.tab))]
  haplo.big.no <- (1:n.haplo)[block.list.tab.sort >= n.haplo.mark.min]
  set.seed(seed[1])
  haplo.gene.site <- sample(haplo.big.no, size = 1)
  haplo.marks <- block.list[block.list[, 1] == paste0("haploblock_", haplo.gene.site), 2]
  
  set.seed(seed[1])
  causal.site <- sample(haplo.marks, n.causal.site)
  
  x.causal.site <- x[, causal.site]
  beta.haplo.eff <- 1 / apply(x.causal.site, 2, sd)
  
  cor.x.causal.site <- cor(x.causal.site)[1, ]
  cor.x.causal.site.minus <- cor.x.causal.site <= 0
  beta.haplo.eff[cor.x.causal.site.minus] <- -beta.haplo.eff[cor.x.causal.site.minus]
  if(direction.eff == "minus"){
    n.minus <- n.causal.site %/% 2
    beta.haplo.eff[(n.causal.site - n.minus + 1):n.causal.site] <- 
      -beta.haplo.eff[(n.causal.site - n.minus + 1):n.causal.site]
  }
  
  beta.haplo.eff <- beta.haplo.eff * weight.haplo
  
  
  
  #### 3.2.2. Effect from another major gene ####
  set.seed(seed[2])
  gene.another <- sample((1:n.marker)[chr != chr[causal.site[1]]], n.gene - 1)
  x.gene.another <- x[, gene.another, drop = F]
  
  beta.gene.another <- 1 / apply(x.gene.another, 2, sd)
  
  
  
  #### 3.2.3. Merge fixed effects ####
  causals <- c(causal.site, gene.another)
  X <- cbind(x.causal.site, x.gene.another)
  beta <- c(beta.haplo.eff, beta.gene.another)
  
  Xbeta <- c(X %*% beta)
  
  chr.causals <- chr[causals]
  pos.causals <- pos[causals]
  marker.causals <- marker[causals]
  map.causals <- map[causals, ]
  
  
  
  ##### 3.3. Random effects #####
  #### 3.3.1. Estimate additive genetic relationship matrix ####
  if(is.null(K.A)){
    K.A <- A.mat(x)
  }  
  
  
  #### 3.3.1. Polygenetic effects ####
  Vu <- c(var(Xbeta)) / prop
  set.seed(seed[3])
  u <- mvrnorm(n = 1, mu = rep(0, n.line), Sigma = Vu * K.A)
  
  Xbetaplusu <- Xbeta + u
  Xbetaplusu.scaled <- Xbetaplusu / sd(Xbetaplusu)
  
  
  
  ##### 3.4. Residuals #####
  Ve <- (1 - h2) / h2 
  set.seed(seed[4])
  e <- rnorm(n = n.line, mean = 0, sd = Ve)
  
  
  
  ##### 3.5. Simulated phenotypic values #####
  (y <- Xbetaplusu.scaled + e)
  
  
  
  ##### 3.6. Save the results #####
  all.sim.res <- list(phenotype = y, causal.info = map.causals,
                      Xbetaplusu = Xbetaplusu.scaled, e = e, weight.haplo = weight.haplo)
  
  file.map <- paste0(dir.iter, "/causal_information.csv")
  fwrite(map.causals, file.map)
  
  file.all <- paste0(dir.iter, "/all_simulation_results.RData")
  save(all.sim.res, file = file.all)
  
  phenotypes[i, ] <- y
}


dir.data.pheno <- paste0("data/phenotype")
file.pheno <- paste0(dir.data.pheno, "/number_of_causals=", n.causal.site,
                     "_direction_of_effect=", direction.eff,
                     "_trial_no=", trial.no, ".csv")
colnames(phenotypes) <- line.names
fwrite(data.frame(phenotypes), file.pheno)
