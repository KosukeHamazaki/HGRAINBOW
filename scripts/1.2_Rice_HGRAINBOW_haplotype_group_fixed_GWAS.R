#############################################################
####### 1.2_Rice_HGRAINBOW_haplotype_group_fixed_GWAS #######
#############################################################



require(RAINBOW)
require(ape)


score_HGF_geneset <- function(M, y, gene.set, ZETA, num.hap = NULL, n.PC = 0,
                              covariate = NULL, map, grouping.method = "phylo", count = TRUE) {
  chr <- map[, 2]
  chr.tab <- table(chr)
  chr.max <- max(chr)
  chr.cum <- cumsum(chr.tab)
  
  if(is.null(covariate)){
    covariate <- matrix(1, nrow = length(y), 1)
  }
  X <- covariate
  
  K.A <- ZETA[[1]]$K
  Z.A <- ZETA[[1]]$Z
  
  if (n.PC > 0) {
    eigen.K.A <- eigen(K.A)
    eig.K.vec <- eigen.K.A$vectors
    PC.part <- Z.A %*% eig.K.vec[, 1:n.PC]
    colnames(PC.part) <- paste0("n.PC_", 1:n.PC)
    X <- cbind(X, PC.part)
  }
  
  EMM.res.0 <- EMM.cpp(y = y, X = covariate, ZETA = ZETA, return.Hinv = TRUE)
  Hinv <- EMM.res.0$Hinv
  
  
  
  gene.names <- as.character(gene.set[, 1])
  mark.id <- as.character(gene.set[, 2])
  gene.name <- as.character(unique(gene.names))
  n.scores <- length(unique(gene.set[, 1]))
  
  scores <- rep(NA, n.scores)
  
  pb <- txtProgressBar(min = 1, max = n.scores, style = 3)
  n.scores2 <- n.scores - n.scores %% 100
  
  start.scorecalc <- Sys.time()
  for (i in 1:n.scores) {
    if (count) {
      setTxtProgressBar(pb, i)
      if (i == (n.scores2/100 + 1) | i == (n.scores2/10 + 
                                           1) | i == (n.scores2/2 + 1)) {
        cat("\n")
        end.scorecalc <- Sys.time()
        jikan.scorecalc <- (end.scorecalc - start.scorecalc) * 
          (n.scores - i + 1)/(i - 1)
        print(paste0((i - 1) * 100/n.scores2, "%...Done. ", 
                     round(jikan.scorecalc, 2), " ", attr(jikan.scorecalc, 
                                                          "units"), " to end.  Sceduled end time : ", 
                     end.scorecalc + jikan.scorecalc))
      }
    }
    
    
    mark.name.now <- mark.id[gene.names == gene.name[i]]
    Mis.range <- match(mark.name.now, map[, 1])
    Mis.range2 <- 1:length(Mis.range)
    weighting.center <- FALSE
    
    Mis.0 <- as.matrix(M[, Mis.range])
    
    
    Mis.fac <- factor(apply(Mis.0, 1, function(x) paste(x, collapse = "")))
    if(length(levels(Mis.fac)) > 1){
      if(is.null(num.hap)){
        Mis <- as.matrix(Mis.0[!duplicated(as.numeric(Mis.fac)), ])
        bango <- as.factor(as.numeric(Mis.fac))
        levels(bango) <- order(unique(bango))
        bango <- as.numeric(as.character(bango))
        num.hap <- nrow(Mis)
      }else{
        if(grouping.method == "kmed"){
          kmed.res <- cluster::pam(Mis.0, k = num.hap)
          Mis <- kmed.res$medoids
          bango <- kmed.res$clustering
        }
        
        if(grouping.method == "phylo"){
          dist.mat <- as.matrix(dist(Mis.0, method = "manhattan"))
          dist.ev.mat <- - 3 * log(1 - (4 * (dist.mat / (2 * ncol(x)))/3)) / 4
          
          hc.res <- hclust(as.dist(dist.ev.mat), method = "average")
          bango <- cutree(hc.res, k = num.hap)
        }
        
        if(!(grouping.method %in% c("kmed", "phylo"))){
          stop("We only support 'kmed' and 'phylo'!")
        }
      }
      
      
      Xi.test <- as.matrix(Matrix::sparseMatrix(i = 1:nrow(M), j = bango, x = rep(1, nrow(M)),
                                                dims = c(nrow(M), num.hap)))
    }else{
      Xi.test <- Mis.0[, 1, drop = F]
    }
    
    
    
    ni <- length(y)
    yi <- matrix(y, ni, 1)
    Xi <- cbind(X, Xi.test[, -ncol(Xi.test)])
    p <- ncol(Xi)
    v1 <- ncol(Xi.test)
    v2 <- ni - p
    
    beta.stat <- try(GWAS_F_test(y = yi, x = Xi, hinv = Hinv, 
                                 v1 = v1, v2 = v2, p = p), silent = TRUE)
    if (class(beta.stat) != "try-error") {
      scores[i] <- -log10(pbeta(beta.stat, v2/2, v1/2))
    }
  }
  
  return(scores)
}
