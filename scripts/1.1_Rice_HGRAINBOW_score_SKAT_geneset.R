#################################################
####### 1.1_Rice_HGPAG_score_SKAT_geneset #######
#################################################



require(SKAT)

score_SKAT_geneset <- function(M, y, gene.set, covariate = NULL, map, count = TRUE) {
  if(is.null(covariate)){
    obj <- SKAT_Null_Model(y ~ 1, out_type="C")
  }else{
    obj <- SKAT_Null_Model(y ~ covariate, out_type="C")    
  } 
  
  chr <- map[, 2]
  chr.tab <- table(chr)
  chr.max <- max(chr)
  chr.cum <- cumsum(chr.tab)
  
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

    scores[i] <- -log10(SKAT(Mis.0, obj)$p.value)
  }
  
  return(scores)
}
