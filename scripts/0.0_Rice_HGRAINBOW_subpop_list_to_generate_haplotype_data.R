#######################################################################
###### 0.0_Rice_HGRAINBOW_subpop_list_to_generate_haplotype_data ######
#######################################################################

###### 1. Settings ######
##### 1.0. Reset workspace ######
#rm(list=ls())



##### 1.1. Setting working directory to the "HGRAINBOW" directory #####
project <- "HGRAINBOW"
setwd(paste0("/media/hamazaki/d4000953-5e56-40ce-97ea-4cfee57fc91a/research/rice/Project/", project))



##### 1.2. Import packages #####
require(data.table)



##### 1.3. Setting some parameters #####
subpop.sel <- c("ind1A", "ind1B")




###### 2. Modification of raw data ######
##### 2.1. Read original files into R #####
info.list <- fread("raw_data/extra/L3024_404K_core_SNP.full_data.name.csv")



##### 2.2. Extract information for data modyfing #####
info.iris <- as.data.frame(info.list)[, 2]
info.iris[247:253] <- paste0("CX", c(2, 3, 4, 5, 6, 8, 9))
#View(cbind(colnames(data.0)[-c(1:2)], info.iris))   ### checking if the lines of genotype and list match

info.subpop <- as.data.frame(info.list)[, 4]

(subpop.tab <- table(info.subpop))
subpop.name <- names(subpop.tab)

subpop.no <- match(subpop.sel, subpop.name)
subpop.which <- match(info.subpop, subpop.sel)
subpop.which[is.na(subpop.which)] <- 0
subpop.which <- as.logical(subpop.which)


subpop.list.now <- data.frame(info.iris[subpop.which], info.iris[subpop.which])
file.name <- paste0("raw_data/extra/", paste(subpop.sel, collapse = "_"), "_name_list.txt")
fwrite(subpop.list.now, file.name, col.names = FALSE, sep = " ")
