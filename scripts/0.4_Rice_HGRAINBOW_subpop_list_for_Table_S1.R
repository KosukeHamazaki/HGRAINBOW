#########################################################
###### 0.4_Rice_HGRAINBOW_subpop_list_for_Table_S1 ######
#########################################################

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
accession.all.list <- read.csv("raw_data/extra/L3024_404K_core_SNP.full_data.name.csv")

accession.list <- accession.all.list[as.character(accession.all.list$SUBPOPULATION) %in% subpop.sel, ] 


accession.list.modi <- data.frame(IRIS_ID = accession.list[, 2], Vairety_name = accession.list[, 1],
                                  Subpopulation = accession.list[, 4],
                                  Country = accession.list[, 5])


write.csv(accession.list.modi, "data/extra/Table_S1_origin.csv")
