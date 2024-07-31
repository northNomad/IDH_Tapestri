wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

library(tapestri.tools)

# -------------------------------------------------------------------------
data_dir <- "../../DataFiles" #Data directory
## Files needed to run this script:
## **All these files must be in the same directory called 'data_dir'
## (1) 'D_Vehicle_23.11.30.dna+protein.h5'
## (2) 'D_IVO_23.11.30.dna+protein.h5'
## (3) 'D_ENA-DNA+Protein_10.25.23.dna+protein.h5'
## (4) 'D_IVO_ENA-DNA+Protein_23.10.25.dna+protein.h5'
##
##  These .h5 files are scDNA results exported from MisisonBio's Tapestri Pipeline.


## What the script outputs:
## (1): 'ls_sce_filtered_mixed.rds'
##       An R object. A list of single cell experiments containing filtered scDNA seq. results with clonal annotation.


## Read Tapestri .h5 files and store in list
h5_veh <- read_h5(file.path(data_dir, "CaseStudy2", "D_Vehicle_23.11.30.dna+protein.h5"))
h5_ivo <- read_h5(file.path(data_dir, "CaseStudy2", "D_IVO_23.11.30.dna+protein.h5"))
h5_ena <- read_h5(file.path(data_dir, "CaseStudy2", "D_ENA-DNA+Protein_10.25.23.dna+protein.h5"))
h5_ivo.ena <- read_h5(file.path(data_dir, "CaseStudy2", "D_IVO_ENA-DNA+Protein_23.10.25.dna+protein.h5"))

ls_h5 <- list(veh=h5_veh, ivo=h5_ivo, ena=h5_ena, ivo.ena=h5_ivo.ena)

## relevant variants
c(
  DNMT3A_R882C = "chr2:25457243:G/A",
  IDH2_R140Q = "chr15:90631934:C/T",
  IDH1_R132H = "chr2:209113112:C/T",
  NPM1c = "chr5:170837543:C/CTCTG",
  KRAS_K117N = "chr12:25378647:T/G", #From SKM1 spikein
  TP53_R248Q = "chr17:7577538:C/T", #From SKM1 spikein
  EZH2_Y641C = "chr7:148508727:T/C", #From SKM1 spikein
  FLT3_ITD = "chr13:28608262:T/TTCATATTCTCTGAAATCAACG"
) -> int_var

pathogenic_var <- int_var[c("DNMT3A_R882C", "IDH1_R132H", "NPM1c", "IDH2_R140Q")]
spikein <- int_var[c("KRAS_K117N", "TP53_R248Q", "EZH2_Y641C")]
## Read tapestri variants
ls_var <- lapply(ls_h5, function(x) get_variants(x, "data.table"))

## The FLT3ITD calls are not always consistent. Be careful.
## There are 6 FLT3ITD variants that passed the standard filter, likely arising from 
## two true variants, one for each patient.
ls_var %>% 
  lapply(function(x) x %>% 
           subset(filtered == "00" & SYMBOL == "FLT3") %>%
           subset(nchar(ALT) > 1) %>% 
           .$id) %>% 
  do.call(c, .) %>% 
  unique() -> FLT3ITD_tapestri

names(FLT3ITD_tapestri) <- paste0("FLT3ITD_", 1:length(FLT3ITD_tapestri))


# -------------------------------------------------------------------------


lapply(ls_h5[c("veh", "ivo", "ivo.ena", "ena")], function(x){
  read_assays_variants(x, 
                       c("AF", "DP", "GQ", "NGT"), 
                       index_variants = get_variants_index(x, c(pathogenic_var, FLT3ITD_tapestri)),
                       format = "SingleCellExperiment")
}) -> ls_sce

## Add in protein data 
for(i in names(ls_sce)){
  protein_raw <- read_protein_counts(ls_h5[[i]], format="matrix",  normalization = "raw")
  protein_clr <- read_protein_counts(ls_h5[[i]], format="matrix", normalization = "clr")
  
  #Store protein names
  names_protein_raw <- gsub(" ", "", rownames(protein_raw))
  names_protein_clr <- gsub(" ", "", rownames(protein_clr))
  
  #Remove attributes
  attributes(protein_raw) <- attributes(protein_raw)[1]
  attributes(protein_clr) <- attributes(protein_clr)[1]
  
  #Make s.e. object
  protein_raw <- SummarizedExperiment(list(Protein_raw=protein_raw), 
                                      rowData=DataFrame(Protein=names_protein_raw))
  protein_clr <- SummarizedExperiment(list(Protein_clr=protein_clr), 
                                      rowData=DataFrame(Protein=names_protein_clr))
  
  altExp(ls_sce[[i]], "Protein_raw") <- protein_raw
  altExp(ls_sce[[i]], "Protein_clr") <- protein_clr
}

# Filter - for IVO & VEH
# Keep only cells where (1) TP53_R248Q (2) KRAS_K117N" (3) "EZH2_Y641C" are all NGT = 0
# These three variants all come from SKM1 spikein
read_assays_variants(ls_h5[["ivo"]], "NGT", 
                     get_variants_index(ls_h5[["ivo"]], spikein),
                     format="list")[[1]] -> NGT_spikein_ivo
index_keep_ivo <- which(apply(NGT_spikein_ivo, 2, function(x) all(x == 0)))

read_assays_variants(ls_h5[["veh"]], "NGT", 
                     get_variants_index(ls_h5[["veh"]], spikein),
                     format="list")[[1]] -> NGT_spikein_veh
index_keep_veh <- which(apply(NGT_spikein_veh, 2, function(x) all(x == 0)))

ls_sce[["ivo"]] <- ls_sce[["ivo"]][, index_keep_ivo]
ls_sce[["veh"]] <- ls_sce[["veh"]][, index_keep_veh]

# Filter - remove NGT where GQ < 30

for(i in names(ls_sce)){
  low_GQ <- assays(ls_sce[[i]])[["GQ"]] < 30
  ngt <- assays(ls_sce[[i]])[["NGT"]]
  ngt[low_GQ] <- 3
  assays(ls_sce[[i]])[["NGT_filter_GQ30"]] <- ngt
}

# Filter - remove NGT where DP < 10

for(i in names(ls_sce)){
  low_DP <- assays(ls_sce[[i]])[["DP"]] < 10
  ngt <- assays(ls_sce[[i]])[["NGT_filter_GQ30"]]
  ngt[low_DP] <- 3
  assays(ls_sce[[i]])[["NGT_filter_GQ30_DP10"]] <- ngt
}

# # Collapse all 6 FLT3ITD variants
# 
for(i in names(ls_sce)){
  assays(ls_sce[[i]])[["NGT_filter_GQ30_DP10"]] %>%
    .[paste0("FLT3ITD_", 1:6), ] %>%
    apply(2, function(x) any(x == 1 | x == 2)) %>%
    which() -> index_FLT3ITD_positive
  
  assays(ls_sce[[i]])[["NGT_filter_GQ30_DP10"]]["FLT3ITD_1", ] <- 0
  assays(ls_sce[[i]])[["NGT_filter_GQ30_DP10"]]["FLT3ITD_1", ][index_FLT3ITD_positive] <- 1
  
  rn <- rownames(ls_sce[[i]])
  rn <- rn[!grepl("FLT3ITD", rn)]
  rn <- c(rn, "FLT3ITD_1")
  ls_sce[[i]] <- ls_sce[[i]][rn, ]
}

# Filter - remove cells where not all variants are genotyped

ls_sce_filtered <- list()
for(i in names(ls_sce)){
  ngt <- assays(ls_sce[[i]])[["NGT_filter_GQ30_DP10"]]
  index <- apply(ngt, 2, function(x) any(x == 3))
  index <- which(!index) #which cells have no NGT==3
  
  ls_sce_filtered[[i]] <- ls_sce[[i]][, index]
}

#Remove doublets:
#(1) Cells carry both DNMT3A_R882C & IDH2_R140Q
#(2) Cells carry both IDH1_R132H & IDH2_R140Q
ls_sce_filtered2 <- list()
for(i in names(ls_sce_filtered)){
  ngt <- assays(ls_sce_filtered[[i]])[["NGT_filter_GQ30_DP10"]]

  mut_IDH2_R140Q <- ngt["IDH2_R140Q", ] == 1 | ngt["IDH2_R140Q", ] == 2
  mut_IDH1_R132H <- ngt["IDH1_R132H", ] == 1 | ngt["IDH1_R132H", ] == 2
  mut_DNMT3A_R882C <- ngt["DNMT3A_R882C", ] == 1 | ngt["DNMT3A_R882C", ] == 2
  
  index_doublet <- (mut_IDH2_R140Q & mut_IDH1_R132H) | (mut_IDH2_R140Q & mut_DNMT3A_R882C)
  index_keep <- !index_doublet
  
  ls_sce_filtered2[[i]] <- ls_sce_filtered[[i]][, index_keep]
}

#How many cells pre and post filter?
lapply(ls_sce, ncol) %>% unlist()
lapply(ls_sce_filtered, ncol) %>% unlist()
lapply(ls_sce_filtered2, ncol) %>% unlist()

#
## Assign clone
#Clone0 - root
#Clone1 - DNMT3A only
#Clone2 - DI1
#Clone3 - DI1N
#Clone4 - DI1NF
#Clone5 - I2
#Clone6 - I2N
#Clone7 - I2NF
lapply(ls_sce_filtered2, function(x){
  assays(x)[["NGT_filter_GQ30_DP10"]] %>% 
    apply(c(1, 2), function(x) gsub(2, 1, x)) %>% 
    as.data.table() %>% 
    lapply(paste, collapse="") %>% 
    unlist()
  
}) -> ls_ngt

ls_clone <- list()

for(i in names(ls_ngt)){
  x <- ls_ngt[[i]]
  
  #Initiate matrix to store clone
  y <- matrix(ncol=1, nrow=length(ls_ngt[[i]]), dimnames = list(NULL, "Clone"))
  
  #get indices
  index_0 <- which(x == "00000") #root
  index_1 <- which(x == "10000") #DNMT3A only
  index_2 <- which(x == "11000" | x == "01000") #DI1
  index_3 <- which(x == "11100" | x == "10100" | x == "01100") #DI1N
  
  # index_4 <- which(x == "11101" | x == "10001" | x == "01001" | x == "011") #DI1NF
  index_4 <- str_sub(x, 5, 5) == "1" 
  index_4_2 <- str_sub_all(x, c(1, 2), c(1, 2)) %>% 
    lapply(function(x) sum(as.numeric(x))) %>% 
    unlist()
  index_4_2 <- index_4_2 > 0
  index_4 <- index_4 & index_4_2
  
  index_5 <- which(x == "00010") #I2
  index_6 <- which(x == "00110") #I2N
  index_7 <- which(x == "00111" | x == "00011") #I2NF
  
  y[index_0, "Clone"] <- 0
  y[index_1, "Clone"] <- 1
  y[index_2, "Clone"] <- 2
  y[index_3, "Clone"] <- 3
  y[index_4, "Clone"] <- 4
  y[index_5, "Clone"] <- 5
  y[index_6, "Clone"] <- 6
  y[index_7, "Clone"] <- 7
  
  ls_clone[[i]] <- y
}

#Add clone info
for(i in names(ls_clone)){
  colData(ls_sce_filtered2[[i]]) <- DataFrame(ls_clone[[i]])
}

#Remove ambiguous cells
for(i in names(ls_sce_filtered2)){
  index_rm <- is.na(colData(ls_sce_filtered2[[i]])$Clone)
  ls_sce_filtered2[[i]] <- ls_sce_filtered2[[i]][, !index_rm]
}

#Save
write_rds(ls_sce_filtered2, "ls_sce_filtered_mixed.rds")

# > sessionInfo()
# R version 4.4.0 (2024-04-24 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] tapestri.tools_0.0.9000     lubridate_1.9.3             forcats_1.0.0               stringr_1.5.1              
# [5] purrr_1.0.2                 readr_2.1.5                 tidyr_1.3.1                 tibble_3.2.1               
# [9] tidyverse_2.0.0             SingleCellExperiment_1.26.0 SummarizedExperiment_1.34.0 Biobase_2.64.0             
# [13] MatrixGenerics_1.16.0       matrixStats_1.3.0           rtracklayer_1.64.0          rhdf5_2.48.0               
# [17] rentrez_1.2.3               ggridges_0.5.6              ggplot2_3.5.1               GenomicRanges_1.56.0       
# [21] GenomeInfoDb_1.40.1         IRanges_2.38.0              S4Vectors_0.42.0            BiocGenerics_0.50.0        
# [25] magrittr_2.0.3              maftools_2.20.0             dplyr_1.1.4                 data.table_1.15.4          
# [29] compositions_2.0-8         
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.1         Biostrings_2.72.0        bitops_1.0-7             RCurl_1.98-1.14         
# [5] tensorA_0.36.2.1         GenomicAlignments_1.40.0 XML_3.99-0.16.1          timechange_0.3.0        
# [9] lifecycle_1.0.4          survival_3.5-8           compiler_4.4.0           rlang_1.1.3             
# [13] tools_4.4.0              utf8_1.2.4               yaml_2.3.8               S4Arrays_1.4.1          
# [17] curl_5.2.1               DelayedArray_0.30.1      RColorBrewer_1.1-3       abind_1.4-5             
# [21] BiocParallel_1.38.0      withr_3.0.0              grid_4.4.0               fansi_1.0.6             
# [25] colorspace_2.1-0         Rhdf5lib_1.26.0          scales_1.3.0             MASS_7.3-60.2           
# [29] cli_3.6.2                crayon_1.5.2             generics_0.1.3           rstudioapi_0.16.0       
# [33] robustbase_0.99-2        httr_1.4.7               tzdb_0.4.0               rjson_0.2.21            
# [37] bayesm_3.1-6             DNAcopy_1.78.0           zlibbioc_1.50.0          splines_4.4.0           
# [41] parallel_4.4.0           XVector_0.44.0           restfulr_0.0.15          vctrs_0.6.5             
# [45] Matrix_1.7-0             jsonlite_1.8.8           hms_1.1.3                glue_1.7.0              
# [49] DEoptimR_1.1-3           codetools_0.2-20         stringi_1.8.4            gtable_0.3.5            
# [53] BiocIO_1.14.0            UCSC.utils_1.0.0         munsell_0.5.1            pillar_1.9.0            
# [57] rhdf5filters_1.16.0      GenomeInfoDbData_1.2.12  R6_2.5.1                 lattice_0.22-6          
# [61] Rsamtools_2.20.0         Rcpp_1.0.12              SparseArray_1.4.8        pkgconfig_2.0.3        
