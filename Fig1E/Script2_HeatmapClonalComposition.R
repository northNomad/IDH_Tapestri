wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

library(tidyverse)
library(data.table)
library(tapestri.tools)
library(patchwork)
library(ComplexHeatmap)
library(circlize)

ls_sce_filtered <- read_rds("ls_sce_filtered.rds")

lapply(ls_sce_filtered, function(x){
  colData(x) %>% 
    as.data.table() %>% 
    .$Clone_annotated %>% 
    table() %>% 
    as.data.table() -> x
  colnames(x) <- c("Clone", "N")
  x
}) %>% 
  bind_rows(.id = "Sample") -> d_clone

d_clone %>% 
  pivot_wider(Sample, names_from = "Clone", values_from = "N") -> d2

d2[is.na(d2)] <- 0
d2$Clone1 <- 0
d2 <- as.data.frame(d2)

#Order rows and columns
row_order <- c("s1", "s2", "s4", "s6", "s12", "s14", "s8", "s10")
col_order <- c("Sample", paste0("Clone", 0:9))

rownames(d2) <- d2$Sample
d2[row_order, col_order] %>% 
  .[, -1] %>% 
  t() %>% 
  scale() %>%  
  t() -> d2_scale

d2 <- d2[row_order, col_order]
d_number <- d2[, -1]

d_percent <- d2[, -1]
d_percent$sum <- rowSums(d_percent)
d_percent[, 1:10] <- d_percent[, 1:10] / d_percent[, "sum"]
d_percent <- round(d_percent[, 1:10]*100, 2) 

d_label <- matrix(NA, nrow=nrow(d_percent), ncol=ncol(d_percent))
for(i in 1:nrow(d_percent)){
  for(j in 1:ncol(d_percent)){
    label <- paste0(d_number[i, j], "\n", "(", d_percent[i, j], ")")
    d_label[i, j] <- label
  }
}

d2_scale <- t(d2_scale)
d_label <- t(d_label)
col_fun <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
Heatmap(d2_scale, cluster_rows = FALSE, cluster_columns = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(d_label[i, j], x, y, gp = gpar(fontsize = 10))
        },
        column_labels = c("Primary sample", "Baseline",
                          "Vehicle", "Ivosidenib", "Venetoclax", "Ivo. + Ven.",
                          "Azacitidine", "Ivo. + Aza."),
        column_names_rot = -30, 
        row_names_side = "left", 
        column_names_side = "top",
        name = "Z-Scaled proportion",
        col = col_fun
        ) -> p_hm_ClonesAcrossSamples_v3
p_hm_ClonesAcrossSamples_v3

svg(filename = "p_hm_ClonesAcrossSamples.svg", width = 8, height = 5)
p_hm_ClonesAcrossSamples_v3
dev.off()



# -------------------------------------------------------------------------
# 
# sessionInfo()
# R version 4.4.0 (2024-04-24 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
# [3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.utf8    
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
#   [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] circlize_0.4.16             ComplexHeatmap_2.20.0       patchwork_1.2.0            
# [4] tapestri.tools_0.0.9000     SingleCellExperiment_1.26.0 SummarizedExperiment_1.34.0
# [7] Biobase_2.64.0              MatrixGenerics_1.16.0       matrixStats_1.3.0          
# [10] rtracklayer_1.64.0          rhdf5_2.48.0                rentrez_1.2.3              
# [13] ggridges_0.5.6              GenomicRanges_1.56.0        GenomeInfoDb_1.40.1        
# [16] IRanges_2.38.0              S4Vectors_0.42.0            BiocGenerics_0.50.0        
# [19] magrittr_2.0.3              maftools_2.20.0             compositions_2.0-8         
# [22] data.table_1.15.4           lubridate_1.9.3             forcats_1.0.0              
# [25] stringr_1.5.1               dplyr_1.1.4                 purrr_1.0.2                
# [28] readr_2.1.5                 tidyr_1.3.1                 tibble_3.2.1               
# [31] ggplot2_3.5.1               tidyverse_2.0.0            
# 
# loaded via a namespace (and not attached):
#   [1] bitops_1.0-7             rlang_1.1.3              clue_0.3-65              GetoptLong_1.0.5        
# [5] compiler_4.4.0           png_0.1-8                vctrs_0.6.5              shape_1.4.6.1           
# [9] pkgconfig_2.0.3          crayon_1.5.2             fastmap_1.2.0            XVector_0.44.0          
# [13] utf8_1.2.4               Rsamtools_2.20.0         rmarkdown_2.27           tzdb_0.4.0              
# [17] UCSC.utils_1.0.0         xfun_0.44                zlibbioc_1.50.0          jsonlite_1.8.8          
# [21] rhdf5filters_1.16.0      DelayedArray_0.30.1      Rhdf5lib_1.26.0          BiocParallel_1.38.0     
# [25] cluster_2.1.6            parallel_4.4.0           R6_2.5.1                 stringi_1.8.4           
# [29] RColorBrewer_1.1-3       pagedown_0.20            DNAcopy_1.78.0           pkgload_1.3.4           
# [33] Rcpp_1.0.12              iterators_1.0.14         knitr_1.46               Matrix_1.7-0            
# [37] splines_4.4.0            timechange_0.3.0         tidyselect_1.2.1         rstudioapi_0.16.0       
# [41] abind_1.4-5              yaml_2.3.8               doParallel_1.0.17        codetools_0.2-20        
# [45] curl_5.2.1               lattice_0.22-6           withr_3.0.0              evaluate_0.23           
# [49] survival_3.5-8           bayesm_3.1-6             Biostrings_2.72.0        pillar_1.9.0            
# [53] tensorA_0.36.2.1         foreach_1.5.2            generics_0.1.3           RCurl_1.98-1.14         
# [57] hms_1.1.3                munsell_0.5.1            scales_1.3.0             glue_1.7.0              
# [61] tools_4.4.0              BiocIO_1.14.0            robustbase_0.99-2        GenomicAlignments_1.40.0
# [65] XML_3.99-0.16.1          colorspace_2.1-0         GenomeInfoDbData_1.2.12  restfulr_0.0.15         
# [69] cli_3.6.2                fansi_1.0.6              S4Arrays_1.4.1           gtable_0.3.5            
# [73] DEoptimR_1.1-3           digest_0.6.35            SparseArray_1.4.8        rjson_0.2.21            
# [77] htmltools_0.5.8.1        lifecycle_1.0.4          httr_1.4.7               GlobalOptions_0.1.2     
# [81] MASS_7.3-60.2           
