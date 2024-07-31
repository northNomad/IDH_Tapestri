wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

library(tidyverse)
library(data.table)
library(tapestri.tools)
library(patchwork)
library(ComplexHeatmap)

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
        name = "Z-Scaled proportion"
        ) -> p_hm_ClonesAcrossSamples_v2
p_hm_ClonesAcrossSamples_v2

svg(filename = "p_hm_ClonesAcrossSamples.svg", width = 8, height = 5)
p_hm_ClonesAcrossSamples_v2
dev.off()


# -------------------------------------------------------------------------
# > sessionInfo()
# R version 4.1.3 (2022-03-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
# [3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.1252    
# 
# attached base packages:
#   [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ComplexHeatmap_2.8.0        patchwork_1.1.1             tapestri.tools_0.0.9000     SingleCellExperiment_1.14.1
# [5] SummarizedExperiment_1.24.0 Biobase_2.54.0              MatrixGenerics_1.6.0        matrixStats_0.63.0         
# [9] rtracklayer_1.52.1          rhdf5_2.36.0                rentrez_1.2.3               ggridges_0.5.3             
# [13] GenomicRanges_1.46.1        GenomeInfoDb_1.30.1         IRanges_2.28.0              S4Vectors_0.32.4           
# [17] BiocGenerics_0.40.0         magrittr_2.0.3              maftools_2.8.05             compositions_2.0-4         
# [21] data.table_1.14.6           forcats_0.5.1               stringr_1.5.0               dplyr_1.0.10               
# [25] purrr_0.3.4                 readr_2.1.1                 tidyr_1.2.1                 tibble_3.1.8               
# [29] ggplot2_3.4.0               tidyverse_1.3.1            
# 
# loaded via a namespace (and not attached):
#   [1] colorspace_2.0-3         rjson_0.2.20             ellipsis_0.3.2           circlize_0.4.15         
# [5] XVector_0.34.0           GlobalOptions_0.1.2      fs_1.5.2                 clue_0.3-60             
# [9] rstudioapi_0.13          fansi_1.0.3              lubridate_1.8.0          xml2_1.3.3              
# [13] codetools_0.2-18         splines_4.1.3            doParallel_1.0.17        robustbase_0.95-0       
# [17] jsonlite_1.8.4           Rsamtools_2.8.0          Cairo_1.5-15             broom_1.0.3             
# [21] cluster_2.1.2            dbplyr_2.2.1             png_0.1-8                compiler_4.1.3          
# [25] httr_1.4.4               backports_1.4.1          assertthat_0.2.1         Matrix_1.4-0            
# [29] cli_3.4.1                tools_4.1.3              gtable_0.3.1             glue_1.6.2              
# [33] GenomeInfoDbData_1.2.7   Rcpp_1.0.9               cellranger_1.1.0         vctrs_0.5.1             
# [37] Biostrings_2.62.0        rhdf5filters_1.4.0       iterators_1.0.14         tensorA_0.36.2          
# [41] rvest_1.0.2              lifecycle_1.0.3          restfulr_0.0.15          XML_3.99-0.13           
# [45] DEoptimR_1.0-11          zlibbioc_1.40.0          MASS_7.3-55              scales_1.2.1            
# [49] hms_1.1.2                parallel_4.1.3           RColorBrewer_1.1-3       yaml_2.3.6              
# [53] stringi_1.7.8            BiocIO_1.2.0             foreach_1.5.2            BiocParallel_1.26.2     
# [57] shape_1.4.6              rlang_1.0.6              pkgconfig_2.0.3          bitops_1.0-7            
# [61] lattice_0.20-45          Rhdf5lib_1.14.2          GenomicAlignments_1.28.0 tidyselect_1.2.0        
# [65] R6_2.5.1                 generics_0.1.3           DelayedArray_0.20.0      DBI_1.1.3               
# [69] pillar_1.8.1             haven_2.4.3              withr_2.5.0              survival_3.2-13         
# [73] RCurl_1.98-1.9           bayesm_3.1-4             modelr_0.1.8             crayon_1.5.2            
# [77] utf8_1.2.2               tzdb_0.2.0               GetoptLong_1.0.5         readxl_1.3.1            
# [81] reprex_2.0.1             digest_0.6.31            munsell_0.5.0           


