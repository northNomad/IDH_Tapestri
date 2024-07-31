wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

library(tidyverse)
library(data.table)
library(patchwork)
library(tapestri.tools)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggtext)
library(EnhancedVolcano)
library(ggsci)

# -------------------------------------------------------------------------
## What is needed to run the script: 
## (1): 'ls_sce_filtered.rds' 
##       List of filtered tapestri data in single cell experiment format.
##       It is generated from 'Script1_QC.R'
ls_sce_filtered <- read_rds("../Fig1E/ls_sce_filtered.rds")

## What the script outputs:
## (1): 'p_volcano_ivo_veh.png' 
##       This is Figure 2A
## (2): 'hm_protein_v_clone.svg'
##       This is Figure 2C
## (3): 'p_CD117_CD69_ridge.png'
##       This is Figure 2B

# -------------------------------------------------------------------------
map_treatment <- c(s1="Primary_sample", s4="VEH", s6="IVO")

# Keep only the samples with protein data
protein_samples <- c("s1", "s4", "s6")
ls_sce_filtered <- ls_sce_filtered[protein_samples]

# -------------------------------------------------------------------------
## Normalize protein by CLR
ls_sce_filtered %>% 
  lapply(function(x){
    prot <- rowData(altExp(x, "Protein_clr"))$Protein
    prot <- gsub("-", "_", prot)
    d <- assays(altExp(x, "Protein_clr"))[[1]]
    rownames(d) <- prot
    d <- as.data.table(t(d))
    d$Clone <- colData(x)$Clone_annotated
    d
  }) %>% 
  bind_rows(.id="Sample") -> d_clr

d_clr$clone <- unlist(lapply(ls_sce_filtered, function(x) colData(x)$Clone_annotated))

proteins <- colnames(d_clr)[2:44]

## Total effect of ivo on immunophenotype (pseudobulk)

ls_result <- list()
d_clr_sub <- subset(d_clr, Sample != "s1")
for(i in proteins){
  f <- as.formula(paste0(i, "~ Sample"))
  lm <- lm(f, d_clr_sub)
  ls_result[[i]] <- broom::tidy(lm)
}

ls_result %>% 
  bind_rows(.id="protein") %>% 
  mutate(fdr=p.adjust(p.value, method="BH")) -> d_result

d_result$term <- gsub("\\(Intercept\\)", "Intercept", d_result$term)

d_result$significant <- ifelse(d_result$fdr < .05, "yes", "no")

d_result %>% 
  subset(term == "Samples6") %>% 
  ggplot(aes(estimate, -log10(fdr), color=significant)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label=protein)) +
  scale_color_manual(values=c(yes="firebrick", no="grey70"), name="FDR<.05",
                     guide="none") +
  geom_hline(yintercept = -log10(0.05)) +
  labs(x="", y="-log10 ( FDR )") +
  theme_bw(15) +
  theme(plot.title.position = "plot") -> p_volcano_ivo_veh
p_volcano_ivo_veh
ggsave("p_volcano_ivo_veh.png", p_volcano_ivo_veh, width = 5.5, height = 4.5, units = "in",
       dpi="retina")

# -------------------------------------------------------------------------

## Calculate mean normalized expression for each clone across treatment
d_clr %>% 
  pivot_longer(2:44, names_to="Protein", values_to="CLR") %>% 
  group_by(Sample, Clone, Protein) %>% 
  dplyr::summarise(mu = mean(CLR)) %>% 
  pivot_wider(id_cols=c("Sample", "Protein"), names_from="Clone", values_from="mu") -> dmu

dmu <- dmu[, c("Sample", "Protein", paste0("Clone", 3:9))]

dmu_ivo <- dmu %>% subset(Sample == "s6")
dmu_veh <- dmu %>% subset(Sample == "s4")

## Show differentially expressed proteins in heatmap
hm_proteins <- d_result %>% 
  subset(term == "Samples6") %>% 
  subset(significant == "yes") %>% 
  .$protein

m_fc <- subset(dmu_ivo, Protein %in% hm_proteins)[, -c(1, 2)] - subset(dmu_veh, Protein %in% hm_proteins)[, -c(1, 2)]
m_fc <- t(m_fc)
colnames(m_fc) <- subset(dmu_ivo, Protein %in% hm_proteins)$Protein
m_fc <- as.matrix(m_fc)
# 
# d_result %>% 
#   subset(term == "Samples6") %>% 
#   subset(protein %in% colnames(m_fc)) %>% 
#   dplyr::select(protein, estimate) -> fc_bulk
# 
# m_fc_bulk <- matrix(fc_bulk$estimate, nrow=1)
# colnames(m_fc_bulk) <- fc_bulk$protein
# rownames(m_fc_bulk) <- "Pseudobulk"
# 
# m_fc <- rbind(m_fc_bulk, m_fc)
# fc_bulk %>% 
#   arrange(-estimate) %>% 
#   .$protein -> protein_order

m_fc %>% 
  Heatmap(cluster_rows = FALSE, 
        cluster_columns = TRUE, 
        row_names_side = "left") -> hm_protein_v_clone

svg("hm_protein_v_clone.svg", height = 3.5, width = 4)
hm_protein_v_clone
dev.off()

# -------------------------------------------------------------------------

d_clr_sub %>% 
  pivot_longer(cols=c(CD117, CD69), names_to="Protein", values_to="CLR") %>% 
  ggplot(aes(CLR, Sample)) +
  geom_density_ridges(aes(fill=Sample)) +
  facet_wrap(~Protein) +
  scale_y_discrete(label=c("VEH", "IVO")) +
  labs(y=NULL) +
  scale_fill_npg(guide="none") +
  theme_classic(15) -> p_CD117_CD69_ridge
p_CD117_CD69_ridge
ggsave("p_CD117_CD69_ridge.png", p_CD117_CD69_ridge, 
       width = 5.3, height = 3.4, units = "in", dpi="retina")

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
#   [1] ggsci_3.1.0                 ggtext_0.1.2                EnhancedVolcano_1.22.0     
# [4] ggrepel_0.9.5               RColorBrewer_1.1-3          ComplexHeatmap_2.20.0      
# [7] tapestri.tools_0.0.9000     SingleCellExperiment_1.26.0 SummarizedExperiment_1.34.0
# [10] Biobase_2.64.0              MatrixGenerics_1.16.0       matrixStats_1.3.0          
# [13] rtracklayer_1.64.0          rhdf5_2.48.0                rentrez_1.2.3              
# [16] ggridges_0.5.6              GenomicRanges_1.56.0        GenomeInfoDb_1.40.1        
# [19] IRanges_2.38.0              S4Vectors_0.42.0            BiocGenerics_0.50.0        
# [22] magrittr_2.0.3              maftools_2.20.0             compositions_2.0-8         
# [25] patchwork_1.2.0             data.table_1.15.4           lubridate_1.9.3            
# [28] forcats_1.0.0               stringr_1.5.1               dplyr_1.1.4                
# [31] purrr_1.0.2                 readr_2.1.5                 tidyr_1.3.1                
# [34] tibble_3.2.1                ggplot2_3.5.1               tidyverse_2.0.0            
# 
# loaded via a namespace (and not attached):
#   [1] tensorA_0.36.2.1         rstudioapi_0.16.0        jsonlite_1.8.8           shape_1.4.6.1           
# [5] farver_2.1.2             ragg_1.3.2               GlobalOptions_0.1.2      BiocIO_1.14.0           
# [9] zlibbioc_1.50.0          vctrs_0.6.5              Rsamtools_2.20.0         RCurl_1.98-1.14         
# [13] S4Arrays_1.4.1           curl_5.2.1               broom_1.0.6              cellranger_1.1.0        
# [17] Rhdf5lib_1.26.0          SparseArray_1.4.8        StanHeaders_2.32.8       plyr_1.8.9              
# [21] GenomicAlignments_1.40.0 commonmark_1.9.1         lifecycle_1.0.4          iterators_1.0.14        
# [25] pkgconfig_2.0.3          Matrix_1.7-0             R6_2.5.1                 GenomeInfoDbData_1.2.12 
# [29] clue_0.3-65              digest_0.6.35            colorspace_2.1-0         textshaping_0.4.0       
# [33] labeling_0.4.3           fansi_1.0.6              timechange_0.3.0         httr_1.4.7              
# [37] abind_1.4-5              compiler_4.4.0           withr_3.0.0              doParallel_1.0.17       
# [41] backports_1.5.0          inline_0.3.19            BiocParallel_1.38.0      QuickJSR_1.1.3          
# [45] pkgbuild_1.4.4           MASS_7.3-60.2            bayesm_3.1-6             DelayedArray_0.30.1     
# [49] rjson_0.2.21             loo_2.7.0                DNAcopy_1.78.0           tools_4.4.0             
# [53] glue_1.7.0               restfulr_0.0.15          rhdf5filters_1.16.0      gridtext_0.1.5          
# [57] cluster_2.1.6            generics_0.1.3           gtable_0.3.5             tzdb_0.4.0              
# [61] hms_1.1.3                xml2_1.3.6               utf8_1.2.4               XVector_0.44.0          
# [65] foreach_1.5.2            pillar_1.9.0             markdown_1.13            robustbase_0.99-2       
# [69] circlize_0.4.16          splines_4.4.0            lattice_0.22-6           survival_3.5-8          
# [73] tidyselect_1.2.1         Biostrings_2.72.0        gridExtra_2.3            V8_4.4.2                
# [77] xfun_0.44                DEoptimR_1.1-3           rstan_2.32.6             stringi_1.8.4           
# [81] UCSC.utils_1.0.0         yaml_2.3.8               codetools_0.2-20         cli_3.6.2               
# [85] RcppParallel_5.1.7       systemfonts_1.1.0        munsell_0.5.1            Rcpp_1.0.12             
# [89] readxl_1.4.3             png_0.1-8                XML_3.99-0.16.1          parallel_4.4.0          
# [93] bitops_1.0-7             scales_1.3.0             crayon_1.5.2             GetoptLong_1.0.5        
# [97] rlang_1.1.3             