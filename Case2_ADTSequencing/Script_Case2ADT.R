wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

library(tapestri.tools)
library(rethinking)
library(ggrepel)
library(ggdensity)
library(cmdstanr)
library(patchwork)
library(ggsci)
library(ggtext)
library(ComplexHeatmap)
set_cmdstan_path()

# -------------------------------------------------------------------------
## What is needed to run the script: 
## (1): 'ls_sce_filtered_mixed.rds' 
##       List of filtered tapestri data in single cell experiment format.
##       It is generated from 'Script_Case2QC.R'
ls_sce_filtered <- read_rds("../Case2_QC/ls_sce_filtered_mixed.rds")

## What the script outputs:
## (1): 'p_pseudobulk_volcano.png' 
##       This is Figure 4F
## (2): 'p_hm_idh1_bulk.svg'
## (3): 'p_hm_idh1_clone3.svg'
## (4): 'p_hm_idh1_clone4.svg'
## (5): 'p_hm_idh2_bulk.svg'
## (6): 'p_hm_idh2_clone5.svg'
## (7): 'p_hm_idh2_clone6.svg'
## (8): 'p_hm_idh2_clone7.svg'
##       These are Figure 4G.

ls_sce_filtered %>% 
  lapply(function(x){
    prot <- rowData(altExp(x, "Protein_clr"))$Protein
    prot <- gsub("-", "_", prot)
    d <- assays(altExp(x, "Protein_clr"))[[1]]
    rownames(d) <- prot
    d <- as.data.table(t(d))
    d$Clone <- colData(x)$Clone
    d
  }) %>% 
  bind_rows(.id="Sample") -> d_clr

##Annotate sample
d_clr$Patient <- plyr::mapvalues(d_clr$Clone, 0:7,
                                 rep(c("WT", "IDH1", "IDH2"), times=c(1, 4, 3))
)

d_clr$Clone <- paste0("Clone", d_clr$Clone)

d_clr$Sample <- factor(d_clr$Sample, levels=c("veh", "ivo", "ena", "ivo.ena"))

# -------------------------------------------------------------------------
proteins <- colnames(d_clr)[2:44]

d_clr_idh1 <- subset(d_clr, Patient == "IDH1")
d_clr_idh2 <- subset(d_clr, Patient == "IDH2")

ls_result_idh1 <- list()
ls_result_idh2 <- list()

## pseudobulk for IDH1
for(i in proteins){
  f <- as.formula(paste0(i, "~ Sample"))
  lm <- lm(f, d_clr_idh1)
  ls_result_idh1[[i]] <- broom::tidy(lm)
}

## pseudobulk for IDH2
for(i in proteins){
  f <- as.formula(paste0(i, "~ Sample"))
  lm <- lm(f, d_clr_idh2)
  ls_result_idh2[[i]] <- broom::tidy(lm)
}

ls_result_idh1 <- bind_rows(ls_result_idh1, .id="Protein")
ls_result_idh2 <- bind_rows(ls_result_idh2, .id="Protein")


ls_result_idh1$FDR <- p.adjust(ls_result_idh1$p.value, method = "BH")
ls_result_idh2$FDR <- p.adjust(ls_result_idh2$p.value, method = "BH")

ls_result_idh1$Significant <- ifelse(ls_result_idh1$FDR < .05 & 
                                       abs(ls_result_idh1$estimate) > .2, 
                                     "Yes", "No")
ls_result_idh2$Significant <- ifelse(ls_result_idh2$FDR < .05 & 
                                       abs(ls_result_idh2$estimate) > .2
                                     , "Yes", "No")


ls_result_idh1 %>% 
  subset(term == "Sampleivo") %>% 
  ggplot(aes(estimate, -log10(FDR))) +
  geom_point(aes(color=Significant), size=3) + 
  geom_text_repel(aes(label=Protein)) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = .2) +
  geom_vline(xintercept = -.2) +
  scale_color_manual(values=c(Yes="firebrick", No="grey70"), guide="none") +
  labs(x="Diff. in mean normalized counts", y="-log10(FDR)", title="IVO v. VEH") +
  xlim(c(-2, 2)) +
  theme_bw(13) +
  theme(plot.title.position = "plot") -> vol_idh1_ivo_veh

ls_result_idh1 %>% 
  subset(term == "Sampleena") %>% 
  ggplot(aes(estimate, -log10(FDR))) +
  geom_point(aes(color=Significant), size=3) + 
  geom_text_repel(aes(label=Protein)) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = .2) +
  geom_vline(xintercept = -.2) +
  scale_color_manual(values=c(Yes="firebrick", No="grey70"), guide="none") +
  labs(x="Diff. in mean normalized counts", y="-log10(FDR)", title="ENA v. VEH") +
  xlim(c(-2, 2)) +
  theme_bw(13) +
  theme(plot.title.position = "plot") -> vol_idh1_ena_veh

ls_result_idh1 %>% 
  subset(term == "Sampleivo.ena") %>% 
  ggplot(aes(estimate, -log10(FDR))) +
  geom_point(aes(color=Significant), size=3) + 
  geom_text_repel(aes(label=Protein)) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = .2) +
  geom_vline(xintercept = -.2) +
  scale_color_manual(values=c(Yes="firebrick", No="grey70"), guide="none") +
  labs(x="Diff. in mean normalized counts", y="-log10(FDR)", title="IVO+ENA v. VEH") +
  xlim(c(-2, 2)) +
  theme_bw(13) +
  theme(plot.title.position = "plot") -> vol_idh1_ivo.ena_veh


# -------------------------------------------------------------------------

ls_result_idh2 %>% 
  subset(term == "Sampleivo") %>% 
  ggplot(aes(estimate, -log10(FDR))) +
  geom_point(aes(color=Significant), size=3) + 
  geom_text_repel(aes(label=Protein)) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = .2) +
  geom_vline(xintercept = -.2) +
  scale_color_manual(values=c(Yes="firebrick", No="grey70"), guide="none") +
  labs(x="Diff. in mean normalized counts", y="-log10(FDR)", title="IVO v. VEH") +
  xlim(c(-2, 2)) +
  theme_bw(13) +
  theme(plot.title.position = "plot") -> vol_idh2_ivo_veh

ls_result_idh2 %>% 
  subset(term == "Sampleena") %>% 
  ggplot(aes(estimate, -log10(FDR))) +
  geom_point(aes(color=Significant), size=3) + 
  geom_text_repel(aes(label=Protein)) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = .2) +
  geom_vline(xintercept = -.2) +
  scale_color_manual(values=c(Yes="firebrick", No="grey70"), guide="none") +
  labs(x="Diff. in mean normalized counts", y="-log10(FDR)", title="ENA v. VEH") +
  xlim(c(-2, 2)) +
  theme_bw(13) +
  theme(plot.title.position = "plot") -> vol_idh2_ena_veh

ls_result_idh2 %>% 
  subset(term == "Sampleivo.ena") %>% 
  ggplot(aes(estimate, -log10(FDR))) +
  geom_point(aes(color=Significant), size=3) + 
  geom_text_repel(aes(label=Protein)) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = .2) +
  geom_vline(xintercept = -.2) +
  scale_color_manual(values=c(Yes="firebrick", No="grey70"), guide="none") +
  labs(x="Diff. in mean normalized counts", y="-log10(FDR)", title="IVO+ENA v. VEH") +
  xlim(c(-2, 2)) +
  theme_bw(13) +
  theme(plot.title.position = "plot") -> vol_idh2_ivo.ena_veh


# -------------------------------------------------------------------------

(vol_idh1_ivo_veh + vol_idh1_ena_veh + vol_idh1_ivo.ena_veh) / 
  (vol_idh2_ivo_veh + vol_idh2_ena_veh + vol_idh2_ivo.ena_veh) -> p_pseudobulk_volcano

ggsave("p_pseudobulk_volcano.png", p_pseudobulk_volcano, width = 11, height = 5,
       units = "in", dpi="retina")

# -------------------------------------------------------------------------
## Make Heatmap legend same scale
col_fun = circlize::colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))

## Find significantly DEP for each patient
dep_idh1 <- ls_result_idh1 %>% subset(FDR < .05) %>% .$Protein %>% unique()
dep_idh2 <- ls_result_idh2 %>% subset(FDR < .05) %>% .$Protein %>% unique()

## Calculate by mean for each clone
d_clr %>% 
  pivot_longer(2:44, names_to="Protein", values_to="CLR") %>% 
  group_by(Sample, Clone, Protein) %>% 
  dplyr::summarise(mu=mean(CLR)) %>% 
  pivot_wider(id_cols=c(Sample, Clone), names_from="Protein", values_from="mu") -> d_mu_wide

## Order IDh1 by FC with IVO
## Calculate mean without grouping by clone for top heatmap
d_clr_idh1 %>% 
  pivot_longer(2:44, names_to="Protein", values_to="CLR") %>% 
  group_by(Sample, Protein) %>% 
  dplyr::summarise(mu=mean(CLR)) -> d_mu_bulk_idh1

d_clr_idh2 %>% 
  pivot_longer(2:44, names_to="Protein", values_to="CLR") %>% 
  group_by(Sample, Protein) %>% 
  dplyr::summarise(mu=mean(CLR)) -> d_mu_bulk_idh2

## Make them wider for heatmap
d_mu_bulk_idh1 %>% 
  pivot_wider(id_cols="Sample", names_from="Protein", values_from="mu") -> d_mu_bulk_idh1
d_mu_bulk_idh2 %>% 
  pivot_wider(id_cols="Sample", names_from="Protein", values_from="mu") -> d_mu_bulk_idh2

## Pseudobulk

d_mu_bulk_idh1[, -1] %>% 
  .[, dep_idh1] -> hm_fc_bulk_idh1
hm_fc_bulk_idh1[2, ] <- hm_fc_bulk_idh1[2, ] - hm_fc_bulk_idh1[1, ]
hm_fc_bulk_idh1[3, ] <- hm_fc_bulk_idh1[3, ] - hm_fc_bulk_idh1[1, ]
hm_fc_bulk_idh1[4, ] <- hm_fc_bulk_idh1[4, ] - hm_fc_bulk_idh1[1, ]

## Make column order the same as pseudobullk
hm_fc_bulk_idh1[2, ] %>% 
  t() %>% 
  as.data.frame() %>% 
  subset(abs(V1) >= .2) %>% 
  arrange(-V1) %>% 
  rownames() -> col_order_fc_slim_idh1

hm_fc_bulk_idh1[-1, col_order_fc_slim_idh1] %>% 
  Heatmap(cluster_rows = FALSE,
          cluster_columns = FALSE,
          row_labels = c("IVO", "ENA", "IVO+ENA"),
          show_column_names = TRUE,
          column_names_side = "top",
          column_title = "mIDH1 Pseudobulk",
          col=col_fun
  ) -> p_hm_idh1_bulk

## Clone 3
d_mu_wide %>% 
  subset(Clone == "Clone3") %>% 
  .[, col_order_fc_slim_idh1] -> hm_fc_clone3

hm_fc_clone3[2, ] <- hm_fc_clone3[2, ] - hm_fc_clone3[1, ]
hm_fc_clone3[3, ] <- hm_fc_clone3[3, ] - hm_fc_clone3[1, ]
hm_fc_clone3[1, ] <- 0

hm_fc_clone3[, col_order_fc_slim_idh1] %>% 
  Heatmap(cluster_rows = FALSE,
          cluster_columns = FALSE,
          row_labels = c("IVO", "ENA", "IVO+ENA"),
          show_column_names = FALSE,
          column_title = "Clone3",
          col=col_fun
  ) -> p_hm_idh1_clone3

## Clone4
d_mu_wide %>% 
  subset(Clone == "Clone4") %>% 
  .[, col_order_fc_slim_idh1] -> hm_fc_clone4

hm_fc_clone4[2, ] <- hm_fc_clone4[2, ] - hm_fc_clone4[1, ]
hm_fc_clone4[3, ] <- hm_fc_clone4[3, ] - hm_fc_clone4[1, ]
hm_fc_clone4[4, ] <- hm_fc_clone4[4, ] - hm_fc_clone4[1, ]

hm_fc_clone4[-1, col_order_fc_slim_idh1] %>% 
  Heatmap(cluster_rows = FALSE,
          cluster_columns = FALSE,
          row_labels = c("IVO", "ENA", "IVO+ENA"),
          show_column_names = FALSE,
          column_title = "Clone4",
          col=col_fun
  ) -> p_hm_idh1_clone4

## Order IDh2 by FC with ena
d_mu_bulk_idh2[, -1] %>% 
  .[, dep_idh2] -> hm_fc_bulk_idh2
hm_fc_bulk_idh2[2, ] <- hm_fc_bulk_idh2[2, ] - hm_fc_bulk_idh2[1, ]
hm_fc_bulk_idh2[3, ] <- hm_fc_bulk_idh2[3, ] - hm_fc_bulk_idh2[1, ]
hm_fc_bulk_idh2[4, ] <- hm_fc_bulk_idh2[4, ] - hm_fc_bulk_idh2[1, ]

hm_fc_bulk_idh2[3, ] %>% 
  t() %>% 
  as.data.frame() %>% 
  subset(abs(V1) >= .2) %>% 
  arrange(-V1) %>% 
  rownames() -> col_order_fc_slim_idh2

hm_fc_bulk_idh2[-1, col_order_fc_slim_idh2] %>% 
  Heatmap(cluster_rows = FALSE,
          cluster_columns = FALSE,
          row_labels = c("IVO", "ENA", "IVO+ENA"),
          show_column_names = TRUE,
          column_names_side = "top",
          column_title = "mIDH2 Pseudobulk",
          col=col_fun
  ) -> p_hm_idh2_bulk

## Clone5
d_mu_wide %>% 
  subset(Clone == "Clone5") %>% 
  .[, col_order_fc_slim_idh2] -> hm_fc_clone5

hm_fc_clone5[2, ] <- hm_fc_clone5[2, ] - hm_fc_clone5[1, ]
hm_fc_clone5[3, ] <- hm_fc_clone5[3, ] - hm_fc_clone5[1, ]
hm_fc_clone5[4, ] <- hm_fc_clone5[4, ] - hm_fc_clone5[1, ]

hm_fc_clone5[-1, col_order_fc_slim_idh2] %>% 
  Heatmap(cluster_rows = FALSE,
          cluster_columns = FALSE,
          row_labels = c("IVO", "ENA", "IVO+ENA"),
          show_column_names = FALSE,
          column_title = "Clone5",
          col=col_fun
  ) -> p_hm_idh2_clone5


## Clone6
d_mu_wide %>% 
  subset(Clone == "Clone6") %>% 
  .[, col_order_fc_slim_idh2] -> hm_fc_clone6

hm_fc_clone6[2, ] <- hm_fc_clone6[2, ] - hm_fc_clone6[1, ]
hm_fc_clone6[3, ] <- hm_fc_clone6[3, ] - hm_fc_clone6[1, ]
hm_fc_clone6[4, ] <- hm_fc_clone6[4, ] - hm_fc_clone6[1, ]

hm_fc_clone6[-1, col_order_fc_slim_idh2] %>% 
  Heatmap(cluster_rows = FALSE,
          cluster_columns = FALSE,
          row_labels = c("IVO", "ENA", "IVO+ENA"),
          show_column_names = FALSE,
          column_title = "Clone6",          
          col=col_fun
          ) -> p_hm_idh2_clone6

## Clone7
d_mu_wide %>% 
  subset(Clone == "Clone7") %>% 
  .[, col_order_fc_slim_idh2] -> hm_fc_clone7

hm_fc_clone7[2, ] <- hm_fc_clone7[2, ] - hm_fc_clone7[1, ]
hm_fc_clone7[3, ] <- hm_fc_clone7[3, ] - hm_fc_clone7[1, ]
hm_fc_clone7[4, ] <- hm_fc_clone7[4, ] - hm_fc_clone7[1, ]

hm_fc_clone7[-1, col_order_fc_slim_idh2] %>% 
  Heatmap(cluster_rows = FALSE,
          cluster_columns = FALSE,
          row_labels = c("IVO", "ENA", "IVO+ENA"),
          show_column_names = FALSE,
          column_title = "Clone7",
          col=col_fun
  ) -> p_hm_idh2_clone7
#
svg("p_hm_idh1_bulk.svg", width = 6.16, height = 1.89)
p_hm_idh1_bulk 
dev.off()
#
svg("p_hm_idh1_clone3.svg", width = 6.16, height = 1.16)
p_hm_idh1_clone3
dev.off()
#
svg("p_hm_idh1_clone4.svg", width = 6.16, height = 1.16)
p_hm_idh1_clone4
dev.off()

#
svg("p_hm_idh2_bulk.svg", width = 6.16, height = 1.89)
p_hm_idh2_bulk 
dev.off()
#
svg("p_hm_idh2_clone5.svg", width = 6.16, height = 1.16)
p_hm_idh2_clone5
dev.off()
#
svg("p_hm_idh2_clone6.svg", width = 6.16, height = 1.16)
p_hm_idh2_clone6
dev.off()
#
svg("p_hm_idh2_clone7.svg", width = 6.16, height = 1.16)
p_hm_idh2_clone7
dev.off()

# 
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
#   [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ComplexHeatmap_2.20.0       ggtext_0.1.2                ggsci_3.1.0                 patchwork_1.2.0            
# [5] ggdensity_1.0.0             ggrepel_0.9.5               rethinking_2.40             posterior_1.5.0            
# [9] cmdstanr_0.8.0              tapestri.tools_0.0.9000     lubridate_1.9.3             forcats_1.0.0              
# [13] stringr_1.5.1               purrr_1.0.2                 readr_2.1.5                 tidyr_1.3.1                
# [17] tibble_3.2.1                tidyverse_2.0.0             SingleCellExperiment_1.26.0 SummarizedExperiment_1.34.0
# [21] Biobase_2.64.0              MatrixGenerics_1.16.0       matrixStats_1.3.0           rtracklayer_1.64.0         
# [25] rhdf5_2.48.0                rentrez_1.2.3               ggridges_0.5.6              ggplot2_3.5.1              
# [29] GenomicRanges_1.56.0        GenomeInfoDb_1.40.1         IRanges_2.38.0              S4Vectors_0.42.0           
# [33] BiocGenerics_0.50.0         magrittr_2.0.3              maftools_2.20.0             dplyr_1.1.4                
# [37] data.table_1.15.4           compositions_2.0-8         
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3       tensorA_0.36.2.1         rstudioapi_0.16.0        jsonlite_1.8.8          
# [5] shape_1.4.6.1            GlobalOptions_0.1.2      BiocIO_1.14.0            zlibbioc_1.50.0         
# [9] vctrs_0.6.5              Rsamtools_2.20.0         RCurl_1.98-1.14          S4Arrays_1.4.1          
# [13] distributional_0.4.0     curl_5.2.1               broom_1.0.6              Rhdf5lib_1.26.0         
# [17] SparseArray_1.4.8        plyr_1.8.9               GenomicAlignments_1.40.0 lifecycle_1.0.4         
# [21] iterators_1.0.14         pkgconfig_2.0.3          Matrix_1.7-0             R6_2.5.1                
# [25] GenomeInfoDbData_1.2.12  clue_0.3-65              digest_0.6.35            colorspace_2.1-0        
# [29] ps_1.7.6                 fansi_1.0.6              timechange_0.3.0         httr_1.4.7              
# [33] abind_1.4-5              compiler_4.4.0           withr_3.0.0              doParallel_1.0.17       
# [37] backports_1.5.0          BiocParallel_1.38.0      MASS_7.3-60.2            bayesm_3.1-6            
# [41] DelayedArray_0.30.1      rjson_0.2.21             DNAcopy_1.78.0           loo_2.7.0               
# [45] tools_4.4.0              glue_1.7.0               restfulr_0.0.15          rhdf5filters_1.16.0     
# [49] gridtext_0.1.5           checkmate_2.3.1          cluster_2.1.6            generics_0.1.3          
# [53] gtable_0.3.5             tzdb_0.4.0               hms_1.1.3                xml2_1.3.6              
# [57] utf8_1.2.4               XVector_0.44.0           foreach_1.5.2            pillar_1.9.0            
# [61] robustbase_0.99-2        circlize_0.4.16          splines_4.4.0            lattice_0.22-6          
# [65] survival_3.5-8           tidyselect_1.2.1         Biostrings_2.72.0        knitr_1.46              
# [69] xfun_0.44                DEoptimR_1.1-3           stringi_1.8.4            UCSC.utils_1.0.0        
# [73] yaml_2.3.8               codetools_0.2-20         cli_3.6.2                munsell_0.5.1           
# [77] processx_3.8.4           Rcpp_1.0.12              coda_0.19-4.1            png_0.1-8               
# [81] XML_3.99-0.16.1          bitops_1.0-7             mvtnorm_1.2-5            scales_1.3.0            
# [85] crayon_1.5.2             GetoptLong_1.0.5         rlang_1.1.3             
