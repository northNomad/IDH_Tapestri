wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

library(tidyverse)
library(data.table)
library(tapestri.tools)
library(patchwork)
library(ggsci)
library(ggridges)
library(circlize)
library(rethinking)
library(cmdstanr)
library(ComplexHeatmap)
set_cmdstan_path()

## What is needed to run the script: 
## (1): 'ls_sce_filtered_mixed.rds' 
##       List of filtered tapestri data in single cell experiment format.
##       It is generated from 'Script_Case2QC.R'
ls_sce_filtered2 <- read_rds("../Case2_QC/ls_sce_filtered_mixed.rds")
##
## (2): 'StanModel_MultilevelMultiNomialClone2.stan' 
##       Stan model. Needs to be in same directory as this script. 
##       Specifies the multilevel multinomial regression of clonal composition.
##
## (3): 'df_posterior_nCells_8Chains_5e4Iter.rds'
##       Posterior samples of number of harvested human cells across timepoint.
d_n <- read_rds("../Case2_BayesCellCountModel/df_posterior_nCells_8Chains_5e4Iter.rds")


## What the script outputs:
## (1): 'stanfit_multilevel_multinomial_mixed.rds'
##      Fitted stan object saved as an R object.
##
## (2): 'Matrix_PosteriorSamples_MultiNomialClone_Multilevel_mixed.rds'
##      Posterior samples stored as an R object.
##
## (3): 'Posterior_FrequencyPerClone_mixed.rds'
##       Posterior of frequency per clone across treatments stored as an R object.
##
## (4): 'p_hm_ClonesAcrossSamples_mixed.svg'
##       This is Figure 4C
## (5): 'p_mPDX_l2fc_ridge.png'
##       This is Figure 4D


#Count clones
lapply(ls_sce_filtered2, function(x){
  as.data.table(table(colData(x)$Clone))
}) %>% 
  bind_rows(.id="Treatment") %>% 
  mutate(V1=paste0("Clone", V1)) %>% 
  pivot_wider(id_cols="Treatment", names_from="V1", values_from="N") -> d_CloneWide

d_CloneWide[is.na(d_CloneWide)] <- 0


stan_data <- list(
  N = nrow(d_CloneWide),
  Clone0 = d_CloneWide$Clone0,
  Clone1 = d_CloneWide$Clone1,
  Clone2 = d_CloneWide$Clone2,
  Clone3 = d_CloneWide$Clone3,
  Clone4 = d_CloneWide$Clone4,
  Clone5 = d_CloneWide$Clone5,
  Clone6 = d_CloneWide$Clone6,
  Clone7 = d_CloneWide$Clone7,
  
  Treatment = plyr::mapvalues(d_CloneWide$Treatment,
                              c("veh", "ivo", "ena", "ivo.ena"),
                              1:4),
  N_Treatment = length(unique(d_CloneWide$Treatment))
)

## Compile stan model
model <- cmdstanr::cmdstan_model("StanModel_MultilevelMultiNomialClone2.stan")

## Sample from posterior
model$sample(data = stan_data,
             parallel_chains = 8,
             chains = 8,
             iter_warmup = 1000,
             iter_sampling = 20000, 
             seed = 123) -> posterior

## Take a look at coefficient table
posterior$print(max_rows = 1e3) #All Rhat = 1

## Save entire posterior
write_rds(posterior, "stanfit_multilevel_multinomial_mixed.rds")

# Draw samples ------------------------------------------------------------

## Posterior samples
post_samples <- posterior$draws(variables = paste0("n_Clone", 0:7), format = "draws_matrix")

## Save posterior samples of number of cells across clones and treatment into .rds
write_rds(post_samples, "Matrix_PosteriorSamples_MultiNomialClone_Multilevel_mixed.rds")

# -------------------------------------------------------------------------

## Calculate frequency of each clone from posterior
freq_post <- matrix(NA,
                    nrow=nrow(post_samples), ncol=ncol(post_samples),
                    dimnames=list(NULL, colnames(post_samples))
)

for(i in 0:7){ #clone
  for(j in 1:4){ #treatment
    col_which <- paste0("n_Clone", i, "[", j, "]")
    n <- post_samples[, col_which] 
    
    #select all clones of that treatment
    col_all <- paste0("n_Clone", 0:7, "[", j, "]")
    sum <- rowSums(post_samples[, col_all]) 
    
    #Calculate frequency 
    percent <- n / sum
    
    #store in "freq_prior" matrix
    freq_post[, col_which] <- percent
  }
}

write_rds(freq_post, "Posterior_FrequencyPerClone_mixed.rds")

# -------------------------------------------------------------------------
## Draw heatmap
d_hm <- as.data.frame(d_CloneWide)
rownames(d_hm) <- d_hm$Treatment
d_hm <- t(d_hm[, -1])
d_hm <- d_hm[paste0("Clone", 0:7), ]
d_hm <- d_hm[, c("veh", "ivo", "ena", "ivo.ena")]

d_percent <- as.data.frame(d_hm)
d_percent %>% 
  lapply(function(x) round(x*100/sum(x), 2)) %>% 
  do.call(cbind, .) -> d_percent

d_label <- matrix(NA, nrow=nrow(d_percent), ncol=ncol(d_percent))
for(i in 1:nrow(d_percent)){
  for(j in 1:ncol(d_percent)){
    label <- paste0(d_hm[i, j], "\n", "(", d_percent[i, j], ")")
    d_label[i, j] <- label
  }
}
col_fun <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
d_hm %>% 
  scale() %>% 
  Heatmap(cluster_rows=FALSE, 
          cluster_columns=FALSE,
          #
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(d_label[i, j], x, y, gp = gpar(fontsize = 10))
          },
          row_split = c(1, 2, 2, 2, 2, 3, 3, 3),
          name = "Z-scaled proportion",
          column_labels = c("Vehicle", "Ivosidenib", "Enasidenib", "Ivo. + Ena."),
          column_names_rot = -30, 
          row_names_side = "left", 
          column_names_side = "top",
          col = col_fun
          ) -> p_hm_ClonesAcrossSamples_mixed_v2
p_hm_ClonesAcrossSamples_mixed_v2

svg(filename = "p_hm_ClonesAcrossSamples_mixed.svg", width = 5.35, height = 4.8)
p_hm_ClonesAcrossSamples_mixed_v2
dev.off()


# Calculating absolute cell numbers ---------------------------------------

d_n <- d_n[, c("end_veh", "end_ivo", "end_ena", "end_ivo.ena")]
colnames(d_n) <- gsub("end_", "", colnames(d_n))
d_n <- as.matrix(d_n)

#downsample
N <- 5e4 #Number of samples
d_n <- d_n[sample(1:nrow(d_n), N), ]
freq_post <-  freq_post [sample(1:nrow(freq_post ), N), ]

#Initiate matrix to store posterior
m_abs <- matrix(NA, nrow=N, ncol=8*4) #8Clones * 4 treatments

c(paste0(paste0("Clone", 0:7), "_veh"),
  paste0(paste0("Clone", 0:7), "_ivo"),
  paste0(paste0("Clone", 0:7), "_ena"),
  paste0(paste0("Clone", 0:7), "_ivo.ena")) -> colnames(m_abs)

for(i in paste0("Clone", 0:7)){
  for(j in c("veh", "ivo", "ena", "ivo.ena")){
    col_which <- paste0(i, "_", j)
    
    k <- plyr::mapvalues(j, c("veh", "ivo", "ena", "ivo.ena"), 1:4)
    col_freq <- paste0("n_", i, "[", k, "]")
    
    n <- d_n[, j]
    p <- freq_post[, col_freq]
    
    m_abs[, col_which] <- n*p
    
  }
}

# Calculate L2FC ----------------------------------------------

#3 arms v. VEH
#combination v. ivo
#combination v. ena
#ivo. v. ena

m_l2fc <- matrix(NA, nrow=N, ncol=8*6) #8clones * 6comparisons

c(paste0(paste0("Clone", 0:7), "_", "ivo", "_", "veh"),
  paste0(paste0("Clone", 0:7), "_", "ena", "_", "veh"),
  paste0(paste0("Clone", 0:7), "_", "ivo.ena", "_", "veh"),
  paste0(paste0("Clone", 0:7), "_", "ivo.ena", "_", "ivo"),
  paste0(paste0("Clone", 0:7), "_", "ivo.ena", "_", "ena"),
  paste0(paste0("Clone", 0:7), "_", "ivo", "_", "ena")
  ) -> colnames(m_l2fc)

for(i in colnames(m_l2fc)){
  clone <- str_split(i, "_") %>% lapply(function(x) x[1]) %>% unlist()
  t1 <- str_split(i, "_") %>% lapply(function(x) x[2]) %>% unlist()
  t2 <- str_split(i, "_") %>% lapply(function(x) x[3]) %>% unlist()
  
  col_t1 <- paste0(clone, "_", t1)
  col_t2 <- paste0(clone, "_", t2)
  
  l2fc <- log2((m_abs[, col_t1] + 1) / (m_abs[, col_t2] + 1))
  m_l2fc[, i] <- l2fc
}


m_l2fc %>% 
  as.data.table() %>% 
  pivot_longer(1:ncol(.), names_to="id", values_to="l2fc") %>% 
  mutate(Clone=str_sub(id, 1, 6),
         Treatment=str_sub(id, 8, -1L)
         ) -> m_l2fc_long

###
#colors
root1 [label = 'root - PM160345', color = ]
a [label = '<I>DNMT3A</I>@^{R882C}', color = ]
b [label = '<I>IDH1</I>@^{R132H}', color = ]
c [label = '<I>NPM1</I>c@^{+}', color = ]
d [label = '<I>FLT3</I>-ITD@^{ }', color = ]

root2 [label = 'root - PM160053', color = ]
a2 [label = '<I>IDH2</I>@^{R140Q}', color = ]
b2 [label = '<I>NPM1</I>c@^{+}', color = ]
c2 [label = '<I>FLT3</I>-ITD@^{ }', color = ]

colors <- c('#FFF7F3', '#FDE0DD', '#FCC5C0', '#FA9FB5', '#F768A1',
            '#DEEBF7', '#C6DBEF', '#9ECAE1')

colors <- setNames(colors, paste0("Clone", 0:7))

m_l2fc_long_sub <- m_l2fc_long %>% 
  subset(Treatment %in% c("ivo_veh", "ena_veh", "ivo.ena_veh")) %>% 
  mutate(Treatment=plyr::mapvalues(Treatment,
                                   c("ena_veh", "ivo.ena_veh", "ivo_veh"),
                                   c("ENA", "IVO+ENA", "IVO")
                                   )
         ) %>% 
  mutate(Treatment=factor(Treatment, levels=c("IVO", "ENA", "IVO+ENA")),
         Clone=factor(Clone, levels=paste0("Clone", 7:0))) 


m_l2fc_long_sub %>% 
  subset(Clone != "Clone0") %>%
  subset(Clone != "Clone1") %>%
  subset(Clone != "Clone2") %>%
  ggplot(aes(l2fc, Clone, fill=Clone)) +
  geom_density_ridges2(lwd=.6, panel_scaling = FALSE) +
  geom_vline(xintercept=0, color="#5c1010") +
  labs(x="Log2 fold change v. Vehicle", y=NULL,
       title="Posterior distribution of change in cell number") +
  facet_grid(~Treatment) +
  scale_fill_manual(guide="none", values=colors) +
  theme_minimal(14) +
  theme(plot.title.position = "plot") -> p_mPDX_l2fc_ridge
p_mPDX_l2fc_ridge
ggsave("p_mPDX_l2fc_ridge.png", p_mPDX_l2fc_ridge, 
       width=6.2, height = 3.4, units = "in", dpi="retina"
       )


# -------------------------------------------------------------------------
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
#   [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ComplexHeatmap_2.20.0       rethinking_2.40             posterior_1.5.0            
# [4] cmdstanr_0.8.0              ggsci_3.1.0                 patchwork_1.2.0            
# [7] tapestri.tools_0.0.9000     SingleCellExperiment_1.26.0 SummarizedExperiment_1.34.0
# [10] Biobase_2.64.0              MatrixGenerics_1.16.0       matrixStats_1.3.0          
# [13] rtracklayer_1.64.0          rhdf5_2.48.0                rentrez_1.2.3              
# [16] ggridges_0.5.6              GenomicRanges_1.56.0        GenomeInfoDb_1.40.1        
# [19] IRanges_2.38.0              S4Vectors_0.42.0            BiocGenerics_0.50.0        
# [22] magrittr_2.0.3              maftools_2.20.0             compositions_2.0-8         
# [25] data.table_1.15.4           lubridate_1.9.3             forcats_1.0.0              
# [28] stringr_1.5.1               dplyr_1.1.4                 purrr_1.0.2                
# [31] readr_2.1.5                 tidyr_1.3.1                 tibble_3.2.1               
# [34] ggplot2_3.5.1               tidyverse_2.0.0             circlize_0.4.16            
# 
# loaded via a namespace (and not attached):
#   [1] bitops_1.0-7             rlang_1.1.3              clue_0.3-65              GetoptLong_1.0.5        
# [5] compiler_4.4.0           loo_2.7.0                png_0.1-8                vctrs_0.6.5             
# [9] pkgconfig_2.0.3          shape_1.4.6.1            crayon_1.5.2             backports_1.5.0         
# [13] XVector_0.44.0           utf8_1.2.4               Rsamtools_2.20.0         tzdb_0.4.0              
# [17] UCSC.utils_1.0.0         ps_1.7.6                 xfun_0.44                zlibbioc_1.50.0         
# [21] jsonlite_1.8.8           rhdf5filters_1.16.0      DelayedArray_0.30.1      Rhdf5lib_1.26.0         
# [25] BiocParallel_1.38.0      cluster_2.1.6            R6_2.5.1                 stringi_1.8.4           
# [29] RColorBrewer_1.1-3       DNAcopy_1.78.0           iterators_1.0.14         Rcpp_1.0.12             
# [33] knitr_1.46               Matrix_1.7-0             splines_4.4.0            timechange_0.3.0        
# [37] tidyselect_1.2.1         rstudioapi_0.16.0        abind_1.4-5              yaml_2.3.8              
# [41] doParallel_1.0.17        codetools_0.2-20         curl_5.2.1               processx_3.8.4          
# [45] lattice_0.22-6           plyr_1.8.9               withr_3.0.0              coda_0.19-4.1           
# [49] survival_3.5-8           bayesm_3.1-6             Biostrings_2.72.0        pillar_1.9.0            
# [53] tensorA_0.36.2.1         foreach_1.5.2            checkmate_2.3.1          distributional_0.4.0    
# [57] generics_0.1.3           RCurl_1.98-1.14          hms_1.1.3                munsell_0.5.1           
# [61] scales_1.3.0             glue_1.7.0               tools_4.4.0              BiocIO_1.14.0           
# [65] robustbase_0.99-2        GenomicAlignments_1.40.0 mvtnorm_1.2-5            XML_3.99-0.16.1         
# [69] colorspace_2.1-0         GenomeInfoDbData_1.2.12  restfulr_0.0.15          cli_3.6.2               
# [73] fansi_1.0.6              S4Arrays_1.4.1           gtable_0.3.5             DEoptimR_1.1-3          
# [77] digest_0.6.35            SparseArray_1.4.8        rjson_0.2.21             lifecycle_1.0.4         
# [81] httr_1.4.7               GlobalOptions_0.1.2      MASS_7.3-60.2     
