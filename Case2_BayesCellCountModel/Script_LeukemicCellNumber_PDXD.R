wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

library(cmdstanr)
library(readxl)
library(parallel)
library(tidyverse)
library(rethinking)
set_cmdstan_path()

set.seed(123)

# -------------------------------------------------------------------------
## What is needed to run the script: 
## (1): 'CaseStudy2_CellCount.xlsx' 
##       This excel file contains the number of human cells we harvested in the experiment (in the sheet 'CellCount').
## (2): 'StanModel_CellNumber2.stan' 
##       A stan model file that specifies the hierarchical model of leukemic cell
##       number in an animal using a gamma poisson likelihood.    

## What the script outputs:
## (1): 'CellNumber_PosteriorSamples_Case2.rds'
##       This is the posterior samples of the number of human cells in animals.
## (2): 'df_posterior_nCells_8Chains_5e4Iter.rds'
##       Also the posterior samples, but a data.frame for easy wrangling.

# Setup data --------------------------------------------------------------
data_dir <- "../../DataFiles"
fp <- file.path(data_dir, "CaseStudy2", "CaseStudy2_CellCount.xlsx")

d <- readxl::read_xlsx(fp, sheet = "CellCount")

## Treatment:
# 1 = Ivosidenib
# 2 = Enasidenib
# 3 = Ivo + Ena
# 4 = Vehicle

d$Treatment <- plyr::mapvalues(d$Treatment,
                               c("IVO", "ENA", "IVO+ENA", "VEH"), 
                               1:4)

d$Treatment
# Organize data -----------------------------------------------------------

data <- list(
  #Cell numbers
  C_bl = d$n_hCD45hCD33_FrozenBL,
  C_end = d$n_hCD45hCD33_Endpoint.BM,
  Treatment = as.numeric(d$Treatment) #Vehicle is 4, coef set to zero (reference)
)


# Run stan model ----------------------------------------------------------

model <- cmdstanr::cmdstan_model("StanModel_CellNumber2.stan") #compile
#Fit the model. 8 Chains, 5e4 warmup, 5e4 iterations per chain.
post <- model$sample(data=data,
                     iter_warmup=5e4,
                     chains=8,
                     iter_sampling=5e4,
                     parallel_chains=8,
                     max_treedepth=500,
                     seed=123)


post$print(max_rows = 34) #Model fits. All Rhat satisfactory.

#Draw samples from posterior
samples <- post$draws()

#Save R object.
write_rds(samples, "CellNumber_PosteriorSamples_Case2.rds")
samples <- read_rds("CellNumber_PosteriorSamples_Case2.rds")
post <- samples

## Treatment:
# 1 = Ivosidenib
# 2 = Enasidenib
# 3 = Ivo + Ena
# 4 = Vehicle

var <- attributes(post)$dimnames$variable

# Calculate cell number ---------------------------------------------------

## create matrix
matrix(nrow = 5e4 * 8, ncol = 8, 
       dimnames = list(NULL, c("bl_ivo", "bl_ena", "bl_ivo.ena", "bl_veh", 
                               "end_ivo", "end_ena", "end_ivo.ena", "end_veh")
       )
) -> d_count

## calculate posterior cell number for baseline
for(i in 1:4){
  mu0 <- as.numeric(post[, , "mu0_bl"])
  bT <- as.numeric(post[, , paste0("bT_adj_bl[", i, "]")])
  thi <- as.numeric(post[, , "thi_bl"])
  
  n <- rnbinom(nrow(d_count), mu = exp(mu0 + bT), size = thi)
  d_count[, i] <- n
}

## calculate posterior cell number for endpoint
for(i in 1:4){
  mu0 <- as.numeric(post[, , "mu0_end"])
  bT <- as.numeric(post[, , paste0("bT_adj_end[", i, "]")])
  thi <- as.numeric(post[, , "thi_end"])
  
  n <- rnbinom(nrow(d_count), mu = exp(mu0 + bT), size = thi)
  d_count[, i + 4] <- n
}

## Add chain and iteration info
d_count <- as.data.table(d_count)
d_count$Chain <- rep(1:8, each=5e4)
d_count$Iteration <- rep(1:5e4, 8)

##Save
write_rds(d_count, "df_posterior_nCells_8Chains_5e4Iter.rds")

# session info ------------------------------------------------------------

# > sessionInfo()
# R version 4.1.3 (2022-03-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    
# 
# attached base packages:
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggpubr_0.4.0        patchwork_1.1.1     ggsci_2.9           CytoExploreR_1.1.0  openCyto_2.4.0      flowWorkspace_4.4.0 flowCore_2.4.0     
# [8] rethinking_2.21     rstan_2.26.22       StanHeaders_2.26.26 forcats_0.5.1       stringr_1.5.0       dplyr_1.0.10        purrr_0.3.4        
# [15] readr_2.1.1         tidyr_1.2.1         tibble_3.1.8        ggplot2_3.4.0       tidyverse_1.3.1     readxl_1.3.1        cmdstanr_0.5.0     
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.2.2           reticulate_1.25      R.utils_2.12.2       ks_1.13.5            tidyselect_1.2.0     htmlwidgets_1.5.4   
# [7] grid_4.1.3           Rtsne_0.16           rainbow_3.6          aws.signature_0.6.0  munsell_0.5.0        codetools_0.2-18    
# [13] umap_0.2.8.0         withr_2.5.0          colorspace_2.0-3     flowViz_1.56.0       Biobase_2.54.0       knitr_1.39          
# [19] rstudioapi_0.13      stats4_4.1.3         flowClust_3.30.0     robustbase_0.95-0    ggsignif_0.6.3       mnormt_2.1.0        
# [25] farver_2.1.1         changepoint_2.2.3    coda_0.19-4          vctrs_0.5.1          generics_0.1.3       xfun_0.31           
# [31] R6_2.5.1             clue_0.3-60          rsvd_1.0.5           bitops_1.0-7         cachem_1.0.6         assertthat_0.2.1    
# [37] promises_1.2.0.1     scales_1.2.1         gtable_0.3.1         processx_3.6.1       RProtoBufLib_2.4.0   rlang_1.0.6         
# [43] splines_4.1.3        rstatix_0.7.0        hexbin_1.28.2        broom_1.0.3          checkmate_2.1.0      inline_0.3.19       
# [49] reshape2_1.4.4       abind_1.4-5          modelr_0.1.8         backports_1.4.1      httpuv_1.6.6         IDPmisc_1.1.20      
# [55] RBGL_1.68.0          tensorA_0.36.2       tools_4.1.3          ellipsis_0.3.2       jquerylib_0.1.4      posterior_1.4.1.9000
# [61] RColorBrewer_1.1-3   BiocGenerics_0.40.0  Rcpp_1.0.9           plyr_1.8.8           visNetwork_2.1.0     base64enc_0.1-3     
# [67] zlibbioc_1.40.0      RCurl_1.98-1.9       ps_1.7.1             prettyunits_1.1.1    openssl_2.0.5        S4Vectors_0.32.4    
# [73] deSolve_1.32         zoo_1.8-10           haven_2.4.3          cluster_2.1.2        fs_1.5.2             fda_6.0.4           
# [79] magrittr_2.0.3       ncdfFlow_2.38.0      data.table_1.14.6    RSpectra_0.16-1      hdrcde_3.4           reprex_2.0.1        
# [85] mvtnorm_1.1-3        matrixStats_0.63.0   hms_1.1.2            mime_0.12            evaluate_0.15        xtable_1.8-4        
# [91] XML_3.99-0.13        jpeg_0.1-9           mclust_5.4.10        gridExtra_2.3        shape_1.4.6          compiler_4.1.3      
# [97] ellipse_0.4.3        flowStats_4.4.0      KernSmooth_2.23-20   V8_4.2.0             crayon_1.5.2         R.oo_1.25.0         
# [103] fds_1.8              htmltools_0.5.4      corpcor_1.6.10       pcaPP_2.0-1          later_1.3.0          tzdb_0.2.0          
# [109] rrcov_1.7-0          RcppParallel_5.1.5   lubridate_1.8.0      aws.s3_0.3.21        DBI_1.1.3            dbplyr_2.2.1        
# [115] MASS_7.3-55          car_3.1-0            Matrix_1.4-0         cli_3.4.1            flowAI_1.22.0        R.methodsS3_1.8.2   
# [121] pkgconfig_2.0.3      xml2_1.3.3           bslib_0.4.1          rvest_1.0.2          distributional_0.3.0 callr_3.7.3         
# [127] digest_0.6.31        pracma_2.3.8         graph_1.72.0         rmarkdown_2.14       cellranger_1.1.0     curl_4.3.2          
# [133] shiny_1.7.3          gtools_3.9.2.2       lifecycle_1.0.3      jsonlite_1.8.4       carData_3.0-5        askpass_1.1         
# [139] fansi_1.0.3          pillar_1.8.1         lattice_0.20-45      loo_2.5.1            fastmap_1.1.0        httr_1.4.4          
# [145] DEoptimR_1.0-11      pkgbuild_1.4.0       glue_1.6.2           png_0.1-8            rhandsontable_0.3.8  Rgraphviz_2.36.0    
# [151] stringi_1.7.8        sass_0.4.4           latticeExtra_0.6-29  cytolib_2.4.0        EmbedSOM_2.1.1      
