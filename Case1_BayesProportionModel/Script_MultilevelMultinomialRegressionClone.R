wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

library(tidyverse)
library(data.table)
library(tapestri.tools)
library(patchwork)
library(ggsci)
library(ggridges)
library(rethinking)
library(cmdstanr)
set_cmdstan_path()

# -------------------------------------------------------------------------
## What this script does:
## Fits the stan multilevel multinomial regression model to determine clonal 
## composition across treatments. Then it uses the posterior samples of cell numbers
## in animals to calculate the number of cells per clone across treatments.
## It also performs useful visualization.


## What is needed to run the script: 
## (1): 'ls_sce_filtered.rds' 
##       List of filtered tapestri data in single cell experiment format.
##       It is generated from 'Script1_QC.R'
fp_ls_sce_filtered <- "../Fig1E/ls_sce_filtered.rds"
##
## (2): 'StanModel_MultilevelMultiNomialClone.stan' 
##       Stan model. Needs to be in same directory as this script. 
##       Specifies the multilevel multinomial regression of clonal composition.
##
## (3): 'CellNumber_PosteriorSamples.rds'
##       Posterior samples of number of harvested human cells across timepoint.
##       See 'CellCountModel'.
fp_CellNumber_PosteriorSamples <- "../Case1_BayesCellCountModel/CellNumber_PosteriorSamples.rds"


## What the script outputs:
## (1): 'stanfit_multilevel_multinomial.rds'
##      Fitted stan object saved as an R object.
##
## (2): 'Matrix_PosteriorSamples_MultiNomialClone_Multilevel.rds'
##      Posterior samples stored as an R object.
##
## (3): 'Posterior_FrequencyPerClone.rds'
##       Posterior of frequency per clone across treatments stored as an R object.
##
## (4): 'Posterior_CellNumber_Clone_Treatment.rds'
##       Posterior of absolute cell number per clone across treatments stored as an R object.
##
## (5): 'p_posterior_clone_L2FC.png'
##       Density plot showing changes in cell number of each clone against vehicle.
##       This is Figure 1G
##
## (6): 'p_PosteriorShannonDiversityIndex.png'
##       Mean + errorbar showing posterior distribution of shannon diversity index.
##       This is Figure 1F
##
## Brief description of the stan model:
## Multilevel Multivariate model - Each of clone 0(root)~9 is modeled using a poisson likelihood
## with a multilevel varying intercept for each treatment. 

# -------------------------------------------------------------------------

## Colors 
color <- RColorBrewer::brewer.pal(9, "YlOrRd")
color <- colorRampPalette(color)(10)
color

# -------------------------------------------------------------------------

c(s1 = "Primary sample", 
  s12 = "VEN", 
  s14 = "IVO+VEN",
  s4 = "VEH", 
  s6 = "IVO",
  s8 = "AZA", 
  s10 = "IVO+AZA", 
  s2 = "Baseline") -> map_treatment

ls_sce_filtered <- read_rds(fp_ls_sce_filtered)

# d_Clone <- read_rds("d_Clone.rds")
lapply(ls_sce_filtered, function(x){
  colData(x)["Clone_annotated"][, 1] %>% 
    table() %>% 
    as.data.frame() -> d
  
  colnames(d) <- c("Clone", "Count")
  d$Percent <- with(d, Count*100 / sum(Count))
  d
}) %>% 
  bind_rows(.id = "Sample") -> d_Clone

d_Clone$Clone <- factor(d_Clone$Clone, levels = paste0("Clone", 0:9))
d_Clone$Treatment <- plyr::mapvalues(d_Clone$Sample, names(map_treatment), map_treatment)

d_Clone$Treatment <- factor(d_Clone$Treatment,
                            levels = c("Primary sample", "Baseline", "VEH", "IVO", "VEN", "IVO+VEN", "AZA", "IVO+AZA"))

## Turn the data into wide format
d_Clone %>% 
  pivot_wider(id_cols=c("Sample", "Treatment"), 
              names_from = "Clone", 
              values_from = "Count") -> d_CloneWide

d_CloneWide[is.na(d_CloneWide)] <- 0 #set NA to zero. 
d_CloneWide$Clone1 <- 0 #No clone1 detected in any sample

# -------------------------------------------------------------------------

## Prepare data to feed into stan model
## There are 8 treatment levels: 1-8
## (1) Primary sample
## (2) Baseline
## (3) VEH
## (4) IVO
## (5) VEN
## (6) IVO+VEN
## (7) AZA
## (8) IVO+AZA
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
  Clone8 = d_CloneWide$Clone8,
  Clone9 = d_CloneWide$Clone9,
  
  Treatment = plyr::mapvalues(d_CloneWide$Treatment,
                              c("Primary sample", "Baseline", "VEH", "IVO",
                                "VEN", "IVO+VEN", "AZA", "IVO+AZA"),
                              1:8),
  N_Treatment = length(unique(d_CloneWide$Treatment))
)

## Compile stan model
model <- cmdstanr::cmdstan_model("StanModel_MultilevelMultiNomialClone.stan")

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
write_rds(posterior, "stanfit_multilevel_multinomial.rds")

# Draw samples ------------------------------------------------------------

## Posterior samples
## n_Clone + (0~9 clones) + 1~8 treatments
post_samples <- posterior$draws(variables = paste0("n_Clone", 0:9), format = "draws_matrix")

## Save posterior samples of number of cells across clones and treatment into .rds
write_rds(post_samples, "Matrix_PosteriorSamples_MultiNomialClone_Multilevel.rds")

# -------------------------------------------------------------------------

## Calculate frequency of each clone from posterior
freq_post <- matrix(NA,
                    nrow=nrow(post_samples), ncol=ncol(post_samples),
                    dimnames=list(NULL, colnames(post_samples))
)

for(i in 0:9){ #clone
  for(j in 1:8){ #treatment
    col_which <- paste0("n_Clone", i, "[", j, "]")
    n <- post_samples[, col_which] 
    
    #select all clones of that treatment
    col_all <- paste0("n_Clone", 0:9, "[", j, "]")
    sum <- rowSums(post_samples[, col_all]) 
    
    #Calculate frequency 
    percent <- n / sum
    
    #store in "freq_prior" matrix
    freq_post[, col_which] <- percent
  }
}

write_rds(freq_post, "Posterior_FrequencyPerClone.rds")


## Calculate absolute number of cells
post <- read_rds(fp_CellNumber_PosteriorSamples) #Posterior samples
var <- attributes(post)$dimnames$variable

index_n <- grep("n_", var, value = TRUE)
index_n <- index_n[index_n != "n_inj"]

ls_count <- list()
for(i in index_n){
  ls_count[[i]] <- as.numeric(post[,,i]) 
  ls_count[[i]] <- sample(ls_count[[i]], 1e4)
}

ls_count %>% 
  bind_cols() %>% 
  pivot_longer(1:ncol(.), names_to = "Group", values_to = "Count") -> ls_count

ls_count$Group %>% 
  str_split("_") %>% 
  lapply(function(x) x[2]) %>% 
  unlist() %>% 
  str_sub(1, nchar(.) -3) -> ls_count$Timepoint
ls_count$Timepoint <- factor(ls_count$Timepoint, levels = c("bl", "w4", "end"))

ls_count$Group %>% 
  str_sub(nchar(.) - 1, nchar(.) - 1) -> ls_count$Treatment

ls_count$Treatment <- plyr::mapvalues(ls_count$Treatment,
                                      1:6, 
                                      c("IVO+VEN", "VEN", "IVO+AZA", "AZA", "IVO", "Vehicle")
                                      )

ls_count$Treatment <- factor(ls_count$Treatment, 
                             levels = c("Vehicle", "IVO", "AZA", "IVO+AZA", "VEN", "IVO+VEN")
                             )

ls_count %>% 
  subset(Timepoint == "end") %>% 
  mutate(iter = rep(1:1e4, each=6)) %>% 
  pivot_wider("iter", names_from="Treatment", values_from="Count") -> post_bm

#Calculate number of cells of each clone at end point
CellNumber_post <- matrix(NA, 
                          nrow=nrow(post_bm), 
                          ncol=6*10) #6 arms * 10 clones 

index_map <- c("Vehicle" = 3, #Read as: column vehicle in post_bm is 3 in stan multivariate
               "IVO" = 4,
               "VEN" = 5,
               "IVO+VEN" = 6,
               "AZA" = 7,
               "IVO+AZA" = 8
               ) 
#Set column names
c(
  paste(names(index_map)[1], paste0("Clone", 0:9), sep = "_"),
  paste(names(index_map)[2], paste0("Clone", 0:9), sep = "_"),
  paste(names(index_map)[3], paste0("Clone", 0:9), sep = "_"),
  paste(names(index_map)[4], paste0("Clone", 0:9), sep = "_"),
  paste(names(index_map)[5], paste0("Clone", 0:9), sep = "_"),
  paste(names(index_map)[6], paste0("Clone", 0:9), sep = "_")
) -> colnames(CellNumber_post)

freq_post_sub <- freq_post[sample(1:nrow(freq_post), 1e4), ] #Keep just 1e4 samples

for(i in names(index_map)){ #Using post_bm as reference 
  for(j in 0:9){
    col_which <- paste0(i, "_Clone", j) #column in CellNumber_post
    col_freq <- paste0("n_Clone", j, "[", index_map[i], "]") #column in freq_post_sub
    
    n <- unlist(freq_post_sub[, col_freq]) * unlist(post_bm[, i]) #frequency * number of cells
    CellNumber_post[, col_which] <- n
  }
}

write_rds(CellNumber_post, "Posterior_CellNumber_Clone_Treatment.rds")

# -------------------------------------------------------------------------

CellNumber_post %>% 
  as.data.table() %>% 
  pivot_longer(1:ncol(.), names_to="Group", values_to="CellNumber") %>% 
  mutate(Treatment = Group %>% str_split("_") %>% lapply(function(x) x[1]) %>% unlist(),
         Clone = Group %>% str_split("_") %>% lapply(function(x) x[2]) %>% unlist()
         ) -> d_CellNumber_post

#Calculate contrast versus vehicle (L2FC of cell number per clone
CellNumber_L2FC <- matrix(NA,
                          nrow=nrow(CellNumber_post), 
                          ncol=5*10) #5 arms * 10 clones 

colnames_CellNumber_L2FC <- colnames(CellNumber_post)
colnames_CellNumber_L2FC <- colnames_CellNumber_L2FC[!grepl("Vehicle", colnames_CellNumber_L2FC)]
colnames(CellNumber_L2FC) <- colnames_CellNumber_L2FC #L2FC v. Vehicle

for(i in names(index_map)[-1]){
  for(j in 0:9){
    col_which <- paste0(i, "_Clone", j) #Column of the final l2fc matrix. As well as
                                        #Coulmn name of treatment i cell number of clone j
    col_veh <- paste0("Vehicle_Clone", j) #Vehicle cell number of clone j
    
    n_which <- unlist(CellNumber_post[, col_which]) + 1 #add pseudocount
    n_veh <- unlist(CellNumber_post[, col_veh]) + 1 #add pseudocount
    
    l2fc <- log2(n_which / n_veh) #Calculate l2fc
    CellNumber_L2FC[, col_which] <- l2fc #Store in matrix
  }
}

CellNumber_L2FC %>% 
  as.data.table() %>% 
  pivot_longer(1:ncol(.), names_to="Group", values_to="L2FC_v_VEH") %>% 
  mutate(Treatment = Group %>% str_split("_") %>% lapply(function(x) x[1]) %>% unlist(),
         Clone = Group %>% str_split("_") %>% lapply(function(x) x[2]) %>% unlist()
  ) %>% 
  mutate(Treatment = factor(Treatment, 
                            levels = c("IVO", "VEN", "IVO+VEN", "AZA", "IVO+AZA"))
  ) %>% 
  subset(Clone != "Clone0") %>%
  subset(Clone != "Clone1") %>%
  subset(Clone != "Clone2") -> d_new

d_new$Clone <- factor(d_new$Clone, levels=paste0("Clone", 9:3))

ggplot(d_new, aes(L2FC_v_VEH, Clone, fill=Clone)) +
  geom_density_ridges(lwd=0.6) +
  facet_wrap(~Treatment, nrow=1) +
  geom_vline(xintercept=0, color="#5c1010") +
  labs(x="Log2 fold change v. Vehicle", y=NULL,
       title="Posterior distribution of change in cell number") +
  scale_fill_manual(values=setNames(color, paste0("Clone", 0:9)),
                    guide="none") +
  theme_minimal(14) +
  theme(plot.title.position = "plot") -> p_posterior_clone_L2FC

p_posterior_clone_L2FC
ggsave("p_posterior_clone_L2FC.png", p_posterior_clone_L2FC,
       width=14, height=4.5, units="in", dpi="retina")

# -------------------------------------------------------------------------

#Calculate shannon diversity index

ls_sdi <- list()
for(i in 1:8){ #8 treatments
  cols <- paste0("n_Clone", 0:9, "[", i, "]")
  m <- freq_post[, cols]
  
  sdi <- apply(m, 1, function(x){
    x <- x[x != 0]
    -sum(x*log(x))
  }) 
  sdi <- sdi[!is.na(sdi)]
  ls_sdi[[i]] <- sdi
}

d_sdi <- data.frame(H = unlist(ls_sdi), 
                    Treatment = factor(rep(1:8, times=unlist(lapply(ls_sdi, length))))
                    )

d_sdi %>% 
  group_by(Treatment) %>% 
  dplyr::summarise(mu_H = mean(H),
                   CI_low = quantile(H, c(.025, .975))[1],
                   CI_high = quantile(H, c(.025, .975))[2]
                   ) -> d_sdi_plot

d_sdi_plot$Treatment <- factor(d_sdi_plot$Treatment)

ls_sce_filtered %>% 
  lapply(function(x){
    clone <- unlist(colData(x)["Clone_annotated"])
    clone <- table(clone) / ncol(x)
    
    sdi <- sum(clone*log(clone)) * -1
    sdi
  }) %>% 
  bind_rows(.id = Treatment) -> sdi_observed
  
sdi_observed %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(Sample = rownames(.)) %>% 
  mutate(Treatment = factor(c(1, 3, 4, 5, 6, 7, 8, 2))) %>% 
  subset(Sample != "s15" & Sample != "s16") %>%  
  dplyr::rename("H"="V1") -> sdi_observed


ggplot(d_sdi_plot, aes(Treatment, mu_H)) +
  geom_point() +
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=.5) +
  scale_x_discrete(labels = c("Primary sample", "Baseline", "VEH", "IVO", 
                              "VEN", "IVO+VEN", "AZA", "IVO+AZA")) +
  labs(y=NULL, title="Shannon diversity index", 
       subtitle="Posterior mean + 95% credible interval")  +
  theme_classic(11) +
  theme(plot.title.position="plot", axis.text.x=element_text(angle=90)) -> p_PosteriorShannonDiversityIndex
p_PosteriorShannonDiversityIndex
ggsave("p_PosteriorShannonDiversityIndex.png",
       p_PosteriorShannonDiversityIndex,
       width = 2.83, height = 4.43, units = "in", dpi = "retina")


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
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] rethinking_2.21             cmdstanr_0.5.0              rstan_2.26.22               StanHeaders_2.26.26        
# [5] ggsci_2.9                   patchwork_1.1.1             tapestri.tools_0.0.9000     SingleCellExperiment_1.14.1
# [9] SummarizedExperiment_1.24.0 Biobase_2.54.0              MatrixGenerics_1.6.0        matrixStats_0.63.0         
# [13] rtracklayer_1.52.1          rhdf5_2.36.0                rentrez_1.2.3               ggridges_0.5.3             
# [17] GenomicRanges_1.46.1        GenomeInfoDb_1.30.1         IRanges_2.28.0              S4Vectors_0.32.4           
# [21] BiocGenerics_0.40.0         magrittr_2.0.3              maftools_2.8.05             compositions_2.0-4         
# [25] data.table_1.14.6           forcats_0.5.1               stringr_1.5.0               dplyr_1.0.10               
# [29] purrr_0.3.4                 readr_2.1.1                 tidyr_1.2.1                 tibble_3.1.8               
# [33] ggplot2_3.4.0               tidyverse_1.3.1            
# 
# loaded via a namespace (and not attached):
#   [1] colorspace_2.0-3         rjson_0.2.20             ellipsis_0.3.2           XVector_0.34.0          
# [5] fs_1.5.2                 rstudioapi_0.13          farver_2.1.1             mvtnorm_1.1-3           
# [9] fansi_1.0.3              lubridate_1.8.0          xml2_1.3.3               codetools_0.2-18        
# [13] splines_4.1.3            robustbase_0.95-0        knitr_1.39               jsonlite_1.8.4          
# [17] Rsamtools_2.8.0          broom_1.0.3              dbplyr_2.2.1             compiler_4.1.3          
# [21] httr_1.4.4               backports_1.4.1          assertthat_0.2.1         Matrix_1.4-0            
# [25] cli_3.4.1                prettyunits_1.1.1        tools_4.1.3              coda_0.19-4             
# [29] gtable_0.3.1             glue_1.6.2               GenomeInfoDbData_1.2.7   posterior_1.4.1.9000    
# [33] V8_4.2.0                 Rcpp_1.0.9               cellranger_1.1.0         vctrs_0.5.1             
# [37] Biostrings_2.62.0        rhdf5filters_1.4.0       tensorA_0.36.2           xfun_0.31               
# [41] ps_1.7.1                 rvest_1.0.2              lifecycle_1.0.3          restfulr_0.0.15         
# [45] XML_3.99-0.13            DEoptimR_1.0-11          zlibbioc_1.40.0          MASS_7.3-55             
# [49] scales_1.2.1             ragg_1.2.2               hms_1.1.2                inline_0.3.19           
# [53] RColorBrewer_1.1-3       curl_4.3.2               yaml_2.3.6               gridExtra_2.3           
# [57] loo_2.5.1                stringi_1.7.8            BiocIO_1.2.0             checkmate_2.1.0         
# [61] pkgbuild_1.4.0           BiocParallel_1.26.2      shape_1.4.6              systemfonts_1.0.4       
# [65] rlang_1.0.6              pkgconfig_2.0.3          bitops_1.0-7             distributional_0.3.0    
# [69] lattice_0.20-45          Rhdf5lib_1.14.2          labeling_0.4.2           GenomicAlignments_1.28.0
# [73] processx_3.6.1           tidyselect_1.2.0         plyr_1.8.8               R6_2.5.1                
# [77] generics_0.1.3           DelayedArray_0.20.0      DBI_1.1.3                pillar_1.8.1            
# [81] haven_2.4.3              withr_2.5.0              survival_3.2-13          abind_1.4-5             
# [85] RCurl_1.98-1.9           bayesm_3.1-4             modelr_0.1.8             crayon_1.5.2            
# [89] utf8_1.2.2               tzdb_0.2.0               grid_4.1.3               readxl_1.3.1            
# [93] callr_3.7.3              reprex_2.0.1             textshaping_0.3.6        RcppParallel_5.1.5      
# [97] munsell_0.5.0           


