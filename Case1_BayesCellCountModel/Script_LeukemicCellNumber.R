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
## (1): 'CaseStudy1_FlowData.xlsx' 
##       This excel file contains the number of human cells we harvested in the experiment.
## (2): 'StanModel_CellNumber.stan' 
##       A stan model file that specifies the hierarchical model of leukemic cell
##       number in an animal using a gamma poisson likelihood.    

## What the script outputs:
## (1): 'CellNumber_PriorPredictiveCheck.svg' 
##       This is the prior predictive distribution for the stan model
## (2): 'CellNumber_PosteriorSamples.rds'
##       This is the posterior samples of the number of human cells in animals.

# Setup data --------------------------------------------------------------
dir_name <- "../../DataFiles"
fp <- file.path(dir_name, "CaseStudy1", "CaseStudy1_FlowData.xlsx")

d <- readxl::read_xlsx(fp) #Change file path here
d <- d[c(1:10, 26:30, 21:25, 16:20, 11, 12, 15, 13, 14), ] #Change ordering of treatment groups

d$Treatment <- as.numeric(plyr::mapvalues(d$Treatment, 
                                          c("VEH", "IVO", "AZA", "IVO+AZA", "VEN", "IVO+VEN"), 
                                          6:1))
d$Mice <- paste0("Mouse", 1:30)
d$MiceNumber <- 1:30

map_Treatment <- setNames(d$Treatment, d$Mice)

##Calculate number of frozen BM cells. The dilution factor is 75. We took 1 out of 200uL for flow.
d$n_hCD45hCD33_Endpoint.BM <- d$LeukemicBurdenEvents_per_uL_Endpoint.BM*75*199

#Round number of cells to integer
d$n_hCD45hCD33_FrozenBL <- round(d$n_hCD45hCD33_FrozenBL)
d$n_hCD45hCD33_FrozenW4 <- round(d$n_hCD45hCD33_FrozenW4)
d$n_hCD45hCD33_Endpoint.BM <- round(d$n_hCD45hCD33_Endpoint.BM)

# Organize data -----------------------------------------------------------

data <- list(
  #Cell numbers
  C_bl = d$n_hCD45hCD33_FrozenBL,
  C_w4 = d$n_hCD45hCD33_FrozenW4,
  C_end = d$n_hCD45hCD33_Endpoint.BM[1:28],
  Treatment = d$Treatment #Vehicle is 6, coef set to zero (reference)
)

# Prior predictive check --------------------------------------------------

set.seed(123)
N <- 20
#Grand mean, should be between log(0~5e6)
mu0 <- runif(N, 0, 17.7) 

#Drug coefficients
mu_T <- rnorm(N, 0, 7)
sigma_T <- abs(rnorm(N, 0, 2))
bT <- rnorm(N, mu_T, sigma_T) 

mu <- exp(mu0 + bT)
thi <- abs(rnorm(N, 0, 1)) #Overdispersion

svg("CellNumber_PriorPredictiveCheck.svg", width = 5, height = 3.5)
plot(NULL, xlim=c(0, 50e6), ylim=c(0, 1e-6),
     xlab="Number of cells", ylab="Density", main="Prior predictive check")
for(i in 1:N){
  nCell <- rnbinom(1e5, mu=mu[i], size=thi[i])
  lines(density(nCell),
        lwd=3, col=col.alpha(2, .5))
}
dev.off()

# Run stan model ----------------------------------------------------------

model <- cmdstanr::cmdstan_model("StanModel_CellNumber.stan") #compile
#Fit the model. 8 Chains, 1000 warmup, 1e4 iterations per chain.
post <- model$sample(data=data,
                     iter_warmup=1e3, 
                     chains=8,
                     iter_sampling=1e4,
                     parallel_chains=8, 
                     max_treedepth=20,
                     seed=123)

post$print(max_rows = 46) #Model fits. All Rhat satisfactory.

#Draw samples from posterior
samples <- post$draws()

#Save R object.
write_rds(samples, "CellNumber_PosteriorSamples.rds")


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
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] rethinking_2.40 posterior_1.5.0 lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2    
# [8] readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0 readxl_1.4.3    cmdstanr_0.8.0 
# 
# loaded via a namespace (and not attached):
#   [1] tensorA_0.36.2.1     utf8_1.2.4           generics_0.1.3       shape_1.4.6.1        lattice_0.22-6      
# [6] stringi_1.8.4        hms_1.1.3            magrittr_2.0.3       timechange_0.3.0     grid_4.4.0          
# [11] mvtnorm_1.2-5        cellranger_1.1.0     plyr_1.8.9           jsonlite_1.8.8       processx_3.8.4      
# [16] backports_1.5.0      ps_1.7.6             fansi_1.0.6          scales_1.3.0         abind_1.4-5         
# [21] cli_3.6.2            rlang_1.1.3          munsell_0.5.1        withr_3.0.0          tools_4.4.0         
# [26] tzdb_0.4.0           coda_0.19-4.1        checkmate_2.3.1      colorspace_2.1-0     vctrs_0.6.5         
# [31] R6_2.5.1             matrixStats_1.3.0    lifecycle_1.0.4      MASS_7.3-60.2        pkgconfig_2.0.3     
# [36] pillar_1.9.0         gtable_0.3.5         loo_2.7.0            glue_1.7.0           Rcpp_1.0.12         
# [41] xfun_0.44            tidyselect_1.2.1     rstudioapi_0.16.0    knitr_1.46           compiler_4.4.0      
# [46] distributional_0.4.0