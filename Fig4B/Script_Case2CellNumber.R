wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

library(tidyverse)
library(patchwork)

# -------------------------------------------------------------------------
## What is needed to run the script: 
## (1): 'CaseStudy2_CellCount.xlsx' 
##       This excel file contains the number of human cells we harvested in the experiment.

## What the script outputs:
## (1): 'p_case2_rawcount.png' 
##       This is Figure 4B

# Setup data --------------------------------------------------------------
data_dir <- "../../DataFiles"
fp <- file.path(data_dir, "CaseStudy2", "CaseStudy2_CellCount.xlsx")
d_raw <- readxl::read_xlsx(fp, sheet = "CellCount")

d_raw$Treatment <- factor(d_raw$Treatment, levels=c("VEH", "IVO", "ENA", "IVO+ENA"))

## bl
d_raw %>% 
  ggplot(aes(Treatment, n_hCD45hCD33_FrozenBL)) +
  geom_jitter(width=.1, shape=21, size=2) +
  stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. - 0.25, yend=..y..))+
  stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. + 0.25, yend=..y..)) +
  labs(x=NULL, y=NULL, title="Baseline (RFM)") +
  scale_y_log10(breaks=c(10, 10000, 1e5, 1.5e6), 
                labels=c(10, 10000, "100K", "1.5M"), 
                limits=c(.1, 1.5e6)) +
  theme_classic(13) +
  theme(plot.title.position = "plot") -> p_bl


## w4
d_raw %>% 
  ggplot(aes(Treatment, n_hCD45hCD33_FrozenW4)) +
  geom_jitter(width=.1, shape=21, size=2) +
  stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. - 0.25, yend=..y..))+
  stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. + 0.25, yend=..y..)) +
  labs(x=NULL, y=NULL, title="Week4 (LFM)") +
  scale_y_log10(breaks=c(10, 10000, 1e5, 1.5e6), 
                labels=c(10, 10000, "100K", "1.5M"), 
                limits=c(.1, 1.5e6)) +
  theme_classic(13) +
  theme(plot.title.position = "plot") -> p_w4

## endpoint
d_raw %>% 
  ggplot(aes(Treatment, n_hCD45hCD33_Endpoint.BM)) +
  geom_jitter(width=.1, shape=21, size=2) +
  stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. - 0.25, yend=..y..))+
  stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. + 0.25, yend=..y..)) +
  labs(x=NULL, y=NULL, title="Endpoint (Whole BM)") +
  scale_y_log10(breaks=c(10, 1e5, 1e6, 2e6, 3e6, 4e6),
                labels=c(10, "100K", "1M", "2M", "3M", "4M"),
                limits=c(5000, 5e6)) +
  theme_classic(13) +
  theme(plot.title.position = "plot") -> p_end


p_case2_rawcount <- p_bl + p_w4 + p_end
p_case2_rawcount
ggsave("p_case2_rawcount.png",
       p_case2_rawcount,
       width = 9, height = 3.8, units = "in", dpi = "retina")

# > sessionInfo()
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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] patchwork_1.2.0 lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2    
# [7] readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.5      compiler_4.4.0    crayon_1.5.2      tidyselect_1.2.1  Rcpp_1.0.12      
# [6] textshaping_0.4.0 systemfonts_1.1.0 scales_1.3.0      readxl_1.4.3      R6_2.5.1         
# [11] plyr_1.8.9        generics_0.1.3    munsell_0.5.1     pillar_1.9.0      tzdb_0.4.0       
# [16] rlang_1.1.3       utf8_1.2.4        stringi_1.8.4     timechange_0.3.0  cli_3.6.2        
# [21] withr_3.0.0       magrittr_2.0.3    grid_4.4.0        rstudioapi_0.16.0 hms_1.1.3        
# [26] lifecycle_1.0.4   vctrs_0.6.5       glue_1.7.0        farver_2.1.2      cellranger_1.1.0 
# [31] ragg_1.3.2        fansi_1.0.6       colorspace_2.1-0  tools_4.4.0       pkgconfig_2.0.3  