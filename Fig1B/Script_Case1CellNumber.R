wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

library(tidyverse)
library(patchwork)

# -------------------------------------------------------------------------
## What is needed to run the script: 
## (1): 'CaseStudy1_FlowData.xlsx' 
##       This excel file contains the number of human cells we harvested in the experiment.


## What the script outputs:
## (1): 'p_case1_rawcount.png' 
##       This is Figure 1B

# Setup data --------------------------------------------------------------
data_dir <- "../../DataFiles" #Change file path here
fp <- file.path(data_dir, "CaseStudy1", "CaseStudy1_FlowData.xlsx")

d <- readxl::read_xlsx(fp) 
d <- d[c(1:10, 26:30, 21:25, 16:20, 11, 12, 15, 13, 14), ] #Change ordering of treatment groups

d$Treatment <- as.numeric(plyr::mapvalues(d$Treatment, 
                                          c("VEH", "IVO", "AZA", "IVO+AZA", "VEN", "IVO+VEN"), 
                                          6:1))

##Calculate number of frozen BM cells. The dilution factor is 75. We took 1 out of 200uL for flow.
d$n_hCD45hCD33_Endpoint.BM <- d$LeukemicBurdenEvents_per_uL_Endpoint.BM*75*199

#Round number of cells to integer
d$n_hCD45hCD33_FrozenBL <- round(d$n_hCD45hCD33_FrozenBL)
d$n_hCD45hCD33_FrozenW4 <- round(d$n_hCD45hCD33_FrozenW4)
d$n_hCD45hCD33_Endpoint.BM <- round(d$n_hCD45hCD33_Endpoint.BM)



d %>% 
  pivot_longer(starts_with("n_hCD45"), names_to="timepoint", values_to="n") %>% 
  mutate(timepoint_lab = plyr::mapvalues(timepoint,
                                         c("n_hCD45hCD33_FrozenBL",
                                           "n_hCD45hCD33_FrozenW4",
                                           "n_hCD45hCD33_Endpoint.BM"
                                         ),
                                         c("Baseline (RFM)",
                                           "Week4 (LFM)",
                                           "Endpoint (Whole BM)")
  ),
  timepoint_lab = factor(timepoint_lab, c("Baseline (RFM)",
                                          "Week4 (LFM)",
                                          "Endpoint (Whole BM)")
  ),
  treatment_lab = plyr::mapvalues(Treatment, 
                                  6:1, 
                                  c("VEH", "IVO", "AZA", "IVO+AZA", "VEN", "IVO+VEN")
  ),
  treatment_lab = factor(treatment_lab, 
                         c("VEH", "IVO", "VEN", "IVO+VEN", "AZA", "IVO+AZA")
  )
  ) -> d_raw


## bl
d_raw %>% 
  subset(timepoint %in% c("n_hCD45hCD33_FrozenBL")) %>% 
  ggplot(aes(treatment_lab, n)) +
  geom_jitter(width=.1, shape=21, size=2) +
  stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. - 0.25, yend=..y..))+
  stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. + 0.25, yend=..y..)) +
  labs(x=NULL, y=NULL, title="Baseline (RFM)") +
  scale_y_log10(breaks=c(10, 100, 1000, 5000, 10000, 20000),
                labels=c(10, 100, "1K", "5K", "10K", "20K"),
                limits=c(10, 2e4)) +
  theme_classic(13) +
  theme(axis.text.x = element_text(angle=90),
        plot.title.position = "plot") -> p_num_bl


## w4
d_raw %>% 
  subset(timepoint %in% c("n_hCD45hCD33_FrozenW4")) %>% 
  ggplot(aes(treatment_lab, n)) +
  geom_jitter(width=.1, shape=21, size=2) +
  stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. - 0.25, yend=..y..))+
  stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. + 0.25, yend=..y..)) +
  labs(x=NULL, y=NULL, title="Week4 (LFM)") +
  scale_y_log10(breaks=c(10, 100, 1000, 5000, 10000, 20000),
                labels=c(10, 100, "1K", "5K", "10K", "20K"),
                limits=c(10, 2e4)) +
  theme_classic(13) +
  theme(axis.text.x = element_text(angle=90),
        plot.title.position = "plot") -> p_num_w4

## endpoint
d_raw %>% 
  subset(timepoint %in% c("n_hCD45hCD33_Endpoint.BM")) %>% 
  ggplot(aes(treatment_lab, n)) +
  geom_jitter(width=.1, shape=21, size=2) +
  stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. - 0.25, yend=..y..))+
  stat_summary(fun="mean", geom="segment", mapping=aes(xend=..x.. + 0.25, yend=..y..)) +
  labs(x=NULL, y=NULL, title="Endpoint (Whole BM)") +
  scale_y_log10(breaks=c(10, 100, 1000, 1e4, 1e5, 1e6, 10e6),
                labels=c(10, 100, "1K", "10K", "100K", "1M", "10M"),
  ) +
  theme_classic(13) +
  theme(axis.text.x = element_text(angle=90),
        plot.title.position = "plot") -> p_num_end

p_case1_rawcount <- p_num_bl + p_num_w4 + p_num_end
p_case1_rawcount
ggsave("p_case1_rawcount.png",
       p_case1_rawcount,
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
# [6] scales_1.3.0      readxl_1.4.3      R6_2.5.1          plyr_1.8.9        generics_0.1.3   
# [11] munsell_0.5.1     pillar_1.9.0      tzdb_0.4.0        rlang_1.1.3       utf8_1.2.4       
# [16] stringi_1.8.4     timechange_0.3.0  cli_3.6.2         withr_3.0.0       magrittr_2.0.3   
# [21] grid_4.4.0        rstudioapi_0.16.0 hms_1.1.3         lifecycle_1.0.4   vctrs_0.6.5      
# [26] glue_1.7.0        farver_2.1.2      cellranger_1.1.0  fansi_1.0.6       colorspace_2.1-0 
# [31] tools_4.4.0       pkgconfig_2.0.3  
