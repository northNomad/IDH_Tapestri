wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

library(tidyverse)

# -------------------------------------------------------------------------
## What is needed to run the script: 
## (1): 'CaseStudy1_FlowData.xlsx' 
##       This excel file contains the number of human cells we harvested in the experiment.


## What the script outputs:
## (1): 'p_case1_LeukemicBurden_Kinetics.png'
## (2): 'p_case1_LeukemicBurden_Kinetics.normalized.png'
##       These are Figure 1C


data_dir <- "../../DataFiles"
fp <- file.path(data_dir, "CaseStudy1", "CaseStudy1_FlowData.xlsx")
d <- readxl::read_xlsx(fp)


#convert to percentage scale
d$perc_LeukemicBurden_Endpoint.BM <- d$perc_LeukemicBurden_Endpoint.BM*100

#Calculate number of frozen BM cells
d$n_hCD45hCD33_Endpoint.BM <- d$LeukemicBurdenEvents_per_uL_Endpoint.BM*75*199


#select columns for plotting
plot.cols <- c("n_hCD45hCD33_FrozenBL", "n_hCD45hCD33_FrozenW4", "n_hCD45hCD33_Endpoint.BM",#frozen cells
               "Cage", "EarMark","Treatment", "MiceNumber", #metadata 
               "perc_LeukemicBurden_Endpoint.BM", "perc_LeukemicBurden_BL.RFM", "perc_LeukemicBurden_W4.LFM"
               )

d <- d[, plot.cols]

#refactor
d$Treatment <- factor(d$Treatment, levels = c("VEH", "IVO", "AZA", "IVO+AZA", "VEN", "IVO+VEN"))

# plot - leukemic burden (BM) ---------------------------------------------

d.long <- d %>%  
  pivot_longer(c(perc_LeukemicBurden_BL.RFM, perc_LeukemicBurden_W4.LFM, perc_LeukemicBurden_Endpoint.BM), 
               names_to = "Timepoint", values_to = "LeukemicBurden.BM") 

d.long$Timepoint <-  plyr::mapvalues(d.long$Timepoint,
                                        c("perc_LeukemicBurden_BL.RFM", "perc_LeukemicBurden_W4.LFM", "perc_LeukemicBurden_Endpoint.BM"),
                                        c("BL", "W4", "Endpoint")
                                        )

d.long$Timepoint <- factor(d.long$Timepoint, levels = c("BL", "W4", "Endpoint"))
d.long$Treatment <- factor(d.long$Treatment, levels = c("VEH", "IVO", "VEN", "IVO+VEN", "AZA", "IVO+AZA"))

ggplot(d.long, aes(Timepoint, LeukemicBurden.BM)) +
  geom_point() +
  geom_line(aes(group = EarMark)) +
  facet_grid(~Treatment) +
  labs(x = "Timepoint", y = "hCD45+ hCD33+ (%)", 
       title = "Leukemic burden kinetics (raw)") +
  scale_x_discrete(label = c("Baseline\nRFM", "Week4\nLFM", "Week8\nWhole BM")) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), 
                     breaks = c(0, 1, 10, 50, 100),
                     labels = c(0, 1, 10, 50, 100)
  ) +
  theme_classic(13) +
  theme(strip.background = element_blank(), 
        plot.title.position = "plot") -> p_case1_LeukemicBurden_Kinetics
p_case1_LeukemicBurden_Kinetics
ggsave("p_case1_LeukemicBurden_Kinetics.png", plot = p_case1_LeukemicBurden_Kinetics,
       height = 2.7, width = 11.5, units = "in", dpi = "retina")

# plot - leukemic burden (BM) normalized to baseline ----------------------

d %>% 
  group_by(Treatment, EarMark) %>% 
  dplyr::summarise(W4.norm = perc_LeukemicBurden_W4.LFM / perc_LeukemicBurden_BL.RFM,
                   Endpoint.norm = perc_LeukemicBurden_Endpoint.BM / perc_LeukemicBurden_BL.RFM) %>% 
  mutate(BL = 1) %>% 
  pivot_longer(c(W4.norm, Endpoint.norm, BL), names_to = "Timepoint", values_to = "NormalizedLeukemicBurden") %>% 
  mutate(Treatment = factor(Treatment, levels = c("VEH", "IVO", "VEN", "IVO+VEN", "AZA", "IVO+AZA"))) -> dt.burden.norm

dt.burden.norm$Timepoint <- factor(dt.burden.norm$Timepoint, levels = c("BL", "W4.norm", "Endpoint.norm"))

ggplot(dt.burden.norm, aes(Timepoint, NormalizedLeukemicBurden)) +
  geom_point() +
  geom_line(aes(group = EarMark)) +
  facet_grid(~Treatment) +
  labs(x = "Timepoint", y = "Leukemic burden\n(fold change vs baseline)", 
       title = "Leukemic burden kinetics (normalized to baseline)") +
  scale_x_discrete(label = c("Baseline\nRFM", "Week4\nLFM", "Week8\nWhole BM")) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 50, 100, 300, 600),
                labels = c(0.01, 0.1, 1, 10, 50, 100, 300, 600)) +
  theme_classic(13) +
  theme(strip.background = element_blank(), 
        plot.title.position = "plot") -> p_case1_LeukemicBurden_Kinetics.normalized
p_case1_LeukemicBurden_Kinetics.normalized
ggsave("p_case1_LeukemicBurden_Kinetics.normalized.png", plot = p_case1_LeukemicBurden_Kinetics.normalized,
       height = 2.7, width = 11.5, units = "in", dpi = "retina")

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
#   [1] lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2    
# [6] readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0
# 
# loaded via a namespace (and not attached):
#   [1] vctrs_0.6.5       cli_3.6.2         rlang_1.1.3       stringi_1.8.4    
# [5] generics_0.1.3    glue_1.7.0        colorspace_2.1-0  plyr_1.8.9       
# [9] hms_1.1.3         readxl_1.4.3      scales_1.3.0      fansi_1.0.6      
# [13] grid_4.4.0        cellranger_1.1.0  munsell_0.5.1     tzdb_0.4.0       
# [17] lifecycle_1.0.4   compiler_4.4.0    timechange_0.3.0  Rcpp_1.0.12      
# [21] pkgconfig_2.0.3   rstudioapi_0.16.0 farver_2.1.2      R6_2.5.1         
# [25] tidyselect_1.2.1  utf8_1.2.4        pillar_1.9.0      magrittr_2.0.3   
# [29] tools_4.4.0       withr_3.0.0       gtable_0.3.5 