wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

library(readxl)
library(tidyverse)
library(ggsci)

data_dir <- "../../DataFiles/"

## What is needed to run the script: 
## (1): 'EngraftmentScreen_TidyData.xlsx' 
##       Excel file containing engraftment efficiency of each sample as single PDX.
d <- read_xlsx(file.path(data_dir, "Screening", "EngraftmentScreen_TidyData.xlsx"))

## What the script outputs:
## (1): 'p_engraftment_screen.png'
##      This is Figure S2

d %>% 
  group_by(Sample) %>% 
  dplyr::summarise(mu=mean(LeukemicBurden)) %>% 
  arrange(mu) %>% 
  .$Sample -> sample_order

d$Sample <- factor(d$Sample, sample_order)

ggplot(d, aes(Sample, LeukemicBurden)) +
  geom_jitter(aes(shape=Injection, group=Injection), width=.1, size=2) +
  stat_summary(fun.data="mean_se", geom="pointrange", color="#A50F15") +
  labs(x="Sample", y="% Human CD45<sup>+</sup> Cells in BM", title="Leukemic burden in BM") +
  coord_flip() +
  scale_color_npg(guide="none") +
  scale_shape_manual(values=c(IF=21, IV=22)) +
  scale_y_log10() +
  theme_bw(13) +
  theme(plot.title.position="plot", 
        axis.title.x=element_markdown()) -> p_engraftment_screen
p_engraftment_screen

ggsave("p_engraftment_screen.png", 
       p_engraftment_screen,
       width = 8, height = 4.6, units = "in", dpi="retina")

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
#   [1] ggtext_0.1.2       RColorBrewer_1.1-3 ggsci_3.1.0        lubridate_1.9.3   
# [5] forcats_1.0.0      stringr_1.5.1      dplyr_1.1.4        purrr_1.0.2       
# [9] readr_2.1.5        tidyr_1.3.1        tibble_3.2.1       ggplot2_3.5.1     
# [13] tidyverse_2.0.0    readxl_1.4.3      
# 
# loaded via a namespace (and not attached):
#   [1] tensorA_0.36.2.1     utf8_1.2.4           generics_0.1.3       xml2_1.3.6          
# [5] stringi_1.8.4        hms_1.1.3            magrittr_2.0.3       timechange_0.3.0    
# [9] grid_4.4.0           cellranger_1.1.0     plyr_1.8.9           processx_3.8.4      
# [13] backports_1.5.0      ps_1.7.6             fansi_1.0.6          scales_1.3.0        
# [17] textshaping_0.4.0    abind_1.4-5          cli_3.6.2            rlang_1.1.3         
# [21] cmdstanr_0.8.0       commonmark_1.9.1     munsell_0.5.1        withr_3.0.0         
# [25] tools_4.4.0          tzdb_0.4.0           checkmate_2.3.1      colorspace_2.1-0    
# [29] vctrs_0.6.5          posterior_1.5.0      R6_2.5.1             lifecycle_1.0.4     
# [33] ragg_1.3.2           pkgconfig_2.0.3      pillar_1.9.0         gtable_0.3.5        
# [37] glue_1.7.0           Rcpp_1.0.12          systemfonts_1.1.0    xfun_0.44           
# [41] tidyselect_1.2.1     rstudioapi_0.16.0    knitr_1.46           farver_2.1.2        
# [45] labeling_0.4.3       compiler_4.4.0       markdown_1.13        distributional_0.4.0
# [49] gridtext_0.1.5