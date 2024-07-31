wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

library(tidyverse)
library(ggsci)
library(patchwork)
library(ggpubr)

# -------------------------------------------------------------------------
## What is needed to run the script: 
## (1): 'CaseStudy1_FlowData.xlsx' 
##       This excel file contains the number of human cells we harvested in the experiment.


## What the script outputs:
## (1): 'p_DifferentiationMarkers.png' 
##       This is Figure 2D

# Setup data --------------------------------------------------------------
data_dir <- "../../DataFiles" #Change file path here
fp <- file.path(data_dir, "CaseStudy1", "CaseStudy1_FlowData.xlsx")

d <- readxl::read_xlsx(fp) 


d %>% 
  pivot_longer(starts_with("mfi_Endpoint"), 
               values_to = "MFI", 
               names_to = "Marker") -> d.long.mfi

d.long.mfi$Marker <-  plyr::mapvalues(d.long.mfi$Marker,
                                      c("mfi_Endpoint.BM_SSC", "mfi_Endpoint.BM_CD11b", 
                                        "mfi_Endpoint.BM_CD14", "mfi_Endpoint.BM_CD15"),
                                      c("Side scatter", "CD11b-APC", "CD14-PE", "CD15-PC7")
                                      )

d.long.mfi$Treatment <- factor(d.long.mfi$Treatment, 
                               levels = c("VEH", "IVO", "AZA", 
                                          "IVO+AZA", "VEN", "IVO+VEN")
                               )

## plot
ggplot(d.long.mfi, aes(Treatment, MFI, color = Treatment)) +
  geom_point(position = position_jitter(.05, seed = 123), size=2) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2) +
  facet_wrap(~Marker, scales = "free_y", nrow = 2) +
  labs(x = NULL, y = "Median fluorescent intensity") +
  scale_color_npg(guide = "none") +
  theme_classic(13) +
  theme(plot.title.position = "plot",
        axis.text.x = element_text(angle = 30, vjust = .8)
        ) -> p_DifferentiationMarkers
p_DifferentiationMarkers

ggsave("p_DifferentiationMarkers.png", 
       p_DifferentiationMarkers,
       height = 5.4, width = 6.2, units = "in", dpi = "retina")


