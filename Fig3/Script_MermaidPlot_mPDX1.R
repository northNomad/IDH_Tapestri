wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

library(DiagrammeR)
library(scales)
library(RColorBrewer)
library(DiagrammeRsvg)
library(rsvg)
library(ggsci)

library(scales)
pal_simpsons()(6)
show_col(pal_simpsons()(6))

#F2CB9E, #FEE7D2, #F1F5F8, #D5E4F4, #9DCBE4, #A8C8EB
grViz("
digraph a_nice_graph {


# node [fontname=Helvetica, shape=ellipse, style=filled]
## Define roots
node [shape = plaintext, fontname=Helvetica]

PM305256_root [label='root-PM305256', color='#FFF7F3']
PM246514_root [label='root-PM246514', color='#FFF7F3']
PM160950_root [label='root-PM160950', color='#FFF7F3']
PM160345_root [label='root-PM160345', color='#FFF7F3']
PM165009_root [label='root-PM165009', color='#FFF7F3']
PM325267_root [label='root-PM325267', color='#FFF7F3']

## Define nodes with no/close to no cells (plain text)
node [shape = plaintext, fontname=Helvetica]
PM160345_DNMT3A [label='<I>DNMT3A</I>@^{R882C}<br/>N=0', color='#FFF7F3']
PM246514_DNMT3A [label='<I>DNMT3A</I>@^{R882H}<br/>N=0', color='#FFF7F3']
PM165009_JAK2 [label='<I>JAK2</I>@^{V617F}<br/>N=0', color='#FFF7F3']

node [shape=circle, stroke=black, fontname=Helvetica, style=filled, fixedsize=TRUE]

### PM305256
PM305256_IDH1 [label='<I>IDH1</I>@^{R132H}<br/>N=6039', color='#F2CB9E', width=2]
PM305256_root -> PM305256_IDH1

### PM246514

PM246514_WT1 [label='<I>WT1</I>@^{R462W}<br/>N=4', color='#FEE7D2', width=0.006762468]
PM246514_IDH1 [label='<I>IDH1</I>@^{R132H}<br/>N=550', color='#FEE7D2', width=0.929839391]
PM246514_NPM1 [label='<I>NPM1</I>c@^{+}<br/>N=270', color='#FEE7D2', width=0.456466610]
PM246514_FLT3 [label='<I>FLT3</I>-ITD@^{ }<br/>N=359', color='#FEE7D2', width=0.606931530]

PM246514_root -> PM246514_DNMT3A
PM246514_DNMT3A -> PM246514_WT1
PM246514_WT1 -> PM246514_IDH1
PM246514_IDH1 -> PM246514_NPM1
PM246514_NPM1 -> PM246514_FLT3

### PM160950

PM160950_DNMT3A [label='<I>DNMT3A</I>@^{R882S}<br/>N=7', color='#F1F5F8', width=0.1458333]
PM160950_IDH1 [label='<I>IDH1</I>@^{R132S}<br/>N=61', color='#F1F5F8', width=1.2708333]
PM160950_FLT3 [label='<I>FLT3</I>-ITD@^{ }<br/>N=16', color='#F1F5F8', width=0.3333333]
PM160950_RUNX1 [label='<I>RUNX1</I>@^{R320*}<br/>N=12', color='#F1F5F8', width=0.2500000]

PM160950_root -> PM160950_DNMT3A 
PM160950_DNMT3A -> PM160950_IDH1
PM160950_IDH1 -> PM160950_FLT3
PM160950_FLT3 -> PM160950_RUNX1

### PM160345

PM160345_IDH1 [label='<I>IDH1</I>@^{R132H}<br/>N=12', color='#D5E4F4', width=0.4137931]
PM160345_NPM1 [label='<I>NPM1</I>c@^{+}<br/>N=34', color='#D5E4F4', width=1.1724138]
PM160345_FLT3 [label='<I>FLT3</I>-ITD@^{ }<br/>N=12', color='#D5E4F4', width=0.4137931]

PM160345_root -> PM160345_DNMT3A
PM160345_DNMT3A -> PM160345_IDH1
PM160345_IDH1 -> PM160345_NPM1
PM160345_NPM1 -> PM160345_FLT3

### PM165009
PM165009_EZH2_1 [label='<I>EZH2</I>@^{C552R}<br/>N=1', color='#9DCBE4', width=0.1333333]
PM165009_EZH2_2 [label='<I>EZH2</I>@^{P527H}<br/>N=2', color='#9DCBE4', width=0.2666667]
PM165009_RUNX1 [label='<I>RUNX1</I>@^{R320*}<br/>N=11', color='#9DCBE4', width=1.4666667]
PM165009_IDH1 [label='<I>IDH1</I>@^{R132C}<br/>N=1', color='#9DCBE4', width=0.1333333]

PM165009_root -> PM165009_JAK2
PM165009_JAK2 -> PM165009_EZH2_1
PM165009_EZH2_1 -> PM165009_EZH2_2
PM165009_EZH2_2 -> PM165009_RUNX1
PM165009_EZH2_2 -> PM165009_IDH1

### PM325267

PM325267_IDH1 [label='<I>IDH1</I>@^{R132C}<br/>N=1', color='#A8C8EB', width=0.1818182]
PM325267_NPM1 [label='<I>NPM1</I>c@^{+}<br/>N=7', color='#A8C8EB', width=1.2727273]
PM325267_NRAS [label='<I>NRAS</I>@^{G13D}<br/>N=3', color='#A8C8EB', width=0.5454545]

PM325267_root -> PM325267_IDH1
PM325267_IDH1 -> PM325267_NPM1
PM325267_NPM1 -> PM325267_NRAS

rankdir=LR
}

[1]: 'left'
[2]: 10:20
") -> p_scale_by_sample
p_scale_by_sample

p_scale_by_sample %>%
  DiagrammeRsvg::export_svg() %>%
  charToRaw() %>%
  rsvg_svg("PDXA_p_scale_by_sample.svg")


#### No scaling
grViz("
digraph a_nice_graph {


# node [fontname=Helvetica, shape=ellipse, style=filled]
## Define roots
node [shape = plaintext, fontname=Helvetica]

PM305256_root [label='root-PM305256', color='#FFF7F3']
PM246514_root [label='root-PM246514', color='#FFF7F3']
PM160950_root [label='root-PM160950', color='#FFF7F3']
PM160345_root [label='root-PM160345', color='#FFF7F3']
PM165009_root [label='root-PM165009', color='#FFF7F3']
PM325267_root [label='root-PM325267', color='#FFF7F3']

## Define nodes with no/close to no cells (plain text)


node [shape=ellipse, stroke=black, fontname=Helvetica, style=filled]
PM160345_DNMT3A [label='<I>DNMT3A</I>@^{R882C}<br/>N=0', color='#FFF7F3']
PM246514_DNMT3A [label='<I>DNMT3A</I>@^{R882H}<br/>N=0', color='#FFF7F3']
PM165009_JAK2 [label='<I>JAK2</I>@^{V617F}<br/>N=0', color='#FFF7F3']
### PM305256
PM305256_IDH1 [label='<I>IDH1</I>@^{R132H}<br/>N=6039', color='#F2CB9E']
PM305256_root -> PM305256_IDH1

### PM246514

PM246514_WT1 [label='<I>WT1</I>@^{R462W}<br/>N=4', color='#FEE7D2']
PM246514_IDH1 [label='<I>IDH1</I>@^{R132H}<br/>N=550', color='#FEE7D2']
PM246514_NPM1 [label='<I>NPM1</I>c@^{+}<br/>N=270', color='#FEE7D2']
PM246514_FLT3 [label='<I>FLT3</I>-ITD@^{ }<br/>N=359', color='#FEE7D2']

PM246514_root -> PM246514_DNMT3A
PM246514_DNMT3A -> PM246514_WT1
PM246514_WT1 -> PM246514_IDH1
PM246514_IDH1 -> PM246514_NPM1
PM246514_NPM1 -> PM246514_FLT3

### PM160950

PM160950_DNMT3A [label='<I>DNMT3A</I>@^{R882S}<br/>N=7', color='#F1F5F8']
PM160950_IDH1 [label='<I>IDH1</I>@^{R132S}<br/>N=61', color='#F1F5F8']
PM160950_FLT3 [label='<I>FLT3</I>-ITD@^{ }<br/>N=16', color='#F1F5F8']
PM160950_RUNX1 [label='<I>RUNX1</I>@^{R320*}<br/>N=12', color='#F1F5F8']

PM160950_root -> PM160950_DNMT3A 
PM160950_DNMT3A -> PM160950_IDH1
PM160950_IDH1 -> PM160950_FLT3
PM160950_FLT3 -> PM160950_RUNX1

### PM160345

PM160345_IDH1 [label='<I>IDH1</I>@^{R132H}<br/>N=12', color='#D5E4F4']
PM160345_NPM1 [label='<I>NPM1</I>c@^{+}<br/>N=34', color='#D5E4F4']
PM160345_FLT3 [label='<I>FLT3</I>-ITD@^{ }<br/>N=12', color='#D5E4F4']

PM160345_root -> PM160345_DNMT3A
PM160345_DNMT3A -> PM160345_IDH1
PM160345_IDH1 -> PM160345_NPM1
PM160345_NPM1 -> PM160345_FLT3

### PM165009
PM165009_EZH2_1 [label='<I>EZH2</I>@^{C552R}<br/>N=1', color='#9DCBE4']
PM165009_EZH2_2 [label='<I>EZH2</I>@^{P527H}<br/>N=2', color='#9DCBE4']
PM165009_RUNX1 [label='<I>RUNX1</I>@^{R320*}<br/>N=11', color='#9DCBE4']
PM165009_IDH1 [label='<I>IDH1</I>@^{R132C}<br/>N=1', color='#9DCBE4']

PM165009_root -> PM165009_JAK2
PM165009_JAK2 -> PM165009_EZH2_1
PM165009_EZH2_1 -> PM165009_EZH2_2
PM165009_EZH2_2 -> PM165009_RUNX1
PM165009_EZH2_2 -> PM165009_IDH1

### PM325267

PM325267_IDH1 [label='<I>IDH1</I>@^{R132C}<br/>N=1', color='#A8C8EB']
PM325267_NPM1 [label='<I>NPM1</I>c@^{+}<br/>N=7', color='#A8C8EB']
PM325267_NRAS [label='<I>NRAS</I>@^{G13D}<br/>N=3', color='#A8C8EB']

PM325267_root -> PM325267_IDH1
PM325267_IDH1 -> PM325267_NPM1
PM325267_NPM1 -> PM325267_NRAS

rankdir=LR
}

[1]: 'left'
[2]: 10:20
") -> p_no_scale
p_no_scale

p_no_scale %>%
  DiagrammeRsvg::export_svg() %>%
  charToRaw() %>%
  rsvg_svg("PDXA_p_no_scale.svg")


# -------------------------------------------------------------------------
# > sessionInfo()
# R version 4.3.1 (2023-06-16 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 11 x64 (build 22621)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=English_Canada.utf8  LC_CTYPE=English_Canada.utf8    LC_MONETARY=English_Canada.utf8
# [4] LC_NUMERIC=C                    LC_TIME=English_Canada.utf8    
# 
# time zone: America/Toronto
# tzcode source: internal
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] rsvg_2.6.0                  DiagrammeRsvg_0.1           scales_1.2.1               
# [4] DiagrammeR_1.0.10           ggupset_0.3.0               ggpmisc_0.5.5              
# [7] ggpp_0.5.5                  RColorBrewer_1.1-3          tapestri.tools_0.0.9000    
# [10] lubridate_1.9.2             forcats_1.0.0               stringr_1.5.0              
# [13] purrr_1.0.2                 readr_2.1.4                 tidyr_1.3.0                
# [16] tibble_3.2.1                tidyverse_2.0.0             SingleCellExperiment_1.22.0
# [19] SummarizedExperiment_1.30.2 Biobase_2.60.0              MatrixGenerics_1.12.3      
# [22] matrixStats_1.0.0           rtracklayer_1.60.1          rhdf5_2.44.0               
# [25] rentrez_1.2.3               ggridges_0.5.4              ggplot2_3.4.3              
# [28] GenomicRanges_1.52.0        GenomeInfoDb_1.36.1         IRanges_2.34.1             
# [31] S4Vectors_0.38.1            BiocGenerics_0.46.0         magrittr_2.0.3             
# [34] maftools_2.16.0             dplyr_1.1.2                 data.table_1.14.8          
# [37] compositions_2.0-6          vcfR_1.14.0                
# 
# loaded via a namespace (and not attached):
#   [1] bitops_1.0-7             polynom_1.4-1            permute_0.9-7            rlang_1.1.1             
# [5] compiler_4.3.1           mgcv_1.8-42              vctrs_0.6.3              quantreg_5.97           
# [9] fastmap_1.1.1            memuse_4.2-3             pkgconfig_2.0.3          crayon_1.5.2            
# [13] ellipsis_0.3.2           XVector_0.40.0           labeling_0.4.2           utf8_1.2.3              
# [17] Rsamtools_2.16.0         tzdb_0.4.0               bit_4.0.5                MatrixModels_0.5-2      
# [21] zlibbioc_1.46.0          jsonlite_1.8.7           rhdf5filters_1.12.1      DelayedArray_0.26.7     
# [25] Rhdf5lib_1.22.0          BiocParallel_1.34.2      parallel_4.3.1           cluster_2.1.4           
# [29] R6_2.5.1                 stringi_1.7.12           DNAcopy_1.74.1           Rcpp_1.0.11             
# [33] Matrix_1.6-1             splines_4.3.1            timechange_0.2.0         tidyselect_1.2.0        
# [37] rstudioapi_0.15.0        abind_1.4-5              yaml_2.3.7               vegan_2.6-4             
# [41] codetools_0.2-19         curl_5.0.2               lattice_0.21-8           withr_2.5.0             
# [45] pinfsc50_1.2.0           survival_3.5-5           bayesm_3.1-5             Biostrings_2.68.1       
# [49] pillar_1.9.0             tensorA_0.36.2           generics_0.1.3           vroom_1.6.3             
# [53] RCurl_1.98-1.12          hms_1.1.3                munsell_0.5.0            glue_1.6.2              
# [57] tools_4.3.1              BiocIO_1.10.0            robustbase_0.99-0        SparseM_1.81            
# [61] GenomicAlignments_1.36.0 ggvenn_0.1.10            visNetwork_2.1.2         XML_3.99-0.14           
# [65] grid_4.3.1               ape_5.7-1                colorspace_2.1-0         nlme_3.1-162            
# [69] GenomeInfoDbData_1.2.10  restfulr_0.0.15          cli_3.6.1                fansi_1.0.4             
# [73] S4Arrays_1.0.5           viridisLite_0.4.2        V8_4.3.3                 gtable_0.3.3            
# [77] DEoptimR_1.1-1           digest_0.6.33            htmlwidgets_1.6.2        farver_2.1.1            
# [81] rjson_0.2.21             htmltools_0.5.6          lifecycle_1.0.3          bit64_4.0.5             
# [85] MASS_7.3-60   