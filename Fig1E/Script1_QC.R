wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

library(vcfR)
library(tapestri.tools)
library(RColorBrewer)
library(tidyverse)
library(ggvenn)
library(ggpmisc)


# -------------------------------------------------------------------------
data_dir <- "../../DataFiles/CaseStudy1" #Data directory
## Files needed to run this script:
## **All these files must be in the same directory called 'data_dir'
## (1) 'CMIA0005_tumorOnly_human_2.vep.vcf' 
##
##      This .vcf file contains results from ENSEMBL VEP using data from targeted bulk sequencing  
##      on bone marrow cells collected from untreated animals
##
## (2) '1-AML_23_05_01.dna+protein.h5'
## (3) '4-VEH-END_23_05_24.dna+protein.h5'
## (4) '6-IVO-END_23_05_24.dna+protein.h5'
## (5) '12-VENE-END_23_05_24.dna.h5'
## (6) '14-IVOVEN-END_23_05_24.dna.h5'
## (7) '8-AZA-END_DNA_23_07_28.dna.h5'
## (8) '10-IVOAZA-END_DNA_23_07_28.dna.h5'
## (9) '2-BASELINE_DNA_23_07_28.dna.h5'
##
##  These .h5 files are scDNA results exported from MisisonBio's Tapestri Pipeline.
##
## (12) 'AML-v2RevA.bed
##      
##      This .bed file contains the regions covered by the 127 Tapestri AML amplicons in hg19 coordinate.
##

## What the script outputs:
## (1): 'p_venn_tapestri_bulk.svg' 
##       Venn diagram of variants called by scDNA and bulk seq.
##       This is Figure S1A
## (2): 'p_vaf_s4_v_bulk.png'
##       Concordance of scDNA seq and bulk seq (average scDNA vaf v. bulk vaf) of all 33 variants called by bulk.
##       This is Figure S1B
## (3): 'ls_sce_filtered.rds'
##       An R object. A list of single cell experiments containing filtered scDNA seq. results with clonal annotation.

##Useful colors
color <- RColorBrewer::brewer.pal(9, "YlOrRd")
color <- colorRampPalette(color)(10)
color

#vcf file from bulk sequencing
vcf <- vcfR::read.vcfR(file = file.path(data_dir, "CMIA0005_tumorOnly_human_2.vep.vcf")) 

#GRanges of Tapestri Amplicons 
amplicon <- tapestri.tools::read_amplicons_bed(
  bed = file.path(data_dir, "AML-v2RevA.bed"),
  format = "GRanges"
  ) 

#Turn vcf into GRanges
GRanges(seqnames=vcf@fix[, "CHROM"], 
        ranges=IRanges(as.numeric(vcf@fix[, "POS"]),
                       width=nchar(vcf@fix[, "ALT"])
                       ),
        mcols=data.frame(
          REF=vcf@fix[, "REF"],
          ALT=vcf@fix[, "ALT"],
          INFO=vcf@fix[, "INFO"],
          FORMAT=vcf@gt[, "CMIA0005"]
          )
        ) -> gr_vcf

# Add DP and AF column
mcols(gr_vcf)[["mcols.FORMAT"]] %>% 
  str_split(":") %>% 
  lapply(function(x) x[3]) %>% 
  unlist() %>% 
  as.numeric() -> gr_vcf$mcols.AF

mcols(gr_vcf)[["mcols.FORMAT"]] %>% 
  str_split(":") %>% 
  lapply(function(x) x[4]) %>% 
  unlist() %>% 
  as.numeric() -> gr_vcf$mcols.DP


gr_vcf <- subsetByOverlaps(gr_vcf, amplicon) #Keep variants covered by amplicons
sum(grepl(",", gr_vcf$mcols.ALT)) #Three variants are ambiguous in intronic regions. Remove those.
gr_vcf <- gr_vcf[!grepl(",", gr_vcf$mcols.ALT)]

gr_vcf #33 variants in .vcf (out of original 100455), 14 are coding

#Turn vcf into tapestri id format"
paste0(
  seqnames(gr_vcf), ":", start(gr_vcf), ":", 
  mcols(gr_vcf)[, "mcols.REF"], "/", mcols(gr_vcf)[, "mcols.ALT"]
) -> vcf_id

gr_vcf$id <- vcf_id

## Read Tapestri .h5 files
s1 <- read_h5(file.path(data_dir, "1-AML_23_05_01.dna+protein.h5"))
s4 <- read_h5(file.path(data_dir, "4-VEH-END_23_05_24.dna+protein.h5"))
s6 <- read_h5(file.path(data_dir, "6-IVO-END_23_05_24.dna+protein.h5"))
s12 <- read_h5(file.path(data_dir, "12-VENE-END_23_05_24.dna.h5"))
s14 <- read_h5(file.path(data_dir, "14-IVOVEN-END_23_05_24.dna.h5"))
s8 <- read_h5(file.path(data_dir, "8-AZA-END_DNA_23_07_28.dna.h5"))
s10 <- read_h5(file.path(data_dir, "10-IVOAZA-END_DNA_23_07_28.dna.h5"))
s2 <- read_h5(file.path(data_dir, "2-BASELINE_DNA_23_07_28.dna.h5"))

# Treatment map - vector mapping sample and treatment
c(s1 = "Primary sample", 
  s12 = "VEN", 
  s14 = "IVO+VEN",
  s4 = "VEH", 
  s6 = "IVO",
  s8 = "AZA", 
  s10 = "IVO+AZA", 
  s2 = "Baseline") -> map_treatment

# Store h5f in list
ls_h5 <- list(s1 = s1, 
              s4 = s4, 
              s6 = s6, 
              s12 = s12, 
              s14 = s14, 
              s8 = s8,
              s10 = s10,
              s2 = s2)

# -------------------------------------------------------------------------

#Intersect of primary single cell + bulk seq.

#Use Tapestri FLT3ITD call in bulk
gr_vcf[gr_vcf$id == "chr13:28608336:T/TAGGGTCATATTCTC"]$id <- "chr13:28608262:./ATCTGTAGCTGGCTAGGG" 

var1 <- get_variants(ls_h5[["s4"]], "GRanges") #Use VEH PDX for comparision

subsetByOverlaps(var1, gr_vcf)
subset(var1, id %in% gr_vcf$id) %>% #32 out of 33
  .$id -> intersect_id

names(intersect_id) <- paste0("Variant", 1:length(intersect_id))


#VENN DIMENSION: width700, height450 to save svg directly
svg("p_venn_tapestri_bulk.svg", width = 7, height = 4.5)
ggvenn::ggvenn(list("scDNA variants (unfiltered)" = var1$id, 
                    "Bulk targeted seq."=gr_vcf$id), 
               fill_color = ggsci::pal_npg()(2)
               )
dev.off()


read_assays_variants(h5f = ls_h5[["s4"]],
                     included_assays = "AF",
                     index_variants = get_variants_index(ls_h5[["s4"]], intersect_id),
                     format = "list")[[1]] -> m_vaf
data.frame(
  Variant = names(intersect_id),
  id = intersect_id,
  mean_vaf_s4 = apply(m_vaf, 1, mean)
  ) -> d_intersect

d_intersect <- left_join(as.data.frame(gr_vcf), d_intersect, by = "id")

# For the TET2 undetected variant, set percent mutated to 0
d_intersect[is.na(d_intersect$mean_vaf_s4), "mean_vaf_s4"] <- 0

d_intersect %>% 
  ggplot(aes(mcols.AF, mean_vaf_s4/100)) +
  geom_point() +
  ggpmisc::stat_poly_line(formula = y~x) +
  ggpmisc::stat_poly_eq(formula = y~x) +
  labs(x="PDX bulk VAF", y="Avg. PDX scDNA VAF", caption = "All variants (n=33)") +
  theme_classic(13) -> p_vaf_s4_v_bulk
p_vaf_s4_v_bulk

ggsave("p_vaf_s4_v_bulk.png", p_vaf_s4_v_bulk, width=5, height=5, dpi="retina", units="in")

# -------------------------------------------------------------------------

c(
  DNMT3A_R882H = "chr2:25457242:C/T",
  IDH1_R132H = "chr2:209113112:C/T",
  NPM1c = "chr5:170837543:C/CTCTG",
  WT1_R462W = "chr11:32413566:G/A",
  WT1_R434C = "chr11:32414251:G/A",
  FLT3_ITD = "chr13:28608262:./ATCTGTAGCTGGCTAGGG",
  KIT_H40Qfs = "chr4:55561718:ACCAT/A", #vaf oicr 2.958580 (chr4:55561719-55561722, CCAT del, varsome pathogenic
  DNMT3A_D531del = "chr2:25467481:CCGT/C", #vaf oicr 1.880878 (chr2:25467482-25467484, CGT del), varsome likely pathogenic
  EZH2_E404del = "chr7:148514997:CTCT/C", #Not called in oicr. (chr7:148514998-148515000, TCT del), varsome uncertain
  TET2_P29R = "chr4:106155185:C/G", #vaf oicr = 44.794189, varsome benign
  TET2_H924R = "chr4:106157870:A/G", #vaf oicr = 46.906637, varsome benign
  TET2_H1778R = "chr4:106197000:A/G", #vaf oicr = 45.254075, varsome benign
  TET2_Q1548del = "chr4:106196299:CCAG/C" #Not called in oicr. varsome uncertain
) -> int_var #13 out of the 14 protein coding variants


# A closer look at the 8 pathogenic + 1 unknown variants ----------------------------------------

int_var[c("DNMT3A_R882H", "WT1_R462W", "IDH1_R132H", "NPM1c", "TET2_Q1548del",
          "FLT3_ITD", "WT1_R434C", "DNMT3A_D531del", "KIT_H40Qfs")] -> pathogenic_var


# Filter cells ------------------------------------------------------------
var <- pathogenic_var
spikein <- "chr17:7577517:A/C" #TP53_I255S from PM150437
AF_TP53_co <- 10 #vaf cutoff to filter out spikein cells

# -------------------------------------------------------------------------

ls_sce <- list() 

for(i in names(ls_h5)){
  ls_sce[[i]] <- read_assays_variants(ls_h5[[i]],
                                      c("AF", "NGT", "GQ", "DP"),
                                      index_variants = get_variants_index(ls_h5[[i]], var), 
                                      format = "SingleCellExperiment"
                                      )
}

## Add in protein data for some samples
for(i in c("s1", "s4", "s6")){
  protein_raw <- read_protein_counts(ls_h5[[i]], format="matrix", normalization = "raw")
  protein_clr <- read_protein_counts(ls_h5[[i]], format="matrix", normalization = "clr")
  
  #Store protein names
  names_protein_raw <- gsub(" ", "", rownames(protein_raw))
  names_protein_clr <- gsub(" ", "", rownames(protein_clr))
  
  #Remove attributes
  attributes(protein_raw) <- attributes(protein_raw)[1]
  attributes(protein_clr) <- attributes(protein_clr)[1]

  #Make s.e. object
  protein_raw <- SummarizedExperiment(protein_raw, rowData = DataFrame(Protein = names_protein_raw))
  protein_clr <- SummarizedExperiment(protein_clr, rowData = DataFrame(Protein = names_protein_clr))
  
  altExp(ls_sce[[i]], "Protein_raw") <- protein_raw
  altExp(ls_sce[[i]], "Protein_clr") <- protein_clr
}

# Filter - keep only TP53 I255S WT cells for AZA & IVO+AZA
read_assays_variants(ls_h5[["s8"]], "AF", 
                     get_variants_index(ls_h5[["s8"]], spikein),
                     format="list")[[1]][1, ] -> AF_TP53_s8

read_assays_variants(ls_h5[["s10"]], "AF", 
                     get_variants_index(ls_h5[["s10"]], spikein),
                     format="list")[[1]][1, ] -> AF_TP53_s10

index_keep_s8 <- which(AF_TP53_s8 <= AF_TP53_co)
index_keep_s10 <- which(AF_TP53_s10 <= AF_TP53_co)

ls_sce[["s8"]] <- ls_sce[["s8"]][, index_keep_s8]
ls_sce[["s10"]] <- ls_sce[["s10"]][, index_keep_s10]

# Filter - remove NGT where GQ < 30

for(i in names(ls_sce)){
  low_GQ <- assays(ls_sce[[i]])[["GQ"]] < 30
  ngt <- assays(ls_sce[[i]])[["NGT"]]
  ngt[low_GQ] <- 3
  assays(ls_sce[[i]])[["NGT_filter_GQ30"]] <- ngt
}

# Filter - remove NGT where DP < 10

for(i in names(ls_sce)){
  low_DP <- assays(ls_sce[[i]])[["DP"]] < 10
  ngt <- assays(ls_sce[[i]])[["NGT_filter_GQ30"]]
  ngt[low_DP] <- 3
  assays(ls_sce[[i]])[["NGT_filter_GQ30_DP10"]] <- ngt
}

# Filter - remove cells where not all variants are genotyped

ls_sce_filtered <- list()
for(i in names(ls_sce)){
  ngt <- assays(ls_sce[[i]])[["NGT_filter_GQ30_DP10"]] 
  index <- apply(ngt, 2, function(x) any(x == 3))
  index <- which(!index) #which cells have no NGT==3
  
  ls_sce_filtered[[i]] <- ls_sce[[i]][, index]
}

#How many cells pre and post filter?
lapply(ls_sce, ncol) %>% unlist()
lapply(ls_sce_filtered, ncol) %>% unlist()


# Count mutant cells ------------------------------------------------------

count_cells_sce <- function(sce, ngt, percent_mutated=TRUE, percent_genotyped=TRUE){
  dt.ngt <- assays(sce)[[ngt]]
  NGT0 <- apply(dt.ngt, 1, function(x) sum(x == 0))
  NGT1 <- apply(dt.ngt, 1, function(x) sum(x == 1))
  NGT2 <- apply(dt.ngt, 1, function(x) sum(x == 2))
  NGT3 <- apply(dt.ngt, 1, function(x) sum(x == 3))
  m <- data.frame(Variant = rownames(sce),
                  NGT0 = NGT0, 
                  NGT1 = NGT1,
                  NGT2 = NGT2, 
                  NGT3 = NGT3)
  if (percent_mutated == TRUE) {
    m$percent_mutated <- with(m, (NGT1 + NGT2) * 100/(NGT0 + NGT1 + NGT2))
  }
  if (percent_genotyped == TRUE) {
    m$percent_genotyped <- with(m, (NGT0 + NGT1 + NGT2) * 100/(NGT0 + NGT1 + NGT2 + NGT3))
  }
  return(m)
}

# d_nCell <- lapply(ls_sce_filtered, function(x) count_cells_sce(x, "NGT_final"))
d_nCell <- lapply(ls_sce_filtered, function(x) count_cells_sce(x, "NGT_filter_GQ30_DP10"))

d_nCell %>% 
  bind_rows(.id = "Sample") -> d_nCell

d_nCell$Treatment <- plyr::mapvalues(d_nCell$Sample,
                                     names(map_treatment),
                                     map_treatment)

d_nCell$Treatment <- factor(d_nCell$Treatment,
                            levels = c("Primary sample", "Baseline", "VEH", "IVO", "VEN", "IVO+VEN", "AZA", "IVO+AZA")
                            )
# -------------------------------------------------------------------------
# Clonal architecture

c("DNMT3A_R882H", "WT1_R462W", "IDH1_R132H", "NPM1c", "FLT3_ITD", "WT1_R434C", 
  "TET2_Q1548del", "DNMT3A_D531del", "KIT_H40Qfs") -> variant_ordered

#Possible clones under infinite sites assumption
c(D = "1,0,0,0,0,0,0,0,0", #D
  DW = "1,1,0,0,0,0,0,0,0", #DW1
  DWI = "1,1,1,0,0,0,0,0,0", #DW1I
  DWIN = "1,1,1,1,0,0,0,0,0", #DW1IN
  DWINF = "1,1,1,1,1,0,0,0,0", #DW1INF
  DWINFW = "1,1,1,1,1,1,0,0,0", #DW1INFW2
  DWINFT = "1,1,1,1,1,1,1,0,0", #DW1INFW2T
  DWINFTD = "1,1,1,1,1,1,0,1,0", #DW1INFW2D
  DWINFTDK = "1,1,1,1,1,1,0,0,1" #DW1INFW2K
  ) -> clones_assumed

#Annotate clone of each cell. 
#Two annotations per cell. One with infinite sites assumption, one relaxed.
for(i in names(ls_sce_filtered)){
  
  #Read NGT of 9 variants
  ngt <- assays(ls_sce_filtered[[i]])[["NGT_filter_GQ30_DP10"]][variant_ordered, ] 
  
  #Paste NGT to define clone
  v_relaxed <- apply(ngt, 2, function(x) paste0(x, collapse = ","))
  v_relaxed <- gsub("2", "1", v_relaxed) #No need to keep homozygous mutant call. Turn it into 1.
  
  #Which cells fit infinite sites assumption model?
  v_assumed <- plyr::mapvalues(v_relaxed, clones_assumed,names(clones_assumed))
  index_unexplained <- v_relaxed %in% clones_assumed
  index_unexplained <- !index_unexplained
  v_assumed[index_unexplained] <- "Unexplained"
    
  #Which cells are doublets?
  #(1) Any cell carrying more than one of "TET2_Q1548del", "DNMT3A_D531del", "KIT_H40Qfs"
  v_relaxed %>% 
    str_split(",") %>% 
    lapply(function(x) sum(as.numeric(x[7:9])) > 1) %>% 
    unlist() -> is.doublet
  
  colData(ls_sce_filtered[[i]]) <- DataFrame(Clone_relaxed = v_relaxed,
                                             Clone_assumed = v_assumed, 
                                             Doublet = is.doublet)
}

#Remove doublets 
for(i in names(ls_sce_filtered)){
  index_doublet <- colData(ls_sce_filtered[[i]])["Doublet"][, 1]
  ls_sce_filtered[[i]] <- ls_sce_filtered[[i]][, !index_doublet] #Keep only non-doublets
}

#Final annotation of clone
for(i in names(ls_sce_filtered)){
  Clone_annotated <- vector(mode = "character", length = ncol(ls_sce_filtered[[i]]))
  Clone_relaxed <- colData(ls_sce_filtered[[i]])["Clone_relaxed"][, 1]
  
  # Clone9 - final mutation is KIT_H40Qfs
  index_clone9 <- Clone_relaxed %>% str_split(",") %>% lapply(function(x) x[9] == 1) %>% unlist()
  Clone_annotated[index_clone9] <- "Clone9"
  
  # Clone8 - final mutation is DNMT3A_D531del
  index_clone8 <- Clone_relaxed %>% str_split(",") %>% lapply(function(x) x[8] == 1) %>% unlist()
  Clone_annotated[index_clone8] <- "Clone8"
  
  # Clone7 - final mutation is DNMT3A_D531del
  index_clone7 <- Clone_relaxed %>% str_split(",") %>% lapply(function(x) x[7] == 1) %>% unlist()
  Clone_annotated[index_clone7] <- "Clone7"
  
  # Clone6 - final mutation is WT1_R434C
  index_clone6 <- Clone_relaxed %>% 
    str_split(",") %>% 
    lapply(function(x){
      all(
        x[6] == 1,
        x[7] == 0,
        x[8] == 0,
        x[9] == 0
        )
      }
    ) %>% unlist()
  Clone_annotated[index_clone6] <- "Clone6"
  
  # Clone5 - final mutation is FLT3_ITD
  index_clone5 <- Clone_relaxed %>% 
    str_split(",") %>% 
    lapply(function(x){
      all(
        x[5] == 1,
        x[6] == 0,
        x[7] == 0,
        x[8] == 0,
        x[9] == 0
      )
    }
    ) %>% unlist()
  Clone_annotated[index_clone5] <- "Clone5"
  
  # Clone4 - final mutation is NPM1c
  index_clone4 <- Clone_relaxed %>% 
    str_split(",") %>% 
    lapply(function(x){
      all(
        x[4] == 1,
        x[5] == 0,
        x[6] == 0,
        x[7] == 0,
        x[8] == 0,
        x[9] == 0
      )
    }
    ) %>% unlist()
  Clone_annotated[index_clone4] <- "Clone4"
  
  # Clone3 - final mutation is IDH1_R132H
  index_clone3 <- Clone_relaxed %>% 
    str_split(",") %>% 
    lapply(function(x){
      all(
        x[3] == 1,
        x[4] == 0,
        x[5] == 0,
        x[6] == 0,
        x[7] == 0,
        x[8] == 0,
        x[9] == 0
      )
    }
    ) %>% unlist()
  Clone_annotated[index_clone3] <- "Clone3"
  
  # Clone2 - final mutation is WT1_R462W
  index_clone2 <- Clone_relaxed %>% 
    str_split(",") %>% 
    lapply(function(x){
      all(
        x[2] == 1,
        x[3] == 0,
        x[4] == 0,
        x[5] == 0,
        x[6] == 0,
        x[7] == 0,
        x[8] == 0,
        x[9] == 0
      )
    }
    ) %>% unlist()
  Clone_annotated[index_clone2] <- "Clone2"
  
  # Clone1 - final mutation is WT1_R462W
  index_clone1 <- Clone_relaxed %>% 
    str_split(",") %>% 
    lapply(function(x){
      all(
        x[1] == 1,
        x[2] == 0,
        x[3] == 0,
        x[4] == 0,
        x[5] == 0,
        x[6] == 0,
        x[7] == 0,
        x[8] == 0,
        x[9] == 0
      )
    }
    ) %>% unlist()
  Clone_annotated[index_clone1] <- "Clone1"
  
  # Clone0 - WT for all
  index_clone0 <- Clone_relaxed %>% 
    str_split(",") %>% 
    lapply(function(x){
      all(
        x[1] == 0,
        x[2] == 0,
        x[3] == 0,
        x[4] == 0,
        x[5] == 0,
        x[6] == 0,
        x[7] == 0,
        x[8] == 0,
        x[9] == 0
      )
    }
    ) %>% unlist()
  Clone_annotated[index_clone0] <- "Clone0"
  
  #Change ColData
  colData(ls_sce_filtered[[i]]) <- DataFrame(Clone_relaxed = Clone_relaxed, Clone_annotated = Clone_annotated)
}

#Check if there are any unannotated cells
lapply(ls_sce_filtered, function(x) sum(colData(x)["Clone_annotated"][, 1] == "")) %>% unlist()

write_rds(ls_sce_filtered, "ls_sce_filtered.rds")


# -------------------------------------------------------------------------
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
#   [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggvenn_0.1.10               ggpmisc_0.4.7               ggpp_0.4.4                  chameleon_0.2-2            
# [5] ggpointdensity_0.1.0        RColorBrewer_1.1-3          tapestri.tools_0.0.9000     forcats_0.5.1              
# [9] stringr_1.5.0               purrr_0.3.4                 readr_2.1.1                 tidyr_1.2.1                
# [13] tibble_3.1.8                tidyverse_1.3.1             SingleCellExperiment_1.14.1 SummarizedExperiment_1.24.0
# [17] Biobase_2.54.0              MatrixGenerics_1.6.0        matrixStats_0.63.0          rtracklayer_1.52.1         
# [21] rhdf5_2.36.0                rentrez_1.2.3               ggridges_0.5.3              ggplot2_3.4.0              
# [25] GenomicRanges_1.46.1        GenomeInfoDb_1.30.1         IRanges_2.28.0              S4Vectors_0.32.4           
# [29] BiocGenerics_0.40.0         magrittr_2.0.3              dplyr_1.0.10                data.table_1.14.6          
# [33] compositions_2.0-4          maftools_2.8.05             vcfR_1.12.0                
# 
# loaded via a namespace (and not attached):
#   [1] colorspace_2.0-3         rjson_0.2.20             ellipsis_0.3.2           XVector_0.34.0           fs_1.5.2                
# [6] rstudioapi_0.13          farver_2.1.1             MatrixModels_0.5-0       fansi_1.0.3              lubridate_1.8.0         
# [11] xml2_1.3.3               splines_4.1.3            memuse_4.2-1             robustbase_0.95-0        knitr_1.39              
# [16] polynom_1.4-1            jsonlite_1.8.4           Rsamtools_2.8.0          broom_1.0.3              cluster_2.1.2           
# [21] dbplyr_2.2.1             compiler_4.1.3           httr_1.4.4               backports_1.4.1          assertthat_0.2.1        
# [26] Matrix_1.4-0             cli_3.4.1                quantreg_5.93            tools_4.1.3              gtable_0.3.1            
# [31] glue_1.6.2               GenomeInfoDbData_1.2.7   posterior_1.4.1.9000     Rcpp_1.0.9               cellranger_1.1.0        
# [36] vctrs_0.5.1              Biostrings_2.62.0        rhdf5filters_1.4.0       ape_5.6-2                nlme_3.1-155            
# [41] pinfsc50_1.2.0           tensorA_0.36.2           xfun_0.31                rvest_1.0.2              lifecycle_1.0.3         
# [46] restfulr_0.0.15          XML_3.99-0.13            DEoptimR_1.0-11          zlibbioc_1.40.0          MASS_7.3-55             
# [51] scales_1.2.1             ragg_1.2.2               hms_1.1.2                parallel_4.1.3           SparseM_1.81            
# [56] yaml_2.3.6               gridExtra_2.3            UpSetR_1.4.0             stringi_1.7.8            BiocIO_1.2.0            
# [61] checkmate_2.1.0          permute_0.9-7            BiocParallel_1.26.2      cmdstanr_0.5.0           systemfonts_1.0.4       
# [66] rlang_1.0.6              pkgconfig_2.0.3          bitops_1.0-7             distributional_0.3.0     lattice_0.20-45         
# [71] Rhdf5lib_1.14.2          labeling_0.4.2           GenomicAlignments_1.28.0 tidyselect_1.2.0         ggsci_2.9               
# [76] plyr_1.8.8               R6_2.5.1                 generics_0.1.3           DelayedArray_0.20.0      DBI_1.1.3               
# [81] pillar_1.8.1             haven_2.4.3              withr_2.5.0              mgcv_1.8-39              abind_1.4-5             
# [86] survival_3.2-13          RCurl_1.98-1.9           bayesm_3.1-4             modelr_0.1.8             crayon_1.5.2            
# [91] utf8_1.2.2               tzdb_0.2.0               readxl_1.3.1             vegan_2.6-2              reprex_2.0.1            
# [96] textshaping_0.3.6        munsell_0.5.0            viridisLite_0.4.1 


