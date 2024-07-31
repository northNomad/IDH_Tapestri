wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

library(tapestri.tools)
library(ggsci)
library(ComplexHeatmap)
library(maftools)
library(RColorBrewer)
library(ggupset)
library(maftools)


# -------------------------------------------------------------------------
data_dir <- "../../DataFiles/MixedPDX/" #Data directory
## Files needed to run this script:
## **All these files must be in the same directory called 'data_dir'
## (1) 'mPDX_ID_A.dna+protein.h5'
## (3) 'mpdxC-10_MB_90.dna+protein.h5'
## (4) 'mpdxD-18_MB_33_38.dna+protein.h5'
##
##  These .h5 files are scDNA results exported from MisisonBio's Tapestri Pipeline.


## What the script outputs:
## (1): Cell numbers from each clone that are used to make the mermaid plots.
##      This is not saved.

## (2): 'p_hm_supplementary_complexIDH1.svg' 
## (3): 'p_hm_supplementary_TP53Branching.svg'
## (4): 'p_hm_supplementary_IsozymeSwitching.svg'
##       These are Figure S3


sample_color <- ggsci::pal_npg()(9)
names(sample_color) <- c("PM325267", "PM248808", "PM165009", "PM305256", "PM246514", "PM160950",
                         "PM160345", "PM150437", "PM160053")

pdxa <- read_h5(file.path(data_dir, "mPDX_ID_A.dna+protein.h5"))
var_a <- get_variants(pdxa)


#New one where columns are split by samples
c("DNMT3A_R882H", "DNMT3A_R882C", "DNMT3A_R882S",
  "IDH1_R132H", "IDH1_R132C", "IDH1_R132S", 
  "NPM1c", "FLT3ITD_1", "JAK2_V617F", "EZH2_C552R", "EZH2_P527H",
  "RUNX1_R320Ter", "WT1_R462W", "NRAS_G13D") -> row_order_pdxa

c(
  TET2_L1721W = "chr4:106196829:T/G",
  DNMT3A_R882C = "chr2:25457243:G/A",
  NPM1c = "chr5:170837543:C/CTCTG",
  IDH1_R132H = "chr2:209113112:C/T",
  NRAS_G13D = "chr1:115258744:C/T",
  IDH1_R132C = "chr2:209113113:G/A",
  DNMT3A_R882S = "chr2:25457243:G/T",
  DNMT3A_R882H = "chr2:25457242:C/T",
  JAK2_V617F = "chr9:5073770:G/T",
  IDH1_R132S = "chr2:209113113:G/T",
  WT1_R462W = "chr11:32413566:G/A",
  EZH2_C552R = "chr7:148512024:A/G",
  EZH2_P527H = "chr7:148512098:G/T",
  RUNX1_R320Ter = "chr21:36171607:G/A"
) -> int_a



FLT3ITD_pdxa <- get_flt3itd(pdxa)
FLT3ITD_pdxa <- setNames(FLT3ITD_pdxa$id, paste0("FLT3ITD_", 1:nrow(FLT3ITD_pdxa)))

sce_a <- read_assays_variants(pdxa,
                              c("AF", "NGT", "GQ", "DP"), 
                              index_variants=get_variants_index(pdxa, c(int_a, FLT3ITD_pdxa)),
                              format="SingleCellExperiment"
) 
ls_sce <- list(a=sce_a)


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

# # Collapse all FLT3ITD variants
# 
for(i in names(ls_sce)){
  assays(ls_sce[[i]])[["NGT_filter_GQ30_DP10"]] %>%
    .[paste0("FLT3ITD_", 1:6), ] %>%
    apply(2, function(x) any(x == 1 | x == 2)) %>%
    which() -> index_FLT3ITD_positive
  
  assays(ls_sce[[i]])[["NGT_filter_GQ30_DP10"]]["FLT3ITD_1", ] <- 0
  assays(ls_sce[[i]])[["NGT_filter_GQ30_DP10"]]["FLT3ITD_1", ][index_FLT3ITD_positive] <- 1
  
  rn <- rownames(ls_sce[[i]])
  rn <- rn[!grepl("FLT3ITD", rn)]
  rn <- c(rn, "FLT3ITD_1")
  ls_sce[[i]] <- ls_sce[[i]][rn, ]
}

m_ngt_a <- assays(ls_sce$a)[["NGT_filter_GQ30_DP10"]] 

m_ngt_a %>% 
  apply(c(1, 2), function(x) gsub(2, 1, x)) %>%
  apply(c(1, 2), as.numeric) -> m_ngt_a

m_ngt_a <- m_ngt_a[c("TET2_L1721W", "DNMT3A_R882H", "DNMT3A_R882C", "DNMT3A_R882S",
                     "IDH1_R132H", "IDH1_R132C", "IDH1_R132S", 
                     "NPM1c", "FLT3ITD_1", "JAK2_V617F", "EZH2_C552R", "EZH2_P527H",
                     "RUNX1_R320Ter", "WT1_R462W", "NRAS_G13D"), ]


m_ngt_a <- t(m_ngt_a)
m_ngt_a <- as.data.table(m_ngt_a)

#Is it PM160345?
m_ngt_a$is_PM160345 <- ifelse(m_ngt_a$DNMT3A_R882C %in% c(1, 2) &
                                m_ngt_a$DNMT3A_R882H %in% c(0, 3) &
                                m_ngt_a$DNMT3A_R882S %in% c(0, 3) &
                                m_ngt_a$IDH1_R132C %in% c(0, 3) &
                                m_ngt_a$IDH1_R132S %in% c(0, 3) &
                                m_ngt_a$JAK2_V617F %in% c(0, 3) &
                                m_ngt_a$EZH2_C552R %in% c(0, 3) &
                                m_ngt_a$EZH2_P527H %in% c(0, 3) &
                                m_ngt_a$RUNX1_R320Ter %in% c(0, 3) &
                                m_ngt_a$NRAS_G13D %in% c(0, 3) &
                                m_ngt_a$WT1_R462W %in% c(0, 3) 
                                , 1, 0)

#Is it PM165009
m_ngt_a$is_PM165009 <- ifelse((m_ngt_a$JAK2_V617F %in% c(1, 2) |
                                 m_ngt_a$EZH2_C552R %in% c(1, 2) |
                                 m_ngt_a$EZH2_P527H %in% c(1, 2) 
                               ) & 
                                m_ngt_a$IDH1_R132H %in% c(0, 3) &
                                m_ngt_a$IDH1_R132S %in% c(0, 3) &
                                m_ngt_a$NPM1c %in% c(0, 3) &
                                m_ngt_a$NRAS_G13D %in% c(0, 3) &
                                m_ngt_a$WT1_R462W %in% c(0, 3)  &
                                m_ngt_a$DNMT3A_R882H %in% c(0, 3) &
                                m_ngt_a$DNMT3A_R882S %in% c(0, 3) &
                                m_ngt_a$DNMT3A_R882C %in% c(0, 3) &
                                m_ngt_a$TET2_L1721W %in% c(0, 3) 
                              , 1, 0)
#Is it PM325267
m_ngt_a$is_PM325267 <- ifelse((m_ngt_a$TET2_L1721W %in% c(1, 2) |
                                 m_ngt_a$NRAS_G13D %in% c(1, 2)) &
                                m_ngt_a$IDH1_R132H %in% c(0, 3) &
                                m_ngt_a$DNMT3A_R882H %in% c(0, 3) &
                                m_ngt_a$DNMT3A_R882S %in% c(0, 3) &
                                m_ngt_a$DNMT3A_R882C %in% c(0, 3) &
                                m_ngt_a$JAK2_V617F %in% c(0, 3) &
                                m_ngt_a$EZH2_C552R %in% c(0, 3) &
                                m_ngt_a$EZH2_P527H %in% c(0, 3) &
                                m_ngt_a$RUNX1_R320Ter %in% c(0, 3)
                              , 1, 0)
#Is it PM246514
m_ngt_a$is_PM246514 <- ifelse((m_ngt_a$DNMT3A_R882H %in% c(1, 2) |
                                 m_ngt_a$WT1_R462W %in% c(1, 2)
                               ) &
                                m_ngt_a$DNMT3A_R882C %in% c(0, 3) &
                                m_ngt_a$DNMT3A_R882S %in% c(0, 3) &
                                m_ngt_a$IDH1_R132C %in% c(0, 3) &
                                m_ngt_a$IDH1_R132S %in% c(0, 3) &
                                m_ngt_a$TET2_L1721W %in% c(0, 3) &
                                m_ngt_a$NRAS_G13D %in% c(0, 3) &
                                m_ngt_a$JAK2_V617F %in% c(0, 3) &
                                m_ngt_a$EZH2_C552R %in% c(0, 3) &
                                m_ngt_a$EZH2_P527H %in% c(0, 3) 
                              , 1, 0)
#Is it PM160950
m_ngt_a$is_PM160950 <- ifelse(m_ngt_a$DNMT3A_R882S %in% c(1, 2) &
                                m_ngt_a$DNMT3A_R882H %in% c(0, 3) &
                                m_ngt_a$DNMT3A_R882C %in% c(0, 3) &
                                m_ngt_a$IDH1_R132H %in% c(0, 3) &
                                m_ngt_a$IDH1_R132C %in% c(0, 3) &
                                m_ngt_a$TET2_L1721W %in% c(0, 3) &
                                m_ngt_a$NRAS_G13D %in% c(0, 3) &
                                m_ngt_a$JAK2_V617F %in% c(0, 3) &
                                m_ngt_a$EZH2_C552R %in% c(0, 3) &
                                m_ngt_a$EZH2_P527H %in% c(0, 3) &
                                m_ngt_a$NPM1c %in% c(0, 3) 
                              , 1, 0)
#Is it PM305256
ifelse(m_ngt_a$DNMT3A_R882S %in% c(0, 3) &
         m_ngt_a$DNMT3A_R882H %in% c(0, 3) &
         m_ngt_a$DNMT3A_R882C %in% c(0, 3) &
         m_ngt_a$NPM1c %in% c(0, 3) &
         m_ngt_a$WT1_R462W %in% c(0, 3) &
         m_ngt_a$IDH1_R132C %in% c(0, 3) &
         m_ngt_a$IDH1_R132S %in% c(0, 3) &
         m_ngt_a$NRAS_G13D %in% c(0, 3) &
         m_ngt_a$FLT3ITD_1 %in% c(0, 3) &
         m_ngt_a$JAK2_V617F %in% c(0, 3) &
         m_ngt_a$EZH2_C552R %in% c(0, 3) &
         m_ngt_a$EZH2_P527H %in% c(0, 3) &
         m_ngt_a$RUNX1_R320Ter %in% c(0, 3) &
         m_ngt_a$IDH1_R132H %in% c(1, 2)
       , 1, 0) -> m_ngt_a$is_PM305256

####
m_ngt_a$is_multiplet <- ifelse(apply(m_ngt_a[, 16:21], 1, sum) > 1, 1, 0)
m_ngt_a$is_unknown <- ifelse(apply(m_ngt_a[, 16:21], 1, sum) == 0, 1, 0)

sum(m_ngt_a$is_multiplet)
sum(m_ngt_a$is_unknown)

m_ngt_a <- m_ngt_a %>% subset(is_multiplet == 0 & is_unknown == 0) 



m_ngt_a %>% 
  .[, 16:21] %>% 
  apply(1, function(x) {
    names(which(x == 1))
  }) %>% 
  unlist() -> m_ngt_a$sample

m_ngt_a %>% 
  group_by(sample) %>% 
  group_map(function(x, y){
    apply(x, 2, function(i) sum(i == 1))
  }) -> n_pdxa

n_pdxa %>% 
  lapply(function(x) x[x != 0])

## How many cells per clone 
## PM160345
m_ngt_a %>% 
  subset(sample == "is_PM160345") %>% 
  mutate(clone=paste0(DNMT3A_R882C, IDH1_R132H, NPM1c, FLT3ITD_1)) %>% 
  .$clone %>% 
  sapply(function(x){
    z <- str_split(x, "", simplify = TRUE) 
    max(which(z == "1"))
    }
  ) %>% 
  table()
  
## PM246514
m_ngt_a %>% 
  subset(sample == "is_PM246514") %>% 
  mutate(clone=paste0(DNMT3A_R882H, WT1_R462W, IDH1_R132H, NPM1c, FLT3ITD_1)) %>% 
  .$clone %>% 
  sapply(function(x){
    z <- str_split(x, "", simplify = TRUE) 
    max(which(z == "1"))
  }
  ) %>% 
  table()

## PM160950
m_ngt_a %>% 
  subset(sample == "is_PM160950") %>% 
  mutate(clone=paste0(DNMT3A_R882S, IDH1_R132S, FLT3ITD_1, RUNX1_R320Ter)) %>% 
  .$clone %>% 
  sapply(function(x){
    z <- str_split(x, "", simplify = TRUE) 
    max(which(z == "1"))
  }
  ) %>% 
  table()

## PM165009
m_ngt_a %>% 
  subset(sample == "is_PM165009") %>% 
  mutate(clone=paste0(JAK2_V617F, EZH2_C552R, EZH2_P527H, RUNX1_R320Ter, IDH1_R132C)) %>% 
  .$clone %>% 
  sapply(function(x){
    z <- str_split(x, "", simplify = TRUE) 
    max(which(z == "1"))
  }
  ) %>% 
  table()

## PM325267
m_ngt_a %>% 
  subset(sample == "is_PM325267") %>% 
  mutate(clone=paste0(IDH1_R132C, NPM1c, NRAS_G13D)) %>% 
  .$clone %>% 
  sapply(function(x){
    z <- str_split(x, "", simplify = TRUE) 
    max(which(z == "1"))
  }
  ) %>% 
  table()

### Calculate size for mermaid plot
## PM305256 - 2
##PM160345
r_PM160345 <- c(0, 12, 34, 12)
r_PM160345 <- r_PM160345 * 2 / sum(r_PM160345)
r_PM160345

##PM246514
r_PM246514 <- c(0, 4, 550, 270, 359)
r_PM246514 <- r_PM246514 * 2 / sum(r_PM246514)
r_PM246514

##PM160950
r_PM160950 <- c(7, 61, 16, 12)
r_PM160950 <- r_PM160950 * 2 / sum(r_PM160950)
r_PM160950

##PM165009
r_PM165009 <- c(0, 1, 2, 11, 1)
r_PM165009 <- r_PM165009 * 2 / sum(r_PM165009)
r_PM165009

##PM325267
r_PM325267 <- c(1, 7, 3)
r_PM325267 <- r_PM325267 * 2 / sum(r_PM325267)
r_PM325267



# mPDX2 --------------------------------------------------------------------

### mPDX2
pdxc <- read_h5(file.path(data_dir, "mpdxC-10_MB_90.dna+protein.h5"))
var_c <- get_variants(pdxc)
var_c %>% subset(filtered == "00")

c(
  DNMT3A_R882C = "chr2:25457243:G/A",
  NPM1c = "chr5:170837543:C/CTCTG",
  IDH1_R132H = "chr2:209113112:C/T",
  TP53_I255S = "chr17:7577517:A/C"
  
) -> int_c

FLT3ITD_pdxc <- get_flt3itd(pdxc)
FLT3ITD_pdxc <- setNames(FLT3ITD_pdxc$id, paste0("FLT3ITD_", 1:nrow(FLT3ITD_pdxc)))

sce_c <- read_assays_variants(pdxc,
                              c("AF", "NGT", "GQ", "DP"), 
                              index_variants=get_variants_index(pdxc, c(int_c, FLT3ITD_pdxc)),
                              format="SingleCellExperiment"
) 
ls_sce <- list(c=sce_c)

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

# # Collapse all FLT3ITD variants
# 
for(i in names(ls_sce)){
  assays(ls_sce[[i]])[["NGT_filter_GQ30_DP10"]] %>%
    .[paste0("FLT3ITD_", 1:3), ] %>%
    apply(2, function(x) any(x == 1 | x == 2)) %>%
    which() -> index_FLT3ITD_positive
  
  assays(ls_sce[[i]])[["NGT_filter_GQ30_DP10"]]["FLT3ITD_1", ] <- 0
  assays(ls_sce[[i]])[["NGT_filter_GQ30_DP10"]]["FLT3ITD_1", ][index_FLT3ITD_positive] <- 1
  
  rn <- rownames(ls_sce[[i]])
  rn <- rn[!grepl("FLT3ITD", rn)]
  rn <- c(rn, "FLT3ITD_1")
  ls_sce[[i]] <- ls_sce[[i]][rn, ]
}

m_ngt_c <- assays(ls_sce$c)[["NGT_filter_GQ30_DP10"]] 

m_ngt_c %>% 
  apply(c(1, 2), function(x) gsub(2, 1, x)) %>%
  apply(c(1, 2), as.numeric) -> m_ngt_c


m_ngt_c <- m_ngt_c[c("DNMT3A_R882C", "IDH1_R132H", "NPM1c", "FLT3ITD_1", "TP53_I255S"), ]
m_ngt_c <- as.data.table(t(m_ngt_c))

m_ngt_c$is_PM150437 <- ifelse(m_ngt_c$TP53_I255S %in% c(1, 2) &
                                m_ngt_c$DNMT3A_R882C %in% c(0, 3) &
                                m_ngt_c$IDH1_R132H %in% c(0, 3) &
                                m_ngt_c$NPM1c %in% c(0, 3) &
                                m_ngt_c$FLT3ITD_1 %in% c(0, 3)
                                , 1, 0)
m_ngt_c$is_PM160345 <- ifelse((m_ngt_c$DNMT3A_R882C %in% c(1, 2) |
                                 m_ngt_c$IDH1_R132H %in% c(1, 2) |
                                 m_ngt_c$NPM1c %in% c(1, 2) |
                                 m_ngt_c$FLT3ITD_1 %in% c(1, 2)) &
                                m_ngt_c$TP53_I255S %in% c(0, 3)
                                , 1, 0)

m_ngt_c$is_multiplet <- ifelse(apply(m_ngt_c[, 6:7], 1, sum) > 1, 1, 0)
m_ngt_c$is_unknown <- ifelse(apply(m_ngt_c[, 6:7], 1, sum) == 0, 1, 0)

sum(m_ngt_c$is_multiplet)
sum(m_ngt_c$is_unknown)

m_ngt_c <- m_ngt_c %>% subset(is_multiplet == 0 & is_unknown == 0) 

m_ngt_c %>% 
  .[, 6:7] %>% 
  apply(1, function(x) {
    names(which(x == 1))
  }) %>% 
  unlist() -> m_ngt_c$sample

table(m_ngt_c$sample)

## How many cells per clone 
## PM160345
m_ngt_c %>% 
  subset(sample == "is_PM160345") %>% 
  mutate(clone=paste0(DNMT3A_R882C, IDH1_R132H, NPM1c, FLT3ITD_1)) %>% 
  .$clone %>% 
  sapply(function(x){
    z <- str_split(x, "", simplify = TRUE) 
    max(which(z == "1"))
  }
  ) %>% 
  table()


# mPDX3 --------------------------------------------------------------------

pdxd <- read_h5(file.path(data_dir, "mpdxD-18_MB_33_38.dna+protein.h5"))

var_d <- get_variants(pdxd)


c(
  DNMT3A_R882C = "chr2:25457243:G/A",
  NPM1c = "chr5:170837543:C/CTCTG",
  IDH1_R132H = "chr2:209113112:C/T",
  IDH2_R140Q = "chr15:90631934:C/T"
) -> int_d

FLT3ITD_pdxd <- get_flt3itd(pdxd)
FLT3ITD_pdxd <- setNames(FLT3ITD_pdxd$id, paste0("FLT3ITD_", 1:nrow(FLT3ITD_pdxd)))

sce_d <- read_assays_variants(pdxd,
                              c("AF", "NGT", "GQ", "DP"), 
                              index_variants=get_variants_index(pdxd, c(int_d, FLT3ITD_pdxd)),
                              format="SingleCellExperiment"
) 
ls_sce <- list(d=sce_d)

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

# # Collapse all FLT3ITD variants
# 
for(i in names(ls_sce)){
  assays(ls_sce[[i]])[["NGT_filter_GQ30_DP10"]] %>%
    .[paste0("FLT3ITD_", 1:18), ] %>%
    apply(2, function(x) any(x == 1 | x == 2)) %>%
    which() -> index_FLT3ITD_positive
  
  assays(ls_sce[[i]])[["NGT_filter_GQ30_DP10"]]["FLT3ITD_1", ] <- 0
  assays(ls_sce[[i]])[["NGT_filter_GQ30_DP10"]]["FLT3ITD_1", ][index_FLT3ITD_positive] <- 1
  
  rn <- rownames(ls_sce[[i]])
  rn <- rn[!grepl("FLT3ITD", rn)]
  rn <- c(rn, "FLT3ITD_1")
  ls_sce[[i]] <- ls_sce[[i]][rn, ]
}


m_ngt_d <- assays(ls_sce$d)[["NGT_filter_GQ30_DP10"]] 

m_ngt_d %>% 
  apply(c(1, 2), function(x) gsub(2, 1, x)) %>%
  apply(c(1, 2), as.numeric) -> m_ngt_d

m_ngt_d <- m_ngt_d[c("DNMT3A_R882C", "IDH1_R132H", "NPM1c", 
                     "FLT3ITD_1", "IDH2_R140Q"), ]
m_ngt_d <- as.data.table(t(m_ngt_d))

m_ngt_d$is_PM160345 <- ifelse((m_ngt_d$DNMT3A_R882C %in% c(1, 2) |
                                m_ngt_d$IDH1_R132H %in% c(1, 2)) &
                                m_ngt_d$IDH2_R140Q %in% c(0, 3), 
                              1, 0
                                )

m_ngt_d$is_PM160053 <- ifelse(m_ngt_d$IDH2_R140Q %in% c(1, 2) &
                                m_ngt_d$DNMT3A_R882C %in% c(0, 3) &
                                m_ngt_d$IDH1_R132H %in% c(0, 3), 
                              1, 0
                              )

m_ngt_d$is_multiplet <- ifelse(apply(m_ngt_d[, 6:7], 1, sum) > 1, 1, 0)
m_ngt_d$is_unknown <- ifelse(apply(m_ngt_d[, 6:7], 1, sum) == 0, 1, 0)

sum(m_ngt_d$is_multiplet)
sum(m_ngt_d$is_unknown)

m_ngt_d <- m_ngt_d %>% subset(is_multiplet == 0 & is_unknown == 0) 

m_ngt_d %>% 
  .[, 6:7] %>% 
  apply(1, function(x) {
    names(which(x == 1))
  }) %>% 
  unlist() -> m_ngt_d$sample

table(m_ngt_d$sample)

## How many cells per clone 
## PM160345
m_ngt_d %>% 
  subset(sample == "is_PM160345") %>% 
  mutate(clone=paste0(DNMT3A_R882C, IDH1_R132H, NPM1c, FLT3ITD_1)) %>% 
  .$clone %>% 
  sapply(function(x){
    z <- str_split(x, "", simplify = TRUE) 
    max(which(z == "1"))
  }
  ) %>% 
  table()

## PM160053
m_ngt_d %>% 
  subset(sample == "is_PM160053") %>% 
  mutate(clone=paste0(IDH2_R140Q, NPM1c, FLT3ITD_1)) %>% 
  .$clone %>% 
  sapply(function(x){
    z <- str_split(x, "", simplify = TRUE) 
    max(which(z == "1"))
  }
  ) %>% 
  table()


# NGT heatmap for supplementary -------------------------------------------


### Complex IDH1
var_PM305256 <- c("IDH1_R132H")
var_PM246514 <- c("DNMT3A_R882H", "WT1_R462W", "IDH1_R132H", "NPM1c", "FLT3ITD_1")
var_PM160950 <- c("DNMT3A_R882S", "IDH1_R132S", "FLT3ITD_1", "RUNX1_R320Ter") 
var_PM160345 <- c("DNMT3A_R882C", "IDH1_R132H", "NPM1c", "FLT3ITD_1")  
var_PM165009 <- c("JAK2_V617F", "EZH2_C552R", "EZH2_P527H", "RUNX1_R320Ter", "IDH1_R132C")
var_PM325267 <- c("IDH1_R132C", "NPM1c", "NRAS_G13D")

#Italic font
pdxa_row_labels <- c(
  DNMT3A_R882H="*DNMT3A*^R882H", 
  DNMT3A_R882C="*DNMT3A*^R882C", 
  DNMT3A_R882S="*DNMT3A*^R882S",
  IDH1_R132H="*IDH1*^R132H",
  IDH1_R132C="*IDH1*^R132C",
  IDH1_R132S="*IDH1*^R132S",
  NPM1c="*NPM1c*^+",
  FLT3ITD_1="*FLT3*-ITD",
  JAK2_V617F="*JAK2*^V617F",
  EZH2_C552R="*JAK2*^C552R",
  EZH2_P527H="*JAK2*^P527H",
  RUNX1_R320Ter="*RUNX1*^R320*",
  WT1_R462W="*WT1*^R467W",
  NRAS_G13D="*NRAS*^G13D"
)

m_ngt_a$sample %>% table

## PM305256
m_ngt_a %>% 
  subset(is_PM305256 == 1) %>% 
  .[, 2:15] %>% 
  t() %>% 
  .[unique(c(var_PM305256, row_order_pdxa)), ] %>% 
  Heatmap(cluster_rows = FALSE, 
          use_raster = FALSE,
          col=c("#FEE0D2", "#A50F15", "grey70"),
          show_heatmap_legend = FALSE,
          column_title="PM305256 (6,039 cells)",
          row_split=c(rep("Mutant", length(var_PM305256)), rep("Wildtype", 14 - length(var_PM305256))),
          row_labels = gt_render(pdxa_row_labels[unique(c(var_PM305256, row_order_pdxa))])
          ) -> hm_PM305256

## PM246514
m_ngt_a %>% 
  subset(is_PM246514 == 1) %>% 
  .[, 2:15] %>% 
  t() %>% 
  .[unique(c(var_PM246514, row_order_pdxa)), ] %>%
  Heatmap(cluster_rows = FALSE,
          use_raster = FALSE,
          col=c("#FEE0D2", "#A50F15", "grey70"),
          show_heatmap_legend = FALSE,
          column_title="PM246514 (1,183 cells)",
          row_split=c(rep("Mutant", length(var_PM246514)), rep("Wildtype", 14 - length(var_PM246514))),
          row_labels = gt_render(pdxa_row_labels[unique(c(var_PM246514, row_order_pdxa))])
  ) -> hm_PM246514

## PM160950
m_ngt_a %>% 
  subset(is_PM160950 == 1) %>% 
  .[, 2:15] %>% 
  t() %>% 
  .[unique(c(var_PM160950, row_order_pdxa)), ] %>%
  Heatmap(cluster_rows = FALSE, 
          use_raster = FALSE,
          col=c("#FEE0D2", "#A50F15", "grey70"),
          show_heatmap_legend = FALSE,
          column_title="PM160950 (96 cells)",
          row_split=c(rep("Mutant", length(var_PM160950)), rep("Wildtype", 14 - length(var_PM160950))),
          row_labels = gt_render(pdxa_row_labels[unique(c(var_PM160950, row_order_pdxa))])
  ) -> hm_PM160950


## PM160345
m_ngt_a %>% 
  subset(is_PM160345 == 1) %>% 
  .[, 2:15] %>% 
  t() %>% 
  .[unique(c(var_PM160345, row_order_pdxa)), ] %>%
  Heatmap(cluster_rows = FALSE, 
          use_raster = FALSE,
          col=c("#FEE0D2", "#A50F15", "grey70"),
          show_heatmap_legend = FALSE,
          column_title="PM160345 (58 cells)",
          row_split=c(rep("Mutant", length(var_PM160345)), rep("Wildtype", 14 - length(var_PM160345))),
          row_labels = gt_render(pdxa_row_labels[unique(c(var_PM160345, row_order_pdxa))])
  ) -> hm_PM160345

## PM165009
m_ngt_a %>% 
  subset(is_PM165009 == 1) %>% 
  .[, 2:15] %>% 
  t() %>% 
  .[unique(c(var_PM165009, row_order_pdxa)), ] %>%
  Heatmap(cluster_rows = FALSE, 
          use_raster = FALSE,
          col=c("#FEE0D2", "#A50F15", "grey70"),
          show_heatmap_legend = FALSE,
          column_title="PM165009 (15 cells)",
          row_split=c(rep("Mutant", length(var_PM165009)), rep("Wildtype", 14 - length(var_PM165009))),
          row_labels = gt_render(pdxa_row_labels[unique(c(var_PM165009, row_order_pdxa))])
  ) -> hm_PM165009

## PM325267
m_ngt_a %>% 
  subset(is_PM325267 == 1) %>% 
  .[, 2:15] %>% 
  t() %>% 
  .[unique(c(var_PM325267, row_order_pdxa)), ] %>%
  Heatmap(cluster_rows = FALSE, 
          use_raster = FALSE,
          col=c("#FEE0D2", "#A50F15", "grey70"),
          show_heatmap_legend = FALSE,
          column_title="PM325267 (11 cells)",
          row_split=c(rep("Mutant", length(var_PM325267)), rep("Wildtype", 14 - length(var_PM325267))),
          row_labels = gt_render(pdxa_row_labels[unique(c(var_PM325267, row_order_pdxa))])
  ) -> hm_PM325267


grob_list_pdxa <- list(hm_PM305256=hm_PM305256,
                       hm_PM246514=hm_PM246514,
                       hm_PM160950=hm_PM160950,
                       hm_PM160345=hm_PM160345,
                       hm_PM165009=hm_PM165009,
                       hm_PM325267=hm_PM325267
                       )
grob_list_pdxa <- lapply(grob_list_pdxa, function(x) grid.grabExpr(draw(x)))

dev.off()
p_hm_supplementary_complexIDH1 <- gridExtra::arrangeGrob(grobs=grob_list_pdxa, nrow=2)


### Save
svg("p_hm_supplementary_complexIDH1.svg", width = 15, height = 10)
grid.draw(p_hm_supplementary_complexIDH1)
dev.off()


# -------------------------------------------------------------------------



### TP53 Branching
var_PM160345 <- c("DNMT3A_R882C", "IDH1_R132H", "NPM1c", "FLT3ITD_1")  
var_PM150437 <- c("TP53_I255S")

row_order_pdxc <- c("DNMT3A_R882C", "IDH1_R132H", "NPM1c", "FLT3ITD_1", "TP53_I255S")

#Italic font
pdxc_row_labels <- c(
  DNMT3A_R882C="*DNMT3A*^R882C", 
  IDH1_R132H="*IDH1*^R132H",
  NPM1c="*NPM1c*^+",
  FLT3ITD_1="*FLT3*-ITD",
  TP53_I255S="*TP53*^I255S"
)


m_ngt_c$sample %>% table

## PM160345
m_ngt_c %>% 
  subset(is_PM160345 == 1) %>% 
  .[, 1:5] %>% 
  t() %>% 
  .[unique(c(var_PM160345, row_order_pdxc)), ] %>% 
  Heatmap(cluster_rows = FALSE, 
          use_raster = FALSE,
          col=c("#FEE0D2", "#A50F15", "grey70"),
          show_heatmap_legend = FALSE,
          column_title="PM160345 (3,720 cells)",
          row_split=c(rep("Mutant", length(var_PM160345)), rep("Wildtype", 5 - length(var_PM160345))),
          row_labels = gt_render(pdxc_row_labels[unique(c(var_PM160345, row_order_pdxc))])
  ) -> hm_pdxc_PM160345

## PM150437
m_ngt_c %>% 
  subset(is_PM150437 == 1) %>% 
  .[, 1:5] %>% 
  t() %>% 
  .[unique(c(var_PM150437, row_order_pdxc)), ] %>% 
  Heatmap(cluster_rows = FALSE, 
          use_raster = FALSE,
          col=c("#FEE0D2", "#A50F15", "grey70"),
          show_heatmap_legend = FALSE,
          column_title="PM150437 (3,156 cells)",
          row_split=c(rep("Mutant", length(var_PM150437)), rep("Wildtype", 5 - length(var_PM150437))),
          row_labels = gt_render(pdxc_row_labels[unique(c(var_PM150437, row_order_pdxc))])
  ) -> hm_pdxc_PM150437


grob_list_pdxc <- list(hm_pdxc_PM160345=hm_pdxc_PM160345,
                       hm_pdxc_PM150437=hm_pdxc_PM150437)
  
grob_list_pdxc <- lapply(grob_list_pdxc, function(x) grid.grabExpr(draw(x)))

dev.off()
p_hm_supplementary_TP53Branching <- gridExtra::arrangeGrob(grobs=grob_list_pdxc, nrow=1)


### Save
svg("p_hm_supplementary_TP53Branching.svg", width = 10, height = 5)
grid.draw(p_hm_supplementary_TP53Branching)
dev.off()


# -------------------------------------------------------------------------


### Isozyme Switching

var_PM160345 <- c("DNMT3A_R882C", "IDH1_R132H", "NPM1c", "FLT3ITD_1")  
var_PM160053 <- c("IDH2_R140Q", "NPM1c", "FLT3ITD_1")

row_order_pdxd <- c("DNMT3A_R882C", "IDH1_R132H", "NPM1c", "FLT3ITD_1", "IDH2_R140Q")

#Italic font
pdxd_row_labels <- c(
  DNMT3A_R882C="*DNMT3A*^R882C", 
  IDH1_R132H="*IDH1*^R132H",
  NPM1c="*NPM1c*^+",
  FLT3ITD_1="*FLT3*-ITD",
  IDH2_R140Q="*IDH2*^R140Q"
)


m_ngt_d$sample %>% table


## PM160345
m_ngt_d %>% 
  subset(is_PM160345 == 1) %>% 
  .[, 1:5] %>% 
  t() %>% 
  .[unique(c(var_PM160345, row_order_pdxd)), ] %>% 
  Heatmap(cluster_rows = FALSE, 
          use_raster = FALSE,
          col=c("#FEE0D2", "#A50F15", "grey70"),
          show_heatmap_legend = FALSE,
          column_title="PM160345 (3,046 cells)",
          row_split=c(rep("Mutant", length(var_PM160345)), rep("Wildtype", 5 - length(var_PM160345))),
          row_labels = gt_render(pdxd_row_labels[unique(c(var_PM160345, row_order_pdxd))])
  ) -> hm_pdxd_PM160345

## PM160053
m_ngt_d %>% 
  subset(is_PM160053 == 1) %>% 
  .[, 1:5] %>% 
  t() %>% 
  .[unique(c(var_PM160053, row_order_pdxd)), ] %>% 
  Heatmap(cluster_rows = FALSE, 
          use_raster = FALSE,
          col=c("#FEE0D2", "#A50F15", "grey70"),
          show_heatmap_legend = FALSE,
          column_title="PM160053 (818 cells)",
          row_split=c(rep("Mutant", length(var_PM160053)), rep("Wildtype", 5 - length(var_PM160053))),
          row_labels = gt_render(pdxd_row_labels[unique(c(var_PM160053, row_order_pdxd))])
  ) -> hm_pdxd_PM160053

grob_list_pdxd <- list(hm_pdxd_PM160345=hm_pdxd_PM160345,
                       hm_pdxd_PM160053=hm_pdxd_PM160053)

grob_list_pdxd <- lapply(grob_list_pdxd, function(x) grid.grabExpr(draw(x)))


p_hm_supplementary_IsozymeSwitching <- gridExtra::arrangeGrob(grobs=grob_list_pdxd, nrow=1)


### Save
svg("p_hm_supplementary_IsozymeSwitching.svg", width = 10, height = 5)
grid.draw(p_hm_supplementary_IsozymeSwitching)
dev.off()

