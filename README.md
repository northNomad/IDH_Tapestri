**Patient samples used in this study**

Human AML samples from peripheral blood or BM biopsies were obtained with written consent according to procedures approved by the University Health Network ethics committee and in accordance with the Declaration of Helsinki. The samples were kept in liquid nitrogen at the Princess Margaret Leukemia Tissue Bank. Samples containing >80% blasts with a mutant IDH1 or IDH2 VAF of >40% that were collected prior to any exposure to IDH inhibitors were included. More information on patient characteristics can be found in Supplementary Table S1.

**Generation of IDH1-mutated patient derived xenograft model**

All animal experiments were done in accordance with institutional guidelines approved by the University Health Network Animal Care Committee. Four to six-week-old female NSG mice (NOD/SCID/IL2rg-/-) (RRID:IMSR_JAX:005557) were sub-lethally irradiated with 255 rads. The next day, we transplanted 1 million (unless otherwise noted) CD3 depleted (STEMCELL Technologies, Cat# 17851) primary cells into each animal by tail vein or intrafemoral injection. The mice were kept under Enrofloxacin (Selleckchem, Cat# S3059) (50 ug/mL; drinking water) for 14 days after irradiation.

**Drug administration for animal experiments**

Ivosidenib (gift from Servier Pharmaceuticals) and enasidenib (Selleckchem, Cat# S8205) stock solutions were prepared by re-suspending the drug in vehicle (0.5% methylcellulose + 0.2% Tween80) at a final concentration of 50 mg/mL and 7.5 mg/mL, respectively. Venetoclax (Selleckchem, Cat# S8048) stock solution was prepared by re-suspending the drug in 10% ethanol, 30% polyethylene glycol 400, and 60% PHOSAL 50 PG, at a final concentration of 10 mg/mL. Azacitidine (Sigma Aldrich, Cat# 320-67-2) stock solution was prepared at 1 mg/mL by dissolving the drug in 0.9% NaCl. The drugs were prepared fresh every week and kept at four degrees for storage. Ivosidenib and enasidenib were administered orally twice a day at a dose of 300 and 45 mg/kg, respectively. Venetoclax was administered orally once a day at a dose of 100 mg/kg. Azacitidine was administered by subcutaneous injection twice a week at a dose of 5 mg/kg (on the first and third day of the 5-days-cycle).

**Bone marrow aspiration**

The mice were injected with Meloxicam (5 mg/kg) and anesthetized with isofluorane before the surgery. One draw of bone marrow aspirate from the (left or right) femur was collected through the kneecap with a 28G insulin syringe (fisher scientific, Cat# 14-826-79). The bone marrow aspirates were kept in 100 μL of heparin to prevent clotting. The red cells in the aspirates were lysed with 1 mL of ACK lysis buffer (BD Biosciences, Cat# 555899) for 5 minutes at room temperature. The samples were passed through a 40-micron cell strainer and washed with 3 mL of PBS. Finally, the harvested cells were re-suspended in 1 mL of CryoStor® CS10 (STEMCELL Technologies, Cat# 07930) for cryo-preservation.

**Endpoint sample collection**

The mice were euthanized using carbon dioxide after eight weeks of treatment. The tibia, femur, and iliac crest of both sides were crushed in a motor and pestle, then flushed with excess PBS to collect as many cells as possible. The red cells were lysed with 1 mL of ACK lysis buffer (BD Biosciences, Cat# 555899) for 5 minutes at room temperature. The samples were passed through a 40-micron cell strainer and washed with 3 mL of PBS. Finally, the harvested cells were re-suspended in 1mL of CryoStor® CS10 (STEMCELL Technologies, Cat# 07930) for cryo-preservation.

**Flow cytometry analysis**

From the 1 mL of cryo-preserved cells, 5 μL was taken and added to 95 μL of antibody master mix for a final staining volume of 100 μL. The antibody master mix comprised the following components: (1) PE-conjugated anti human CD14 antibody (RRID:AB_830678) (2) PE/Cy7-conjugated anti human CD15 antibody (RRID:AB_2561669) (3) APC-conjugated anti human CD11b antibody (RRID:AB_1210558) (4) APC/Fire750-conjugated anti mouse CD45.1 antibody (RRID:AB_2629805) (5) BV421-conjugated anti human CD45 antibody (RRID:AB_10900423) (6) BV711-conjugated anti human CD33 antibody (RRID:AB_2565774 ) (7) FITC-conjugated anti mouse TER-119 antibody (RRID: 
AB_313706) (8) Helix NP Green viability dye (Biolegend, Cat# 425303). All antibodies were diluted to a final concentration of 1:100, while Helix NP green (BioLegend, Cat# 425303) was diluted to a final concentration of 1:10,000 in FACS buffer (PBS with 10% FBS and 1% sodium azide). The cells were stained for 30 minutes at 4 degrees Celsius, and the data was acquired using a Cytoflex flow cytometer (RRID:SCR_019627). The leukemic burden is defined as the percentage of mCD45.1-TER-119-hCD45+hCD33+ cells, gated on all live cells. The differentiation response is assessed by the expression of myeloid differentiation markers CD11b, CD14, and CD15, as well as granularity (side-scatter), gated on live mCD45.1-hCD45+ cells.

**Variant calling using targeted deep sequencing**

Bone marrow cells from untreated animals were harvested for targeted deep sequencing using IDT’s Target-Seq SMART-AML panel. Adapter trimmed reads were aligned to a mm10 + hg19 hybrid reference genome using bowtie2 (version 2.4.5) (RRID:SCR_016368). Only reads that mapped to the hg19 portion of the hybrid genome was kept. The reads were sorted by chromosomal positions using samtools (version 1.14) (RRID:SCR_002105) and duplicated reads originating from PCR artifacts were marked using Picard (version 2.10.9) (RRID:SCR_006525). The base quality scores were re-calibrated using BaseRecalibrator (GATK, version 4.2.5.0). The variants were called using Mutect2 tumor only mode (GATK, version 4.2.5.0) (RRID:SCR_001876).

**Single-cell proteogenomic sequencing**

Except for the primary sample, all samples were thawed, pooled, and sorted for live human cells at the Princess Margaret Flow Cytometry Facility in the morning and immediately transferred to the Princess Margaret Genomic Center for library preparation. In the first case study, we spiked-in cells carrying a homozygous TP53p.I255S mutation (covered by Tapestri amplicon) in the AZA and IVO+AZA arms to reach sufficient cells (20,000) for scDNA sequencing (started with 10,000 and 2,500 cells after pooling, respectively). In the second case study, we spiked-in SKM1 cells carrying a homozygous TP53p.R248Q mutation (covered by the same Tapestri amplicon as TP53p.I255S) in vehicle and ivosidenib treated samples to offset excessive cell loss during the wash steps in the protein sequencing workflow. The variants were called using MissionBio’s Tapestri pipeline v2. Protein data was normalized by taking the centered log-ratio across each cell.

**Antibody sequencing analysis and hypothesis testing**

For pseduobulk analysis of differentially expressed surface proteins, p-values were derived from two-sided student’s t-test of normalized protein counts. P-values were corrected using the Benjamini–Hochberg method to control for multiple comparisons.

**Quality control of single-cell genotype data**

First, genotype calls with a Phred quality score <30 or a sequencing depth <10 are removed. We then remove any cell where not all variants of interest are genotyped. For the AZA and IVO+AZA treated samples from case study 1, an additional filter was used to remove spike-in cells carrying the TP53p.I255S mutation. For the vehicle and ivosidenib treated samples from case study 2, an additional filter was used to remove spike-in cells carrying the TP53p.R248Q mutation. In mixed PDX models, cells carrying combination of mutations that are known to be mutually exclusive to patients are marked as doublets and removed. 

**Attaching cells to phylogenetic tree**

Each cell that passed quality control is attached to one clone of the corresponding phylogenetic tree using its most descendant mutation. Under the infinite sites assumption, we assume that any wildtype-call from ancestral variants arises due to dropout of the mutant allele. Using case study 1 as an example, a cell that carries WT1p.R439C but is wildtype for NPM1 and/or FLT3-ITD (and wildtype at loci TET2p.Q1548del¬, DNMT3Ap.D531del, KITp.H40Qfs) will be attached to Clone 6.

**Bayesian modelling of human cell number in xenotransplanted animals**

We model the number of human cells N_ij collected at timepoint i under treatment j as following a gamma-poisson distribution with mean μ_ij and dispersion ϕ_i that accounts for experimental variability of each timepoint i. We use a logarithmic link log(μ_ij )=μ_0j+bT[j], where μ_0j is the grand mean and bT[j] is the coefficient for treatment j. The treatments effects bT is represented using a multilevel model with prior bT~Normal(bT_μ,bT_σ ) and hyperpriors bT_μ~Normal(0,7) and bT_σ~Normal(0,2). We assume a set of weakly informative priors μ_0j~Uniform(0,17.7) and ϕ_i~Uniform(0,1) that cover a wide range of possible outcomes within expectation.

**Bayesian modelling of clonal proportion across treatments**

We fit a multilevel multinomial regression model to estimate the proportion P_ij  of each clone i under treatment j. First, we model the number of sequenced cells C_ij  from clone i under treatment j as following a Poisson distribution with mean λ_ij. We use a logarithmic link log(λ_ij )=bT[i,j], where bT models the varying effects of treatments with prior bT~Normal(bT_μ,bT_σ ), and hyperpriors bT_μ~Normal(0,10) and bT_σ~Gamma(3,1). We then derive P_ij=C_ij/(∑_(i=0)^k▒C_ij ), where k is the total number of mutant clones in the phylogenetic tree.
The Posterior samples of P_ij were used to calculate the Shannon diversity index. The posteriors of clonal proportion and leukemic cell number for each treatment were used to derive the log2-fold change in absolute cell number of each clone compared to the vehicle control. A pseudocount of 1 cell was added to avoid division by zero.

**Model checking**

Bayesian models were fitted using Stan’s implementation of Hamiltonian Monte Carlo (version 2.18) (RRID:SCR_018459) through R (cmdstanr, version 0.6.1). We ran eight separate chains for all models. We inspected the traceplots of each chain to ensure proper exploration of the posterior. We made sure that all chains converged by checking that the Rhat for every parameter is equal to 1. 
