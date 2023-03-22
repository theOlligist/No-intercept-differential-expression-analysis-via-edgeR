#Differential expression during 2018
rm(list = ls())
source("~/Dropbox/R_general_working_directory/Dictionary_of_Colors.R")

library(tidyverse); library(edgeR); library(reshape2); library(WGCNA); library(patchwork);library(jishonoiro);library(vegan);library(corrplot)
setwd("~/Dropbox/SMP_R_working_directory/MetaTranscriptomics/2018_assembly/")

one41= c("#FF5200", "#0D2B52","#A6FF47")
pal1 = c("#8c6510", "#0D2B52","#FF7399")
color_pal = c("black", "light grey", "steel blue")

# Load and process the raw count table ------------------------------------------
#filter to the dominant diatom taxa from 2018 with enough reads to meaningfully do so (Figure S#-diatoms). I left off Rhizosolenia because there just weren't enough reads.
Diatom_raw_df = read.delim("raw-count-data-metaT-SMP2018.txt") %>% 
  filter(str_detect(Taxonomy, "Thalassiosirales") | str_detect(Taxonomy, "Bacillariales")) %>% 
  mutate(tax = case_when(str_detect(Taxonomy, "Thalassiosirales") ~"Thalassiosira", 
                         str_detect(Taxonomy, "Bacillariales") ~"Pseudo_nitz")) %>% 
  select(tax, KO, everything()) %>% 
  unite("ID", c(tax,KO), sep = "-")

Diatom_raw_df %>% 
  filter(str_detect(ID, "Thalas")) %>% nrow()
nrow(Diatom_raw_df) # there are approx 26,000 transcripts
#There's only a little over 600 genes for rizosolenia but I'll roll with it.

#And count the number of KO terms. I plan to use the KO IDs as the rownames;
#therefore I'll will aggregate them for each sample. Flattening the tax diversity to FX diversity.
#KO_total = Diatom_raw_df %>% distinct(KO) %>% nrow() #count the number of distinct KO
# I can expect a max o ~5622

#because I have filtered down this dataset to be diatoms, The subject of this dataframe is the functions themselves, not the taxa.
#However, there is duplication in KO terms because this is still a mixed population of diatoms expressing similar genes in some cases.
# I will aggregate the counts from each KO in each sample. This is an important step in order to assign the KO's the rownames.

#Aggregate Table
agg_df = Diatom_raw_df %>%
  pivot_longer(cols = starts_with("SMP"), names_to = "sample", values_to = "value") %>% 
  #dcast(ID~sample) %>% column_to_rownames("ID") %>% colSums()
  mutate(sample = str_replace(sample, "SMPier\\.", ""),
         sample = str_replace(sample, regex("\\.\\d{2}\\.\\d{2}\\.2018_S\\d+"), ""),
         sample = factor(sample, levels = c("Day1.A", "Day1.B", "Day1.C", 
                                            "Day3.A", "Day3.B", "Day3.C", 
                                            "Day5.A", "Day5.B", "Day5.C", 
                                            "Day9.A", "Day9.B", "Day9.C", 
                                            "Day11.A", "Day11.B", "Day11.C"))) %>% 
  group_by(ID, sample) %>%
  summarise(count = sum(value, na.rm = TRUE)) %>%
  dcast(ID~sample)

#Assign rownames as KO
rownames(agg_df) = agg_df$ID

#Take a peek at the first three rows and first six columns
agg_df[1:3, 1:16]
#sanity check. These should be the same
colSums(agg_df[-1]) == colSums(Diatom_raw_df[-1:-2])

# # Normalize the agg_df dataframe with edgeR -----------------------------
#(1) create a dge object and (2) calculate normalizing factors and then (3) return the cpm as a dataframe
#Ingredients for the dge object: dataframe (counts only!), sample_list, & metadata.

#Make sure these line up with the order of the columns of the agg_df
# Create a sample list 
sample_list = c("Day1","Day3","Day5","Day9","Day11")
# Create the metadata
metadata = factor(c(rep("day1",3),
                    rep("day3",3),
                    rep("day5",3),
                    rep("day9",3),
                    rep("day11",3)),
                  levels = str_to_lower(sample_list)) # levels is important for purposefully setting the sample order.


# Use this information to create the DGEList object
dge_obj = DGEList(counts = as.data.frame.matrix(agg_df[-1]), # Counts is the actual columns of the dataframe
                  genes = agg_df[1], # This is the name of the column(similar to genes)
                  group = metadata)

#calculate the normalizing factors using TMM
dge_obj = calcNormFactors(dge_obj, method = "TMM")


# Differential Expression analysis (edgeR) -------------------------------------------------------------
# the goal here is to produce a vector of KOID that are differentially expressed over the course of the samples.
# Filter out low count rows using a conditional filter--that is, to keep rows that meet a certian condition.
keep = filterByExpr(dge_obj)
dge_obj = dge_obj[keep, ,keep.lib.sizes=FALSE]
# setup the design matrix without an intercept as day 1!
design = model.matrix(~0+group, data = dge_obj$samples)

#Set the column names for the design to match the
colnames(design) = levels(dge_obj$samples$group)

# using makeContrasts() to specify the pariwise comparisons
conts = makeContrasts(day3-day1, day5-day1, day9-day1, day11-day1, levels = str_to_lower(sample_list)) #for contrasting all samples to prebloom state

# estimate common and tagwise(pairwise) dispersion accrding the study design matrix
disp = estimateGLMCommonDisp(dge_obj, design = design)
disp = estimateGLMTagwiseDisp(disp, design = design)

# Determine differentially expressed genes
fit = glmQLFit(disp, design = design)
DEGs = glmQLFTest(fit, contrast = conts)
DEG_table = DEGs$table %>% data.frame()

# adjust p values due to so many comparissons to decrease chance of type 1 error.
DEG_table$P.adjusted = p.adjust(DEG_table$PValue, method = "fdr")
DEG_table$KO = rownames(DEG_table)
DEG_table_Pfilt = DEG_table %>% 
  filter(P.adjusted < 0.001) %>% 
  separate(KO, into = c("tax","KO"), sep = "-")

# End