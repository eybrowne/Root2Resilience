
#a sequence table with chimeras removed, taxonomy table, OTU table, and metadata sheet is needed for this code
#This code uses DESeq differential analysis with a heatmap and complex upset plot output. 




# Clean-up the memory and start a new session
#############################################################
rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################

library("plyr")
library("grid")
library("phyloseq") 
library("DESeq2")
library("ggplot2")
library("vegan")
library ("ape")
library("dplyr")
library("tidyverse")
library("forcats")
library("PMCMRplus")
library ("UpSetR")
library("agricolae")
library("tibble")
library("UpSetR")
library("ComplexUpset")
library("tidyr")
library("magrittr")
library("writexl")
library("agricolae")
library("devtools")
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis") 
library("pairwiseAdonis")

seqtab.nochim <- readRDS("seqtab.nochim.rds")
taxid <-readRDS("taxid.rds")
ids <-readRDS("ids.rds")
#attach metadata
metadata <- read.csv("metadata.csv", header=T, row.names = 1)
metadata2 <- read.csv("metadata2.csv", header=T, row.names = 1) #metadata file with outliers removed

#create phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(metadata), 
               tax_table(taxid))

#change colnames of ASV table into ASV1 to ASVXXXXX and store the DNA sequences of our ASVs in the refseq slot of the phyloseq object
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps


#Convert the metadata to sample_data format
sample_metadata <- sample_data(metadata)

#Add sample data to  phyloseq object
ps <- merge_phyloseq(ps, sample_metadata)

#Check if the metadata has been added
sample_data(ps)

#if taxa are in columns
otu_data <- otu_table(ps)

# Transpose the OTU table if taxa are in columns
if (!taxa_are_rows(ps)) {
  otu_data <- t(otu_data)  # Transpose the matrix
}

# Assign the transposed OTU table back to the phyloseq object
otu_table(ps) <- otu_table(otu_data, taxa_are_rows = TRUE)

# Verify if the change has been applied
taxa_are_rows(ps)  # Should now return TRUE
head(otu_table(ps))

sample_names(ps)
rank_names(ps)
ITS_tax_table<-tax_table(ps)

#check dimenstions of phyloseq object
dim(tax_table(ps))

#filter out zero read samples 
otu_filtered <- otu_table(ps)[rowSums(otu_table(ps)) > 0, ]
physeq_filtered <- phyloseq(otu_table(otu_filtered, taxa_are_rows(ps)), sample_data(ps))

##ratio filtered reads/total reads
ratio <- sum(sample_sums(physeq_filtered))/sum(sample_sums(ps))*100
ratio

####################################################################################
# Rarefy at an even sequencing depth (200,000 for our NextSeq data) and "freeze" these objects for downstream analyses
#####################################################################################

#ASVs : ignore the warnings, the object will be saved right after
physeq_200K <- rarefy_even_depth(ps, 200000)

saveRDS(physeq_200K, file ="R2R_ITS_rare_ASV_200K.rds")

physeq_200K <- readRDS("R2R_ITS_rare_ASV_200K.rds")

#check if taxa are rows
taxa_are_rows(physeq_200K)
#if not, convert to rows
otu_table(physeq_200K) <- otu_table(physeq_200K, taxa_are_rows = TRUE)
# Step 2: Convert the metadata to sample_data format
sample_metadata1 <- sample_data(metadata1)

# Step 3: Add sample data to your phyloseq object
physeq <- merge_phyloseq(ps, sample_metadata)
physeq_200K <- merge_phyloseq(physeq_200K, sample_metadata)
sample_names(physeq_200K)

#prune samples with no reads and those which were identified as outliers
samples_to_remove <- c("NX_8_JC_D_UDI418_S407", "NX_8_JC_L_S397")

# Prune samples by keeping the ones that are NOT in the samples_to_remove
physeq_200K_pruned <- prune_samples(!(sample_names(physeq_200K) %in% samples_to_remove), physeq_200K)

#save the pruned phyloseq object
saveRDS(physeq_200K_pruned, file ="R2R_ITS_200K_pruned.rds")

physeq_200K_pruned <- readRDS("R2R_ITS_200K_pruned.rds")

#check remaining samples 
nsamples(physeq_200K_pruned)
sample_names(physeq_200K_pruned)

metadata2 <- read.csv("metadata2.csv", header=T, row.names = 1) #metadata file with pruned samples removed


###########################################################################################
##B-diversity across protocols
#Constrained ordination: constrained for Description2
#Rarefied data

JH33_CAP <- ordinate(physeq_200K_pruned, "CAP", "bray", ~ Description2)
p <-plot_ordination(physeq_200K_pruned, JH33_CAP, color = "Method")
p1 <-plot_ordination(physeq_200K_pruned, JH33_CAP, color = "Description2")
p
p1

#ANOVA on the axis
aov <- anova(JH33_CAP , permutations=5000)
aov

#assign shapes to Soil type and color to Sample type
p = plot_ordination(physeq_200K_pruned, JH33_CAP , color = "Description2", shape = "Description")
p = p + geom_point (size = 5, alpha = 0.80)
p = p + scale_colour_manual(values = c("#0072B2","#56B4E9","#E69F00","#D55E00","#009E73", "#0072B2","#D55E00"))
p + ggtitle("R2R Pilot Trial ITS data, Bray distance-rarefied samples")

design <- metadata2
row.names(design)
colnames(design)

#Protocol effect (rhizosphere only)
physeq_rhizo <- subset_samples(physeq_200K_pruned, Description == "Barke")
design_rhizo <- design[rownames(otu_table(physeq_rhizo)), ]
BC <- phyloseq::distance(physeq_rhizo, method = "bray")

#run adonis 
set.seed(123)  # For reproducibility
ad <- adonis2(BC ~ Description2, data = design_rhizo, permutations = 5000)
ad

#Pairwise MANOVA
pair.mod<-pairwise.adonis(BC,factors=design_rhizo$Description2)
pair.mod


#assign shapes to Soil type and color to Sample type

p = plot_ordination(JH33_rhizo, JH33_CAP_rhizo , color = "Description2", shape = "Description", label = "Sample")
p = p + geom_point (size = 5, alpha = 0.80)
p = p + scale_colour_manual(values = c("#0072B2","#56B4E9","#E69F00","#D55E00","#009E73", "#0072B2","#D55E00"))
p + ggtitle("JH33 CAP 16S data, Bray distance-rarefied samples")


#assign shapes to Soil type and color to Sample type
p2 =plot_ordination(physeq_200K_pruned, JH33_CAP , color = "Description2", shape = "Description", label = "Sample")
p2 = p2 + geom_point (size = 5, alpha = 0.80)
p2 = p2 + scale_colour_manual(values = c("#0072B2","#56B4E9","#E69F00","#D55E00", "#0072B2","#D55E00"))
p2 + ggtitle("JH33 CAP 16S data, Bray distance-rarefied samples")


#Adonis on axis 
design_all <- design[colnames(otu_table(physeq_200K_pruned)), ]
BC_all <- phyloseq::distance(physeq_200K_pruned, "bray")
ad <-adonis2(BC_all ~ Description2, data= design_all , permutations = 5000)
ad


##Alpha diversity ########################
###############################################################################

g = plot_richness(physeq_200K_pruned , x = "Description2", color = "Method", measures=c("Observed", "Shannon"))
g = g + scale_colour_manual(values = c("#0072B2","#56B4E9","#E69F00","#D55E00","#009E73", "#0072B2","#D55E00"))
g

##################################################################################
Protocol_alpha <-  estimate_richness(physeq_200K_pruned, measures = c("Observed", "Shannon", "Chao1")) 
Protocol_alpha 

Protocol_otu_table <-otu_table(physeq_200K_pruned)


#Data frame Genotype_Description 
colData = metadata[rownames(Protocol_otu_table),]
rownames(colData)
colData

#Description 
design_Protocol  <- as.data.frame(colData[, 4]) 
rownames(design_Protocol) <- rownames(colData) 
colnames(design_Protocol) <- c("Description2") 
design_Protocol 


###############################
#### Alpha diversity OBSERVED 
##############################
#Observed ASVs 
Protocol_alpha_Observed <- as.data.frame(Protocol_alpha[ ,1]) 
rownames(Protocol_alpha_Observed) <- rownames(Protocol_alpha) 
colnames(Protocol_alpha_Observed) <- c("Observed") 

#Combine the dataset sample description and Observed OTUs 
Protocol_alpha_Observed_TD <- cbind(design_Protocol , Protocol_alpha_Observed) 
Protocol_alpha_Observed_TD <- as.data.frame(Protocol_alpha_Observed_TD)
Protocol_alpha_Observed_TD$Description2
#drop DB Bulk due to no reads


Protocol_alpha_Observed_TD$Description2 <- droplevels(Protocol_alpha_Observed_TD$Description2)
levels(Protocol_alpha_Observed_TD$Description2)

Protocol_alpha_Observed_TD$Description2 <- ordered(Protocol_alpha_Observed_TD$Description2, levels=c("Bulk_EO", "Barke_WET","Barke_DRY+WET", "Barke_ALL","Barke_DB"))

#Plotting
with(Protocol_alpha_Observed_TD, boxplot(Observed~Description2, xlab = "Description2", ylab = "Number of taxa", ylim=c(0,600),main = "Taxa observed richness", col=c("#D55E00","#D55E00","#E69F00","#56B4E9","#0072B2")))
#jitter points included 
with(Protocol_alpha_Observed_TD, stripchart(Observed~Description2, xlab = "Description2", ylab = "Number of taxa", vertical=T, add=T, method ="jitter", col=c('black'), pch=16))

#ANOVA 
Observed_OTUs_stats <- aov(Observed ~ Description2, data = Protocol_alpha_Observed_TD) 
summary(Observed_OTUs_stats)
#Tukey HSD
HSD.test(Observed_OTUs_stats, trt = c("Description2"), console = TRUE)

###########################
#### Alpha diversity CHAO1 
##########################
#Chao1 ASVs 
Protocol_alpha_Chao1 <- as.data.frame(Protocol_alpha[ ,2]) 
rownames(Protocol_alpha_Chao1) <- rownames(Protocol_alpha) 
colnames(Protocol_alpha_Chao1) <- c("Chao1") 
Protocol_alpha_Chao1_TD <- cbind(design_Protocol , Protocol_alpha_Chao1) 
Protocol_alpha_Chao1_TD <- as.data.frame(Protocol_alpha_Chao1_TD) 
Protocol_alpha_Chao1_TD$Description2
Protocol_alpha_Chao1_TD$Description2 <- droplevels(Protocol_alpha_Chao1_TD$Description2)
levels(Protocol_alpha_Chao1_TD$Description2)
#levels in defined order
Protocol_alpha_Chao1_TD$Description2 <- ordered(Protocol_alpha_Chao1_TD$Description2, levels=c("Bulk_EO", "Barke_WET","Barke_DRY+WET", "Barke_ALL","Barke_DB"))
#Plotting
with(Protocol_alpha_Chao1_TD, boxplot(Chao1~Description2, xlab = "Protocol", ylab = "Number of ASVs", ylim=c(0,600), main = "ASVs Chao richness", col=c("#D55E00","#D55E00","#E69F00","#56B4E9","#0072B2")))
with(Protocol_alpha_Chao1_TD, stripchart(Chao1~Description2, xlab = "Protocol", ylab = "Number of ASVs", vertical=T, add=T, method ="jitter", col=c('black'), pch=16))

#ANOVA 
Chao1_OTUs_stats <- aov(Chao1 ~ Description2, data = Protocol_alpha_Chao1_TD) 
summary(Chao1_OTUs_stats) 
HSD.test(Chao1_OTUs_stats, trt = c("Description2"), console = TRUE)

#############################
#### Alpha diversity SHANNON ###Assessing fungal abundances is tricky due to differences in ITS copies between taxa
#############################
#Shannon ASvs 
Protocol_alpha_Shannon <- as.data.frame(Protocol_alpha[ ,4]) 
rownames(Protocol_alpha_Shannon) <- rownames(Protocol_alpha) 
colnames(Protocol_alpha_Shannon) <- c("Shannon") 


#Combine the dataset sample description and Shannon OTUs 
Protocol_alpha_Shannon_TD <- cbind(design_Protocol , Protocol_alpha_Shannon) 
Protocol_alpha_Shannon_TD <- as.data.frame(Protocol_alpha_Shannon_TD) 
Protocol_alpha_Shannon_TD$Description2
Protocol_alpha_Shannon_TD$Description2 <- droplevels(Protocol_alpha_Shannon_TD$Description2)
levels(Protocol_alpha_Chao1_TD$Description2)

#Order the levels according to a defined order 
Protocol_alpha_Shannon_TD$Description2 <- ordered(Protocol_alpha_Shannon_TD$Description2, levels=c("Bulk_EO", "Barke_WET","Barke_DRY+WET", "Barke_ALL","Barke_DB"))

#Plotting
with(Protocol_alpha_Shannon_TD, boxplot(Shannon~Description2, xlab = "Protocol", ylab = "Number of ASVs", ylim=c(0,6), main = "ASVs Shannon richness",col=c("#D55E00","#D55E00","#E69F00","#56B4E9","#0072B2")))
with(Protocol_alpha_Shannon_TD, stripchart(Shannon~Description2, xlab = "Protocol", ylab = "Number of ASVs", vertical=T, add=T, method ="jitter", col=c('black'), pch=16))


#ANOVA 
Shannon_OTUs_stats <- aov(Shannon ~ Description2, data = Protocol_alpha_Shannon_TD) 
summary(Shannon_OTUs_stats) 
HSD.test(Shannon_OTUs_stats, trt = c("Description2"), console = TRUE)

#########################################################################################

##Code to compare the rhizosphere samples (Barke/soil batch)

Allsoils_phyloseq_rhizo <- subset_samples(physeq_200K_pruned, Description =="BARKE")
Allsoils_phyloseq_rhizo 

##############################################################
## DESeq
###############################################################

#formatting 
abund_table <- otu_table(physeq_200K_pruned)

meta_table <- metadata2[rownames(abund_table),]

abund_table<-abund_table[rownames(abund_table) %in% rownames(meta_table),]

abund_table<-abund_table[,colSums(abund_table)>0]

#We will convert our table to DESeqDataSet object
countData = round(as(abund_table, "matrix"), digits = 0)
#we need to flip it so the sample names are rows 
countData1 = t(countData)

#run DeSeq
dds <- DESeqDataSetFromMatrix(countData1, meta_table, as.formula(~ Description2))
data_deseq_test = DESeq(dds)


#extract the results
res = results(data_deseq_test, cooksCutoff = FALSE)
res_tax = cbind(as.data.frame(res), as.matrix(countData1[rownames(res), ]), OTU = rownames(res))
plot.point.size = 2
label=F
tax.display = NULL
tax.aggregate = "OTU"

sig = 0.05
fold = 2
res_tax_sig = subset(res_tax, padj < sig & fold < abs(log2FoldChange))

res_tax_sig <- res_tax_sig[order(res_tax_sig$padj),]

levels(data_deseq_test$Description2)

saveRDS(data_deseq_test, file="PilotTrial_cds_test_ASVs.rds")
JH33_cds_test <- data_deseq_test

###Preparing data for Upset visualisation ##########

############ EO_BUlk as control of all methods 25/06/2025 ################
####### EO_Bulk vs DB_ALL #############
rhizo_EO_DB_ALL <- results(JH33_cds_test, contrast = c("Description2","Barke_ALL","Bulk_EO")) 

#inspect a result file
rhizo_EO_DB_ALL 

# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
rhizo_EO_DB_ALL_FDR005 <- rhizo_EO_DB_ALL[(rownames(rhizo_EO_DB_ALL)[which(rhizo_EO_DB_ALL$padj <0.05)]), ]

#enriched in rhizosphere (negative fold change, HID is at the second term of comparison)
rhizo_EO_DB_ALL_enriched <-  rhizo_EO_DB_ALL[(rownames(rhizo_EO_DB_ALL)[which(rhizo_EO_DB_ALL$log2FoldChange > 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
rhizo_EO_DB_ALL_enriched_FDR005 <- intersect(rownames(rhizo_EO_DB_ALL_FDR005), rownames(rhizo_EO_DB_ALL_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(rhizo_EO_DB_ALL_enriched_FDR005)
#Enriched in rhizo_EO_DB_ALL 


############ EO_BUlk vs DB ################
rhizo_EO_DB <- results(JH33_cds_test, contrast = c("Description2","Barke_DB","Bulk_EO")) 

#inspect a result file
rhizo_EO_DB 

# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
rhizo_EO_DB_FDR005 <- rhizo_EO_DB[(rownames(rhizo_EO_DB)[which(rhizo_EO_DB$padj <0.05)]), ]

#enriched in rhizosphere (negative fold change, HID is at the second term of comparison)
rhizo_EO_DB_enriched <-  rhizo_EO_DB[(rownames(rhizo_EO_DB)[which(rhizo_EO_DB$log2FoldChange > 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
rhizo_EO_DB_enriched_FDR005 <- intersect(rownames(rhizo_EO_DB_FDR005), rownames(rhizo_EO_DB_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(rhizo_EO_DB_enriched_FDR005)
#Enriched in rhizo_EO_DB_ALL

############ EO_BUlk vs EO_WET ################
rhizo_EO_wet <- results(JH33_cds_test, contrast = c("Description2","Barke_WET","Bulk_EO")) 

#inspect a result file
rhizo_EO_wet 

# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
rhizo_EO_wet_FDR005 <- rhizo_EO_wet[(rownames(rhizo_EO_wet)[which(rhizo_EO_wet$padj <0.05)]), ]

#enriched in rhizosphere (negative fold change, HID is at the second term of comparison)
rhizo_EO_wet_enriched <-  rhizo_EO_wet[(rownames(rhizo_EO_wet)[which(rhizo_EO_wet$log2FoldChange > 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
rhizo_EO_wet_enriched_FDR005 <- intersect(rownames(rhizo_EO_wet_FDR005), rownames(rhizo_EO_wet_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(rhizo_EO_wet_enriched_FDR005)

#Enriched in EO_wet      


############ EO_BUlk vs EO_DRY+WET ################
rhizo_EO_drywet <- results(JH33_cds_test, contrast = c("Description2","Barke_DRY+WET","Bulk_EO")) 

#inspect a result file
rhizo_EO_drywet 

# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
rhizo_EO_drywet_FDR005 <- rhizo_EO_drywet[(rownames(rhizo_EO_drywet)[which(rhizo_EO_drywet$padj <0.05)]), ]

#enriched in rhizosphere (negative fold change, HID is at the second term of comparison)
rhizo_EO_drywet_enriched <-  rhizo_EO_drywet[(rownames(rhizo_EO_drywet)[which(rhizo_EO_drywet$log2FoldChange > 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
rhizo_EO_drywet_enriched_FDR005 <- intersect(rownames(rhizo_EO_drywet_FDR005), rownames(rhizo_EO_drywet_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(rhizo_EO_drywet_enriched_FDR005)
#Enriched in EO_drywet      


#####################################################################################################
##Comparison between protocols
###########################################################################################################
DB_vs_DBall <- results(JH33_cds_test, contrast = c("Description2","Barke_DB","Barke_ALL")) 

#inspect a result file
DB_vs_DBall 
mcols(DB_vs_DBall  , use.names=TRUE)
# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
DB_vs_DBall_FDR005 <- DB_vs_DBall[(rownames(DB_vs_DBall)[which(DB_vs_DBall$padj <0.05)]), ]

#enriched in rhizosphere (negative fold change, HID is at the second term of comparison)
DB_vs_DBall_enriched <-  DB_vs_DBall[(rownames(DB_vs_DBall)[which(DB_vs_DBall$log2FoldChange > 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
DB_vs_DBall_enriched_FDR005 <- intersect(rownames(DB_vs_DBall_FDR005), rownames(DB_vs_DBall_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(DB_vs_DBall_enriched_FDR005)
#Enriched in DB rhizo 


DB_vs_DBall_enriched_FDR005_rhizo <- intersect(DB_vs_DBall_enriched_FDR005, rhizo_EO_DB_enriched_FDR005)
length(DB_vs_DBall_enriched_FDR005_rhizo)


#remove bulk soil whose abundances are markedly different from plant-associated ones
otu_table_JH33 <- otu_table(JH33_plants)

# Transpose the OTU table
otu_table_transposed <- t(otu_table_JH33)

# Create a new OTU table object from the transposed data
new_otu_table <- otu_table(otu_table_transposed, taxa_are_rows = TRUE)

#  Create a new phyloseq object
JH33_plants <- phyloseq(new_otu_table, sample_data(JH33_plants), tax_table(JH33_plants))

# Verify the new dimensions
print(dim(otu_table(JH33_plants)))


JH33_plants_Bulk <- subset_samples(JH33_plants, Description !="Bulk")
JH33_plants_otu_table <-otu_table(JH33_plants_Bulk)

taxa_names(JH33_plants_Bulk)


###############################################################################
EO_vs_EOdrywet <- results(JH33_cds_test, contrast = c("Description2","Barke_WET","Barke_DRY+WET")) 

#inspect a result file
EO_vs_EOdrywet 
mcols(EO_vs_EOdrywet  , use.names=TRUE)
# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
EO_vs_EOdrywet_FDR005 <- EO_vs_EOdrywet[(rownames(EO_vs_EOdrywet)[which(EO_vs_EOdrywet$padj <0.05)]), ]

#enriched in rhizosphere (negative fold change, HID is at the second term of comparison)
EO_vs_EOdrywet_enriched <-  EO_vs_EOdrywet[(rownames(EO_vs_EOdrywet)[which(EO_vs_EOdrywet$log2FoldChange > 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
EO_vs_EOdrywet_enriched_FDR005 <- intersect(rownames(EO_vs_EOdrywet_FDR005), rownames(EO_vs_EOdrywet_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(EO_vs_EOdrywet_enriched_FDR005)

EO_vs_EOdrywet_enriched_FDR005_rhizo <- intersect(EO_vs_EOdrywet_enriched_FDR005,rhizo_EO_wet_enriched_FDR005)
length(EO_vs_EOdrywet_enriched_FDR005_rhizo)


#Subset for taxa DE in DBs and EOs protocols

JH33_EO_vs_EOdrywet <- prune_taxa(EO_vs_EOdrywet_enriched_FDR005_rhizo, JH33_plants_Bulk)
JH33_EO_vs_EOdrywet

All_EO_vs_EOdrywet <- median(sample_sums(JH33_EO_vs_EOdrywet))/median(sample_sums(JH33_plants_Bulk)) * 100
All_EO_vs_EOdrywet


#####################################################################################
DB_vs_ALL <- results(JH33_cds_test, contrast = c("Description2","Barke_DB","Barke_ALL")) 

#inspect a result file
DB_vs_ALL 
mcols(DB_vs_EOwet  , use.names=TRUE)
# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
DB_vs_ALL_FDR005 <- DB_vs_ALL[(rownames(DB_vs_ALL)[which(DB_vs_ALL$padj <0.05)]), ]

#enriched in rhizosphere (negative fold change, HID is at the second term of comparison)
DB_vs_ALL_enriched <-  DB_vs_ALL[(rownames(DB_vs_ALL)[which(DB_vs_ALL$log2FoldChange < 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
DB_vs_ALL_enriched_FDR005 <- intersect(rownames(DB_vs_ALL_FDR005), rownames(DB_vs_ALL_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(DB_vs_ALL_enriched_FDR005)

# db vs DB_ALL  #  DB_ALL vs db

DB_vs_ALL_enriched_FDR005_rhizo <- intersect(DB_vs_ALL_enriched_FDR005,rhizo_EO_DB_ALL_enriched_FDR005)
length(DB_vs_ALL_enriched_FDR005_rhizo)


JH33_DB_vs_ALL <- prune_taxa(DB_vs_ALL_enriched_FDR005_rhizo, JH33_plants)
JH33_DB_vs_ALL

All_DB_vs_ALL <- median(sample_sums(JH33_DB_vs_ALL))/median(sample_sums(JH33_plants)) * 100
All_DB_vs_ALL

###################################################################################
DB_vs_EOwet <- results(JH33_cds_test, contrast = c("Description2","Barke_DB","Barke_WET")) 

#inspect a result file
DB_vs_EOwet 
mcols(DB_vs_EOwet  , use.names=TRUE)
# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
DB_vs_EOwet_FDR005 <- DB_vs_EOwet[(rownames(DB_vs_EOwet)[which(DB_vs_EOwet$padj <0.05)]), ]

#enriched in rhizosphere (negative fold change, HID is at the second term of comparison)
DB_vs_EOwet_enriched <-  DB_vs_EOwet[(rownames(DB_vs_EOwet)[which(DB_vs_EOwet$log2FoldChange < 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
DB_vs_EOwet_enriched_FDR005 <- intersect(rownames(DB_vs_EOwet_FDR005), rownames(DB_vs_EOwet_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(DB_vs_EOwet_enriched_FDR005)


DB_vs_EOwet_enriched_FDR005_rhizo <- intersect(DB_vs_EOwet_enriched_FDR005,rhizo_EO_wet_enriched_FDR005)
length(DB_vs_EOwet_enriched_FDR005_rhizo)




JH33_DB_vs_EOwet <- prune_taxa(DB_vs_EOwet_enriched_FDR005_rhizo, JH33_plants)
JH33_DB_vs_EOwet

All_DB_vs_EOwet <- median(sample_sums(JH33_DB_vs_EOwet))/median(sample_sums(JH33_plants)) * 100
All_DB_vs_EOwet

############################################################################################################

DB_vs_EOdrywet <- results(JH33_cds_test, contrast = c("Description2","Barke_DB","Barke_DRY+WET")) 

#inspect a result file
DB_vs_EOdrywet 
mcols(DB_vs_EOdrywet  , use.names=TRUE)
# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
DB_vs_EOdrywet_FDR005 <- DB_vs_EOdrywet[(rownames(DB_vs_EOdrywet)[which(DB_vs_EOdrywet$padj <0.05)]), ]

#enriched in rhizosphere (negative fold change, HID is at the second term of comparison)
DB_vs_EOdrywet_enriched <-  DB_vs_EOdrywet[(rownames(DB_vs_EOdrywet)[which(DB_vs_EOdrywet$log2FoldChange > 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
DB_vs_EOdrywet_enriched_FDR005 <- intersect(rownames(DB_vs_EOdrywet_FDR005), rownames(DB_vs_EOdrywet_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(DB_vs_EOdrywet_enriched_FDR005)


DB_vs_EOdrywet_enriched_FDR005_rhizo <- intersect(DB_vs_EOdrywet_enriched_FDR005,rhizo_EO_DB_enriched_FDR005)
length(DB_vs_EOdrywet_enriched_FDR005_rhizo)



JH33_DB_vs_EOdrywet <- prune_taxa(DB_vs_EOdrywet_enriched_FDR005_rhizo, JH33_plants)
JH33_DB_vs_EOdrywet

All_DB_vs_EOdrywet <- median(sample_sums(JH33_DB_vs_EOdrywet))/median(sample_sums(JH33_plants)) * 100
All_DB_vs_EOdrywet


#############################################################################################################
DB_all_vs_EOWET  <- results(JH33_cds_test, contrast = c("Description2","Barke_ALL","Barke_WET")) 
DB_all_vs_EOWET 
mcols(DB_all_vs_EOWET  , use.names=TRUE)
# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
DB_all_vs_EOWET_FDR005 <- DB_all_vs_EOWET[(rownames(DB_all_vs_EOWET)[which(DB_all_vs_EOWET$padj <0.05)]), ]

#enriched in rhizosphere (negative fold change, HID is at the second term of comparison)
DB_all_vs_EOWET_enriched <-  DB_all_vs_EOWET[(rownames(DB_all_vs_EOWET)[which(DB_all_vs_EOWET$log2FoldChange > 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
DB_all_vs_EOWET_enriched_FDR005 <- intersect(rownames(DB_all_vs_EOWET_FDR005), rownames(DB_all_vs_EOWET_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(DB_all_vs_EOWET_enriched_FDR005)


DB_all_vs_EOWET_enriched_FDR005_rhizo <- intersect(DB_all_vs_EOWET_enriched_FDR005,rhizo_EO_DB_ALL_enriched_FDR005)
length(DB_all_vs_EOWET_enriched_FDR005_rhizo)



JH33_DBall_vs_EOWET <- prune_taxa(DB_all_vs_EOWET_enriched_FDR005_rhizo, JH33_plants)
JH33_DBall_vs_EOWET

All_DBall_vs_EOWET <- median(sample_sums(JH33_DBall_vs_EOWET))/median(sample_sums(JH33_plants)) * 100
All_DBall_vs_EOWET

###############################################################################################
DB_all_vs_EO_drywet  <- results(JH33_cds_test, contrast = c("Description2","Barke_ALL","Barke_DRY+WET")) 
DB_all_vs_EO_drywet 
mcols(DB_all_vs_EO_drywet  , use.names=TRUE)
# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
DB_all_vs_EO_drywet_FDR005 <- DB_all_vs_EO_drywet[(rownames(DB_all_vs_EO_drywet)[which(DB_all_vs_EO_drywet$padj <0.05)]), ]

#enriched in rhizosphere (negative fold change, HID is at the second term of comparison)
DB_all_vs_EO_drywet_enriched <-  DB_all_vs_EO_drywet[(rownames(DB_all_vs_EO_drywet)[which(DB_all_vs_EO_drywet$log2FoldChange > 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
DB_all_vs_EO_drywet_enriched_FDR005 <- intersect(rownames(DB_all_vs_EO_drywet_FDR005), rownames(DB_all_vs_EO_drywet_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(DB_all_vs_EO_drywet_enriched_FDR005)


DB_all_vs_EO_drywet_enriched_FDR005_rhizo <- intersect(DB_all_vs_EO_drywet_enriched_FDR005,rhizo_EO_DB_ALL_enriched_FDR005)
length(DB_all_vs_EO_drywet_enriched_FDR005_rhizo)


JH33_DB_all_vs_EO_drywet <- prune_taxa(DB_all_vs_EO_drywet_enriched_FDR005_rhizo, JH33_plants)
JH33_DB_all_vs_EO_drywet

All_DB_all_vs_EO_drywet <- median(sample_sums(JH33_DB_all_vs_EO_drywet))/median(sample_sums(JH33_plants)) * 100
All_DB_all_vs_EO_drywet

###################################################################################
###################################################################################
#Data visualisation using UpSetR
###################################################################################

#Prepare the data for UpSetR visualisation

DB_enriched_vs_BulkEO_FDR_005<-rhizo_EO_DB_FDR005[rhizo_EO_DB_enriched_FDR005,]
DBall_enriched_vs_BulkEO_FDR_005<-rhizo_EO_DB_ALL_FDR005[rhizo_EO_DB_ALL_enriched_FDR005,]
EOwet_enriched_vs_BulkEO_FDR_005<-rhizo_EO_wet_FDR005[rhizo_EO_wet_enriched_FDR005,]
EOdrywet_enriched_vs_BulkEO_FDR_005<-rhizo_EO_drywet_FDR005[rhizo_EO_drywet_enriched_FDR005,]
#Extract baseMean (1st column)

#########I changed this to _EO_FDR_005_df so we know it is using the EO_Bulk control 25/06/2025 ###########
DB_enriched_vs_EO_FDR_005_df <- as.data.frame(DB_enriched_vs_BulkEO_FDR_005[ ,1])
DBall_enriched_vs_DBall_FDR_005_df <- as.data.frame(DBall_enriched_vs_BulkEO_FDR_005[ ,1])
EOwet_enriched_vs_EOwet_FDR_005_df <- as.data.frame(EOwet_enriched_vs_BulkEO_FDR_005[ ,1])
EOdrywet_enriched_vs_EOdrywet_FDR_005_df <- as.data.frame(EOdrywet_enriched_vs_BulkEO_FDR_005[ ,1])

#rename rows
rownames(DB_enriched_vs_EO_FDR_005_df)<-rhizo_EO_DB_enriched_FDR005
rownames(DBall_enriched_vs_DBall_FDR_005_df)<-rhizo_EO_DB_ALL_enriched_FDR005
rownames(EOwet_enriched_vs_EOwet_FDR_005_df)<-rhizo_EO_wet_enriched_FDR005
rownames(EOdrywet_enriched_vs_EOdrywet_FDR_005_df)<-rhizo_EO_drywet_enriched_FDR005

#rename columns
colnames(DB_enriched_vs_EO_FDR_005_df)<-c("counts_rhizo_DB")
colnames(DBall_enriched_vs_DBall_FDR_005_df)<-c("counts_rhizo_DBall")
colnames(EOwet_enriched_vs_EOwet_FDR_005_df)<-c("counts_rhizo_EOwet")
colnames(EOdrywet_enriched_vs_EOdrywet_FDR_005_df)<-c("counts_rhizo_EOdrywet")

DB_enriched_vs_EO_FDR_005_df[DB_enriched_vs_EO_FDR_005_df > 1] <- 1
DBall_enriched_vs_DBall_FDR_005_df[DBall_enriched_vs_DBall_FDR_005_df > 1] <- 1
EOwet_enriched_vs_EOwet_FDR_005_df[EOwet_enriched_vs_EOwet_FDR_005_df > 1] <- 1
EOdrywet_enriched_vs_EOdrywet_FDR_005_df[EOdrywet_enriched_vs_EOdrywet_FDR_005_df > 1] <- 1

dim(DB_enriched_vs_EO_FDR_005_df)

dim(DBall_enriched_vs_DBall_FDR_005_df)

dim(EOwet_enriched_vs_EOwet_FDR_005_df)

dim(EOdrywet_enriched_vs_EOdrywet_FDR_005_df)
 

#define a list of unique OTUs
Rhizo_list <- unique(c(rownames(DB_enriched_vs_EO_FDR_005_df), 
                       rownames(DBall_enriched_vs_DBall_FDR_005_df), 
                       rownames(EOwet_enriched_vs_EOwet_FDR_005_df), 
                       rownames(EOdrywet_enriched_vs_EOdrywet_FDR_005_df)))
length(Rhizo_list)


#Pellet
DB_enriched_merging<- as.data.frame(DB_enriched_vs_EO_FDR_005_df[Rhizo_list, ])
DBall_enriched_merging <- as.data.frame(DBall_enriched_vs_DBall_FDR_005_df[Rhizo_list, ])
EOwet_enriched_merging<- as.data.frame(EOwet_enriched_vs_EOwet_FDR_005_df[Rhizo_list, ])
EOdrywet_enriched_merging<- as.data.frame(EOdrywet_enriched_vs_EOdrywet_FDR_005_df[Rhizo_list, ])


colnames(DB_enriched_merging) <- c("Rhizo_DB")
colnames(DBall_enriched_merging) <- c("Rhizo_DBall")
colnames(EOwet_enriched_merging) <- c("Rhizo_EOwet")
colnames(EOdrywet_enriched_merging) <- c("Rhizo_EOdrywet")


row.names(DB_enriched_merging) <- as.vector(Rhizo_list)
row.names(DBall_enriched_merging) <- as.vector(Rhizo_list)
row.names(EOwet_enriched_merging) <- as.vector(Rhizo_list)
row.names(EOdrywet_enriched_merging) <- as.vector(Rhizo_list)

#Merge the dataset
enriched_Rhizo <- cbind(DB_enriched_merging, DBall_enriched_merging, EOwet_enriched_merging, EOdrywet_enriched_merging)

enriched_Rhizo <-as.data.frame(enriched_Rhizo)

######### INPUT FOR COMPLEX UPSET PLOT!!!!


#####visualization

##I think the Upset plot cannot be ordered by protocol because its is ordered by default by 'Set size' by default. Try keep.order= T
upset_ordered <- upset(
  enriched_Rhizo,
  sets = c("Rhizo_DB", "Rhizo_DBall", "Rhizo_EOwet", "Rhizo_EOdrywet"),
  sets.bar.color = c("#D55E00", "#E69F00", "#56B4E9", "#0072B2")
)
#order.by = "freq" , empty.intersections = "on")
upset_ordered


# Step 1: Identify your sets (these must be logical or binary columns)
sets <- c("Rhizo_DB", "Rhizo_DBall", "Rhizo_EOwet", "Rhizo_EOdrywet")

# Step 2: Use upset_data() directly on your original dataframe, not the plot object
intersection_data <- upset_data(enriched_Rhizo, intersect = sets)

# View the resulting dataframe
head(intersection_data)





# Extract intersection data
intersection_data <- upset_data(upset_ordered)
# Example: Extract information about 4 intersections
intersection_info <- intersection_data[1:4, ]

# Print the information
print(intersection_info)
##################################################################################################



##Comples Upset just intersection information
# Your merged data
enriched_Rhizo <- as.data.frame(enriched_Rhizo)
enriched_Rhizo$taxon <- rownames(enriched_Rhizo)

# Make sure all columns are numeric 
enriched_Rhizo$rhizo_EO_DB <- as.numeric(enriched_Rhizo$rhizo_EO_DB)
enriched_Rhizo$rhizo_EO_DB_ALL <- as.numeric(enriched_Rhizo$rhizo_EO_DB_ALL)
enriched_Rhizo$rhizo_EOwet <- as.numeric(enriched_Rhizo$rhizo_EOwet)
enriched_Rhizo$rhizo_EOdrywet <- as.numeric(enriched_Rhizo$rhizo_EOdrywet)

# Plot
complex_upset_plot <- upset(
  enriched_Rhizo,
  intersect = c("Rhizo_DB", "Rhizo_DBall", "Rhizo_EOwet", "Rhizo_EOdrywet"),
  name = "Enrichment",
  base_annotations = list(
    'Intersection size' = intersection_size(
      text = list(vjust = -0.5)
    )
  ),
  set_sizes = upset_set_size()
)

print(complex_upset_plot)

####Adding the abundance information
#####################################################################

## Because I changed the rownames for Upset plot, I need to add an additional column to design 'Description 3' with new colnames

design_abundance2 <- metadata3

design_abundance2 <- as.data.frame(design_abundance2)
colnames(design_abundance2)
# set rownames from the first column
rownames(design_abundance2) <- design_abundance2[[1]]

# optionally remove that first column now that it's in the rownames
design_abundance2 <- design_abundance2[, -1, drop = FALSE]

## 2) Abundance info################################


JH33_plants_otu_table <-otu_table(JH33_plants)

# transpose so samples are rows
otu_t <- as.data.frame((JH33_plants_otu_table))

otu_t <- otu_t[rownames(otu_t) %in% rownames(design_abundance2), ]



# Add Description3 grouping variable from metadata
# Re-assign Description3 correctly
otu_t$Description3 <- design_abundance2[rownames(otu_t), "Description3"]


sum(is.na(otu_t$Description3))  # Should be 0 ideally
rownames(otu_t)[is.na(otu_t$Description3)]

# Summarize by Description3, summing over selected numeric columns
otu_grouped <- otu_t %>%
  filter(!is.na(Description3)) %>%
  group_by(Description3) %>%
  summarise(
    .groups = "drop",
    across(where(is.numeric), ~ sum(.x, na.rm = TRUE))
  )

# Now correctly assign as rownames
otu_grouped <- column_to_rownames(otu_grouped, var = "Description3")

# Transpose back to have OTUs as rows
otu_grouped_t <- as.data.frame(t(otu_grouped))

# Add Phylum info and collapse
taxonomy_df <- as.data.frame(tax_table(JH33_plants))
all(rownames(otu_grouped_t) %in% rownames(taxonomy_df))  # should be TRUE

otu_grouped_t$phylum <- taxonomy_df[rownames(otu_grouped_t), "phylum"]

numeric_cols_phy <- sapply(otu_grouped_t, is.numeric)
head (numeric_cols_phy)

###if this returned FALSE, convert to numeric
#otu_grouped_t[, 1:4] <- lapply(otu_grouped_t[, 1:4], as.numeric)

# Recheck which columns are numeric
numeric_cols_phy <- sapply(otu_grouped_t, is.numeric)
head(numeric_cols_phy)


# Summarize by Phylum, summing numeric columns
phylum_grouped <- otu_grouped_t %>%
  group_by(phylum) %>%
  summarise_at(names(otu_grouped_t)[numeric_cols_phy], sum, na.rm = TRUE) %>%
  as.data.frame()

head(phylum_grouped)
######### INPUT FOR COMPLEX UPSET PLOT!!!!

#########################################

##Upset info into Phylum level#################################

##########################################

# Merge enriched_Rhizo with taxonomy table (assuming 'Phylum' column in taxonomy)

# Your merged data
enriched_Rhizo <- as.data.frame(enriched_Rhizo)
enriched_Rhizo$Taxon <- rownames(enriched_Rhizo)

# Make sure all columns are numeric 0/1 for presence/absence
enriched_Rhizo$Rhizo_DB <- as.numeric(enriched_Rhizo$Rhizo_DB)
enriched_Rhizo$Rhizo_DBall <- as.numeric(enriched_Rhizo$Rhizo_DBall)
enriched_Rhizo$Rhizo_EOwet <- as.numeric(enriched_Rhizo$Rhizo_EOwet)
enriched_Rhizo$Rhizo_EOdrywet <- as.numeric(enriched_Rhizo$Rhizo_EOdrywet)

#Tax table into data.frame
taxonomy_df<-tax_table(JH33_plants)

taxonomy_df<-as.data.frame(taxonomy_df)
head(taxonomy_df)
#how many NAs are at the Phylum level 
sum(is.na(taxonomy_df[, "phylum"]))
#[1] 857

taxonomy_df$Taxon <- rownames(taxonomy_df)

# Identify numeric columns
num_cols_enrich <- sapply(enriched_Rhizo, is.numeric)

# Merge with taxonomy: Phylum
enriched_phylum <- enriched_Rhizo %>%
  left_join(taxonomy_df, by = c("Taxon")) %>%
  group_by(phylum) %>%
  summarise_at(vars(names(enriched_Rhizo)[num_cols_enrich]), sum, na.rm = TRUE) %>%
  as.data.frame()

head (enriched_phylum)

###########################################################



#Use 'phylum_grouped for relative abundance
#The same for 'phylum_grouped' I need this as a character
phylum_grouped$phylum <- as.character(phylum_grouped$phylum)

#Used for Carmen's test #phylum_grouped$Abundance <- as.character(phylum_grouped$Abundance)

head(phylum_grouped)
is.character(phylum_grouped$phylum)
# Merge presence/absence with relative abundance
phylum_data <- enriched_phylum %>%
  left_join(phylum_grouped, by = "phylum")

phylum_data<-as.data.frame(phylum_data)
head(phylum_data)
#write.table(phylum_data, file="phylum_data_conplexUpset_0425.txt ", sep="\t")
colnames(phylum_data)
dim(phylum_data)
#5, 9 
phylum_data$phylum[is.na(phylum_data$phylum)] <- "Unclassified"

# Plot

# Custom palette with 16 strong colors
custom_colors <- colorRampPalette(c("#999999", "red", "#0072B2", "#F0E442", "#009E73", "#D55E00"))(16)

###########################

# Custom palette

custom_colors <- rep(c("steelblue", "red", "darkgreen", "orange", "purple", "grey"), 16)

custom_colors <- scales::hue_pal(c = 100, l = 60)(16)


# ComplexUpset plot
upset(
  phylum_data,
  intersect = c("Rhizo_DB", "Rhizo_DBall", "Rhizo_EOwet", "Rhizo_EOdrywet"),
  name = "Enrichment",
  base_annotations = list(
    'Intersection size' = intersection_size(text = list(vjust = -0.5)),
    'Abundance' = (
      ggplot(mapping = aes(x = .data$intersection, y = Rhizo_DB, fill = phylum)) +
        geom_col() +
        theme_minimal() +
        scale_fill_manual(values = custom_colors) +
        labs(y = "Abundance", x = NULL)
    )
  ),
  set_sizes = upset_set_size()
)

###############

##At  ASVs level##################################

# Add Taxon column to intersection matrix
enriched_Rhizo$Taxon <- rownames(enriched_Rhizo)

# Extract ASV abundances (already transposed)
otu_grouped_t$Taxon <- rownames(otu_grouped_t)
otu_grouped_t$phylum <- taxonomy_df[rownames(otu_grouped_t), "phylum"]

# Merge presence/absence with abundance and Phylum info
asv_data <- enriched_Rhizo %>%
  left_join(otu_grouped_t, by = "Taxon")

# Use one representative abundance column for plotting
# You can replace 'DB_rhizo_sample' with a sample group or mean abundance
# Use one representative abundance column for plotting
# You can replace 'DB_rhizo_sample' with a sample group or mean abundance
asv_data$Abundance <- asv_data$Rhizo_DB.y 

# Ensure Phylum is character
asv_data$phylum <- as.character(asv_data$phylum)

asv_data$phylum[is.na(asv_data$phylum)] <- "Unclassified"


library(dplyr)
library(tidyr)

sets <- c("Rhizo_DB.x","Rhizo_DBall.x","Rhizo_EOwet.x","Rhizo_EOdrywet.x")

asv_data_clean <- asv_data %>%
  mutate(across(all_of(sets), ~replace_na(., 0))) %>%  # NAs -> 0
  mutate(across(all_of(sets), ~. > 0))                 # convert to TRUE/FALSE


c16<-c("darkorange4", "#E31A1C", "#009E73", "#6A3D9A", "#D55E00", "black", "palegreen2", "#88CCEE", "#0072B2", "#DDCC77", "#999999", "#999933","gold1" , "darkturquoise", "#CC79A7", "brown")
# Plot at ASV level: bar per ASV, colored by phylum
p <- upset(
  asv_data_clean,
  intersect = c("Rhizo_DB.x", "Rhizo_DBall.x", "Rhizo_EOwet.x", "Rhizo_EOdrywet.x"),
  name = "ASV Enrichment",
  base_annotations = list(
    'Intersection size' = intersection_size(
      aes(y = Abundance, fill = phylum),
      text = list(vjust = -0.5)
    )
  ),
  set_sizes = upset_set_size()
) + scale_fill_manual(values = c16) + labs(fill = "phylum", y = "ASV Abundance")

p


upset(
  asv_data_clean,
  intersect = c("Rhizo_DB.x", "Rhizo_DBall.x", "Rhizo_EOwet.x", "Rhizo_EOdrywet.x"),
  name = "ASV Enrichment",
  base_annotations = list(
    'Intersection size' = intersection_size(
      aes(y = Abundance, fill = phylum),
      text = list(vjust = -0.5)
    )
  ),
  set_sizes = upset_set_size()
) + scale_fill_manual(values = c16) + labs(fill = "phylum", y = "ASV Abundance")


unique(asv_data$phylum)


###Log scale
upsetlog_clean <-upset(
  asv_data_clean,
  intersect = c("Rhizo_DB.x", "Rhizo_DBall.x", "Rhizo_EOwet.x", "Rhizo_EOdrywet.x"),
  name = "ASV Enrichment",
  base_annotations = list(
    'Intersection size' = (
      intersection_size(
        aes(y = log10(Abundance), fill = phylum),
        text = list(vjust = -0.5)
      ) +
        scale_fill_manual(values = c16) +
        labs(fill = "phylum", y = "log10(Phylum Abundance)"))))

upsetlog_clean


c16 <-c("darkorange4", "#E31A1C", "#009E73", "#6A3D9A", "#D55E00", "black", "palegreen2", "#88CCEE", "#0072B2", "#DDCC77", "#999999", "#999933","gold1" , "darkturquoise", "#CC79A7", "brown")

unique(asv_data$phylum)
#replace NA with 'Unclassified'
asv_data$phylum[is.na(asv_data$phylum)] <- "Unclassified"


asv_data$phylum <- factor(asv_data$phylum, levels = names(custom_colors)) 
head(asv_data)


###Log scale
upsetlog <-upset(
  asv_data,
  intersect = c("Rhizo_DB", "Rhizo_DBall", "Rhizo_EOwet", "Rhizo_EOdrywet"),
  name = "ASV Enrichment",
  base_annotations = list(
    'Intersection size' = (
      intersection_size(
        aes(y = log10(Abundance), fill = phylum),
        text = list(vjust = -0.5)
      ) +
        scale_fill_manual(values = custom_colors) +
        labs(fill = "phylum", y = "log10(Phylum Abundance)"))))

upsetlog


############################# End Upset
##Heatmap phylum level
###############################################



# ## agglomerate at the Family taxonomic rank, this Is a new phyloseq object

#JH33_plants is all the rhizosphere sample taxa, not the 'enriched in the rhizopshere', so I will make a new phlyloseq withjust the
#taxa enriched in the rhizosphere


JH33_EO_vs_Bulk <- prune_taxa(rhizo_EO_wet_enriched_FDR005, JH33_plants)
JH33_EO_vs_Bulk

JH33_EOdrywet_vs_Bulk <- prune_taxa(rhizo_EO_drywet_enriched_FDR005, JH33_plants)
JH33_EOdrywet_vs_Bulk


JH33_DB_vs_EO_Bulk <- prune_taxa(rhizo_EO_DB_enriched_FDR005, JH33_plants)
JH33_DB_vs_EO_Bulk


JH33_DBall_vs_EO_Bulk <- prune_taxa(rhizo_EO_DB_ALL_enriched_FDR005, JH33_plants)
#this kept failing, so I checked why and there are no taxa from DB_All enriched in the rhizosphere? 
JH33_taxa <- taxa_names(JH33_plants)
sum(JH33_taxa %in% taxa_names(rhizo_DB_all_enriched_FDR005))

JH33_DBall_vs_EO_Bulk
#I have had to remove DBALL to make this run
JH33_rhizo_enriched<-merge_phyloseq(JH33_EO_vs_Bulk,JH33_EOdrywet_vs_Bulk,JH33_DB_vs_EO_Bulk)#,JH33_DBall_vs_Bulk)

##Here we can use JH33_plants or JH33_rhizo_enriched

FE_phylum <- tax_glom(JH33_rhizo_enriched, taxrank="phylum") 
# Suppose you have the sample IDs in a column called "SampleID"
rownames(design_abundance2) <- design_abundance2$'Illumina code'

# Then reorder:
design_abundance2 <- design_abundance2[sample_names(FE_phylum), , drop = FALSE]

# Now assign
sample_data(FE_phylum) <- design_abundance2

design_metadata <- sample_data(design_abundance2)

# Now assign it in a safe way
FE_phylum <- merge_phyloseq(FE_phylum, design_metadata)

##############################################
##############################################################
## FE_phylum is the phyloseq of enriched rhizosphere reads across protocols
###############################################################
#workaround code Deseq in R derived from JH07 original code
#extract count data : 
Phylum_counts <- otu_table(FE_phylum)
countData = as.data.frame(Phylum_counts)
head(countData)
colnames(Phylum_counts)

Phylum_tax<-tax_table(FE_phylum)

#the design file containing sample information
colData = design_abundance2
rownames(colData) <- colData$`Illumina code`
rownames(colData)
colData


#construct a DESeq dataset combining count data and sample information
#I changed the design to 'Method' for me 
Phylum_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design=~Description2)

#execute the differential count analysis with the function DESeq 
Phylum_cds_test <- DESeq(Phylum_cds, fitType="local", betaPrior = FALSE) 
levels(Phylum_cds_test$Description2)

saveRDS(Phylum_cds_test, file="Phylum_cds_test_ASVs_TEAGASC.rds")

Phylum_cds_test<-Phylum_cds_test_ASVs_TEAGASC


############################################################################### Table S4
EO_vs_EOdrywet_Phylum <- results(Phylum_cds_test, contrast = c("Description2","Barke_WET","Barke_DRY+WET")) 

#inspect a result file
EO_vs_EOdrywet_Phylum 
mcols(EO_vs_EOdrywet_Phylum  , use.names=TRUE)
# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
EO_vs_EOdrywet_Phylum_FDR005 <- EO_vs_EOdrywet_Phylum[(rownames(EO_vs_EOdrywet_Phylum)[which(EO_vs_EOdrywet_Phylum$padj <0.05)]), ]

#enriched in rhizosphere (negative fold change, HID is at the second term of comparison)
EO_vs_EOdrywet_Phylum_enriched <-  EO_vs_EOdrywet_Phylum[(rownames(EO_vs_EOdrywet_Phylum)[which(EO_vs_EOdrywet_Phylum$log2FoldChange > 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
EO_vs_EOdrywet_Phylum_enriched_FDR005 <- intersect(rownames(EO_vs_EOdrywet_Phylum_FDR005), rownames(EO_vs_EOdrywet_Phylum_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(EO_vs_EOdrywet_Phylum_enriched_FDR005)
#0

Phylum_tax[EO_vs_EOdrywet_Phylum_enriched_FDR005]
###################################################################################
DB_vs_EOwet_Phylum <- results(Phylum_cds_test, contrast = c("Description2","Barke_DB","Barke_WET")) 

#inspect a result file
DB_vs_EOwet_Phylum 
mcols(DB_vs_EOwet_Phylum  , use.names=TRUE)
# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
DB_vs_EOwet_Phylum_FDR005 <- DB_vs_EOwet_Phylum[(rownames(DB_vs_EOwet_Phylum)[which(DB_vs_EOwet_Phylum$padj <0.05)]), ]

#enriched in rhizosphere (negative fold change, Barke-wet is at the second term of comparison)
DB_vs_EOwet_Phylum_enriched <-  DB_vs_EOwet_Phylum[(rownames(DB_vs_EOwet_Phylum)[which(DB_vs_EOwet_Phylum$log2FoldChange > 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
DB_vs_EOwet_Phylum_enriched_FDR005 <- intersect(rownames(DB_vs_EOwet_Phylum_FDR005), rownames(DB_vs_EOwet_Phylum_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(DB_vs_EOwet_Phylum_enriched_FDR005)

#1
Phylum_tax[DB_vs_EOwet_Phylum_enriched_FDR005]
#Ascomycota, NA

############################################################################################################

DB_vs_EOdrywet_Phylum <- results(Phylum_cds_test, contrast = c("Description2","Barke_DB","Barke_DRY+WET")) 

#inspect a result file
DB_vs_EOdrywet_Phylum 
mcols(DB_vs_EOdrywet_Phylum  , use.names=TRUE)
# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
DB_vs_EOdrywet_Phylum_FDR005 <- DB_vs_EOdrywet_Phylum[(rownames(DB_vs_EOdrywet_Phylum)[which(DB_vs_EOdrywet_Phylum$padj <0.05)]), ]

#enriched in rhizosphere (negative fold change, HID is at the second term of comparison)
DB_vs_EOdrywet_Phylum_enriched <-  DB_vs_EOdrywet_Phylum[(rownames(DB_vs_EOdrywet_Phylum)[which(DB_vs_EOdrywet_Phylum$log2FoldChange > 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
DB_vs_EOdrywet_Phylum_enriched_FDR005 <- intersect(rownames(DB_vs_EOdrywet_Phylum_FDR005), rownames(DB_vs_EOdrywet_Phylum_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(DB_vs_EOdrywet_Phylum_enriched_FDR005)

#1 
Phylum_tax[DB_vs_EOdrywet_Phylum_enriched_FDR005]
#Ascomycota, NA
#############################################################################################################
DB_all_vs_EOt_Phylum  <- results(Phylum_cds_test, contrast = c("Description2","Barke_ALL","Barke_WET")) 
DB_all_vs_EOt_Phylum 
mcols(DB_all_vs_EOt_Phylum  , use.names=TRUE)
# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
DB_all_vs_EOt_Phylum_FDR005 <- DB_all_vs_EOt_Phylum[(rownames(DB_all_vs_EOt_Phylum)[which(DB_all_vs_EOt_Phylum$padj <0.05)]), ]

#enriched in rhizosphere (negative fold change, HID is at the second term of comparison)
DB_all_vs_EOt_Phylum_enriched <-  DB_all_vs_EOt_Phylum[(rownames(DB_all_vs_EOt_Phylum)[which(DB_all_vs_EOt_Phylum$log2FoldChange > 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
DB_all_vs_EOt_Phylum_enriched_FDR005 <- intersect(rownames(DB_all_vs_EOt_Phylum_FDR005), rownames(DB_all_vs_EOt_Phylum_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(DB_all_vs_EOt_Phylum_enriched_FDR005)


Phylum_tax[DB_all_vs_EOt_Phylum_enriched_FDR005]

###############################################################################################
DB_all_vs_EO_drywet_Phylum  <- results(Phylum_cds_test, contrast = c("Description2","Barke_ALL","Barke_DRY+WET")) 
DB_all_vs_EO_drywet_Phylum 
mcols(DB_all_vs_EO_drywet_Phylum  , use.names=TRUE)
# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
DB_all_vs_EO_drywet_Phylum_FDR005 <- DB_all_vs_EO_drywet_Phylum[(rownames(DB_all_vs_EO_drywet_Phylum)[which(DB_all_vs_EO_drywet_Phylum$padj <0.05)]), ]

#enriched in rhizosphere (negative fold change, HID is at the second term of comparison)
DB_all_vs_EO_drywet_Phylum_enriched <-  DB_all_vs_EO_drywet_Phylum[(rownames(DB_all_vs_EO_drywet_Phylum)[which(DB_all_vs_EO_drywet_Phylum$log2FoldChange < 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
DB_all_vs_EO_drywet_Phylum_enriched_FDR005 <- intersect(rownames(DB_all_vs_EO_drywet_Phylum_FDR005), rownames(DB_all_vs_EO_drywet_Phylum_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(DB_all_vs_EO_drywet_Phylum_enriched_FDR005)

Phylum_tax[DB_all_vs_EO_drywet_Phylum_enriched_FDR005]
#Olpidiomycota, NA
############

# Combine significant OTUs from all contrasts
all_enriched_phyla <- unique(c(
  EO_vs_EOdrywet_Phylum_enriched_FDR005,
  DB_vs_EOwet_Phylum_enriched_FDR005,
  DB_vs_EOdrywet_Phylum_enriched_FDR005,
  DB_all_vs_EOt_Phylum_enriched_FDR005,
  DB_all_vs_EO_drywet_Phylum_enriched_FDR005
))

# Extract normalized counts
normalized_counts <- assay(varianceStabilizingTransformation(Phylum_cds_test, blind=FALSE))[all_enriched_phyla, ]

# Create annotation for samples
annotation_col <- as.data.frame(colData(Phylum_cds_test)[, "Description2", drop=FALSE])

sample_data(FE_phylum)

# Extract Phylum names for enriched OTUs
phylum_labels <- as.character(Phylum_tax[rownames(normalized_counts), "phylum"])
rownames(normalized_counts) <- make.unique(phylum_labels)


# Plot heatmap with row annotation
pheatmap::pheatmap(
  normalized_counts,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  scale = "row",
  color = viridis::viridis(100),
  fontsize_row = 8
)

#####End


