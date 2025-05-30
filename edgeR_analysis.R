#Loading Packages

if (!require("BiocManager", quietly = TRUE)) install.packages("BioCManager")
BiocManager::install("Rsubread") 
library(Rsubread)
library(tibble)
library(edgeR)
library(limma)
library(dplyr)
library(gridExtra)


# Reading in the gene counts table produced using STAR

counts <- read.table("gene_counts.txt", header = TRUE)

# Updating sample names based on the stages of development

colnames(counts) [7] = "MB1"
colnames(counts) [8] = "MB2"
colnames(counts) [9] = "MB3"

colnames(counts) [10] = "TB1"
colnames(counts) [11] = "TB2"
colnames(counts) [12] = "TB3"

colnames(counts) [13] = "EB1"
colnames(counts) [14] = "EB2"
colnames(counts) [15] = "EB3"

count_matrix <- counts[, 7:ncol(counts)]

sample_info <- data.frame(
  sample = colnames(count_matrix),
  group = c("Mature biofilm", "Mature biofilm", "Mature biofilm", "Thin biofilm", "Thin biofilm", "Thin biofilm", "Early biofilm", "Early biofilm", "Early biofilm")
)


DGE <- DGEList(counts = count_matrix, group = sample_info$group)

DGE$genes <- counts[, 1:6]




librarysize = c(9030851, 6912070, 8765210, 5980833, 6089675,6309292, 6734945, 6897699, 7528054)

normfactor = c(1, 1, 1, 1, 1, 1, 1, 1, 1)

sample_info$libsize = librarysize

sample_info$norm.factors = normfactor

o <- order(rowSums(DGE$counts), decreasing=TRUE)

DGE <- DGE[o,]

d <- duplicated(DGE$genes$Geneid)

DGE <- DGE[!d, ]


# Since the data has not been already filtered, I will apply those steps here to this analyis

keep<-filterByExpr(DGE)

table(keep)

DGE<-DGE[keep,,keep.lib.sizes=FALSE]

DGE<- calcNormFactors(DGE,method="TMM") 




plotMDS(DGE, col=c("red", "blue", "green")[DGE$samples$group], main="MDS Plot of Biofilm Samples")
legend("bottomright", legend=levels(DGE$samples$group), col=c("red", "blue", "green"), pch=16)

group <- factor(c("Mature", "Mature", "Mature", 
                  "Thin", "Thin", "Thin", 
                  "Early", "Early", "Early"))


design <- model.matrix(~ 0 + group)

colnames(design) <- levels(group)

DGE <- estimateDisp(DGE, design)


DGE2 <- estimateDisp(DGE, design, robust=TRUE)

plotBCV(DGE2)
plotBCV(DGE)

# Creating a contrast matrix so that I can compare the different conditions:

contrast.matrix <- makeContrasts(
  Mature_vs_Thin = Mature - Thin,
  Mature_vs_Early = Mature - Early,
  Thin_vs_Early = Thin - Early,
  Early_vs_Later = Early - (Thin + Mature) / 2,
  levels = design
)

fit<-glmQLFit(DGE2, design)


#Better for gene specific variability
#Log com is average expression level of a gene across all samples




M.Tfit <- glmQLFTest(fit, contrast=contrast.matrix[,"Mature_vs_Thin"])

DEG <- topTags(M.Tfit, n=Inf)$table

# Using FDR to prevent less false positives from being interpreted as significant.
# To prevent false positives I will be using FDR < 0.05 as the threshold cut-off for deteriming the list of differential expressed genes

DEG_sig <- DEG[DEG$FDR < 0.05, ]

DEG_sig <- select(DEG_sig, 1, 7, 8, 10, 11) 

rownames(DEG_sig) <- NULL

grid.table(head(DEG_sig, 10))


summary(decideTests(M.Tfit,method="fdr",p.value=0.05)) 



M.Efit <- glmQLFTest(fit, contrast=contrast.matrix[,"Mature_vs_Early"])

DEG2 <- topTags(M.Efit, n=Inf)$table

DEG2_sig <- DEG2[DEG2$FDR < 0.05, ]


DEG2_sig <- select(DEG2_sig, 1, 7, 8, 10, 11) 

rownames(DEG2_sig) <- NULL

grid.table(head(DEG2_sig, 10))

summary(decideTests(M.Efit,method="fdr",p.value=0.05)) 




T.Efit <- glmQLFTest(fit, contrast=contrast.matrix[,"Thin_vs_Early"])

DEG3 <- topTags(T.Efit, n=Inf)$table

DEG3_sig <- DEG3[DEG3$FDR < 0.05, ]

DEG3_sig <- select(DEG3_sig, 1, 7, 8, 10, 11) 

rownames(DEG3_sig) <- NULL

grid.table(head(DEG3_sig, 10))


summary(decideTests(T.Efit,method="fdr",p.value=0.05)) 


E.Lfit <- glmQLFTest(fit, contrast=contrast.matrix[,"Early_vs_Later"])

DEG4 <- topTags(E.Lfit, n=Inf)$table

DEG4_sig <- DEG4[DEG4$FDR < 0.05, ]


DEG4_sig <- select(DEG4_sig, 1, 7, 8, 10, 11) 

rownames(DEG4_sig) <- NULL

grid.table(head(DEG4_sig, 10))

summary(decideTests(E.Lfit,method="fdr",p.value=0.05)) 







# Graphing volcano plot
dt1 <- decideTests(M.Tfit)

isDE <- as.logical(dt1)

DEnames <- rownames(M.Tfit)[isDE] 
plotSmear(M.Tfit, de.tags=DEnames, main="Smear Plot of Differential Expression - Mature vs. Thin") 
abline(h=c(-1,1), col="blue")


dt2 <- decideTests(M.Efit)

isDE2 <- as.logical(dt2)

DEnames2 <- rownames(M.Efit)[isDE2] 
plotSmear(M.Efit, de.tags=DEnames2,  main="Smear Plot of Differential Expression - Mature vs. Early") 
abline(h=c(-1,1), col="blue")


dt3 <- decideTests(T.Efit)

isDE3 <- as.logical(dt3)

DEnames3 <- rownames(T.Efit)[isDE3] 
plotSmear(T.Efit, de.tags=DEnames3,  main="Smear Plot of Differential Expression - Thin vs. Early") 
abline(h=c(-1,1), col="blue")





#Genes that are significantly differentially expressed are typically highlighted



# Early vs Thin + Early vs. Mature
# Which are expressed common amongst both of these tables?



#Export lists. 

all_DE_genes <- rbind(DEG_sig, DEG2_sig, DEG3_sig)

rownames(all_DE_genes) <- NULL

write.table(all_DE_final, file = "all_DE_genes.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(DEG4_sig, file = "early-later_DE_genes.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)















#edgeR corrects for sequencing depth automatically (TMM normalization).
#CPM is useful for exploratory data analysis but not for differential expression testing.
#Gene length correction is only needed when comparing expression between genes, not for differential expression analysis.

# library size total number of raw sequencing reads in a sample
#RNA-Seq Sequencing is the Depth Total number of reads per sample, Normalization used to account for sequencing depth across samples. 






