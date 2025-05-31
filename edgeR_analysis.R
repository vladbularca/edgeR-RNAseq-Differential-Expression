#Loading Packages

if (!require("BiocManager", quietly = TRUE)) install.packages("BioCManager")
BiocManager::install("Rsubread") 
library(Rsubread)
library(tibble)
library(edgeR)
library(limma)
library(dplyr)
library(gridExtra)


# Importing the gene counts table produced using STAR

counts <- read.table("gene_counts.txt", header = TRUE)

# Updating column names with sample names based on stages of development

colnames(counts) [7] = "MB1"
colnames(counts) [8] = "MB2"
colnames(counts) [9] = "MB3"

colnames(counts) [10] = "TB1"
colnames(counts) [11] = "TB2"
colnames(counts) [12] = "TB3"

colnames(counts) [13] = "EB1"
colnames(counts) [14] = "EB2"
colnames(counts) [15] = "EB3"

# Preparing the different components needed for a DGEList object. First is a matrix containing raw counts for each sample 

count_matrix <- counts[, 7:ncol(counts)]

# Sample_info contains the sample information and groups to identify which stage each sample falls into

sample_info <- data.frame(
  sample = colnames(count_matrix),
  group = c("Mature biofilm", "Mature biofilm", "Mature biofilm", "Thin biofilm", "Thin biofilm", "Thin biofilm", "Early biofilm", "Early biofilm", "Early biofilm"))

# Creating a DGEL list and adding the genes component from the original counts table

DGE <- DGEList(counts = count_matrix, group = sample_info$group)
DGE$genes <- counts[, 1:6]


# Adding library sizes, setting normalization factors to 1 for each sample and adding these to the sample_info table

librarysize = c(9030851, 6912070, 8765210, 5980833, 6089675,6309292, 6734945, 6897699, 7528054)

normfactor = c(1, 1, 1, 1, 1, 1, 1, 1, 1)
sample_info$libsize = librarysize
sample_info$norm.factors = normfactor

# Sorting the data by total counts, checking to see if there are duplicates and removing them

o <- order(rowSums(DGE$counts), decreasing=TRUE)
DGE <- DGE[o,]
d <- duplicated(DGE$genes$Geneid)
DGE <- DGE[!d, ]


# Filtering out genes that have low expression across all samples, updating the library sizes and then normalizing the data:

keep<-filterByExpr(DGE)
table(keep)
DGE<-DGE[keep,,keep.lib.sizes=FALSE]
DGE<- calcNormFactors(DGE,method="TMM") 


# MDS plot to observe sample clustering and check initial data quality – see Figure 1

plotMDS(DGE, col=c("red", "blue", "green")[DGE$samples$group], main="MDS Plot of Biofilm Samples")
legend("bottomright", legend=levels(DGE$samples$group), col=c("red", "blue", "green"), pch=16)

# Creating a factor variable to define the experimental conditions

group <- factor(c("Mature", "Mature", "Mature", 
                  "Thin", "Thin", "Thin", 
                  "Early", "Early", "Early"))


# Creating a design matrix and labeling it based on the levels from group

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Estimating dispersion, which measures within treatment variation:

DGE2 <- estimateDisp(DGE, design, robust=TRUE)

plotBCV(DGE2)


# Creating a contrast matrix so that different conditions can be compared 

contrast.matrix <- makeContrasts(
  Mature_vs_Thin = Mature - Thin,
  Mature_vs_Early = Mature - Early,
  Thin_vs_Early = Thin - Early,
  Early_vs_Later = Early - (Thin + Mature) / 2,
  levels = design
)

# Early vs. Later to get DE genes between early stage and both later stages combined 

# Fitting the quasi-likelihood negative binomial generalized log-linear model to the DGE2

fit<-glmQLFit(DGE2, design)

# Determining the Differentially expressed (DE) genes for Mature vs. Thin stages

M.Tfit <- glmQLFTest(fit, contrast=contrast.matrix[,"Mature_vs_Thin"])
DEG <- topTags(M.Tfit, n=Inf)$table

# Filtering DE genes based on FDR instead of p-value to ensure less false positives are interpreted as significant. Using FDR < 0.05 as the threshold cut-off to determine the list of differentially expressed genes 

DEG_sig <- DEG[DEG$FDR < 0.05, ]
DEG_sig <- select(DEG_sig, 1, 7, 8, 10, 11) 
rownames(DEG_sig) <- NULL

# Creating a table with top 10 DE genes for Mature vs. Thin

grid.table(head(DEG_sig, 10))
summary(decideTests(M.Tfit,method="fdr",p.value=0.05)) 


# Determining the DE genes for Mature vs. Early stages

M.Efit <- glmQLFTest(fit, contrast=contrast.matrix[,"Mature_vs_Early"])
DEG2 <- topTags(M.Efit, n=Inf)$table
DEG2_sig <- DEG2[DEG2$FDR < 0.05, ]
DEG2_sig <- select(DEG2_sig, 1, 7, 8, 10, 11) 
rownames(DEG2_sig) <- NULL


grid.table(head(DEG2_sig, 10))
summary(decideTests(M.Efit,method="fdr",p.value=0.05)) 


# Determining the DE genes for Thin vs. Early stages

T.Efit <- glmQLFTest(fit, contrast=contrast.matrix[,"Thin_vs_Early"])
DEG3 <- topTags(T.Efit, n=Inf)$table
DEG3_sig <- DEG3[DEG3$FDR < 0.05, ]
DEG3_sig <- select(DEG3_sig, 1, 7, 8, 10, 11) 
rownames(DEG3_sig) <- NULL

grid.table(head(DEG3_sig, 10))
summary(decideTests(T.Efit,method="fdr",p.value=0.05)) 

# Determining the DE genes for Thin vs. both Later stages

E.Lfit <- glmQLFTest(fit, contrast=contrast.matrix[,"Early_vs_Later"])
DEG4 <- topTags(E.Lfit, n=Inf)$table
DEG4_sig <- DEG4[DEG4$FDR < 0.05, ]
DEG4_sig <- select(DEG4_sig, 1, 7, 8, 10, 11) 
rownames(DEG4_sig) <- NULL

grid.table(head(DEG4_sig, 10))
summary(decideTests(E.Lfit,method="fdr",p.value=0.05)) 



# Generating Smear plots for the 3 pairwise comparisons

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


# Exporting a list of all DE genes for pairwise comparisons (all_DE_genes.txt) and a list of all DE genes between the Early stage and both later stages combined (early-later_DE_genes.txt) 


all_DE_genes <- rbind(DEG_sig, DEG2_sig, DEG3_sig)
rownames(all_DE_genes) <- NULL

write.table(all_DE_genes, file = "all_DE_genes.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(DEG4_sig, file = "early-later_DE_genes.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

