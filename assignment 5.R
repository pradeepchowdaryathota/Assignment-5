library(readxl)
library(readr)
library(dplyr)
library(viridis)
library(pheatmap)



# a.Load data
df1<- read_excel("Gene_Expression_Data.xlsx")
df2 <- read_csv("Gene_Information.csv")
df3 <- read.delim("Sample_Information.tsv")

# b. Change Sample Names
mapping <- df3$group
df1_renamed <- c("Probe_ID", mapping)

# c. Split Merged Data and Split the DataFrame into 'tumor' and 'normal'
tumor_columns <- grep("tumor", df1_renamed, value = TRUE)
tumor_df <- df1[, c("Probe_ID", tumor_columns)]

normal_columns <- grep("normal", df1_renamed, value = TRUE)
normal_df <- df1[, c("Probe_ID", normal_columns)]

# d. Compute the average expression for all genes from the 2 data sets from part d
tumor_avg_expression <- rowMeans(tumor_df[, -1], na.rm = TRUE)
normal_avg_expression <- rowMeans(normal_df[, -1], na.rm = TRUE)

# e. Determine the fold change for each Probe between the two groups ((Tumour – Control) / Control)
fold_change <- (tumor_avg_expression - normal_avg_expression) / normal_avg_expression

# f. Use the data from part e and "Gene_Information.csv" to identify all genes fold change magnitude (absolute value) was greater than 5
fold_change_genes <- data.frame(Probe_ID = df2$Probe_ID, Fold_Change = fold_change)
fold_change_genes <- fold_change_genes[abs(fold_change_genes$Fold_Change) > 5, ]

# g. Add a column to the result of part f to include if the gene was higher expressed in "Normal" or "Tumor" samples
fold_change_genes$Higher_Expression <- ifelse(fold_change_genes$Fold_Change > 0, "Tumor", "Normal")

# Merging df2 and fold_change_genes
fold_change_genes <- merge(fold_change_genes, df2[, c("Probe_ID", "Chromosome")], by = "Probe_ID", all.x = TRUE)

#2 a. Perform exploratory data analysis on the genes from part 1g
plot(1:nrow(fold_change_genes), fold_change_genes$Fold_Change, xlab = "Index", ylab = "Fold Change", main = "Scatter Plot of Fold Change Values")

#b.create histogram
chromosome_counts <- table(fold_change_genes$Chromosome)
hist(chromosome_counts, breaks = 1:max(chromosome_counts), col = "yellow", 
     main = "Distribution of DEGs by Chromosome", xlab = "Chromosome", ylab = "Frequency")

# c. Histogram of DEGs by Chromosome Segregated by Sample Type
normal_genes <- subset(fold_change_genes, Higher_Expression == "Normal")
tumor_genes <- subset(fold_change_genes, Higher_Expression == "Tumor")

par(mfrow = c(1, 2), mar = c(5, 4, 4, 2) + 0.1)  # Split the plotting area into two panels
hist(as.factor(normal_genes$Chromosome), col = "blue", main = "Normal", xlab = "Chromosome", ylab = "Frequency", breaks = 20, labels = TRUE)
hist(as.factor(tumor_genes$Chromosome), col = "red", main = "Tumor", xlab = "Chromosome", ylab = "Frequency", breaks = 20, labels = TRUE)

# d. Bar Chart of Upregulated and Downregulated DEGs
upregulated_percentage <- sum(fold_change_genes$Higher_Expression == "Tumor") / nrow(fold_change_genes) * 100
downregulated_percentage <- 100 - upregulated_percentage

barplot(c(upregulated_percentage, downregulated_percentage), names.arg = c("Upregulated", "Downregulated"), col = c("black", "orange"), xlab = "Gene Expression", ylab = "Percentage", main = "Percentage of DEGs Upregulated and Downregulated in Tumor Samples")

#e.
# Remove the "Probe_ID" column from df1_renamed
dfx <- df1_renamed[, -which(names(df1_renamed) == "Probe_ID")]
# Create the heatmap
pheatmap(dfx, cluster_rows = FALSE, cluster_cols = FALSE, main = "Heatmap of Gene Expression by Sample")
#f.
# Create the clustermap
pheatmap(dfx, main = "Clustermap of Gene Expression by Sample")


