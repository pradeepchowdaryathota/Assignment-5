library(readxl)
library(readr)
library(dplyr)
library(viridis)
library(pheatmap)

# a.Load data
df1<- read_excel("Gene_Expression_Data.xlsx")
df2 <- read_csv("Gene_Information.csv")
df3 <- read_tsv("Sample_Information.tsv")

#b
mapping <- setNames(df3$group, df3$sample_names)
colnames(df1)[-1] <- mapping[match(colnames(df1)[-1], names(mapping))]
df1_renamed <- df1

#c
tumor_columns <- grep("tumor", colnames(df1_renamed), value = TRUE)
normal_columns <- grep("normal", colnames(df1_renamed), value = TRUE)

if ("Probe_ID" %in% colnames(df1_renamed)) {
  tumor_columns <- c("Probe_ID", tumor_columns)
  normal_columns <- c("Probe_ID", normal_columns)
}
tumor_df <- df1_renamed[, tumor_columns, drop = FALSE]
normal_df <- df1_renamed[, normal_columns, drop = FALSE]

#d
# Compute the average expression for all genes from the two data sets
tumor_avg_expression <- rowMeans(tumor_df[, -1])
normal_avg_expression <- rowMeans(normal_df[, -1])
# Combine average expression into a single DataFrame
avg_expression <- data.frame(Probe_ID = tumor_df$Probe_ID, Tumor_Avg_Expression = tumor_avg_expression, Normal_Avg_Expression = normal_avg_expression)

#e.Determine the log2 fold change for each Probe between the two groups log2((Tumour – Control) / Control)
log2_fold_change <- log2((tumor_avg_expression - normal_avg_expression) / normal_avg_expression)
# Combine log2 fold change into the DataFrame
fold_change <- data.frame(Probe_ID = tumor_df$Probe_ID, Log2_Fold_Change = log2_fold_change)

#f. Use the data from part e and “Gene_Information.csv” to identify all genes fold change magnitude (absolute value) was greater than 5
# Merge fold change data with gene information
fold_change_genes_info <- merge(fold_change, df2, by = "Probe_ID")
# Filter genes where fold change magnitude (absolute value) is greater than 5
significant_genes <- fold_change_genes_info[abs(fold_change_genes_info$Log2_Fold_Change) > 5, ]
significant_genes

#g. Add a column to the result of part f to include if the gene was higher expressed in “Normal” or “Tumor” samples
significant_genes <- significant_genes %>%
  mutate(HigherIn = ifelse(fold_change_genes_info > 0, "Tumor", "Normal"))


# Create a new column 'Higher_Expression' based on the condition
fold_change_genes_info$Higher_Expression <- ifelse(fold_change_genes_info$Log2_Fold_Change > 0, "Tumor", "Normal")


#2.
# Convert 'Chromosome' to factor
fold_change_genes_info$Chromosome <- factor(fold_change_genes_info$Chromosome, levels=unique(fold_change_genes_info$Chromosome))

# Create the histogram
ggplot(fold_change_genes_info, aes(x = Chromosome)) +
  geom_histogram(stat = "count", fill = "yellow") + 
  labs(title = "Distribution of DEGs by Chromosome", x = "Chromosome", y = "Frequency") +
  theme_bw()

#c
# Create the bar chart with automatic binning
ggplot(fold_change_genes_info, aes(x = Chromosome, fill = Higher_Expression)) +
  geom_bar(stat = "count") +  # Use stat="count" for frequency
  labs(title = "Distribution of DEGs by Chromosome Segregated by Sample Type", x = "Chromosome", y = "Frequency") +
  scale_fill_manual(values = c("yellow", "red"), labels = c("Normal", "Tumor")) +
  theme_bw() +
  facet_wrap(~ Higher_Expression)


#d.Bar Chart of Upregulated and Downregulated DEGs
# Calculate percentages
upregulated_percentage <- (nrow(subset(fold_change_genes_info, Higher_Expression == "Tumor")) / nrow(fold_change_genes_info)) * 100
downregulated_percentage <- 100 - upregulated_percentage
# Create a dataframe for plotting
percentage_data <- data.frame(
  Gene_Expression = c("Upregulated", "Downregulated"),
  Percentage = c(upregulated_percentage, downregulated_percentage)
)
# Create the bar chart
ggplot(data = percentage_data, aes(x = Gene_Expression, y = Percentage, fill = Gene_Expression)) +
  geom_bar(stat = "identity") +
  labs(title = "Percentage of DEGs Upregulated and Downregulated in Tumor Samples",
       x = "Gene Expression",
       y = "Percentage") +
  scale_fill_manual(values = c("black", "orange")) +
  theme_minimal()

#e.Heatmap of Gene Expression by Sample
# Remove the Probe_ID column
dfx <- df1_renamed[, !names(df1_renamed) %in% "Probe_ID"]
# Transpose the dataframe for ggplot
dfx <- t(dfx)
# Create the heatmap
heatmap(as.matrix(dfx),
        scale = "none",  
        col = c("blue", "green", "yellow", "red"),  
        Rowv = NA, Colv = NA,  
        labRow = "Gene", labCol = "Sample",
        main = "Heatmap of Gene Expression by Sample")

#f.
# Create the clustermap
pheatmap(dfx,
         cluster_cols = TRUE,    
         cluster_rows = TRUE,    
         show_colnames = TRUE,  
         show_rownames = TRUE,   
         color = viridis::viridis(256), 
         main = "Clustermap of Gene Expression by Sample")


# g .write a few sentences of analysis
#scatter plot of fold change indicates the significant variation in gene expression between tumor and normal samples exhibiting foldcahnges greater than 5.
# Histogram of DEG by chromosome indicates the non uniform distribution of DEG's across the chromosome implicating the tumors genomic regions
# the higher percentage of upregulated genes in tumor samples compared to down regulated suggests shifts towards gene expression levels in tumor samples.
# the heat map and cluster map provide the understanding of the correaltion and clustering  patterns of gene expression. this analysis comprehensively emphasize the importance of understanding of gene expression in biology.





