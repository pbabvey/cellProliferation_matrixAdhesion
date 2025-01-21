################################################################################
# Part 0, 1, 2 are data preparation steps. The rest of the sections are independent.

# Table of content
# 0. Install & Load Packages (required)
# 1. Data Prep (required)
# 2. Read Geneset Data (required)
# 3. Quality Control (QC)
# 4. ITGB4- trinary saved GO terms
# 5. Enhanced Volcano Plot
# 6. GO Slim
# 7. TopGO
# 8. COGENT
# 9. Pairwise correlation Analysis
# 14. TSNE plot
################################################################################





# ---------- Part 0: Install & Load Packages ----------


""" 
First, let's install some general packages for matrix transformation, visualizations, etc. including:
 - here: The here package automatically sets the working directory to the root of the project when running the script.
 - Matrix: Provides dense and sparse matrix classes and methods for efficient computation.
 - impute: Implements imputation methods for missing values, often used in microarray data.
 - vctrs: Implements foundational tools for creating and manipulating vector types, ensuring type consistency and extensibility in R.
 - gridExtra: Provides tools to arrange and combine multiple grid-based plots, such as ggplot2 plots, in a flexible layout.
 - ggplot2: A popular package for creating elegant and versatile data visualizations.
 - tidyverse: A collection of R packages for data manipulation, visualization, and analysis.
 - doParallel: for parallel processing and faster matrix transormations
 - parallel: A base R package providing tools for parallel computing to improve performance.
 - ggthemes: Provides additional themes, scales, and color palettes for ggplot2 visualizations.
"""

# here: The here package automatically sets the working directory to the root of the project when running the script.
if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", dependencies = TRUE)
# it sets the  working directory
library(here)

getwd()

# Matrix: Provides dense and sparse matrix classes and methods for efficient computation.
if (!requireNamespace("Matrix", quietly = TRUE))
  install.packages("Matrix", dependencies = TRUE) # 'dependencies = TRUE' ensures all dependent packages are also updated
library(Matrix)

# impute: Implements imputation methods for missing values, often used in microarray data.
if (!requireNamespace("impute", quietly = TRUE))
  BiocManager::install("impute")


# vctrs: Implements foundational tools for creating and manipulating vector types, ensuring type consistency and extensibility in R.
if (!requireNamespace("vctrs", quietly = TRUE))
  install.packages("vctrs")
library(vctrs)

# gridExtra: Provides tools to arrange and combine multiple grid-based plots, such as ggplot2 plots, in a flexible layout.
if (!requireNamespace("gridExtra", quietly = TRUE))
  install.packages("gridExtra")
library(gridExtra)

# ggplot2: A popular package for creating elegant and versatile data visualizations.
if (!requireNamespace("ggplot2", quietly = TRUE)) 
  install.packages("ggplot2")
library(ggplot2)


# tidyverse: A collection of R packages for data manipulation, visualization, and analysis.
if (!requireNamespace("tidyverse", quietly = TRUE))
  install.packages("tidyverse")
library(tidyverse)

### Activating parallel computing for faster computation
if (!requireNamespace("doParallel", quietly = TRUE))
  install.packages("doParallel")
library(doParallel, quietly = TRUE)

# parallel: A base R package providing tools for parallel computing to improve performance.
if (!requireNamespace("parallel", quietly = TRUE)) {
  install.packages("parallel")
}
library(parallel)

# Detect the number of cores, and leave two cores free for system processes
no_cores <- detectCores() - 2  

cl <- makeCluster(no_cores)
registerDoParallel(cl)


# ggthemes: Provides additional themes, scales, and color palettes for ggplot2 visualizations.
if (!requireNamespace("ggthemes", quietly = TRUE))
install.packages("ggthemes")
library("ggthemes")

""" 
Second, we install packages we need to collect data from bio databases as well as bio tools for our anlayses, including:
-  devtools: A set of tools to simplify R package development and installation from GitHub and other sources.
 - BiocManager: Manages installation and updates of Bioconductor packages in R.
 - GO.db: A Bioconductor package providing Gene Ontology data for biological processes, molecular functions, and cellular components.
 - org.Hs.eg.db: A Bioconductor package with comprehensive annotation data for Homo sapiens (human).
 - biomaRt: Enables querying of the BioMart database for biological datasets and annotations.
 - AnnotationDbi: A Bioconductor package for database-like annotation access in R.
 - 
"""

# devtools: A set of tools to simplify R package development and installation from GitHub and other sources.
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

# BiocManager: Manages installation and updates of Bioconductor packages in R.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)


# GO.db: A Bioconductor package providing Gene Ontology data for biological processes, molecular functions, and cellular components.
if (!requireNamespace("GO.db", quietly = TRUE))
  BiocManager::install("GO.db")
library(GO.db)

# org.Hs.eg.db: A Bioconductor package with comprehensive annotation data for Homo sapiens (human).
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

# biomaRt: Enables querying of the BioMart database for biological datasets and annotations.
if (!requireNamespace("biomaRt", quietly = TRUE)) 
  BiocManager::install("biomaRt")
library(biomaRt)

# AnnotationDbi: A Bioconductor package for database-like annotation access in R.
if (!requireNamespace("AnnotationDbi", quietly = TRUE))
  BiocManager::install("AnnotationDbi")
library(AnnotationDbi)


""" 
Finally, let's install packages we need to preprocess and curate data, and draw GSVA, WGCNA, and gene enrichment anlayses, including:
 - preprocessCore: Provides methods for preprocessing high-throughput genomic data.
 - DESeq2: A Bioconductor package for differential gene expression analysis using count data. Use BiocManager to install DESeq2
 - GSEABase: Provides classes and methods for gene set enrichment analysis in R.
 - GSVA: Implements Gene Set Variation Analysis to estimate pathway activity in transcriptomics data.
 - WGCNA: Implements Weighted Gene Co-expression Network Analysis for finding co-expression modules.
 - clusterProfiler: Facilitates statistical analysis and visualization of functional profiles for genes and gene clusters.
 - topGO: A Bioconductor package for performing Gene Ontology enrichment analysis.
 - COGENT: acilitates gene co-expression network analysis and visualization to study gene-gene relationships.
 - CorLevelPlot: Used to create correlation level plots for visualizing and comparing correlations between variables.
"""

# preprocessCore: Provides methods for preprocessing high-throughput genomic data.
if (!requireNamespace("preprocessCore", quietly = TRUE))
  BiocManager::install("preprocessCore")
library(preprocessCore)

# DESeq2: A Bioconductor package for differential gene expression analysis using count data. Use BiocManager to install DESeq2
if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")

# GSEABase: Provides classes and methods for gene set enrichment analysis in R.
if (!requireNamespace("GSEABase", quietly = TRUE))
  BiocManager::install("GSEABase")
library(GSEABase)

# GSVA: Implements Gene Set Variation Analysis to estimate pathway activity in transcriptomics data.
if (!requireNamespace("GSVA", quietly = TRUE))
  BiocManager::install("GSVA")
library(GSVA)

# WGCNA: Implements Weighted Gene Co-expression Network Analysis for finding co-expression modules.
if (!requireNamespace("WGCNA", quietly = TRUE))
  install.packages("WGCNA")
library(WGCNA)
enableWGCNAThreads(nThreads = no_cores) 

# clusterProfiler: Facilitates statistical analysis and visualization of functional profiles for genes and gene clusters.
if (!requireNamespace("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")
library(clusterProfiler)


# topGO: A Bioconductor package for performing Gene Ontology enrichment analysis.
if (!requireNamespace("topGO", quietly = TRUE))
  BiocManager::install("topGO")

# COGENT: acilitates gene co-expression network analysis and visualization to study gene-gene relationships.
if (!requireNamespace("COGENT", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) 
    install.packages("devtools")
  devtools::install_github("lbozhilova/COGENT")
}
library("COGENT")


# CorLevelPlot: Used to create correlation level plots for visualizing and comparing correlations between variables.
#if (!requireNamespace("CorLevelPlot", quietly = TRUE)) {
#  if (!requireNamespace("BiocManager", quietly = TRUE)) 
#    install.packages("BiocManager")
  
#  Install CorLevelPlot using BiocManager
#  BiocManager::install("CorLevelPlot")
#}
# library(CorLevelPlot)

# ---------- Part 1: Data Prep ----------

# Fixing the random seed to make results reproducible
set.seed(101019)



#setwd("/Users/pouria/Dropbox/Pouria & Sadaf's shared stuff!/Bioinformatics")


# Define the data path using here()
data_path <- here("data", "combined_data_for_WGCNA.tsv")

# Read the data
data <- read.delim(data_path, header = TRUE)


# Load the combined data
data_path <- "12_31_24_WGCNA/combined_data_for_WGCNA.tsv"
data <- read.delim(data_path, header = TRUE)
# In many regex implementations, a single backslash \ is enough to escape a character, but in R strings, 
# you need to use two backslashes \\ to represent a single backslash in the regex pattern.
data[,'gene_id'] <- gsub("\\..*$", "", data[,'gene_id'] )
data <- data[!duplicated(data$gene_id), ]
rownames(data) <- data$gene_id

dim(data)

# gene name column removed
data1 <- data[,-1] 


"""
We can remove the following code, I think it's for checking that there is just
one single column as rowname with nonnumerical vals
"""
# Assuming df is your dataframe
row_names <- rownames(data)
# Extract the substring before the dot
substrings <- gsub("\\..*$", "", row_names)

# Count occurrences
counts <- table(substrings)

# View the counts
counts
all(counts == 1)


"""
Total Count Normalization (TCN)
Also known as library size normalization, this method involves scaling each sample
by the total number of reads or counts detected in that sample. 
Its a basic form of normalization suitable for some analyses but does not account
for the variability in gene expression levels across different genes or other systematic biases.
"""
# checking sum to see if tpm normalize is applied and the sum of gene experessions for a cell are equal to 10^6
# it seems it has already been applied
col_sums <- colSums(data1)
row_sums <- rowSums(data1)
col_sums
row_sums

# remove the rows with NA value
na_count_per_row <- apply(data1, 1, function(x) sum(is.na(x)))
rows_with_na <- sum(na_count_per_row > 0)
print(rows_with_na)
data2 <- data1[na_count_per_row == 0, ]


# select 10k genes with the highest variances
gene_variances <- apply(data2, 1, var)

#filtered_variances <- gene_variances [gene_variances >= 0 & gene_variances <= 10000]
filtered_variances <- gene_variances [gene_variances > 0]
any(is.na(filtered_variances))
# Plot the histogram
# Assuming 'data_vector' is your numeric vector of data
breaks <- c(-10, seq(0, 90000, by = 10000), 10^8)  # Define breaks from 0 to 1000 and anything above 1000 goes into the last bin
hist(filtered_variances, breaks = breaks, main = "Histogram with Custom Breaks",
     xlab = "Value Range", ylab = "Frequency", right = TRUE)
# Sort variances in descending order and get the names of the top 10,000 genes
top_genes_indices <- order(gene_variances, decreasing = TRUE)[1:10000]
high_var_genes <- rownames(data2)[top_genes_indices]
data3 <- data2[top_genes_indices,]
data_filtered <- data3



"""
there is another method with hard threshold which I haven't used
"""
#threshold <- 5
#high_var_genes <- rownames(data2)[gene_variances > threshold]
#length(high_var_genes)
#data_filtered <- data3[gene_variances > threshold , ]

dim(data_filtered)
class(data_filtered)
#data_t <- t(data)
#data_t <- as.data.frame(data_t)

"""
I don't know what does it mean but I have removed or commented this part
"""
# Create DESeq2 data set object
# DESeq2 step has been skipped since TPM was already a normalization method, and we have to round the TPM value to feed into Deseq2 which is not a accurate method
# Create a DESeqDataSet object
#dds <- DESeqDataSetFromMatrix(countData = data, colData = colData, design = ~ condition)
# Variance stabilizing transformation
#dds_norm <- vst(dds)
#norm.counts <- assay(dds_norm) %>% t()


# ---------- Part 2: Read Geneset Data ----------
# Define the path to the 'Intermediate Data' folder
#data_path_A <- "12_31_24_WGCNA/GeneSets/GOBP_CELL_CYCLE.v2023.1.Hs.gmt"
#data_path_A <- "12_31_24_WGCNA/GeneSets/small/G1_S_TRANSITION_OF_MITOTIC_CELL_CYCLE.v2023.2.Hs.gmt"
data_path_A <- "12_31_24_WGCNA/GeneSets/KEGG_CELL_CYCLE.v2023.2.Hs.gmt"
#data_path_A <- "12_31_24_WGCNA/GeneSets/PID_TAP63_PATHWAY.v2023.2.Hs.gmt"

# Assuming intermediate_Data is still defined and correct
gene_set_A_df <- read.delim(data_path_A, header = TRUE)
gene_set_A <- colnames(gene_set_A_df)
#gene_set_A <- c('PKMYT1', 'MCM7', 'TTK', 'CDC20', 'BUB1B', 'PTTG1', 'CDC6',
#                'MAD2L1', 'BUB1', 'CCNB2', 'CDC25A', 'CCNB1', 'CDC45', 'ORC1',
#                'PLK1', 'CDC25C', 'CCNA2', 'CCNA1', 'CDK1', 'SFN', 'E2F2')




data_path_B <- "12_31_24_WGCNA/GeneSets/GOBP_CELL_MATRIX_ADHESION.v7.5.1.gmt"
#data_path_B <- "12_31_24_WGCNA/GeneSets/small/WP_INTEGRIN_MEDIATED_CELL_ADHESION.v2023.2.Hs.gmt"
# Assuming intermediate_Data is still defined and correct
gene_set_B_df <- read.delim(data_path_B, header = TRUE)
gene_set_B <- colnames(gene_set_B_df)
#gene_set_B <- c('JUP', 'ITGB4', 'L1CAM', 'ITGB8', 'ITGA11', 'THBS3', 'FN1', 
#                'ITGB3', 'MSLN', 'RHOD', 'THBS1', 'COL16A1', 'JAG1', 'ADAMTS13',
#                'SNED1', 'CDKN2A', 'VEGFA', 'TMEM8B', 'PLAU', 'ADAMTS12')


library("biomaRt")
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_info_A <- getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), 
                     filters = 'external_gene_name', 
                     values = gene_set_A, 
                     mart = ensembl)


ensembl_ids_set_A <- gene_info_A[gene_info_A$external_gene_name %in% gene_set_A, "ensembl_gene_id"]

gene_info_B <- getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), 
                     filters = 'external_gene_name', 
                     values = gene_set_B, 
                     mart = ensembl)
ensembl_ids_set_B <-  gene_info_B[gene_info_B$external_gene_name %in% gene_set_B, "ensembl_gene_id"]

ensembl_ids_set_A_filtered <- ensembl_ids_set_A[ensembl_ids_set_A %in% data_filtered$gene_id]
ensembl_ids_set_B_filtered <- ensembl_ids_set_B[ensembl_ids_set_B %in% data_filtered$gene_id]



"""
I have no idea what this code is doing, don't run it.
"""
checkExpressionDF(data)




# ---------- Part 3: Quality Control (QC) ----------
"""
I am not using this filteration for now
NA values are already removed.
"""


# Store gene_ids in a separate object and remove it from data
#gene_ids <- data$gene_id
#data <- data[,-1]  # Removes the first column, which is 'gene_id' since it has characters not just numbers

# data should be numeric now, but let's make sure by converting it explicitly
data_filtered <- as.data.frame(lapply(data_filtered, as.numeric))

# Check for NA values
if(any(is.na(data_filtered))) {
  stop("Data contains NA values. Please remove or impute these before proceeding.")
}




# the data should now be ready for WGCNA
#WGCNA typically requires data in a "samples-in-rows and genes-in-columns" format 

# Transpose the data so that genes are columns and samples are rows
data_filtered_t <- t(data_filtered)
# Identify outliers among genes and samples
gsg <- goodSamplesGenes(data_filtered_t)
summary(gsg)
gsg$allOK

#If gsg$allOK was false, then:
#To check if the dimension of data and length of gsg are matched
dim(data_filtered_t)
length(gsg$goodGenes)



# Now apply the filtering
data_filtered_2 <- data_filtered_t[gsg$goodSamples, gsg$goodGenes]

# Hierarchical clustering without outlier samples using filtered, transposed data
# Hierarchical clustering to detect outlier samples
# Note: Use data2 as it is already in the "samples-in-rows and genes-in-columns" format
htree <- hclust(dist(data_filtered_2), method = "average")
plot(htree)

# Exclude outlier samples
# Note: Make sure to have a list of samples to be excluded based on the hclust results
samples_to_be_excluded <- c('X9490a06f.43ec.42d6.86ec.4b5480224489', 'adc9bef6.520d.45d5.b36f.a66f8eccc69c')  # Replace with actual sample IDs to exclude
data3 <- data2[!(rownames(data2) %in% samples_to_be_excluded), ]



# Define the full path for the file to be saved
file_path <- file.path(intermediate_Data, "data3.rds")

# Save the processed data
saveRDS(data3, file = file_path)





# ---------- Part 4. ITGB4- trinary saved GO terms ----------
# ITGB4 Ensemble ID: 'ENSG00000132470'
rownames(data_filtered) <- sub("\\..*$", "", rownames(data_filtered))
data_filtered['ENSG00000132470',]

data_filtered <- as.data.frame(data_filtered)
my_vector <- as.numeric(data_filtered['ENSG00000132470',])

# Create a histogram
hist(my_vector, main = "Distribution Plot of ITGB4 Expression Values", xlab = "Values", ylab = "Frequency", col = "blue")

# If you want a density plot instead, you can use:
plot(density(my_vector), main = "Density Plot of ITGB4 Expression Values", xlab = "Values", ylab = "Density", col = "red")


high_threshold <- quantile(my_vector , probs = 0.80)
low_threshold <- quantile(my_vector , probs = 0.20)

high_ITGB4 <- data_filtered[, data_filtered['ENSG00000132470', ] >= high_threshold]
high_ITGB4_means <- rowMeans(high_ITGB4)
length(high_ITGB4_means)

mid_ITGB4 <- data_filtered[, high_threshold > data_filtered['ENSG00000132470', ] & data_filtered['ENSG00000132470', ] > low_threshold]
mid_ITGB4_means <- rowMeans(mid_ITGB4)
length(mid_ITGB4_means)

low_ITGB4 <- data_filtered[,  low_threshold > data_filtered['ENSG00000132470', ] ]
low_ITGB4_means <- rowMeans(low_ITGB4)
length(low_ITGB4_means)

"""
Approach one: p-values
"""
# Assuming your dataframe is named 'data_frame'
# Adjust the column indices as per your actual dataframe

# Initialize a matrix or dataframe to store the results
results <- data.frame(gene = rownames(data_filtered), p_value_higher = NA, p_value_lower = NA)

# Loop through each gene
for (i in 1:nrow(data_filtered)) {
  # Extract gene expression data for this gene
  #low_expression <- data_frame[i, 1:80]  # columns for low ITGB4 expression
  #high_expression <- data_frame[i, 81:160]  # columns for high ITGB4 expression
  
  if (i %% 100 == 0) {
    print(i)
  }
  
  # Perform t-test to compare high vs. low
  test_result <- t.test(high_ITGB4[i,] , low_ITGB4[i,], alternative = "greater")
  
  # Store p-value for the test of higher expression in high ITGB4
  results$p_value_higher[i] <- test_result$p.value
  
  # Perform t-test for lower expression (reverse groups)
  test_result_lower <- t.test(low_ITGB4[i,], high_ITGB4[i,], alternative = "greater")
  
  # Store p-value for the test of lower expression in high ITGB4
  results$p_value_lower[i] <- test_result_lower$p.value
}

# Adjust p-values for multiple testing if necessary, e.g., using Benjamini-Hochberg method
results$p_value_higher_adjusted <- p.adjust(results$p_value_higher, method = "BH")
results$p_value_lower_adjusted <- p.adjust(results$p_value_lower, method = "BH")

# Print results
head(results)


# Sort the dataframe by 'score' in descending order
results_higher_adjusted <- results[order(results$p_value_higher_adjusted, decreasing = FALSE), ]

# View the first few rows of the dataframe to confirm
head(results_higher_adjusted )

upregulated_genes <- head(results_higher_adjusted$gene, 1000)

library(clusterProfiler)
library(org.Hs.eg.db)
ego <- enrichGO(gene = upregulated_genes,
                OrgDb = "org.Hs.eg.db", # Choose the appropriate OrgDb for your organism
                keyType = "ENSEMBL", # Ensure this matches your gene ID type
                ont = "BP", # For Biological Processes; use "MF" for Molecular Functions, "CC" for Cellular Components
                pAdjustMethod = "BH", # Method for adjusting p values
                qvalueCutoff = 0.05, # Set your significance level
                readable = TRUE) # To convert IDs to readable gene names

# View results
summary(ego)

# Convert the ego object to a data frame
ego_df <- as.data.frame(ego)

write.csv(ego_df, file = "12_31_24_WGCNA/pos_corr_ITBG4_ego_results.csv", row.names = TRUE)



# Sort the dataframe by 'score' in descending order
results_lower_adjusted <- results[order(results$p_value_lower_adjusted, decreasing = FALSE), ]

# View the first few rows of the dataframe to confirm
head(results_lower_adjusted )

downregulated_genes <- head(results_lower_adjusted$gene, 1000)

ego <- enrichGO(gene = downregulated_genes,
                OrgDb = "org.Hs.eg.db", # Choose the appropriate OrgDb for your organism
                keyType = "ENSEMBL", # Ensure this matches your gene ID type
                ont = "BP", # For Biological Processes; use "MF" for Molecular Functions, "CC" for Cellular Components
                pAdjustMethod = "BH", # Method for adjusting p values
                qvalueCutoff = 0.05, # Set your significance level
                readable = TRUE) # To convert IDs to readable gene names

# View results
summary(ego)

# Convert the ego object to a data frame
ego_df <- as.data.frame(ego)

write.csv(ego_df, file = "12_31_24_WGCNA/neg_corr_ITBG4_ego_results.csv", row.names = FALSE)




# Select a subset of results for plotting
selected_terms <- ego[ego$p.adjust < 0.005 & ego$Count > 10, ]




"""
Approach two: using mean of gene expressions per gene
"""
# Step 1: Create the data frame
ITGB4_data_frame <- data.frame(low = low_ITGB4_means, mid = mid_ITGB4_means, high = high_ITGB4_means)

# Step 2: Calculate variance for each row
ITGB4_data_frame$std <- apply(ITGB4_data_frame , 1, sd)

# Step 3: Calculate the difference between 'high' and 'low'
ITGB4_data_frame$difference <- ITGB4_data_frame$high - ITGB4_data_frame$low

# Step 4: corr with ITGB4
ITGB4_data_frame$correlated <- ITGB4_data_frame$difference / (ITGB4_data_frame$std+10^-15)

# View the first few rows of the dataframe to confirm
head(ITGB4_data_frame)



# Sort the dataframe by 'score' in descending order
sorted_ITGB4_data_frame_desc <- ITGB4_data_frame[order(ITGB4_data_frame$correlated, decreasing = TRUE), ]

# View the first few rows of the dataframe to confirm
head(sorted_ITGB4_data_frame_desc )

upregulated_genes <- head(rownames(sorted_ITGB4_data_frame_desc), 1000)

ego <- enrichGO(gene = upregulated_genes,
                OrgDb = "org.Hs.eg.db", # Choose the appropriate OrgDb for your organism
                keyType = "ENSEMBL", # Ensure this matches your gene ID type
                ont = "BP", # For Biological Processes; use "MF" for Molecular Functions, "CC" for Cellular Components
                pAdjustMethod = "BH", # Method for adjusting p values
                qvalueCutoff = 0.05, # Set your significance level
                readable = TRUE) # To convert IDs to readable gene names

# View results
summary(ego)

# Convert the ego object to a data frame
ego_df_pos <- as.data.frame(ego)

# Save the data frame to a text file
#write.table(ego_df, file = "12_31_24_WGCNA/ego_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(ego_df_pos, file = "12_31_24_WGCNA/pos_corr_ITBG4_ego_results.csv", row.names = FALSE)
# get the data

downregulated_genes <- tail(rownames(sorted_ITGB4_data_frame_desc), 1000)

ego <- enrichGO(gene = downregulated_genes,
                OrgDb = "org.Hs.eg.db", # Choose the appropriate OrgDb for your organism
                keyType = "ENSEMBL", # Ensure this matches your gene ID type
                ont = "BP", # For Biological Processes; use "MF" for Molecular Functions, "CC" for Cellular Components
                pAdjustMethod = "BH", # Method for adjusting p values
                qvalueCutoff = 0.05, # Set your significance level
                readable = TRUE) # To convert IDs to readable gene names

ego_df_neg <- as.data.frame(ego)

# View results
summary(ego)


# Save the data frame to a text file
#write.table(ego_df, file = "12_31_24_WGCNA/ego_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(ego_df_neg, file = "12_31_24_WGCNA/neg_corr_ITBG4_ego_results.csv", row.names = FALSE)


ego_df_pos_20 <- head(ego_df_pos, 20)

pos_process_list <- c("cell-matrix adhesion", "focal adhesion assembly", "regulation of cell shape",
                      "regulation of cell morphogenesis", "protein localization to plasma membrane",
                      "extracellular matrix organization", "vacuolar transport")

selected_rows <- ego_df_pos[ego_df_pos$Description %in% pos_process_list, ]
#selected_rows$minus_log10_p_value <- -log10(selected_rows$p.adjust)
selected_rows$minus_log10_p_value <- -log10(selected_rows$pvalue)
# Creating the plot
plot <- ggplot(selected_rows, aes(x = reorder(Description, minus_log10_p_value), y = minus_log10_p_value)) +
  geom_bar(stat = "identity", fill = "#ffb3b3", color= "#ff8080", width = 0.7) +  # Bars filled in red
  coord_flip() +  # Flip coordinates to make it horizontal
  geom_text(aes(label = Description), 
            hjust = 1, color = "black", size = 5, nudge_y = -1)+ 
  labs(x = "GO terms", y = expression(-log[10](italic("p-value"))), 
       title = "Gene Set Upregulation Linked to Elevated ITGB4") +
  theme_minimal() +  # Minimal theme
  theme(
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14), 
    axis.title.y = element_blank(),       # Remove y-axis label
    #axis.text.y = element_text(angle = 45, hjust = 1, size = 10),  # Adjust inclination and size of text
    axis.ticks.y = element_blank(), 
    axis.text.y = element_blank(),        # Completely remove y-axis text
    #axis.title.y = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 20, hjust = 0.5)
  ) 

# Print the plot
print(plot)




ego_df_neg_20 <- head(ego_df_neg, 20)

neg_process_list <- c( "nucleoside triphosphate metabolic process", "DNA replication", 
                      "chromosome segregation", "cell cycle checkpoint signaling", "sister chromatid segregation", 
                      "negative regulation of cell cycle phase transition", "nuclear division") #"ribosome biogenesis",


selected_rows <- ego_df_neg[ego_df_neg$Description %in% neg_process_list, ]
#selected_rows$minus_log10_p_value <- -log10(selected_rows$p.adjust)
selected_rows$minus_log10_p_value <- -log10(selected_rows$pvalue)

# Creating the plot
plot <- ggplot(selected_rows, aes(x = reorder(Description, minus_log10_p_value), y = minus_log10_p_value)) +
  geom_bar(stat = "identity", fill = "#b3d1ff", color = "#1a75ff", width = 0.7) +
  coord_flip() +  # Flip coordinates to make it horizontal
  geom_text(aes(label = Description), 
            hjust = 1, color = "black", size = 5, nudge_y = -0.25)+ 
  labs(x = "GO terms", y = expression(-log[10](italic("p-value"))), 
       title = "Gene Set Downregulation Linked to Elevated ITGB4") +
  theme_minimal() +  # Minimal theme
  theme(
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14), 
    axis.title.y = element_blank(),       # Remove y-axis label
    #axis.text.y = element_text(angle = 45, hjust = 1, size = 10),  # Adjust inclination and size of text
    axis.ticks.y = element_blank(), 
    axis.text.y = element_blank(),        # Completely remove y-axis text
    #axis.title.y = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 20, hjust = 0.5)
  ) 

print(plot)
  
  
library(ggplot2)

# Sample data
data <- data.frame(
  GeneSet = c("blood coagulation", 
              "endodermal cell differentiation", 
              "extracellular matrix organization", 
              "leukocyte cell-cell adhesion", 
              "extracellular matrix disassembly"),
  Value = c(7, 5, 6, 3, 4)
)

# Plot with x-axis on top
plot <- ggplot(data, aes(x = Value, y = reorder(GeneSet, Value))) +
  geom_bar(stat = "identity", fill = "grey70", color = "black") +
  labs(title = "MAM Cluster 8", x = "-log₁₀ adj.p-value", y = "") +
  scale_x_continuous(position = "top") + # Place x-axis on top
  theme_minimal() +
  theme(
    axis.title.y = element_blank(), # Remove y-axis label if desired
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5) # Center the title
  )

print(plot)


write.csv(ego_df_pos, file = "12_31_24_WGCNA/neg_corr_ITBG4_bar_plot.csv", row.names = FALSE)







# ---------- Part 5. Enhanced Volcano plot ----------
"""
You need to run the previous section for this part
"""
if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")

if (!requireNamespace("EnhancedVolcano", quietly = TRUE))
  BiocManager::install("EnhancedVolcano")

library(DESeq2)
library(EnhancedVolcano)

all(rownames(high_ITGB4) == rownames(low_ITGB4)) 


combined_df <- cbind(high_ITGB4, low_ITGB4)
combined_df[] <- round(combined_df)
# Ensure column names are unique and descriptive
colnames(combined_df) <- c(paste0("Control_", seq_len(ncol(high_ITGB4))),
                           paste0("Treatment_", seq_len(ncol(low_ITGB4))))

# Create a sample information data frame
sample_info <- data.frame(
  condition = rep(c("Treatment", "Control" ), each = 86),
  row.names = colnames(combined_df)
)

# Create a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = combined_df,
                              colData = sample_info,
                              design = ~ condition)

dds <- DESeq(dds)



results_deseq <- results(dds)
# Optionally, shrink log2 fold changes for better visualization and interpretation
results_deseq <- lfcShrink(dds, coef = 2, type = "apeglm")

# View results
head(results_deseq)



library("biomaRt")
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_names <- getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), 
                     filters = 'ensembl_gene_id', 
                     values = rownames(results_deseq), 
                     mart = ensembl)
filtered_results_deseq <-  results_deseq[rownames(results_deseq) %in% gene_names$ensembl_gene_id, ]




sum(duplicated(gene_names$external_gene_name))


# If there are duplicates, you might need to create unique names by appending a counter or the ID
gene_names$unique_gene_name <- make.unique(as.character(gene_names$external_gene_name))


# Create a named vector for easy lookup
name_lookup <- setNames(gene_names$unique_gene_name, gene_names$ensembl_gene_id)

# Replace rownames of filtered_df1
rownames(filtered_results_deseq) <- name_lookup[rownames(filtered_results_deseq)]
sum(is.na(rownames(filtered_results_deseq)))


rows_to_remove <- c('ITGB4', 'PAGE2')
filtered_results_deseq <- filtered_results_deseq[!rownames(filtered_results_deseq) %in% rows_to_remove, ]

EnhancedVolcano(filtered_results_deseq,
                lab = rownames(filtered_results_deseq),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlab = 'Log2 Fold Change',
                ylab = '-Log10 p-value',
                title = 'GO Terms Enrichment Volcano Plot')



# ---------- Part 6. GO Slim ----------
"""
You need part ITGB4 for this to give you upregulated and downregulated genes
"""



library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
# Assuming 'ego' is your enrichGO output
ego <- enrichGO(gene = downregulated_genes,
                OrgDb = org.Hs.eg.db,
                keyType = "ENSEMBL",
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

ego_df <- as.data.frame(ego)

# Assuming 'ego' is your enrichGO result and it has GO IDs in the row names
go_ids <- ego_df$geneID



# Define URLs for the required files
go_url <- "http://purl.obolibrary.org/obo/go/go-basic.obo"
goslim_url <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"

# Define paths for downloading
go_path <- "12_31_24_WGCNA/go-basic.obo"
goslim_path <- "12_31_24_WGCNA/goslim_generic.obo"

# Download files
download.file(go_url, go_path)
download.file(goslim_url, goslim_path)

# Example GO terms
go_terms <- ego_df$ID

# Write these to a file
write(go_terms, "12_31_24_WGCNA/go_terms.txt")

# Assuming you have a vector of GO IDs
go_ids <- c("GO:0072332")




# Use bitr function to map from GO to GO Slim
slim_mapper <- bitr(go_ids, fromType = "GO", toType = "GOSLIM_GO", OrgDb = "org.Hs.eg.db")

# Check the results
head(slim_mapper)








data_path_C <- "12_31_24_WGCNA/GeneSets/NABA_CORE_MATRISOME.v2023.2.Hs.gmt"

# Assuming intermediate_Data is still defined and correct
gene_set_C_df <- read.delim(data_path_C, header = TRUE)
gene_set_C <- colnames(gene_set_C_df)


#BiocManager::install("biomaRt")
library("biomaRt")
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_info_C <- getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), 
                     filters = 'external_gene_name', 
                     values = gene_set_C, 
                     mart = ensembl)


ensembl_ids_set_C <- gene_info_C[gene_info_C$external_gene_name %in% gene_set_C, "ensembl_gene_id"]


ensembl_ids_set_C_filtered <- ensembl_ids_set_C[ensembl_ids_set_C %in% rownames(data_filtered_original)]



# Select only rows corresponding to the selected genes
high_df <- high_ITBG4[ensembl_ids_set_C_filtered, ]


# Select only rows corresponding to the selected genes
mid_df <- mid_ITBG4[ensembl_ids_set_C_filtered, ]


# Select only rows corresponding to the selected genes
low_df <- low_ITBG4[ensembl_ids_set_C_filtered, ]





# Calculate the average for each gene across all selected columns (patients)
high_gene_averages <- rowMeans(high_df , na.rm = TRUE)
# The resulting vector contains the average for each selected gene
high_gene_averages

# Calculate the average for each gene across all selected columns (patients)
mid_gene_averages <- rowMeans(mid_df , na.rm = TRUE)
# The resulting vector contains the average for each selected gene
mid_gene_averages

# Calculate the average for each gene across all selected columns (patients)
low_gene_averages <- rowMeans(low_df , na.rm = TRUE)
# The resulting vector contains the average for each selected gene
low_gene_averages

# Merge into a dataframe
combined_gene_averages <- data.frame(high_ITGB = high_gene_averages, mid_ITGB = mid_gene_averages,
                                     low_ITGB = low_gene_averages)

# Display the resulting dataframe
dim(combined_gene_averages)

rownames(combined_gene_averages) <- convertEnsemblToGeneNames(rownames(combined_gene_averages))


# Set the color scale extremes based on your data range
color_min <- -0.61
color_max <- 0.61

# Create a color palette
blue_to_red <- colorRampPalette(c("blue", "white", "red"))(n = 999)

# Define the breaks to cover the range of your data
break_points <- seq(color_min, color_max, length.out = 1000)

#install.packages("gplots")
library(gplots)
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
} 
library(pheatmap)

# Visualize the heatmap with defined color scaling
pheatmap(as.matrix(combined_gene_averages), 
        height = 5,
        width = 20,
        color = blue_to_red,
          #col = blue_to_red, 
          #Rowv = NA, 
          #Colv = NA,
          #breaks = break_points,
          #trace = 'none', 
          main = "Impact of ITGB4 Expression on Matrisome",
          #cex.main = 0.6,
          xlab = "",
          ylab = "NABA_CORE_MATRISOME",
          #cexRow = 0.3,  # Adjust row label size
          #cexCol = 0.8,  # Adjust column label size
          fontsize_col = 12,
        fontsize_row = 3,
        cluster_cols = FALSE, 
          labRow = rownames(combined_gene_averages),  # Ensure row labels are correct
          labCol = colnames(combined_gene_averages),
          #margins=c(5,5), 
          #lhei=c(2,40), 
          #lwid=c(2,1)
        )



# ---------- Part 7. TopGO ----------
library(topGO)
library(ALL)
data(ALL)
data(geneList)

topDiffGenes(geneList)

sum(topDiffGenes(geneList))



# ---------- Part 8. COGENT ----------

buildPearson <- function(df, quant=0.90){
  check <- checkExpressionDF(df)
  A <- cor(t(df[,colnames(df)!="Name"]), use="pairwise.complete.obs", method="pearson")
  threshold <- quantile(A[upper.tri(A)], quant, na.rm=TRUE)
  A <- 1*(A>=threshold); diag(A) <- 0
  colnames(A) <- rownames(A) <- df$Name
  return(A)
}

'''
buildKendall<- function(df, quant=0.90){
  check <- checkExpressionDF(df)
  A <- cor(t(df[,colnames(df)!="Name"]), use="pairwise.complete.obs", method="kendall")
  threshold <- quantile(A[upper.tri(A)], quant, na.rm=TRUE)
  A <- 1*(A>=threshold); diag(A) <- 0
  colnames(A) <- rownames(A) <- df$Name
  return(A)

'''


#rm(data3)
#gc()
#cor <- WGCNA::cor
pearsonA <- buildPearson(data_filtered)

#kendallA <- buildKendall(data_filtered)

#PKcomparison <- getEdgeSimilarity(list(pearsonA, kendallA), align=FALSE)
#PKcomparison$nodeCount





adj_matrix <- pearsonA


# Extract the relevant parts of the adjacency matrix for both gene sets
adj_A <- adj_matrix[ensembl_ids_set_A_filtered, ensembl_ids_set_A_filtered]
adj_B <- adj_matrix[ensembl_ids_set_B_filtered, ensembl_ids_set_B_filtered]

# Calculate mean co-expression within each gene set
mean_coexp_A <- mean(adj_A[upper.tri(adj_A)])
mean_coexp_B <- mean(adj_B[upper.tri(adj_B)])

# Calculate co-expression between gene sets
adj_AB <- adj_matrix[ensembl_ids_set_A_filtered, ensembl_ids_set_B_filtered]
mean_coexp_AB <- mean(adj_AB)

# Output results
mean_coexp_A
mean_coexp_B
mean_coexp_AB


# ---------- Part 9. Pairwise correlation Analysis ----------
rownames(data_filtered) <- data_filtered_original$gene_id

common_elements <- intersect(ensembl_ids_set_A_filtered , ensembl_ids_set_B_filtered)
ensembl_ids_set_B_filtered <- setdiff(ensembl_ids_set_B_filtered, common_elements)

data_A <- data_filtered [ensembl_ids_set_A_filtered , ]
data_B <- data_filtered [ensembl_ids_set_B_filtered , ]

itgb4 <- 'ENSG00000132470'
tp63 <- 'ENSG00000073282'
cdkn2a <- 'ENSG00000147889'
cdk4 <- 'ENSG00000135446'
ccnd1 <- 'ENSG00000110092'
p53 <- 'ENSG00000141510'
mdm2 <- 'ENSG00000135679'
my_genes <- c(itgb4, tp63, cdkn2a, cdk4, ccnd1, p53, mdm2)
my_genes <- c('ENSG00000181631', 'ENSG00000163623')
my_data <- data_filtered [my_genes , ]
my_data
my_corr <-cor(t(my_data), t(my_data))
my_corr
# Compute pairwise correlations
pairwise_correlations <- cor(t(data_A), t(data_B))

# Create a red to blue color palette
blue_to_red <- colorRampPalette(c("grey", "black"))(n = 299)

# Visualize the correlation matrix
heatmap(as.matrix(pairwise_correlations), col = blue_to_red )


# Calculate pairwise correlations
cor_matrix <- cor(t(data_A), t(data_B))

# Define thresholds for high positive and negative correlations
high_pos_threshold <- 0.0 # Example threshold for high positive correlation
high_neg_threshold <- -0.0 # Example threshold for high negative correlation

# Find gene pairs with high positive correlation
high_pos_cor_genes <- which(cor_matrix > high_pos_threshold, arr.ind = TRUE)

# Find gene pairs with high negative correlation
high_neg_cor_genes <- which(cor_matrix < high_neg_threshold, arr.ind = TRUE)

# Extract gene names for high positive correlations
genes_high_pos_cor_A <- rownames(cor_matrix)[high_pos_cor_genes[, "row"]]
genes_high_pos_cor_B <- colnames(cor_matrix)[high_pos_cor_genes[, "col"]]
genes_high_pos_cor_A 
genes_high_pos_cor_B 
# Extract gene names for high negative correlations
genes_high_neg_cor_A <- rownames(cor_matrix)[high_neg_cor_genes[, "row"]]
genes_high_neg_cor_B <- colnames(cor_matrix)[high_neg_cor_genes[, "col"]]
genes_high_neg_cor_A 
genes_high_neg_cor_B 


library(biomaRt)
# Function to convert Ensembl IDs to External Gene Names, maintaining structure
convertEnsemblToGeneNames <- function(ensembl_ids) {
  # Connect to Ensembl database
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Query for External Gene Names
  gene_info <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
                     filters = 'ensembl_gene_id', 
                     values = ensembl_ids, 
                     mart = ensembl)
  
  # Create a named vector for easy lookup
  names(gene_info$external_gene_name) <- gene_info$ensembl_gene_id
  
  # Replace each Ensembl ID with the corresponding External Gene Name
  gene_names <- sapply(ensembl_ids, function(id) gene_info$external_gene_name[which(gene_info$ensembl_gene_id == id)])
  
  return(gene_names)
}

convertEnsemblToGeneNames(c('ENSG00000181631', 'ENSG00000163623'))
# Example Usage
#ensembl_ids <- genes_high_pos_cor_A  # Replace with your Ensembl IDs
gene_names_high_pos_A <- convertEnsemblToGeneNames(genes_high_pos_cor_A)
print(table(gene_names_high_pos_A))

gene_names_high_pos_B <- convertEnsemblToGeneNames(genes_high_pos_cor_B)
print(table(gene_names_high_pos_B))

df_high_corr <- data.frame(gene_A = gene_names_high_A, gene_B = gene_names_high_B)
selected_rows <- df_high_corr [df_high_corr$gene_B  == 'ITGB4', ]

gene_names_high_neg_A <- convertEnsemblToGeneNames(genes_high_neg_cor_A)
print(table(gene_names_high_neg_A))

gene_names_high_neg_B <- convertEnsemblToGeneNames(genes_high_neg_cor_B)
print(table(gene_names_high_neg_B))

df_high_neg_corr <- data.frame(gene_A = gene_names_high_neg_A, gene_B = gene_names_high_neg_B)
selected_rows <- df_high_neg_corr [df_high_neg_corr$gene_B  == 'ITGB4', ]
selected_rows



# Threshold for selecting genes
threshold <- 0.5  # Adjust this according to your needs


# Subset the gene sets based on the selection
#selected_genes_A <- gene_set_A[high_cor_genes_A]
#selected_genes_B <- gene_set_B[high_cor_genes_B]

# Output the selected genes
#selected_genes_A
#selected_genes_B

high_pos_threshold <- 0.0 # Example threshold for high positive correlation
high_neg_threshold <- -0.0 # Example threshold for high negative correlation

# Select genes from A with high positive or negative correlation with any gene in B
high_cor_genes_A <- apply(cor_matrix, 1, function(x) any(x >= high_pos_threshold  | x <= high_neg_threshold ))

# Select genes from B with high positive or negative correlation with any gene in A
high_cor_genes_B <- apply(cor_matrix, 2, function(x) any(x >= high_pos_threshold  | x <= high_neg_threshold ))

sub_cor_table <- cor_matrix[high_cor_genes_A , high_cor_genes_B ]
dim(sub_cor_table)

rownames(sub_cor_table) <- convertEnsemblToGeneNames(rownames(sub_cor_table))
colnames(sub_cor_table) <- convertEnsemblToGeneNames(colnames(sub_cor_table))

# Set the color scale extremes based on your data range
color_min <- -0.5
color_max <- 0.5

# Create a color palette
blue_to_red <- colorRampPalette(c("blue", "white", "red"))(n = 999)

# Define the breaks to cover the range of your data
break_points <- seq(color_min, color_max, length.out = 1000)

#install.packages("gplots")
library(gplots)

# Visualize the heatmap with defined color scaling
heatmap.2(as.matrix(sub_cor_table), 
        col = blue_to_red, 
        breaks = break_points,
        trace = 'none', 
        main = "Heatmap with Controlled Color Scaling",
        xlab = "Gene Set B",
        ylab = "Gene Set A",
        cexRow = 0.9,  # Adjust row label size
        cexCol = 0.9,  # Adjust column label size
        labRow = rownames(sub_cor_table),  # Ensure row labels are correct
        labCol = colnames(sub_cor_table))


# Plot the histogram
all_values <- as.vector(as.matrix(sub_cor_table))

hist(all_values, main = "Histogram of All Values in the Dataframe", xlab = "Values", ylab = "Frequency")


###############################################

data_A <- data_filtered [ensembl_ids_set_A_filtered , ]
data_B <- data_filtered [ensembl_ids_set_B_filtered , ]
data_B <- data_B[names(data_A)] # Reordering columns of df2 to match df1
combined_df <- rbind(data_A, data_B)
# Calculate pairwise correlations
cor_matrix <- cor(t(combined_df ), t(combined_df ))

combined_gene_set <- c(high_cor_genes_A , high_cor_genes_B )

# Define the point where gene set A ends and gene set B begins
division_point <- length(high_cor_genes_A)

dim(cor_matrix)
sub_cor_table <- cor_matrix
rownames(sub_cor_table) <- convertEnsemblToGeneNames(rownames(sub_cor_table))
colnames(sub_cor_table) <- convertEnsemblToGeneNames(colnames(sub_cor_table))

# Set the color scale extremes based on your data range
color_min <- -1
color_max <- 1
# Create a color palette
blue_to_red <- colorRampPalette(c("blue", "white", "red"))(n = 999)

# Define the breaks to cover the range of your data
break_points <- seq(color_min, color_max, length.out = 1000)
# Plot the heatmap
heatmap.2(as.matrix(sub_cor_table), 
          breaks = break_points,
          trace = "none",          # Remove trace lines
          Rowv = NA,               # Disable row dendrogram
          Colv = NA,               # Disable column dendrogram
          col = blue_to_red,
          dendrogram = "none",     # No dendrogram
          sepwidth = c(0.05, 0.05), # Width of the separation line
          sepcolor = "black",      # Color of the separation line
          cexRow = 0.6,  # Adjust row label size
          cexCol = 0.6,  # Adjust column label size
          srtRow = 0,              # Text rotation for row names
          #srtCol = 45,             # Text rotation for column names
          margins = c(5, 10))      # Margins around the plot

# Create a vector to mark the separation line
abline_col <- division_point
abline_row <- division_point
# Add lines to indicate the separation
abline(h = abline_row, col = "black", lwd = 2)
#abline(v = abline_col, col = "black", lwd = 2)

# ------------------------ GSVA --------------------

#BiocManager::install("GSVA")
library(GSVA)
ensembl_ids_set_B_filtered  <- c('ENSG00000132470')
geneSets <- list(GeneSetA = ensembl_ids_set_A_filtered , GeneSetB = ensembl_ids_set_B_filtered )

# Perform GSVA
gsvaResults <- gsva(data_filtered, geneSets, method="gsva")

correlationResults <- cor(gsvaResults["GeneSetA",], gsvaResults["GeneSetB",], use="complete.obs")

#plot(gsvaResults["GeneSetA",], gsvaResults["GeneSetB",], xlab="Gene Set A Activity", ylab="Gene Set B Activity", main="Correlation between Gene Sets A and B")

# Calculate correlation
correlation <- round(cor(gsvaResults["GeneSetA",], gsvaResults["GeneSetB",]), 3)

library(ggplot2)
# Plot
#GeneSet_A_scores <- gsvaResults["GeneSetA", ]
#GeneSet_B_scores <- gsvaResults["GeneSetB", ]

# Create a dataframe for plotting
#plot_df <- data.frame(GeneSet_A = GeneSet_A_scores, GeneSet_B = GeneSet_B_scores)
df_transposed <- as.data.frame(t(gsvaResults))

# Calculate extended limits for x and y axes
x_range <- range(df_transposed$GeneSetA, na.rm = TRUE)
y_range <- range(df_transposed$GeneSetB, na.rm = TRUE)
x_lim <- c(x_range[1] - 0.05, x_range[2] + 0.05)
y_lim <- c(y_range[1] - 0.05, y_range[2] + 0.05)

# Fit a linear model
model <- lm(GeneSetB ~ GeneSetA, data = df_transposed)

# Extract R-squared and p-value
r_squared <- summary(model)$r.squared
p_value <- summary(model)$coefficients[2, 4] # p-value for the slope coefficient

# Print the results
print(paste("R-squared:", r_squared))
print(paste("p-value:", p_value))

ggplot(df_transposed, aes(x = GeneSetA, y = GeneSetB)) +
  geom_point() +
  geom_smooth(method = "lm", color = "grey", fill = "#1f78b4") +
  annotate("text", x = max(df_transposed$GeneSetA)+0.03, y = max(df_transposed$GeneSetB)+0.05, 
           label = paste("Correlation:", correlation, "R-Squared:", r_squared), hjust = 1, vjust = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14)) +
  labs(x = "Gene Set A GSVA Score: PID_DELTA_NP63_PATHWAY", y = "Gene Set B GSVA Score: ITGB4")+
  xlim(x_lim[1], x_lim[2]) + # Set extended x-axis limits
  ylim(y_lim[1], y_lim[2])   # Set extended y-axis limits




# --------------------- Heatmap --------------------
# Assuming adj_matrix is your adjacency matrix
library(ggplot2)
install.packages("reshape2")
library(reshape2)
melted_adj_matrix <- melt(as.matrix(adj_A))
ggplot(melted_adj_matrix, aes(Var1, Var2, fill=value)) + 
  geom_tile() +
  scale_fill_gradient2(low="blue", high="red", mid="white", 
                       midpoint=median(melted_adj_matrix$value, na.rm=TRUE)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x="Gene 1", y="Gene 2", fill="Coexpression")





# ------------------------------- PCA plot --------------------

# Load necessary library
library(ggplot2)

# Your gene expression matrix (genes x patients)
# expression_data <- ...

# Vector of selected genes
# selected_genes <- c("Gene1", "Gene2", "Gene3", ...)

# Subset the expression matrix to include only selected genes
data_filtered_2 <- data2[sapply(data2 , is.numeric)]
subsampled_data <- t(data_filtered_2)

normalizeData <- function(data) {
  # Apply function to each row (MARGIN = 1)
  normalized_data <- t(apply(data, MARGIN = 1, function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }))
  return(normalized_data)
}

# Assuming your original data is stored in 'data'
subsampled_data <- normalizeData(subsampled_data)

# Perform PCA
pca_result <- prcomp(t(subsampled_data), scale. = TRUE)
pca_result <- prcomp(subsampled_data, scale. = TRUE)

# Create a data frame for plotting
pca_df <- data.frame(PC1 = pca_result$x[,1], PC2 = pca_result$x[,2], Patient = colnames(subsampled_data))
pca_df <- data.frame(PC1 = pca_result$x[,1], PC2 = pca_result$x[,2], Patient = rownames(subsampled_data))

# Plotting PCA
ggplot(pca_df, aes(x = PC1, y = PC2, label = Patient)) +
  geom_point() +
  theme_minimal() +
  ggtitle("PCA of Selected Genes") +
  xlab("Principal Component 1") +
  ylab("Principal Component 2") #+
  #geom_text(aes(label=Patient), vjust=2, hjust=0.5, check_overlap = TRUE, size=3)

#install.packages("Rtsne")
# ------------------------------- TSNE plot --------------------

library(Rtsne)

#data_filtered[data_filtered$gene_id == 'ENSG00000132470', ]
data_filtered['ENSG00000132470',]

# Assuming your data is in a dataframe called 'genes_data'
# where genes are rows and patients are columns.

# First, transpose the data so that genes are columns and patients are rows.

na_count_per_row <- apply(data_filtered, 1, function(x) sum(is.na(x)))
rows_with_na <- sum(na_count_per_row > 0)
print(rows_with_na)
# remove the rows with NA value
data_filtered <- data_filtered[na_count_per_row == 0, ]

data_filtered_2 <- data_filtered[sapply(data_filtered , is.numeric)]
subsampled_data <- data_filtered_2

data_filtered_unique <- data_filtered[!duplicated(data_filtered), ]


# this code normalize the gene expression with (X-Mu)/Sigma
normalizeRows <- function(df) {
  # Function to normalize a single row
  normalizeRow <- function(row) {
    if (var(row) == 0) {
      return(row)  # Return the row unchanged if variance is 0
    } else {
      return((row - mean(row)) / sqrt(var(row)))
    }
  }
  
  # Apply the normalization function to each row
  normalized_df <- t(apply(df, 1, normalizeRow))
  
  # Setting the row and column names of the normalized dataframe to be the same as the original dataframe
  rownames(normalized_df) <- rownames(df)
  colnames(normalized_df) <- colnames(df)
  
  return(normalized_df)
}


# Assuming your original data is stored in 'data'

data_filtered_unique_normalized <- normalizeRows(data_filtered_unique)
#data3 <- data2
dim(data_filtered_unique_normalized)
gene_variances <- apply(data_filtered_unique_normalized, 1, var)
# Perform t-SNE
#tsne_result <- Rtsne(subset_df[, -ncol(subset_df)], perplexity = 50, theta = 0.5, dims = 2)
tsne_result <- Rtsne(data_filtered_unique_normalized, perplexity = 200, theta = 0.3, dims = 2)



# Assuming 'tsne_result$Y' is your t-SNE result matrix with two dimensions
colors <- rep("grey", nrow(data_filtered_unique))
colors[rownames(data_filtered_unique) %in% ensembl_ids_set_A] <- "red"
colors[rownames(data_filtered_unique) %in% ensembl_ids_set_B] <- "blue"
table(colors)

# Base plot with grey points
plot(tsne_result$Y, col = colors, pch = 19, cex = 0.5, xlab = "t-SNE 1", ylab = "t-SNE 2")

# Add red points on top
red_points <- tsne_result$Y[rownames(data_filtered_unique) %in% ensembl_ids_set_A, ]
points(red_points, col = "red", pch = 19, cex = 0.5)

# Add blue points on top
blue_points <- tsne_result$Y[rownames(data_filtered_unique) %in% ensembl_ids_set_B, ]
points(blue_points, col = "blue", pch = 19, cex = 0.5)

# Adding a legend to the plot
legend("topright", legend=c("Cell Cycle", "Adhesion"), col=c("red", "blue"), pch=19)



title("t-SNE of Genes (by Patients)")












### TSNE for patients
library(Rtsne)

# Transpose your dataframe where rows are genes and columns are patients
transposed_filtered_data <- t(data_filtered_unique)  # Assuming `data_filtered_unique` is your original data matrix

# Run t-SNE on the transposed data
tsne_result <- Rtsne(transposed_filtered_data, perplexity = 40, theta = 0.2, dims = 2, verbose = FALSE)

# Plot the t-SNE results
# Assuming you want to color by some patient metadata or simply use a single color
colors <- rep("blue", nrow(transposed_filtered_data))
# Optionally modify colors based on some metadata criteria if available

plot(tsne_result$Y, col = colors, pch = 19, cex = 0.5, 
     xlab = "t-SNE 1", ylab = "t-SNE 2", main = "Patient t-SNE Visualization")


# ------------------------------- WGCNA Network  --------------------
#install.packages("WGCNA")
library(WGCNA)

subset_df <- data_filtered[high_var_genes, ]
# Convert it to a dataframe
subset_df <- as.data.frame(subset_df)
class(subset_df)

dim(subset_df)






# Transpose the dataframe
transposed_matrix <- t(subset_df)

# Convert the transposed matrix back to a dataframe
transposed_dataframe <- as.data.frame(transposed_matrix)

# Optionally, you can set the row and column names
rownames(transposed_dataframe) <- colnames(subset_df)
colnames(transposed_dataframe) <- rownames(subset_df)

# Choose soft-threshold powers, ###Takes time about 1-2 Hrs
powers <- c(c(1:10), seq(from = 12, to = 50, by = 2))
sft <- pickSoftThreshold(transposed_dataframe,
                         powerVector = powers,
                         networkType = "signed",
                         verbose = 5)

# Plot the scale-free fit index and mean connectivity to help decide the power
sizeGrWindow(9, 5)
#par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, signed R^2",
     type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="blue")
# Choose the power with the highest scale-free topology fit (e.g., power = 6)
#The warnings "Some correlations are NA", indicate that there are missing values in the correlations. 
#This can occur if there are genes with constant or near-constant expression across all samples.
#To remove genes with very low variance, we can use the following R code:

# Calculate variance for each gene (now in columns)

# We opt for a power of 10 as it is the smallest power that gives an R2 value above 0.8. 
# This is a conservative choice and considered a good practice to start with the smallest power. 
# Alternatively, we could opt for higher powers like 12 or 14, which offer even higher R2 values, depending on what we are looking for in our analysis.
# Look at the plot to decide on the 'power' value
# This you have to decide after inspecting the plots


# Optimal soft-threshold power is determined from the plot
# The 'soft threshold' function is used to determine the optimal power for network construction in WGCNA.
# Power 10 shows a high scale-free topology fit index (SFT.R.sq = 0.912), indicating a good fit to a scale-free network.
# However, its mean connectivity (mean.k = 1.39) is relatively high, suggesting a more connected network than might be ideal.
# While Power 10 is a viable choice, considering slightly higher powers where SFT.R.sq remains high but mean connectivity decreases 
# may offer a better balance between scale-free topology and network connectivity.

#------------------------Network Construction and TOM Calculation---------------------------#
# selected 7 based on the 7
soft_power <- 12 # Replace with the selected value

# Choose a network type, for example, unsigned or signed
networkType <- "signed"  # Or "signed" based on your analysis


# Turn the expression data into an adjacency matrix
adjacency <- WGCNA::adjacency(transposed_dataframe, power = soft_power, type = networkType)

dim(adjacency)
mean(adjacency)

adjacency2 <- adjacency
# Assuming 'your_matrix' is your original matrix
# use this for graph visualization
#your_matrix[your_matrix < 0.1] <- 0
#your_matrix[your_matrix >= 0.1] <- 1

library(ggraph)
#install.packages(igraph)

# plot using igraph
library(igraph)

# Assuming adj_matrix is your adjacency matrix
#network <- graph_from_adjacency_matrix(adjacency, mode="undirected", weighted=TRUE)
#plot(network, vertex.size = 10, vertex.label.cex = 0.7, edge.width = E(network)$weight)

#igraph(network, layout = 'grid') + 
#  geom_edge_link() + 
#  geom_node_point() + 
#  theme_graph()

# Turn the adjacency into a TOM - takes time
TOM <- TOMsimilarity(adjacency)

row_sums <- rowSums(TOM, na.rm = TRUE) 
# Create a histogram
hist(row_sums, main = "Distribution Plot of node degrees", xlab = "Values", ylab = "Frequency", col = "blue")


# Convert the TOM into a dissimilarity format
dissTOM = 1 - TOM



geneTree = hclust(as.dist(dissTOM), method = "average")




# Plot the results

#par(mar = c(2, 2, 2, 2) + 0.1)  # Adjust margins (bottom, left, top, right)
plot(geneTree, labels = FALSE, hang = 0.01, main = "Gene Clustering Tree",
     xlab = "", sub = "")


# Dynamic tree cut to identify modules (color labels)
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, 
                            deepSplit = 2, pamRespectsDendro = FALSE, 
                            cutHeight = 0.99, minClusterSize = 20) #
table(dynamicMods)
dynamicColors <- labels2colors(dynamicMods)



names(dynamicMods) <- rownames(subset_df)
subset_df$Module <- as.factor(dynamicMods)

sample_cluster <- subset_df[subset_df$Module == 4, ]



library(clusterProfiler)
library(org.Hs.eg.db)
# Assuming 'geneList' is your vector of genes in a module
ego <- enrichGO(gene = rownames(sample_cluster),
                OrgDb = "org.Hs.eg.db", # Choose the appropriate OrgDb for your organism
                keyType = "ENSEMBL", # Ensure this matches your gene ID type
                ont = "BP", # For Biological Processes; use "MF" for Molecular Functions, "CC" for Cellular Components
                pAdjustMethod = "BH", # Method for adjusting p values
                qvalueCutoff = 0.05, # Set your significance level
                readable = TRUE) # To convert IDs to readable gene names

# View results
summary(ego)

# Convert the ego object to a data frame
ego_df <- as.data.frame(ego)

results <- list()

for(i in 0:33) { #length(table(dynamicMods))
  sample_cluster <- subset_df[subset_df$Module == i, ]
  
  ego <- enrichGO(gene = rownames(sample_cluster),
                  OrgDb = "org.Hs.eg.db", # Choose the appropriate OrgDb for your organism
                  keyType = "ENSEMBL", # Ensure this matches your gene ID type
                  ont = "BP", # For Biological Processes; use "MF" for Molecular Functions, "CC" for Cellular Components
                  pAdjustMethod = "BH", # Method for adjusting p values
                  qvalueCutoff = 0.05, # Set your significance level
                  readable = TRUE) # To convert IDs to readable gene names
  
  ego_df <- as.data.frame(ego)
  print(i)
  if("cell-matrix adhesion" %in% ego_df$Description[1:20]) {
    results[[length(results) + 1]] <- list(cluster = names(table(dynamicMods))[i],
                                           topGO = ego_df$Description[1:20])
    # If it matches, plot the dendrogram with module colors
    print('****')
    #break
  }
}

"synaptic vesicle exocytosis"
"cell-matrix adhesion"

adhesion <- c('JUP', 'ITGB4', 'L1CAM', 'ITGB8', 'ITGA11', 'THBS3', 'FN1', 
                'ITGB3', 'MSLN', 'RHOD', 'THBS1', 'COL16A1', 'JAG1', 'ADAMTS13',
                 'SNED1', 'CDKN2A', 'VEGFA', 'TMEM8B', 'PLAU', 'ADAMTS12')

ego <- enrichGO(gene = adhesion,
                OrgDb = "org.Hs.eg.db", # Choose the appropriate OrgDb for your organism
                keyType = "SYMBOL", # Ensure this matches your gene ID type
                ont = "BP", # For Biological Processes; use "MF" for Molecular Functions, "CC" for Cellular Components
                pAdjustMethod = "BH", # Method for adjusting p values
                qvalueCutoff = 0.05, # Set your significance level
                readable = TRUE) # To convert IDs to readable gene names

# View results
summary(ego)

# Convert the ego object to a data frame
ego_df <- as.data.frame(ego)



#dynamicMods[match(rownames(subset_df), names(dynamicMods))]
#subset_df$Module <- as.factor(subset_df$Module)

#labels <- as.factor(dynamicMods)
#subset_df <- data.frame(subset_df, Module = labels[rownames(subset_df)])


# Assuming 'data' is your dataframe with 17k rows and 400 columns

# Scale the data
scaled_data <- scale(subset_df)


# Assuming 'data' is your dataframe
set.seed(123) # for reproducibility

# Calculate total within-cluster sum of square
wss <- sapply(25:27, function(k){
  kmeans(scaled_data, centers=k, nstart=20)$tot.withinss
})

# Plot elbow plot
plot(25:27, wss, type="b", pch=19, frame=FALSE, xlab="Number of clusters K", ylab="Total within-clusters sum of squares")
# Determine the number of clusters (k)
# This is a simplistic approach; consider more sophisticated methods like the Elbow Method
k <- 35

# Apply k-means clustering
set.seed(123) # Setting a seed for reproducibility
kmeans_result <- kmeans(scaled_data, centers = k, nstart = 25)

# View the results
print(kmeans_result)

# To see the cluster assignment for each row
cluster_assignments <- kmeans_result$cluster

subset_df$Module <- cluster_assignments

library(Rtsne)
tsne_data <- Rtsne(subset_df[, -ncol(subset_df)], perplexity = 40, theta = 0.5,check_duplicates = FALSE)

#tsne_data <- Rtsne(subset_df, perplexity = 40, theta = 0.5,check_duplicates = FALSE)

#pca_result <- prcomp(adjacency, scale. = TRUE)
#tsne_data <- Rtsne(pca_result$x[, 1:50], perplexity = 20, theta = 0.5,check_duplicates = FALSE)

#install.packages("RColorBrewer")
library(RColorBrewer)

# Create a color palette with 12 distinct colors
color_palette <- brewer.pal(k, "Set3")

# Color palette
colors <- rainbow(length(unique(subset_df$Module)))
#colors <- rep("grey", nrow(subset_df))
# Plot
#tsne_data <- Rtsne(t(subset_df[, -ncol(subset_df)]), perplexity = 10, theta = 0.1,check_duplicates = FALSE)
plot(tsne_data$Y, col = color_palette[cluster_assignments], pch = 19, cex = 0.2, asp = 1, xlab = "t-SNE 1", ylab = "t-SNE 2")

plot(tsne_data$Y, col = colors[subset_df$Module], pch = 19, cex = 0.2, asp = 1, xlab = "t-SNE 1", ylab = "t-SNE 2")
legend("topright", legend = levels(subset_df$Module), col = colors, pch = 20)

### UMAP ####

#install.packages("umap")
library(umap)
library(ggplot2)
umap_result <- umap(subset_df[, -ncol(subset_df)])
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(size = .001) +
  theme_minimal() +
  xlab("UMAP 1") + ylab("UMAP 2") +
  ggtitle("UMAP Plot")

# Plot the dendrogram and colors underneath
#pdf(file = "4-module_tree.pdf", width = 8, height = 6);
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE, leaf=FALSE,
                    hang = 0.02,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")


# Assuming you have a dendrogram object 'geneTree'
# And a vector of genes of interest
genes_of_interest <- c("1", "Gene2", "Gene3")

# Convert the dendrogram to a structure easier to manipulate
dendro_data <- as.dendrogram(geneTree)

# Function to modify the labels of genes of interest
modify_labels <- function(dendro){
  if(is.leaf(dendro)){
    label <- attr(dendro, "label")
    if(label %in% genes_of_interest){
      attr(dendro, "label") <- paste0(label, "*")  # Add an asterisk or other marker
    }
  }
  dendro
}

# Apply the label modification function to the dendrogram
modified_dendro <- dendrapply(dendro_data, modify_labels)

# Clear the graphics device
#graphics.off()

# Now use plotDendroAndColors with the modified dendrogram
#pdf(file = "4-module_tree.pdf", width = 8, height = 6);
plotDendroAndColors(dendro = modified_dendro, colors = dynamicColors, 
                    "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")


# Apply the coloring function to the dendrogram
colored_dendro <- dendrapply(dendro_data, color_branches)

# Plot the colored dendrogram
plot(colored_dendro)









if (!requireNamespace("dendextend", quietly = TRUE)) {
  install.packages("dendextend")
}
library(dendextend)

# Convert the hclust object to a dendrogram
dendro <- as.dendrogram(geneTree)

num_modules <- length(unique(dynamicMods))
print(num_modules)
colored_dendro <- color_branches(dendro, k = length(unique(dynamicMods)), groupLabels = dynamicMods)
plot(colored_dendro)
# Color branches based on dynamic modules
colored_dendro <- color_branches(dendro, k = length(unique(dynamicMods)), groupLabels = dynamicMods)

# Plot the colored dendrogram
plot(colored_dendro)

# Add a color bar
# Add the color bar
y_lim <- par("usr")[3]  # Get the current y-axis limits
rect(min(members) - 0.5, y_lim - 1, max(members) + 0.5, y_lim - 0.5, col = dynamicColors[members], border = dynamicColors[members])
abline(h = -0.25, col = "red") # Optional: adds a line for clarity
colors_used <- unique(dynamicColors)
for (color in colors_used) {
  members <- which(dynamicColors == color)
  rect(min(members) - 0.5, -1, max(members) + 0.5, -0.75, col = color, border = color)
}








# ------------------------------- Gene Ontology --------------------




library(clusterProfiler)
library(org.Hs.eg.db)
# Assuming 'geneList' is your vector of genes in a module
ego <- enrichGO(gene = gene_set_A,
                OrgDb = "org.Hs.eg.db", # Choose the appropriate OrgDb for your organism
                keyType = "SYMBOL", # Ensure this matches your gene ID type
                ont = "BP", # For Biological Processes; use "MF" for Molecular Functions, "CC" for Cellular Components
                pAdjustMethod = "BH", # Method for adjusting p values
                qvalueCutoff = 0.05, # Set your significance level
                readable = TRUE) # To convert IDs to readable gene names

# View results
summary(ego)

# Convert the ego object to a data frame
ego_df <- as.data.frame(ego)

# Save the data frame to a text file
#write.table(ego_df, file = "12_31_24_WGCNA/ego_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(ego_df, file = "12_31_24_WGCNA/ego_results.csv", row.names = FALSE)
#### Using `topGO`:
#topGO` is another R package that specializes in testing GO terms while accounting for the hierarchical structure of GO annotations. It's particularly useful for identifying GO terms that are significantly enriched at different levels of the GO hierarchy.


#BiocManager::install("ALL")

library(topGO)
library(ALL)
data(ALL)

## discriminate B-cell from T-cell
classLabel <- as.integer(sapply(ALL$BT, function(x) return(substr(x, 1, 1) == 'T')))

## Differentially expressed genes
geneList <- getPvalues(exprs(ALL), classlabel = classLabel,
                       alternative = "greater", correction = "BY")

hist(geneList, 50) 
# Example data preparation and analysis steps will vary based on your specific dataset and needs

groupGOTerms()

### 3. **Interpretation**
# ------------------------------- Graph Visualization --------------------

#BiocManager::install("ggraph")
library(ggraph)


#install.packages(igraph)
library(igraph)
# Assuming adj_matrix is your adjacency matrix
network <- graph_from_adjacency_matrix(adj_A, mode="undirected", weighted=TRUE)
plot(network, vertex.size=10, vertex.label.cex=0.7, edge.width=E(graph)$weight)
ggraph(network, layout = 'grid') + 
  geom_edge_link() + 
  geom_node_point() + 
  theme_graph()





#install.packages("beepr")
library(beepr)



# Run blockwiseModules and store the output in bwnetModules
bwnetModules <- blockwiseModules(data4,
                                 maxBlockSize = 5000,
                                 TOMType = "signed",
                                 power = soft_power,
                                 mergeCutHeight = 0.25,
                                 numericLabels = FALSE,
                                 randomSeed = 1234)

beepr::beep(sound = 8)  # Different built-in sound


#=============================================================================
# Define the file path for the bwnetModules object
bwnetModules_file_path <- file.path(intermediate_Data, "bwnetModules.rds")

# Save the bwnetModules object using the new file path
saveRDS(bwnetModules, file = bwnetModules_file_path)
#=============================================================================

# Load necessary library
library(WGCNA)


# Keep dissTOM in a separate object
dissTOMObject <- list()
dissTOMObject$dissTOM <- dissTOM
# When you need to use the module colors and dendrogram from blockwiseModules
# Use bwnetModules$colors and bwnetModules$dendrograms

# When you need to use dissTOM
# Use dissTOMObject$dissTOM

# Blockwise module detection

cor <- WGCNA::cor


# Extract module eigengenes
module_eigengenes <- bwnetModules$MEs

# Since we don't have traits, we can skip the module-trait correlation step.

# Continue with other downstream analyses as necessary


# ---------------------------- Visualisation_new GPT suggested method------------------------------#
#This part plots the scale-free topology model fit and the mean connectivity as a function of soft-thresholding power.

#Before plotting, inspect the values of sft$fitIndices[,3] to see if there are any non-positive values.
print(sft$fitIndices[,3])


# Set up the layout for two plots
par(mfrow = c(1, 1), mar = c(5, 4, 4, 5) + 1)
cex1 = 0.9  # Size of the text labels

# Plot for scale-free topology model fit (R²)
plot(sft$fitIndices[,1], sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit (R^2)",
     type="n",
     main = "Scale Independence",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# Adding text labels
text(sft$fitIndices[,1], sft$fitIndices[,2],
     labels=sft$fitIndices[,1],
     cex=0.9, col="red")

# Plot for mean connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity",
     type="n",
     main = "Mean Connectivity",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels=sft$fitIndices[,1], cex=cex1, col="red")



# Visualize the gene dendrogram and the module colors underneath
plotDendroAndColors(bwnetModules$dendrograms[[1]],
                    cbind(bwnetModules$unmergedColors, bwnetModules$colors),
                    c("Unmerged", "Merged"),
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)

str(bwnetModules)

#________________________New Suggested Visualization_______________________#


for(i in 1:length(bwnetModules$dendrograms)) {
  geneTree = bwnetModules$dendrograms[[i]]
  mergedColors = labels2colors(moduleColors)
  
  if(length(mergedColors) == length(geneTree$order)) {
    # If it matches, plot the dendrogram with module colors
    plotDendroAndColors(geneTree, 
                        mergedColors,
                        "Module Colors",
                        dendroLabels = FALSE, 
                        hang = 0.03,
                        addGuide = TRUE, 
                        guideHang = 0.05)
    break
  }
}

# If no matching dendrogram is found, print a message
if(i == length(bwnetModules$dendrograms)) {
  print("No matching dendrogram found for the provided module colors.")
}

length(geneTree$order)

length(mergedColors)

# Find the block number corresponding to your dendrogram
block_number <- which(sapply(bwnetModules$dendrograms, function(x) identical(x, geneTree)))

# Extract colors for the specific block
block_genes <- bwnetModules$blockGenes[[block_number]]
block_colors <- bwnetModules$colors[block_genes]

# Convert these module labels to colors
block_mergedColors <- labels2colors(block_colors)

# Now check if the lengths match
if(length(block_mergedColors) == length(geneTree$order)) {
  plotDendroAndColors(geneTree, 
                      block_mergedColors,
                      "Module Colors",
                      dendroLabels = FALSE, 
                      hang = 0.03,
                      addGuide = TRUE, 
                      guideHang = 0.05)
} else {
  print("Mismatch in length of colors and dendrogram branches for the selected block")
}



#-------------------------------------------------------------------------------
# to check if whoch module we plotted the dendrogram for

# Check the total number of blocks
total_blocks <- length(bwnetModules$dendrograms)
print(paste("Total number of blocks:", total_blocks))

# Find the block number corresponding to your dendrogram
block_number <- which(sapply(bwnetModules$dendrograms, function(x) identical(x, geneTree)))

# Check if block_number is found
if (length(block_number) > 0) {
  print(paste("Block number identified:", block_number))
  
  # Extract the genes and colors for the identified block
  block_genes <- bwnetModules$blockGenes[[block_number]]
  block_colors <- bwnetModules$colors[block_genes]
  
  # Convert these module labels to colors
  block_mergedColors <- labels2colors(block_colors)
  
  # Check if the lengths of colors and dendrogram branches match
  if(length(block_mergedColors) == length(geneTree$order)) {
    # Plot the dendrogram with module colors
    plotDendroAndColors(geneTree, 
                        block_mergedColors,
                        "Module Colors",
                        dendroLabels = FALSE, 
                        hang = 0.03,
                        addGuide = TRUE, 
                        guideHang = 0.05)
  } else {
    print("Mismatch in length of colors and dendrogram branches for the selected block")
  }
} else {
  print("No matching block found for the provided dendrogram.")
}

#-------------------------------------------------------------------------------

# Get the total number of blocks in the bwnetModules object
total_blocks <- length(bwnetModules$dendrograms)

# Loop through each block
for(i in 1:total_blocks) {
  # Extract the dendrogram for the current block
  geneTree <- bwnetModules$dendrograms[[i]]
  
  # Extract the gene indices for the current block
  block_genes <- bwnetModules$blockGenes[[i]]
  
  # Extract the colors assigned to genes in the current block
  block_colors <- bwnetModules$colors[block_genes]
  
  # Convert module labels to colors for plotting
  block_mergedColors <- labels2colors(block_colors)
  
  # Check if the length of colors matches the number of branches in the dendrogram
  if(length(block_mergedColors) == length(geneTree$order)) {
 
    # If it matches, plot the dendrogram with module colors
    plotTitle <- paste("Cluster Dendrogram - Block", i)
    plotDendroAndColors(geneTree, 
                        block_mergedColors,
                        paste("Module Colors"), # Title indicating the block number
                        dendroLabels = FALSE, 
                        hang = 0.03,
                        addGuide = TRUE, 
                        guideHang = 0.05)
  } else {
    print(paste("Mismatch in length of colors and dendrogram branches for Block", i))
  }
}

#-------------------------------------------------------------------------------



# Define the paths to the gene set files
geneset_path1 <- "/Users/sadaf/Library/CloudStorage/Dropbox/#Bioinformatics/2024_01_06_WGCNA/GeneSets/GOBP_CELL_MATRIX_ADHESION.v7.5.1.gmt"
geneset_path2 <- "/Users/sadaf/Library/CloudStorage/Dropbox/#Bioinformatics/2024_01_06_WGCNA/GeneSets/GOBP_CELL_CYCLE.v2023.1.Hs.gmt"

# Load the gene sets
geneset1 <- getGmt(geneset_path1)
geneset2 <- getGmt(geneset_path2)

head(geneset1)
head(geneset2)

# Prepare the data for GSVA or other analysis methods
# Assuming 'data4' is your final expression data used in WGCNA
# Transpose data4 to have genes in rows and samples in columns
data4_gsva <- t(data4)
head(data4_gsva)

# Assuming geneset1 and geneset2 are GeneSetCollection objects
# gene_ids is a vector of gene names from your WGCNA data

# Convert a GeneSetCollection to a binary matrix
geneset_to_matrix <- function(geneset, all_genes) {
  matrix_genes <- lapply(geneset, geneIds)
  binary_matrix <- sapply(matrix_genes, function(gs_genes) as.integer(all_genes %in% gs_genes))
  colnames(binary_matrix) <- names(geneset)
  return(binary_matrix)
}

# Convert gene sets to binary matrices
binary_geneset1 <- geneset_to_matrix(geneset1, gene_ids)
binary_geneset2 <- geneset_to_matrix(geneset2, gene_ids)

# Combine gene set matrices
combined_genesets <- cbind(binary_geneset1, binary_geneset2)

# Ensure the rows of the module eigengenes and gene sets match
aligned_module_eigengenes <- module_eigengenes[rownames(combined_genesets), ]

# Calculate correlation
module_genesets_correlation <- cor(combined_genesets, aligned_module_eigengenes, use = "pairwise.complete.obs")

## Error: Y has zero dimentions

# Check the structure and contents of the binary gene sets
str(binary_geneset1)
str(binary_geneset2)

# Check the structure and contents of the module eigengenes
str(module_eigengenes)

# Check if the gene IDs in the module eigengenes match with gene_ids
any(rownames(module_eigengenes) %in% gene_ids)



##Error:gene ids dont match
# Inspect a few gene identifiers from your gene sets
print(head(geneIds(geneset1)))
print(head(geneIds(geneset2)))

# Inspect a few gene identifiers from your WGCNA data
print(head(gene_ids))

# You might need to convert identifiers to a common format here
# Example: Convert gene symbols in gene_ids to Ensembl IDs (or vice versa)
# This will depend on the formats you have and the desired common format




# Connect to the Ensembl database
#ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

##Error: Ensembl service is currently unavailable
# Specify a mirror, for example, use the US mirror
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")

# Convert gene symbols to Ensembl IDs for geneset1
gene_symbols_geneset1 <- geneIds(geneset1)[[1]]
geneset1_ensembl <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                          filters = 'hgnc_symbol', 
                          values = gene_symbols_geneset1, 
                          mart = ensembl)

# Repeat for geneset2 (assuming you have geneset2 loaded)
gene_symbols_geneset2 <- geneIds(geneset2)[[1]]
geneset2_ensembl <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                          filters = 'hgnc_symbol', 
                          values = gene_symbols_geneset2, 
                          mart = ensembl)



# Function to convert a vector of gene IDs into a binary matrix
# gene_ids_vector: A vector of gene IDs (e.g., Ensembl IDs)
# all_genes: A vector of all gene IDs in the dataset
geneset_to_matrix <- function(gene_ids_vector, all_genes) {
  # Iterate over each gene ID in gene_ids_vector. For each gene ID, check if it is present in all_genes
  # Convert the logical vector (TRUE/FALSE) to a binary vector (1/0)
  binary_matrix <- sapply(gene_ids_vector, function(gene_id) as.integer(all_genes %in% gene_id))
  
  # Return the binary matrix
  return(binary_matrix)
}

# Example usage of the function
# Assuming gene_ids_vector contains gene IDs and all_genes is the list of all genes in your dataset
# binary_geneset1 <- geneset_to_matrix(gene_ids_vector, all_genes)
# Convert Ensembl IDs to binary matrices




binary_geneset1 <- geneset_to_matrix(geneset1_ensembl$ensembl_gene_id, gene_ids)
binary_geneset2 <- geneset_to_matrix(geneset2_ensembl$ensembl_gene_id, gene_ids)

# Combine gene set matrices
combined_genesets <- cbind(binary_geneset1, binary_geneset2)

# Ensure the rows of the module eigengenes and gene sets match
aligned_module_eigengenes <- module_eigengenes[rownames(combined_genesets), ]

# Calculate correlation
module_genesets_correlation <- cor(combined_genesets, aligned_module_eigengenes, use = "pairwise.complete.obs")
##Error y has zero dimention


# Check the dimensions of the combined gene set matrices
print(dim(combined_genesets))

# Check the dimensions of the module eigengenes
print(dim(module_eigengenes))

# Ensure the rownames of combined_genesets are present in module eigengenes
print(any(rownames(combined_genesets) %in% rownames(module_eigengenes)))

# If the above check returns FALSE, it indicates a mismatch in row names
# You might need to align the rownames of combined_genesets with module eigengenes

##Error: not aligned

# Create a named vector of Ensembl IDs
ensembl_ids <- setNames(module_ensembl$ensembl_gene_id, module_ensembl$hgnc_symbol)

# Initialize a copy of module eigengenes with original row names
module_eigengenes_ensembl <- module_eigengenes

# Replace row names with Ensembl IDs where available
for(gene in rownames(module_eigengenes)) {
  if(gene %in% names(ensembl_ids)) {
    rownames(module_eigengenes_ensembl)[rownames(module_eigengenes_ensembl) == gene] <- ensembl_ids[gene]
  }
}

# Proceed with alignment and correlation calculation

# Filter out rows from combined_genesets that are not in module_eigengenes_ensembl
filtered_combined_genesets <- combined_genesets[rownames(combined_genesets) %in% rownames(module_eigengenes_ensembl), ]

# Filter out rows from module_eigengenes_ensembl that are not in filtered_combined_genesets
filtered_module_eigengenes <- module_eigengenes_ensembl[rownames(module_eigengenes_ensembl) %in% rownames(filtered_combined_genesets), ]

# Ensure the order of rownames in both matrices is the same
filtered_combined_genesets <- filtered_combined_genesets[match(rownames(filtered_module_eigengenes), rownames(filtered_combined_genesets)), ]

# Calculate the correlation
module_genesets_correlation <- cor(filtered_combined_genesets, filtered_module_eigengenes, use = "pairwise.complete.obs")

##Erroor X is zero

# Check dimensions of filtered_combined_genesets
print(dim(filtered_combined_genesets))

# Check dimensions of filtered_module_eigengenes
print(dim(filtered_module_eigengenes))

# Check if there are any common gene identifiers
print(any(rownames(filtered_combined_genesets) %in% rownames(filtered_module_eigengenes)))
##Error 

# Check how many gene symbols from geneset1 were successfully converted
print(length(gene_symbols_geneset1))
print(nrow(geneset1_ensembl))

# Repeat for geneset2
print(length(gene_symbols_geneset2))
print(nrow(geneset2_ensembl))

# Check for module eigengenes
print(length(gene_symbols_module))
print(nrow(module_ensembl))



# View the correlation matrix
print(module_genesets_correlation)




#==================================================================
# List all objects in the environment
all_objects <- ls()

# Loop through each object and print its structure
for (obj_name in all_objects) {
  cat("\nStructure of", obj_name, ":\n")
  str(get(obj_name))
}

#===== unused code from data preparation =====
# this code normalize the gene expression with (X-Mu)/Sigma
normalizeRows <- function(df) {
  # Function to normalize a single row
  normalizeRow <- function(row) {
    if (var(row) == 0) {
      return(row)  # Return the row unchanged if variance is 0
    } else {
      return((row - mean(row)) / sqrt(var(row)))
    }
  }
  
  # Apply the normalization function to each row
  normalized_df <- t(apply(df, 1, normalizeRow))
  
  # Setting the row and column names of the normalized dataframe to be the same as the original dataframe
  rownames(normalized_df) <- rownames(df)
  colnames(normalized_df) <- colnames(df)
  
  return(normalized_df)
}


# Assuming your original data is stored in 'data'

data3 <- normalizeRows(data2)
#data3 <- data2
dim(data3)
gene_variances <- apply(data3, 1, var)
threshold <- 0.00
data_filtered <- data3[gene_variances > threshold , ]
data_filtered_original <- data1[gene_variances > threshold , ]
#rownames(data_filtered) <- data_filtered_original$gene_id

dim(data_filtered)
class(data_filtered)
#data_t <- t(data)
#data_t <- as.data.frame(data_t)

