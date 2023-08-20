# Load necessary libraries
load_library("ChAMP")
load_library("doParallel")
load_library("factoextra")
load_library("ggfortify")
load_library("RColorBrewer")
load_library("viridis")
load_library("randomcoloR")
load_library("M3C")
load_library("pheatmap")
load_library("FlowSorted.Blood.EPIC")
load_library("IlluminaHumanMethylationEPICanno.ilm10b3.hg19")
load_library("AnnotationHub")
load_library("ExperimentHub")
load_library("pROC")
load_library("corrr")
load_library("tidyverse")
load_library("dataframes2xls")
load_library("wateRmelon")
load_library("qqman")
load_library("missMethyl")

data_directory=setwd("PATH/directory")

# Load data using ChAMP
load_minfi <- champ_load(directory = data_directory, arraytype = "EPIC", method = "minfi", force = TRUE)

# Perform quality control
qc_results <- champ_QC(beta = load_minfi$beta, pheno = load_minfi$pd$Sample_Group)

# Perform normalization using FunctionalNormalization
normalized_data <- champ_norm(beta = load_minfi$beta, rgSet = load_minfi$rgSet, method = "FunctionalNormalization", arraytype = "EPIC")

# Check sex-age samples
GRset <- mapToGenome(load_minfi$rgSet)
sexTable <- getSex(GRset)
Age<-agep(normalized_data)

# Select significant controls that match in terms of sex and age with affected subjects
load_library("MatchIt")
# Filter out rows with NA in 'Treat' column
filtered_data <- load_minfi$pd[!(is.na(load_minfi$pd$Treat)),]
# Perform matching using MatchIt package
match.it <- matchit(Treat ~ Age + Sex + Slide, data = filtered_data, method = "nearest", ratio = 4)
# Plot the matching results
plot_matching_results(match.it, type = 'jitter', interactive = FALSE)
# Display summary of matching results
matching_summary <- summarize_matching_results(match.it)
# Select the best ratio that provides desired results
best_ratio <- select_best_ratio(matching_summary)
# Extract the corrected dataset after matching
corrected_data <- extract_corrected_data(match.it)

# Estimate cell counts
cell_counts <- estimateCellCounts2(load_minfi$rgSet, compositeCellType = "Blood", processMethod = "auto", probeSelect = "IDOL", cellTypes = ["CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"])

# Perform SVD analysis for batch effects check
svd_results <- champ_SVD(normalized_data, pd = load_minfi$pd)

# Calculate M-values
clean_norm <- select_columns(normalized_data, columns = which_in(tmp_compare, load_minfi$pd$Classification))
clean_mval <- calculate_m_values(clean_norm)

# Create the design matrix using the significant covariates from the svd plot 
design <- create_design_matrix(p, pAge, pSlide, pCellType1, pCellType2, pCellType3, pCellType4, pCellType5)
# Create the contrast matrix
contrast_matrix <- create_contrast_matrix(design)
# Display the contrast matrix
display_contrast_matrix(contrast_matrix)
# Perform linear model fitting
fit <- perform_linear_model_fit(cleanMval, design)
# Perform contrasts fitting
fit2 <- perform_contrasts_fit(fit, contrast_matrix)
# Perform eBayes analysis
fit3 <- perform_eBayes_analysis(fit2)

# Select differentially methylated probes (DMPs) using adjusted p-value and fold change threshold
DMP <- select_DMPs(fit3, coefficient = 1, numProbes = num_rows(cleanMval), adjustMethod = "BH", pValueThreshold = 0.05, foldChangeThreshold = 0.07)
# Select all probes without adjusting p-value
DMP_all <- select_DMPs(fit3, coefficient = 1, numProbes = num_rows(cleanMval), adjustMethod = "BH", pValueThreshold = 1)


# Calculate DMP scores and select top probes
probe_scores <- abs(DMP$logFC) * -log(DMP$adj.P.Val)
sorted_probes <- sort_probes_by_scores(probe_scores)
top_probes <- select_top_probes(sorted_probes, count = 1000)

# Select the best probes using AUC values
selected_probes <- select_probes_by_AUC(top_probes, load_minfi$pd$Classification, normalized_data)

# Calculate correlation matrix and remove highly correlated probes
correlation_matrix <- calculate_correlation_matrix(selected_probes)
low_correlation_probes <- filter_low_correlation_probes(correlation_matrix, threshold = 0.8)

# Refine the final signature
refined_signature <- select_non_correlated_probes(selected_probes, low_correlation_probes)

----------------------------------------------------------------------------------------
#MDS plot of the refined_signature
create_MDS_plot(
  data =refined_signature,
  sampGroups = load_minfi$pd$Classification,
  pch = 19,
  cex = 2.5,
  legendPos = "topleft",
  legendNCol = 1,
  main = "MDS beta values episignature",
  pal = ["#5999FF49", "red"]
)


# Hiearchical Clustering using the refined_signature 
# Define color palette for heatmap
my_palette <- create_color_palette(c("blue", "yellow"), n_colors = 90)
# Define annotation colors for groups
ann_colors = create_annotation_colors(Groups = create_color_map("",""))
# Convert normalized signature matrix to a matrix and remove rows with missing values
signatNormMatrix <- as_matrix(remove_rows_with_missing(refined_signature))
# Get sample groups
Groups <- get_sample_groups(load_minfi$pd$Classification)
# Perform hierarchical clustering using pvclust
fitClust <- perform_hierarchical_clustering(signatNormMatrix, method.hclust="ward.D2", method.dist="euclidean")
# Create a heatmap using pheatmap
create_heatmap(
  data = signatNormMatrix,
  scale = 'none',
  col = my_palette,
  annotation_col = create_annotation_dataframe(Groups),
  annotation_colors = ann_colors,
  fontsize = 15,
  cluster_cols = fitClust$hclust,
  clustering_method = "ward.D2"
)

--------------------------------------------------------------------------------------------
  
# Compute the  differential methylation regions
  
# Annotate CpG sites and determine significant CpG sites
annotated_cpg <- annotate_cpg_sites(data_type = "array", fdr_cutoff = 0.05, m_values = filtered_M_values, design = design, contrast_coefficients = "pwt-paffectedTraining", analysis_type = "differential", annotation = list(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19"), use_contrasts = TRUE, contrast_matrix = contrast_matrix)
# Perform DMR analysis
dmr_output <- perform_dmr_analysis(annotated_cpg, min_cpgs = 5, lambda = 1000, C = 2, beta_cutoff = 0.1)
# Extract DMR information
extracted_dmr <- as.data.frame(extract_dmr_ranges(dmr_output, genome = "hg19"))
extracted_dmr_original <- extract_dmr_ranges(dmr_output, genome = "hg19")
# Assign unique names to DMRs
assign_dmr_names(extracted_dmr)
# Prepare the output in a structured format
output_dmr <- create_dmr_output_list(DMRcateDMR = extracted_dmr, DMRcateGrange = extracted_dmr_original)


# Functional enrichment con missMethyl
gst.region <- goregion(extracted_dmr_original, collection="KEGG", array.type="EPIC", plot.bias=TRUE)
hallmark <- readRDS("Hs.h.all.v7.1.entrez.rds")
gsa.region <- gsaregion(extracted_dmr_original, collection=hallmark, array.type="EPIC")
topGSA(gst.region, n=10)
topGSA(gsa.region, n=10)

-----------------------------------------------------------------------------------------------

## SVM Classifier ##

# Load necessary libraries
library(e1071) # For SVM implementation
library(DMwR)  # For SMOTE implementation

# Load and preprocess data
data <- read.csv("data.csv") # Load your dataset
features <- data[, -ncol(data)] # Extract features
labels <- data$labels # Extract labels

# Split data into training and testing sets
set.seed(123) # Set a seed for reproducibility
train_indices <- sample(1:nrow(data), 0.7 * nrow(data)) # 70% for training
test_indices <- setdiff(1:nrow(data), train_indices) # Remaining for testing

train_data <- features[train_indices, ]
train_labels <- labels[train_indices]
test_data <- features[test_indices, ]
test_labels <- labels[test_indices]

# Apply SMOTE to oversample the minority class
train_data_balanced <- SMOTE(train_data, train_labels, perc.over = 200)

# Train the Support Vector Machine (SVM) classifier
svm_model <- svm(train_labels ~ ., data = train_data_balanced, type="nu-classification", kernel="linear", probability = TRUE)

# Make predictions on the test set
predictions <- predict(svm_model, newdata = test_data)

# Evaluate the performance of the classifier
confusion_matrix <- table(predictions, test_labels)
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)

# Print the confusion matrix and accuracy
print("Confusion Matrix:")
print(confusion_matrix)
print(paste("Accuracy:", accuracy))
-----------------------------------------------------------------------------------------------








