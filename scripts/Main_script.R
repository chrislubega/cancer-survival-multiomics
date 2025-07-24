library(tidyverse)
library(matrixStats)
library(ggfortify)
library(factoextra)
library(survival)
library(survminer)
library(Rtsne)
library(FactoMineR)
library(factoextra)
library(NMF)
library(cluster)

# ----------------------------
# Prepare data
# ----------------------------
expression <- read_excel("data/fpkm.xlsx")
Clinical <- read_excel("data/Clinical.xlsx", sheet = "TCGA-SARC.clinical")
survival <- read_tsv("data/TCGA-SARC.survival.tsv")


geneids <- expression[[1]]
expr_t <- expression[,-1]|>t() |> as.data.frame()
colnames(expr_t) <- geneids
expr_t$sample <- rownames(expr_t)
Data1 <- inner_join(Clinical,expr_t, by="sample")
Data2 <- inner_join(survival,Data1, by ="sample")
Data2 <- Data2%>%
  mutate(across(6:15,as.factor))

# Separate clinical and expression data

clin_data <- Data2[, 1:15]
expr_data <- Data2[, -c(1:15)]


# ----------------------------
# Feature Selection
# ----------------------------

# Filter low variance genes
expr_mat <- as.matrix(expr_data)
gene_vars <- rowVars(t(expr_mat))  # variance across genes

# Keep top N variable genes (e.g. 1000)
top_n <- 1000
top_genes_idx <- order(gene_vars, decreasing = TRUE)[1:top_n]
expr_top <- expr_mat[, top_genes_idx]

# Save for unsupervised analysis
write.csv(expr_top, "data/top_variable_genes.csv", row.names = FALSE)
dim(expr_top)

# ----------------------------
# PCA for Visualization
# ----------------------------

# Combine PCA and clinical for visualization
pca_res <- prcomp(expr_top, scale. = TRUE)
autoplot(pca_res, data = clin_data, colour = "disease_type") +
  theme_minimal() +
  labs(title = "PCA of Top 1000 Variable Genes")


# --- t-SNE Visualization ---
set.seed(42)
tsne_res <- Rtsne(expr_top, perplexity = 30)
tsne_df <- as.data.frame(tsne_res$Y)
colnames(tsne_df) <- c("tSNE1", "tSNE2")

# --- Clustering ---
# Choose number of clusters k
k <- 3  

# KMeans Clustering
set.seed(123)
kmeans_res <- kmeans(expr_top, centers = 3)
cluster_labels <- kmeans_res$cluster

# Alternatively: NMF
#expr_top_t <- t(expr_top)
#nmf_res <- nmf(expr_top_t, rank = k, method = "brunet", nrun = 30)
#cluster_labels1 <- predict(nmf_res)

# Add cluster info to clinical data
clin_data$cluster <- factor(cluster_labels[Data2$sample %in% clin_data$sample])

# --- Visualize Clusters in t-SNE ---
tsne_df$cluster <- clin_data$cluster
ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = cluster)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "t-SNE of Molecular Subtypes (Clusters)")

# Save cluster-annotated expression and clinical data
write.csv(clin_data, "clinical_with_clusters.csv")


# ----------------------------
# Survival Analysis
# ----------------------------

# Fit KM model
km_fit <- survfit(Surv(OS.time, OS) ~ cluster, data = clin_data)

# Plot KM curve
ggsurvplot(
  km_fit, data = clin_data,
  pval = TRUE, risk.table = F,
  title = "Overall Survival by Molecular Subtype",
  xlab = "Time (days)", ylab = "Survival Probability",
  palette = "Dark2"
)

cox_uni <- coxph(Surv(OS.time, OS) ~ cluster, data = clin_data)
summary(cox_uni)

cox_multi <- coxph(Surv(OS.time, OS) ~ cluster + age_at_index.demographic
                   , data = clin_data)
summary(cox_multi)




