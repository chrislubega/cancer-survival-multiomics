library(tidyverse)
library(readxl)
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
library(umap)          # UMAP
library(glmnet)        # Lasso-Cox
library(matrixStats)
# ----------------------------
# Load raw data
# ----------------------------
expression <- read_excel("data/fpkm.xlsx")
mirna <- read_excel("data/mirna.xlsx")
clinical <- read_excel("data/Clinical.xlsx")
survival <- read_tsv("data/TCGA-SARC.survival.tsv")

# ----------------------------
# Prepare gene expression data
# ----------------------------
gene_ids <- expression[[1]]
expr_mat <- expression[, -1] |> t() |> as.matrix()
colnames(expr_mat) <- gene_ids

# Filter low-expressed genes
keep_genes <- colMeans(expr_mat > 1) > 0.5
expr_mat <- expr_mat[, keep_genes]
expr_mat <- expr_mat |> as.data.frame() |>
  rownames_to_column(var = "sample")

# ----------------------------
# Merge clinical and survival data
# ----------------------------
# ----------------------------
metadata<- inner_join(clinical, survival, by = "sample")
clin_data_gene <- inner_join(metadata, expr_mat, by = "sample")
clin_data <- clin_data_gene[,1:15]


clin_data_gene <- clin_data_gene %>%
  column_to_rownames(var = "sample")

# Normalize and convert rownames to column
expr_z <- scale(clin_data_gene[,-c(1:15)])|> as.data.frame()|>
  rownames_to_column(var = "sample")
  

# ----------------------------
# Prepare miRNA expression data
# ----------------------------
mirna_ids <- mirna[[1]]
mirna_mat <- mirna[, -1] |> t() |> as.data.frame()
colnames(mirna_mat) <- mirna_ids

mirna_mat <- mirna_mat %>%
  rownames_to_column(var = "sample")
mirna_mat <- inner_join(metadata, mirna_mat, by = "sample")%>%
  column_to_rownames(var = "sample")

# Normalize and convert rownames to column
mirna_z <- scale(mirna_mat[,-c(1:15)]) |>
  as.data.frame() |>
  rownames_to_column(var = "sample")


# Keep top N variable genes (e.g. 1000)
gene_vars <- colVars(as.matrix(expr_z[,-1]))
top_n <- 1000
top_genes_idx <- order(gene_vars, decreasing = TRUE)[1:top_n]
expr_top <- expr_z[, top_genes_idx] 


# Save for unsupervised analysis
#write.csv(expr_top, "data/top_variable_genes.csv", row.names = FALSE)
dim(expr_top)

# ----------------------------
# PCA for Visualization
# ----------------------------

# Combine PCA and clinical for visualization
pca_res <- prcomp(expr_top, scale. = TRUE)
tiff("results/pca_clusters.tiff", width = 6, height = 5, units = "in", res = 300, compression = "lzw")
autoplot(pca_res, data = clin_data_gene, colour = "disease_type") +
  theme_minimal() +
  theme(legend.position = "bottom")+
  guides(color = guide_legend(nrow = 3)) +
    labs(title = "PCA of Top 1000 Variable Genes")
dev.off()

# --- t-SNE Visualization ---
set.seed(42)
tsne_res <- Rtsne(expr_top, perplexity = 30)
tsne_df <- as.data.frame(tsne_res$Y)
colnames(tsne_df) <- c("tSNE1", "tSNE2")

# UMAP
umap_res <- umap(expr_top)
umap_df <- as.data.frame(umap_res$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")

# --- Clustering ---
# Choose number of clusters k
k <- 3 
set.seed(25)
kmeans_res <- kmeans(expr_top, centers = k, nstart = 25)
sil <- silhouette(kmeans_res$cluster, dist(expr_top))
fviz_silhouette(sil)


# KMeans Clustering
set.seed(123)
kmeans_res <- kmeans(expr_top, centers = 3)
cluster_labels <- kmeans_res$cluster


# Alternatively: NMF
#expr_top_t <- t(expr_top)
#nmf_res <- nmf(expr_top_t, rank = k, method = "brunet", nrun = 30)
#cluster_labels1 <- predict(nmf_res)

# Add cluster info to clinical data
clin_data$cluster <- factor(cluster_labels[expr_z$sample %in% clin_data$sample])

# --- Visualize Clusters in t-SNE ---
tsne_df$cluster <- clin_data$cluster
tiff("results/tsne_clusters.tiff", width = 6, height = 5, units = "in", res = 300, compression = "lzw")
ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = cluster)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "t-SNE of Molecular Subtypes (Clusters)")
dev.off()


umap_df$cluster <- clin_data$cluster
tiff("results/umap_clusters.tiff", width = 6, height = 5, units = "in", res = 300, compression = "lzw")
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "UMAP of Molecular Subtypes")
dev.off()
# Save cluster-annotated expression and clinical data
write.csv(clin_data, "clinical_with_clusters.csv")

# ----------------------------
# Survival Analysis
# ----------------------------

# Fit KM model
km_fit <- survfit(Surv(OS.time, OS) ~ cluster, data = clin_data)

# Plot KM curve
ggsurv <- ggsurvplot(
  km_fit, data = clin_data,
  pval = TRUE, risk.table = F,
  title = "Overall Survival by Molecular Subtype",
  xlab = "Time (days)", ylab = "Survival Probability",
  palette = "Dark2"
)
tiff("results/km_curve.tiff", width = 6, height = 5, units = "in", res = 300, compression = "lzw")
print(ggsurv$plot)
dev.off()

cox_uni <- coxph(Surv(OS.time, OS) ~ cluster, data = clin_data)
summary(cox_uni)

cox_multi <- coxph(Surv(OS.time, OS) ~ cluster + age_at_index.demographic
                   , data = clin_data)
summary(cox_multi)

# ----------------------------
#miRNA Integration
# ----------------------------
#select predictive genes 

y <- Surv(clin_data$OS.time, clin_data$OS)
fit <- cv.glmnet(as.matrix(expr_z[,-1]), y, family = "cox", alpha = 1)
selected <- rownames(coef(fit, s = "lambda.min"))[which(coef(fit, s = "lambda.min") != 0)]
predictive_top <- expr_z[, selected] |> as.data.frame()
predictive_top$sample <- expr_z$sample
top_gene_names <- colnames(predictive_top)[1:10]  
clin_data_predictive <- inner_join(clin_data, predictive_top, by = "sample")

# Prepare survival response
y_mirna <- Surv(clin_data$OS.time, clin_data$OS)

# Fit Lasso-Cox
fit_mirna <- cv.glmnet(as.matrix(mirna_z[,-1]), y_mirna, family = "cox", alpha = 1)

# Extract selected features
selected_mirnas <- rownames(coef(fit_mirna, s = "lambda.min"))[
  which(coef(fit_mirna, s = "lambda.min") != 0)]

# Remove intercept if present
selected_mirnas <- setdiff(selected_mirnas, "(Intercept)")

# Subset normalized miRNA expression
mirna_top <- mirna_mat[, selected_mirnas]|> as.data.frame()
mirna_mat <- mirna_mat %>%
  rownames_to_column(var = "sample")
mirna_top$sample <- mirna_mat$sample
clin_data_mirna <- inner_join(clin_data, mirna_top, by = "sample")


# Univariate Cox on top variable miRNAs and genes

cox_results <- function(features, data, type) {
  results <- lapply(features, function(f) {
    mod <- coxph(Surv(OS.time, OS) ~ data[[f]], data = data)
    tibble(
      feature = f,
      HR = exp(coef(mod)),
      p = summary(mod)$coef[,"Pr(>|z|)"],
      type = type
    )
  })
  bind_rows(results)
}

cox_gene <- cox_results(top_gene_names, clin_data_predictive, type = "Gene")
cox_mirna <- cox_results(selected_mirnas, clin_data_mirna, type = "miRNA")

cox_combined <- bind_rows(cox_gene, cox_mirna) %>%
  mutate(
    HR = as.numeric(HR),  p = as.numeric(p),
    Direction = case_when(
      HR > 1 & p < 0.05 ~ "Risk",
      HR < 1 & p < 0.05 ~ "Protective",
      TRUE ~ "NS" ) ) %>%
  rename(    Feature = feature,
    `Hazard Ratio` = HR,  `P-value` = p,  Type = type)
  
# Save to CSV
write.csv(cox_combined, "results/survival_feature_summary.csv", row.names = FALSE)



