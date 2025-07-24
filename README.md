# TCGA Sarcoma Survival Analysis with Multi-Omics Integration
This project explores molecular subtypes in TCGA Sarcoma (SARC) patients using gene expression data and evaluates their association with survival outcomes. It also integrates microRNA expression, protein expression, and somatic mutation data to build a multi-omics survival model.

-  Identify molecular subtypes of sarcoma using unsupervised clustering of gene expression data
- Evaluate overall survival differences between identified subtypes using Kaplan-Meier and Cox models
- Integrate additional omics layers (miRNA, RPPA protein, somatic mutations) to enrich survival modeling
- Visualize results with dimensionality reduction, survival plots, and feature importance

### Step 1: Data Collection & Cleaning
- Downloaded omics datasets from UCSC Xena
- Matched expression and clinical data by sample IDs
- Filtered to primary tumor samples
### Step 2: Feature Selection & Clustering
- Select top N most variable genes (e.g. top 1000)
- Apply PCA and/or t-SNE for visualization
- Perform unsupervised clustering (e.g. KMeans, NMF) to define molecular subtypes
### Step 3: Survival Analysis
- Merge cluster labels with clinical data
- Performe Kaplan-Meier survival analysis by cluster
- Fit Cox Proportional Hazards models using:
  - Clusters
  - Top features (expression, miRNA, mutation load, etc.)
  - Optional multivariate model with age/sex
### Step 4: Multi-Omics Integration
- Add miRNA expression, mutation burden, and protein levels as model covariates
- Compare survival model performance across modalities
## Key Visualizations
- PCA/t-SNE plots of molecular clusters
- Kaplan-Meier survival curves stratified by cluster
- Forest plots of Cox model hazard ratios
- Heatmaps of top subtype-specific genes or miRNAs

## Research Questions
- Do molecular clusters of sarcoma patients differ in survival outcomes?
- Which genes, miRNAs, or proteins are most predictive of short vs. long-term survival?
- Does mutation burden explain variance in outcomes when combined with expression?




Author 
linkedin.com/in/christopher-lubega-aa5b5a116