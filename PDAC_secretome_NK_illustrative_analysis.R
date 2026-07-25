# ==============================================================================
# PDAC SECRETOME–NK CELL INTEGRATIVE ANALYSIS
# Illustrative analysis script for public deposition
# ==============================================================================
#
# Purpose
# -------
# This script documents the analytical sequence used to integrate:
#
#   1. proteomic secretome profiles from patient-derived PDAC cultures;
#   2. exploratory differential abundance, GSEA and machine-learning analyses;
#   3. single-cell RNA-seq profiles of NK cells from PDAC and adjacent tissue;
#   4. theoretical EV-associated protein filtering;
#   5. high-confidence STRING interaction networks;
#   6. Cytoscape/CytoHubba topological prioritization.
#
# IMPORTANT
# ---------
# This is an ILLUSTRATIVE, NON-EXECUTABLE script prepared for methodological
# transparency. It is not intended to reproduce the analysis without adaptation.
#
# The code preserves the principal object names, thresholds, package choices and
# analytical order used in:
#
#   - Secretoma.Rmd
#   - Single_cell_Pancreas.Rmd
#
# Absolute local paths, exploratory code, unrelated analyses, repeated blocks,
# interactive viewers and intermediate workspace files were removed or replaced
# by placeholders.
#
# The original notebooks remain the authoritative record of the exploratory
# analysis. This file is a concise representation suitable for a GitHub
# repository accompanying the manuscript.
#
# ==============================================================================


# ==============================================================================
# 1. PACKAGES USED IN THE ORIGINAL WORKFLOW
# ==============================================================================

if (FALSE) {

  # Data manipulation and visualization
  library(tidyverse)
  library(readxl)
  library(openxlsx)
  library(broom)
  library(patchwork)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)

  # Statistical analysis and multivariate visualization
  library(pwr)
  library(FactoMineR)
  library(factoextra)
  library(cluster)

  # Functional enrichment
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(ReactomePA)
  library(gprofiler2)

  # Machine learning
  library(binda)
  library(caret)
  library(randomForest)
  library(e1071)
  library(xgboost)
  library(DALEX)
  library(iml)
  library(SHAPforxgboost)
  library(tidymodels)

  # Single-cell RNA-seq
  library(Seurat)
  library(SingleR)
  library(celldex)

  # Interaction networks
  library(httr)
  library(jsonlite)
  library(igraph)
  library(STRINGdb)
  library(RCy3)
}


# ==============================================================================
# 2. INPUT OBJECTS AND FILE ORGANIZATION
# ==============================================================================

if (FALSE) {

  # Local paths used in the original notebooks were replaced by generic paths.
  #
  # Expected inputs:
  #
  # data/
  # ├── proteomics/
  # │   ├── ID_of_PDAC_cells_JI.xlsx
  # │   ├── clinical_metadata.xlsx
  # │   └── Exocarta_pancreatic_cancer_cells.xlsx
  # ├── annotations/
  # │   └── Annotation_RNAseq.rda
  # ├── scRNAseq/
  # │   └── GSE212966_RAW/
  # │       ├── PDAC/
  # │       └── ADJ/
  # └── cytohubba/
  #     ├── Huba_PDCA_VesMet2.xlsx
  #     ├── Huba_PDCA_VesLoca.xlsx
  #     ├── Hubba_ADJ_Meta.xlsx
  #     └── Huba_ADJ_Loca.xlsx

  proteomics_file <- "data/proteomics/ID_of_PDAC_cells_JI.xlsx"
  metadata_file   <- "data/proteomics/clinical_metadata.xlsx"
  exocarta_file   <- "data/proteomics/Exocarta_pancreatic_cancer_cells.xlsx"

  pdac_path <- "data/scRNAseq/GSE212966_RAW/PDAC/"
  adj_path  <- "data/scRNAseq/GSE212966_RAW/ADJ/"

  annotation_file <- "data/annotations/Annotation_RNAseq.rda"
}


# ==============================================================================
# PART I. SECRETOME PROTEOMICS
# ==============================================================================


# ==============================================================================
# 3. IMPORT AND PREPARATION OF THE PROTEOMIC DATASET
# ==============================================================================

if (FALSE) {

  # The original dataset contained 4,304 proteins.
  # Proteins with at least 50% observed values were retained.
  # Remaining missing values were imputed using the protein mean.
  #
  # The resulting working matrix contained 2,604 proteins.

  Base <- read.xlsx(proteomics_file)

  # Clinical groups were derived from the metadata included with the dataset.
  MetaData <- read.xlsx(metadata_file)

  # Locally advanced cultures were excluded from the localized-versus-metastatic
  # comparison. The interactome analysis used localized and primary metastatic
  # cultures, as defined in the original notebook.

  MetaData_f <- MetaData %>%
    filter(!Disease_stage %in% "locally advanced") %>%
    filter(
      Disease_stage %in% c("localized", "metastatic")
    )

  # The original notebook renamed abundance columns using disease stage and
  # sample position.

  Base.Fix_1 <- Base
  names(Base.Fix_1) <- c(
    "Accession",
    paste(MetaData_f$Disease_stage, seq_len(nrow(MetaData_f)), sep = ".")
  )

  # Illustration of the missing-value procedure.
  # Exact column selections depended on the source spreadsheet structure.

  abundance_matrix <- Base.Fix_1[, -1]

  detected_fraction <- rowMeans(!is.na(abundance_matrix))

  abundance_matrix <- abundance_matrix[
    detected_fraction >= 0.50,
    ,
    drop = FALSE
  ]

  abundance_matrix <- t(
    apply(
      abundance_matrix,
      1,
      function(x) {
        x[is.na(x)] <- mean(x, na.rm = TRUE)
        x
      }
    )
  )

  Base.Fix <- data.frame(
    Accession = Base.Fix_1$Accession[
      detected_fraction >= 0.50
    ],
    abundance_matrix,
    check.names = FALSE
  )
}


# ==============================================================================
# 4. LONG-FORM DATA AND UNIVARIATE COMPARISON
# ==============================================================================

if (FALSE) {

  # The original R Markdown converted protein abundances to long format and
  # performed a two-group t.test for each protein.

  Base.long <- Base.Fix %>%
    pivot_longer(
      cols = -Accession,
      names_to = "Sample",
      values_to = "Log.Abunda"
    ) %>%
    mutate(
      Group = case_when(
        str_detect(Sample, regex("metastatic", ignore_case = TRUE)) ~
          "Metastatic",
        str_detect(Sample, regex("localized", ignore_case = TRUE)) ~
          "Localized",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Group))

  Estadistica <- Base.long %>%
    group_by(Accession) %>%
    nest() %>%
    mutate(
      test = map(
        data,
        ~ t.test(Log.Abunda ~ Group, data = .x)
      ),
      tidy_test = map(test, broom::tidy)
    ) %>%
    unnest(tidy_test)

  # Group means and metastatic-minus-localized abundance difference.

  Mean_groups <- Base.long %>%
    group_by(Accession, Group) %>%
    summarise(
      Mean = mean(Log.Abunda, na.rm = TRUE),
      SD = sd(Log.Abunda, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = Group,
      values_from = c(Mean, SD, n)
    ) %>%
    mutate(
      Dif.Abundance = Mean_Metastatic - Mean_Localized,
      Sd.pool = sqrt(
        (
          (n_Metastatic - 1) * SD_Metastatic^2 +
            (n_Localized - 1) * SD_Localized^2
        ) /
          (n_Metastatic + n_Localized - 2)
      ),
      Cohen = Dif.Abundance / Sd.pool
    )

  FULL.TO.1 <- Estadistica %>%
    select(Accession, p.value) %>%
    left_join(Mean_groups, by = "Accession")

  # The original notebook estimated post hoc power using pwr.t.test.

  power <- tidy(
    pwr.t.test(
      n = 10,
      d = abs(FULL.TO.1$Dif.Abundance) / FULL.TO.1$Sd.pool,
      sig.level = 0.05,
      type = "two.sample",
      alternative = "two.sided"
    )
  )

  # Candidate DEPs were selected using nominal p value and standardized
  # effect-size criteria in the exploratory workflow.

  Binding_down <- FULL.TO.1 %>%
    filter(
      p.value < 0.05,
      Cohen <= -1
    )

  Binding_up <- FULL.TO.1 %>%
    filter(
      p.value < 0.05,
      Cohen >= 1
    )

  Binding <- bind_rows(
    Binding_down %>% mutate(Direction = "Lower in metastatic"),
    Binding_up %>% mutate(Direction = "Higher in metastatic")
  )
}


# ==============================================================================
# 5. VOLCANO PLOT AND PCA / K-MEANS VISUALIZATION
# ==============================================================================

if (FALSE) {

  # Volcano plot based on Cohen's effect size and nominal p value.

  Volcano <- ggplot(
    FULL.TO.1,
    aes(
      x = Cohen,
      y = -log10(p.value)
    )
  ) +
    geom_point(alpha = 0.5) +
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed"
    ) +
    labs(
      title = "Metastatic versus localized PDAC secretomes",
      x = "Cohen standardized effect size",
      y = "-log10(p value)"
    ) +
    theme_classic()

  # PCA using the complete filtered protein matrix.

  PCA_input_all <- t(
    as.matrix(
      Base.Fix[, -1]
    )
  )

  PCA_all <- PCA(
    PCA_input_all,
    scale.unit = FALSE,
    graph = FALSE
  )

  PCA_all_plot <- fviz_pca_ind(
    PCA_all,
    habillage = MetaData_f$Disease_stage,
    addEllipses = TRUE,
    repel = TRUE
  )

  # PCA using selected proteins was used as a descriptive visualization of the
  # feature space, not as independent validation.

  DEP_accessions <- Binding$Accession

  PCA_input_DEP <- Base.Fix %>%
    filter(Accession %in% DEP_accessions) %>%
    column_to_rownames("Accession") %>%
    t()

  PCA_DEP <- PCA(
    PCA_input_DEP,
    scale.unit = FALSE,
    graph = FALSE
  )

  hk <- hkmeans(
    PCA_input_DEP,
    k = 2
  )

  PCA_DEP_plot <- fviz_cluster(
    hk,
    data = PCA_input_DEP,
    ellipse.type = "convex"
  ) +
    labs(
      title = "k-means clustering–PCA of selected proteins"
    )
}


# ==============================================================================
# 6. GSEA USING THE METASTATIC-MINUS-LOCALIZED RANKING
# ==============================================================================

if (FALSE) {

  # Protein accessions were mapped to gene symbols before enrichment analysis.
  # The original workflow used a ranked vector based on:
  #
  #   metastatic mean abundance - localized mean abundance

  GSEA_df <- FULL.TO.1 %>%
    left_join(
      protein_annotation,
      by = "Accession"
    ) %>%
    filter(
      !is.na(Gene),
      Gene != ""
    ) %>%
    distinct(
      Gene,
      .keep_all = TRUE
    )

  gene_list <- GSEA_df$Dif.Abundance
  names(gene_list) <- GSEA_df$Gene
  gene_list <- sort(
    gene_list,
    decreasing = TRUE
  )

  # GO Biological Process

  gse_BP <- gseGO(
    geneList = gene_list,
    ont = "BP",
    keyType = "SYMBOL",
    OrgDb = org.Hs.eg.db,
    minGSSize = 3,
    maxGSSize = 800,
    pvalueCutoff = 0.20,
    pAdjustMethod = "BH",
    verbose = FALSE
  )

  # GO Molecular Function

  gse_MF <- gseGO(
    geneList = gene_list,
    ont = "MF",
    keyType = "SYMBOL",
    OrgDb = org.Hs.eg.db,
    minGSSize = 3,
    maxGSSize = 800,
    pvalueCutoff = 0.20,
    pAdjustMethod = "BH",
    verbose = FALSE
  )

  # GO Cellular Component

  gse_CC <- gseGO(
    geneList = gene_list,
    ont = "CC",
    keyType = "SYMBOL",
    OrgDb = org.Hs.eg.db,
    minGSSize = 3,
    maxGSSize = 800,
    pvalueCutoff = 0.20,
    pAdjustMethod = "BH",
    verbose = FALSE
  )

  # KEGG and Reactome were also evaluated in the original notebook.

  kk2 <- gseKEGG(
    geneList = gene_list,
    organism = "hsa",
    pvalueCutoff = 0.20,
    pAdjustMethod = "BH",
    verbose = FALSE
  )

  REACT_gse <- gsePathway(
    geneList = gene_list,
    organism = "human",
    pvalueCutoff = 0.20,
    pAdjustMethod = "BH",
    verbose = FALSE
  )

  # Additional exploratory resources included Disease Ontology, DisGeNET,
  # Network of Cancer Genes, cell/tissue markers, oncogenic signatures and
  # immunologic signatures.
}


# ==============================================================================
# PART II. EXPLORATORY MACHINE LEARNING
# ==============================================================================


# ==============================================================================
# 7. TRAIN / TEST PARTITION
# ==============================================================================

if (FALSE) {

  # The labeled localized/metastatic dataset contained 30 samples.
  # A stratified random split generated:
  #
  #   training set: n = 23
  #   test set:     n = 7
  #
  # This procedure was exploratory and was not a nested or external validation.

  featureMatrix <- t(
    as.matrix(
      Base.Fix[, -1]
    )
  )

  Labels <- factor(
    MetaData_f$Disease_stage,
    levels = c("localized", "metastatic")
  )

  set.seed(123)

  train_index <- createDataPartition(
    Labels,
    p = 0.80,
    list = FALSE
  )

  featureMatrix_Train <- featureMatrix[
    train_index,
    ,
    drop = FALSE
  ]

  featureMatrix_TestUnb <- featureMatrix[
    -train_index,
    ,
    drop = FALSE
  ]

  Labels_train <- Labels[train_index]
  Labels_test  <- Labels[-train_index]
}


# ==============================================================================
# 8. BINARY DISCRIMINANT ANALYSIS FEATURE SELECTION
# ==============================================================================

if (FALSE) {

  # The original notebook used a custom function:
  #
  #   BDA.analysis_Train_test_New()
  #
  # It repeatedly estimated dichotomization thresholds, assessed threshold
  # stability and ranked features using the binda package.

  BDA.analysis_Train_test_New(
    feature.matrix = featureMatrix_Train,
    rowName = Labels_train,
    seed.sr = 123,
    Prob_Unbiased = 0.80
  )

  load("Matrices.BDA.rda")

  # Stable BDA-ranked features were retained.
  # The reported workflow yielded 105 candidate proteins.

  BDA_best <- BDA.Rank.f1_stable %>%
    arrange(desc(abs(t.score))) %>%
    distinct(Picos, .keep_all = TRUE)

  BDA_b <- unique(
    as.character(
      BDA_best$Picos
    )
  )
}


# ==============================================================================
# 9. RANDOM-FOREST FEATURE REFINEMENT
# ==============================================================================

if (FALSE) {

  # Random forest was fitted to the 105 BDA-selected features.

  TRAIN_BDA <- as.data.frame(
    featureMatrix_Train[, BDA_b, drop = FALSE]
  )

  TRAIN_BDA$Labels_train <- Labels_train

  rf_model <- randomForest(
    Labels_train ~ .,
    data = TRAIN_BDA,
    ntree = 100,
    importance = TRUE
  )

  rf_importance <- importance(
    rf_model
  ) %>%
    as.data.frame() %>%
    rownames_to_column("Picos") %>%
    arrange(desc(MeanDecreaseGini))

  # Two RF-derived subsets were evaluated:
  #
  #   - the highest-ranked 4 proteins;
  #   - 95 proteins after removing the lowest-ranked features.

  RF_best_4 <- rf_importance$Picos[1:4]
  RF_best_95 <- rf_importance$Picos[1:95]
}


# ==============================================================================
# 10. PREDICTIVE MODELS
# ==============================================================================

if (FALSE) {

  # Model 1: random forest using 105 BDA-selected features.

  results_RF_105 <- Cross.Val.Caret.Rf_tune_BDA_p2(
    TRAIN = featureMatrix_Train[, BDA_b, drop = FALSE],
    TEST = featureMatrix_TestUnb[, BDA_b, drop = FALSE],
    rowNameTrain = Labels_train,
    rowNameTest = Labels_test
  )

  # Model 2: random forest using the highest-ranked RF subset.

  results_RF_4 <- Cross.Val.Caret.Rf_tune_BDA_p2(
    TRAIN = featureMatrix_Train[, RF_best_4, drop = FALSE],
    TEST = featureMatrix_TestUnb[, RF_best_4, drop = FALSE],
    rowNameTrain = Labels_train,
    rowNameTest = Labels_test
  )

  # Model 3: random forest using the 95-protein subset.

  results_RF_95 <- Cross.Val.Caret.Rf_tune_BDA_p2(
    TRAIN = featureMatrix_Train[, RF_best_95, drop = FALSE],
    TEST = featureMatrix_TestUnb[, RF_best_95, drop = FALSE],
    rowNameTrain = Labels_train,
    rowNameTest = Labels_test
  )

  # Linear, radial and polynomial SVMs were evaluated on the 95-protein subset.

  results_SVM <- Cross.Val.Caret.SVM_BDA_paralell2(
    TRAIN = featureMatrix_Train[, RF_best_95, drop = FALSE],
    TEST = featureMatrix_TestUnb[, RF_best_95, drop = FALSE],
    rowNameTrain = Labels_train,
    rowNameTest = Labels_test
  )

  # XGBoost was also evaluated on the same feature subset.

  results_XGB <- Cross.Val.Caret.XGB_BDA(
    TRAIN = featureMatrix_Train[, RF_best_95, drop = FALSE],
    TEST = featureMatrix_TestUnb[, RF_best_95, drop = FALSE],
    rowNameTrain = Labels_train,
    rowNameTest = Labels_test
  )
}


# ==============================================================================
# 11. SHAP-BASED MODEL INTERPRETATION
# ==============================================================================

if (FALSE) {

  # The original workflow used DALEX to calculate SHAP-style contributions for
  # random-forest models.

  explainer_rf <- explain(
    rf_model,
    data = featureMatrix_Train[, BDA_b, drop = FALSE],
    y = Labels_train
  )

  shap_rf <- predict_parts(
    explainer_rf,
    new_observation = featureMatrix_TestUnb[, BDA_b, drop = FALSE],
    type = "shap"
  )

  plot(
    shap_rf
  )

  # The machine-learning component was used for exploratory feature
  # prioritization. Performance estimates based on the seven-sample test set
  # were not treated as external clinical validation.
}


# ==============================================================================
# 12. OVER-REPRESENTATION ANALYSIS OF THE 95-PROTEIN SET
# ==============================================================================

if (FALSE) {

  BindBDARF <- data.frame(
    Best = RF_best_95
  ) %>%
    separate(
      Best,
      into = c("Gene", "Other"),
      sep = "\\.",
      extra = "merge",
      fill = "right"
    )

  Enri_GO_BP <- enrichGO(
    gene = BindBDARF$Gene,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )

  Enri_GO_MF <- enrichGO(
    gene = BindBDARF$Gene,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )
}


# ==============================================================================
# PART III. SINGLE-CELL RNA-SEQ ANALYSIS
# ==============================================================================


# ==============================================================================
# 13. IMPORT OF GSE212966 PDAC SAMPLES
# ==============================================================================

if (FALSE) {

  # Six PDAC samples were imported from filtered feature-barcode matrices.

  PaTH <- pdac_path

  dirs <- list.dirs(
    path = PaTH,
    recursive = FALSE,
    full.names = FALSE
  )

  for (x in dirs) {

    name <- gsub(
      "_filtered_feature_bc_matrix",
      "",
      x
    )

    cts <- ReadMtx(
      mtx = paste0(PaTH, x, "/matrix.mtx.gz"),
      features = paste0(PaTH, x, "/genes.tsv.gz"),
      cells = paste0(PaTH, x, "/barcodes.tsv.gz")
    )

    assign(
      name,
      CreateSeuratObject(
        counts = cts
      )
    )
  }

  merged_seurat_PBMC <- merge(
    PDAC1,
    y = c(
      PDAC2,
      PDAC3,
      PDAC4,
      PDAC5,
      PDAC6
    ),
    add.cell.ids = c(
      "PDAC1",
      "PDAC2",
      "PDAC3",
      "PDAC4",
      "PDAC5",
      "PDAC6"
    ),
    project = "Dani"
  )
}


# ==============================================================================
# 14. PDAC SINGLE-CELL QUALITY CONTROL AND NORMALIZATION
# ==============================================================================

if (FALSE) {

  merged_seurat_PBMC$sample <- rownames(
    merged_seurat_PBMC@meta.data
  )

  merged_seurat_PBMC@meta.data <- separate(
    merged_seurat_PBMC@meta.data,
    col = "sample",
    into = c(
      "Type",
      "Patient",
      "Barcode"
    ),
    sep = "_"
  )

  merged_seurat_PBMC$mitoPercent <- PercentageFeatureSet(
    merged_seurat_PBMC,
    pattern = "^MT-"
  )

  merged_seurat_PBMC_filt <- subset(
    merged_seurat_PBMC,
    subset =
      nFeature_RNA > 200 &
        nFeature_RNA < 5000 &
        nCount_RNA > 200 &
        nCount_RNA < 20000 &
        mitoPercent < 10
  )

  merged_seurat_PBMC_filt <- NormalizeData(
    merged_seurat_PBMC_filt
  )

  merged_seurat_PBMC_filt <- FindVariableFeatures(
    merged_seurat_PBMC_filt,
    selection.method = "vst",
    nfeatures = 2000
  )

  Variable_genes <- VariableFeatures(
    merged_seurat_PBMC_filt
  )

  merged_seurat_PBMC_filt <- ScaleData(
    object = merged_seurat_PBMC_filt,
    features = Variable_genes
  )

  merged_seurat_PBMC_filt <- RunPCA(
    merged_seurat_PBMC_filt,
    features = VariableFeatures(
      merged_seurat_PBMC_filt
    )
  )

  merged_seurat_PBMC_filt <- FindNeighbors(
    merged_seurat_PBMC_filt,
    dims = 1:15
  )

  merged_seurat_PBMC_filt <- FindClusters(
    merged_seurat_PBMC_filt,
    resolution = 0.5
  )

  merged_seurat_PBMC_filt <- RunUMAP(
    merged_seurat_PBMC_filt,
    dims = 1:15
  )
}


# ==============================================================================
# 15. RPCA INTEGRATION OF PDAC SAMPLES
# ==============================================================================

if (FALSE) {

  merged_seurat_integrated_PBMC <- IntegrateLayers(
    object = merged_seurat_PBMC_filt,
    method = RPCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.rpca",
    verbose = FALSE
  )

  merged_seurat_integrated_PBMC <- FindNeighbors(
    merged_seurat_integrated_PBMC,
    reduction = "integrated.rpca",
    dims = 1:15
  )

  merged_seurat_integrated_PBMC <- FindClusters(
    merged_seurat_integrated_PBMC,
    cluster.name = "rpca_clusters"
  )

  merged_seurat_integrated_PBMC <- RunUMAP(
    merged_seurat_integrated_PBMC,
    reduction = "integrated.rpca",
    dims = 1:15,
    reduction.name = "umap.rpca"
  )

  merged_seurat_integrated_PBMC[["RNA"]] <- JoinLayers(
    merged_seurat_integrated_PBMC[["RNA"]]
  )
}


# ==============================================================================
# 16. SINGLER ANNOTATION AND NK-CELL SELECTION IN PDAC
# ==============================================================================

if (FALSE) {

  hpca <- HumanPrimaryCellAtlasData()

  Seurat_counts <- GetAssayData(
    object = merged_seurat_integrated_PBMC,
    assay = "RNA",
    slot = "counts"
  )

  pred <- SingleR(
    test = Seurat_counts,
    ref = hpca,
    labels = hpca$label.fine
  )

  merged_seurat_integrated_PBMC$singleR.labels <-
    pred$labels[
      match(
        rownames(
          merged_seurat_integrated_PBMC@meta.data
        ),
        rownames(pred)
      )
    ]

  Annotations <- data.frame(
    Labels = unique(
      hpca$label.fine
    )
  ) %>%
    mutate(
      NK_cells = str_detect(
        Labels,
        "NK"
      )
    )

  Annotations_NK_cells <- Annotations %>%
    filter(
      NK_cells == TRUE
    )

  merged_seurat_filt_PBMC_NK <- subset(
    merged_seurat_integrated_PBMC,
    subset =
      singleR.labels %in%
        Annotations_NK_cells$Labels
  )
}


# ==============================================================================
# 17. EXPRESSED PLASMA-MEMBRANE GENES IN PDAC NK CELLS
# ==============================================================================

if (FALSE) {

  counts_matrix_counts_PDAC <-
    merged_seurat_filt_PBMC_NK@assays$RNA$counts

  Counts_Df_NK <- data.frame(
    counts_matrix_counts_PDAC
  )

  # Genes were retained when they had at least 5 counts in at least 2 NK cells.

  keep <- rowSums(
    Counts_Df_NK >= 5
  ) >= 2

  Counts_Df_NK <- Counts_Df_NK[
    keep,
    ,
    drop = FALSE
  ]

  Genes_PDAC_NK <- data.frame(
    Genes = rownames(
      Counts_Df_NK
    )
  )

  load(
    annotation_file
  )

  annot1 <- annot[, c(2, 18)]

  Prot_MP <- annot1 %>%
    mutate(
      MP = str_detect(
        GOComponent,
        "plasma membrane"
      )
    ) %>%
    filter(
      MP == TRUE
    )

  names(Prot_MP)[1] <- "Genes"

  NK_PDAC_MP <- Genes_PDAC_NK %>%
    semi_join(
      Prot_MP,
      by = "Genes"
    )

  save(
    NK_PDAC_MP,
    Counts_Df_NK,
    file = "PDAC_NK_data.rda"
  )
}


# ==============================================================================
# 18. ADJACENT-TISSUE SINGLE-CELL ANALYSIS
# ==============================================================================

if (FALSE) {

  # The adjacent-tissue workflow repeated the same sequence used for PDAC:
  #
  #   ReadMtx -> CreateSeuratObject -> merge six samples -> QC ->
  #   NormalizeData -> FindVariableFeatures -> ScaleData -> RunPCA ->
  #   FindNeighbors -> FindClusters -> RunUMAP -> RPCAIntegration ->
  #   SingleR annotation -> NK-cell subset.

  PaTH <- adj_path

  dirs <- list.dirs(
    path = PaTH,
    recursive = FALSE,
    full.names = FALSE
  )

  for (x in dirs) {

    name <- gsub(
      "_filtered_feature_bc_matrix",
      "",
      x
    )

    cts <- ReadMtx(
      mtx = paste0(PaTH, x, "/matrix.mtx.gz"),
      features = paste0(PaTH, x, "/genes.tsv.gz"),
      cells = paste0(PaTH, x, "/barcodes.tsv.gz")
    )

    assign(
      name,
      CreateSeuratObject(
        counts = cts
      )
    )
  }

  merged_seurat_PFMC <- merge(
    ADJ1,
    y = c(
      ADJ2,
      ADJ3,
      ADJ4,
      ADJ5,
      ADJ6
    ),
    add.cell.ids = c(
      "ADJ1",
      "ADJ2",
      "ADJ3",
      "ADJ4",
      "ADJ5",
      "ADJ6"
    ),
    project = "Dani"
  )

  merged_seurat_PFMC$mitoPercent <- PercentageFeatureSet(
    merged_seurat_PFMC,
    pattern = "^MT-"
  )

  merged_seurat_PFMC_filt <- subset(
    merged_seurat_PFMC,
    subset =
      nFeature_RNA > 200 &
        nFeature_RNA < 5000 &
        nCount_RNA > 200 &
        nCount_RNA < 20000 &
        mitoPercent < 10
  )

  merged_seurat_PFMC_filt <- NormalizeData(
    merged_seurat_PFMC_filt
  )

  merged_seurat_PFMC_filt <- FindVariableFeatures(
    merged_seurat_PFMC_filt,
    selection.method = "vst",
    nfeatures = 2000
  )

  Variable_genes <- VariableFeatures(
    merged_seurat_PFMC_filt
  )

  merged_seurat_PFMC_filt <- ScaleData(
    merged_seurat_PFMC_filt,
    features = Variable_genes
  )

  merged_seurat_PFMC_filt <- RunPCA(
    merged_seurat_PFMC_filt,
    features = VariableFeatures(
      merged_seurat_PFMC_filt
    )
  )

  merged_seurat_integrated_PFMC <- IntegrateLayers(
    object = merged_seurat_PFMC_filt,
    method = RPCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.rpca",
    verbose = FALSE
  )

  merged_seurat_integrated_PFMC <- FindNeighbors(
    merged_seurat_integrated_PFMC,
    reduction = "integrated.rpca",
    dims = 1:15
  )

  merged_seurat_integrated_PFMC <- FindClusters(
    merged_seurat_integrated_PFMC,
    cluster.name = "rpca_clusters"
  )

  merged_seurat_integrated_PFMC <- RunUMAP(
    merged_seurat_integrated_PFMC,
    reduction = "integrated.rpca",
    dims = 1:15,
    reduction.name = "umap.rpca"
  )

  merged_seurat_integrated_PFMC[["RNA"]] <- JoinLayers(
    merged_seurat_integrated_PFMC[["RNA"]]
  )

  Seurat_counts <- GetAssayData(
    object = merged_seurat_integrated_PFMC,
    assay = "RNA",
    slot = "counts"
  )

  pred <- SingleR(
    test = Seurat_counts,
    ref = hpca,
    labels = hpca$label.fine
  )

  merged_seurat_integrated_PFMC$singleR.labels <-
    pred$labels[
      match(
        rownames(
          merged_seurat_integrated_PFMC@meta.data
        ),
        rownames(pred)
      )
    ]

  merged_seurat_filt_PFMC_NK <- subset(
    merged_seurat_integrated_PFMC,
    subset =
      singleR.labels %in%
        Annotations_NK_cells$Labels
  )

  counts_matrix_counts_ADJ <-
    merged_seurat_filt_PFMC_NK@assays$RNA$counts

  Counts_Df_NK_ADJ <- data.frame(
    counts_matrix_counts_ADJ
  )

  keep <- rowSums(
    Counts_Df_NK_ADJ >= 5
  ) >= 2

  Counts_Df_NK_ADJ <- Counts_Df_NK_ADJ[
    keep,
    ,
    drop = FALSE
  ]

  Genes_ADJ_NK <- data.frame(
    Genes = rownames(
      Counts_Df_NK_ADJ
    )
  )

  NK_ADJ_MP <- Genes_ADJ_NK %>%
    semi_join(
      Prot_MP,
      by = "Genes"
    )

  save(
    NK_ADJ_MP,
    Counts_Df_NK_ADJ,
    Prot_MP,
    file = "ADJ_NK_data.rda"
  )
}


# ==============================================================================
# PART IV. THEORETICAL EV-ASSOCIATED SECRETOME PROTEINS
# ==============================================================================


# ==============================================================================
# 19. EXOCARTA FILTER
# ==============================================================================

if (FALSE) {

  Ves_Pancreas <- read.xlsx(
    exocarta_file
  )

  # In the source file, the first row contained the effective column names.

  names(Ves_Pancreas) <- Ves_Pancreas[1, ]
  Ves_Pancreas <- Ves_Pancreas[-1, ]

  Ves_Pancreas <- Ves_Pancreas %>%
    filter(
      Species == "Homo sapiens"
    )

  Ves_Pancreas_genes <- Ves_Pancreas[, c(3, 4)]

  names(
    Ves_Pancreas_genes
  )[1] <- "Gene"
}


# ==============================================================================
# 20. DIRECTIONAL SECRETOME SETS AND EV / MEMBRANE INTERSECTION
# ==============================================================================

if (FALSE) {

  load(
    "Statis_Meta_loca.rda"
  )

  load(
    "PDAC_NK_data.rda"
  )

  load(
    "ADJ_NK_data.rda"
  )

  names(
    Prot_MP
  )[1] <- "Gene"

  # Proteins with positive metastatic-minus-localized differences were assigned
  # to the metastatic secretome set.

  Secretoma_Metas <- GSEA_df %>%
    filter(
      Diff.Abun > 0
    )

  Ves_Secre_Metas <- Secretoma_Metas %>%
    semi_join(
      Ves_Pancreas_genes,
      by = "Gene"
    ) %>%
    semi_join(
      Prot_MP,
      by = "Gene"
    )

  # Proteins with negative differences were assigned to the localized set.

  Secretoma_Loca <- GSEA_df %>%
    filter(
      Diff.Abun < 0
    )

  Ves_Secre_Loca <- Secretoma_Loca %>%
    semi_join(
      Ves_Pancreas_genes,
      by = "Gene"
    ) %>%
    semi_join(
      Prot_MP,
      by = "Gene"
    )

  # These sets represent theoretical EV-associated plasma-membrane candidates.
  # The analysis did not directly demonstrate their localization within EVs.
}


# ==============================================================================
# PART V. STRING INTERACTOME ANALYSIS
# ==============================================================================


# ==============================================================================
# 21. EXAMPLE: METASTATIC EV CANDIDATES × PDAC NK CELLS
# ==============================================================================

if (FALSE) {

  Prot_Ves_Secre_Metas <- Ves_Secre_Metas$Gene
  Prot_NK_PDAC_MP <- NK_PDAC_MP$Genes

  all_genes <- unique(
    c(
      Prot_NK_PDAC_MP,
      Prot_Ves_Secre_Metas
    )
  )

  species_id <- 9606

  gene_list <- paste(
    all_genes,
    collapse = "%0d"
  )

  # STRING REST API query with required_score = 900, equivalent to a combined
  # confidence score of at least 0.9.

  res <- GET(
    paste0(
      "https://string-db.org/api/json/network?identifiers=",
      gene_list,
      "&species=",
      species_id,
      "&required_score=900"
    )
  )

  ppi_data <- fromJSON(
    rawToChar(
      res$content
    )
  )

  # Only interactions with high experimental or curated database evidence were
  # retained.

  ppi_physical <- ppi_data %>%
    filter(
      escore >= 0.9 |
        dscore >= 0.9
    ) %>%
    select(
      preferredName_A,
      preferredName_B,
      escore,
      dscore
    )

  # Only intergroup edges were retained.

  ppi_inter_set <- ppi_physical %>%
    filter(
      (
        preferredName_A %in% Prot_NK_PDAC_MP &
          preferredName_B %in% Prot_Ves_Secre_Metas
      ) |
        (
          preferredName_B %in% Prot_NK_PDAC_MP &
            preferredName_A %in% Prot_Ves_Secre_Metas
        )
    )

  nodos <- data.frame(
    id = unique(
      c(
        ppi_inter_set$preferredName_A,
        ppi_inter_set$preferredName_B
      )
    )
  ) %>%
    mutate(
      group = case_when(
        id %in% Prot_NK_PDAC_MP ~ "NK_PDAC",
        id %in% Prot_Ves_Secre_Metas ~ "VES_metas",
        TRUE ~ "Other"
      )
    )

  edges <- data.frame(
    source = ppi_inter_set$preferredName_A,
    target = ppi_inter_set$preferredName_B,
    score =
      ppi_inter_set$escore +
      ppi_inter_set$dscore
  )

  createNetworkFromDataFrames(
    nodes = nodos,
    edges = edges,
    title = "CM_PPI_INTERSET",
    collection = "Secretion_Interactome2"
  )

  layoutNetwork(
    "force-directed"
  )

  loadTableData(
    nodos,
    data.key.column = "id"
  )
}


# ==============================================================================
# 22. FOUR INTERACTION NETWORKS
# ==============================================================================

if (FALSE) {

  # The preceding STRING/Cytoscape procedure was repeated for four combinations:

  network_sets <- list(

    SM_NK_PDAC = list(
      secretome = Ves_Secre_Metas$Gene,
      nk = NK_PDAC_MP$Genes
    ),

    SL_NK_PDAC = list(
      secretome = Ves_Secre_Loca$Gene,
      nk = NK_PDAC_MP$Genes
    ),

    SM_NK_ADJ = list(
      secretome = Ves_Secre_Metas$Gene,
      nk = NK_ADJ_MP$Genes
    ),

    SL_NK_ADJ = list(
      secretome = Ves_Secre_Loca$Gene,
      nk = NK_ADJ_MP$Genes
    )
  )

  # In the original notebook, the four analyses were written as separate code
  # blocks rather than wrapped in a common function. They are summarized here
  # as a list only to document the design.
}


# ==============================================================================
# PART VI. CYTOHUBBA TOPOLOGICAL ANALYSIS
# ==============================================================================


# ==============================================================================
# 23. IMPORT OF CYTOHUBBA RESULTS
# ==============================================================================

if (FALSE) {

  # CytoHubba calculations were performed in Cytoscape and exported to Excel.
  # The corresponding file was then imported into R.

  hub_data <- read.xlsx(
    "data/cytohubba/Huba_PDCA_VesMet2.xlsx",
    colNames = TRUE
  )

  hubba <- hub_data

  metric_cols <- colnames(
    hubba
  )[-1]

  hubba[metric_cols] <- lapply(
    hubba[metric_cols],
    as.numeric
  )

  hubba <- hubba %>%
    rename(
      node_name = 1
    ) %>%
    mutate(
      grupo = case_when(
        node_name %in% Prot_NK_PDAC_MP ~ "NK_PDAC",
        node_name %in% Prot_Ves_Secre_Metas ~ "VES_metas",
        TRUE ~ "Unknown"
      )
    ) %>%
    filter(
      grupo != "Unknown"
    ) %>%
    mutate(
      node_label = paste0(
        node_name,
        " [",
        grupo,
        "]"
      )
    )
}


# ==============================================================================
# 24. TOPOLOGICAL CATEGORIES
# ==============================================================================

if (FALSE) {

  categorias_metricas <- list(

    "Local Hub" = c(
      "MCC",
      "MNC",
      "DMNC"
    ),

    "Bridge" = c(
      "Betweenness",
      "BottleNeck",
      "Stress"
    ),

    "Global Influencer" = c(
      "Closeness",
      "Radiality"
    )
  )

  get_top_nodes <- function(
      df,
      metrics,
      n = 10
  ) {

    unique(
      unlist(
        lapply(
          metrics,
          function(m) {

            df %>%
              arrange(
                desc(
                  .data[[m]]
                )
              ) %>%
              slice(
                1:n
              ) %>%
              pull(
                1
              )
          }
        )
      )
    )
  }
}


# ==============================================================================
# 25. Z-SCORE HEATMAPS OF CENTRALITY METRICS
# ==============================================================================

if (FALSE) {

  heatmap_list <- list()
  dotplot_list <- list()

  for (
    categoria in names(
      categorias_metricas
    )
  ) {

    met <- categorias_metricas[
      [categoria]
    ]

    top_prots <- get_top_nodes(
      hubba,
      met,
      n = 10
    )

    sub <- hubba %>%
      filter(
        node_name %in% top_prots
      ) %>%
      mutate(
        across(
          all_of(met),
          scale
        )
      ) %>%
      distinct(
        node_name,
        .keep_all = TRUE
      )

    sub_heat <- sub %>%
      select(
        node_label,
        all_of(met)
      ) %>%
      column_to_rownames(
        "node_label"
      )

    pheatmap_obj <- pheatmap(
      sub_heat,
      color = colorRampPalette(
        c(
          "steelblue",
          "white",
          "tomato"
        )
      )(100),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      fontsize_row = 10,
      main = paste(
        "Z-score centralities -",
        categoria
      )
    )

    heatmap_list[
      [categoria]
    ] <- grid::grid.grab()

    sub_long <- sub %>%
      select(
        node_label,
        grupo,
        all_of(met)
      ) %>%
      pivot_longer(
        cols = all_of(met),
        names_to = "metrica",
        values_to = "valor"
      )

    gg <- ggplot(
      sub_long,
      aes(
        x = reorder(
          node_label,
          -valor
        ),
        y = valor,
        color = grupo
      )
    ) +
      geom_point(
        size = 3
      ) +
      facet_wrap(
        ~ metrica,
        scales = "free_y"
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(
          angle = 90,
          hjust = 1
        )
      ) +
      labs(
        title = paste(
          "Dot plot -",
          categoria
        ),
        x = "Protein",
        y = "Z-score",
        color = "Group"
      )

    dotplot_list[
      [categoria]
    ] <- gg
  }
}


# ==============================================================================
# 26. REPEATED CYTOHUBBA IMPORTS
# ==============================================================================

if (FALSE) {

  # The same topological procedure was applied to the four exported files:

  cytohubba_files <- c(
    SM_NK_PDAC = "data/cytohubba/Huba_PDCA_VesMet2.xlsx",
    SL_NK_PDAC = "data/cytohubba/Huba_PDCA_VesLoca.xlsx",
    SM_NK_ADJ  = "data/cytohubba/Hubba_ADJ_Meta.xlsx",
    SL_NK_ADJ  = "data/cytohubba/Huba_ADJ_Loca.xlsx"
  )

  # Each file was annotated with the corresponding secretome and NK source
  # vectors before top-10 selection and Z-score visualization.
}


# ==============================================================================
# 27. OUTPUTS REPORTED IN THE STUDY
# ==============================================================================

if (FALSE) {

  # Main outputs:
  #
  # - filtered secretome protein matrix;
  # - univariate protein-level statistics;
  # - volcano plot;
  # - all-protein and selected-protein PCA visualizations;
  # - GSEA results;
  # - exploratory ML feature rankings and model summaries;
  # - NK-PDAC and NK-ADJ plasma-membrane gene sets;
  # - EV-SM and EV-SL candidate sets;
  # - four high-confidence intergroup STRING networks;
  # - four Cytoscape network visualizations;
  # - CytoHubba local-hub, bridge and global-influencer heatmaps;
  # - comparative table of prioritized proteins.
}


# ==============================================================================
# 28. INTERPRETIVE LIMITS
# ==============================================================================

# The following limitations are integral to interpretation:
#
# 1. The EV-associated protein sets were inferred through ExoCarta and
#    GO Cellular Component annotation; EV localization was not directly
#    demonstrated in the analyzed secretomes.
#
# 2. STRING associations do not necessarily represent direct binding at an
#    EV–NK interface.
#
# 3. NK-cell transcript detection was used as a proxy for potential protein
#    availability.
#
# 4. CytoHubba centrality identifies topological priority, not causal or
#    functional importance.
#
# 5. The machine-learning analysis was exploratory because of the small
#    training and test cohorts.
#
# 6. PRKACA-centered EV transfer and NK-cell reprogramming remain hypotheses
#    requiring experimental validation.
#
# ==============================================================================

