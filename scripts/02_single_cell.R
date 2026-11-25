# NK-cell annotation and membrane-gene selection from GSE212966
# Download/arrange the public matrices before running this script (see README).

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(Matrix)
  library(dplyr)
})

annotation_file <- "data/annotations/Annotation_RNAseq.rda"
if (!file.exists(annotation_file)) stop("Missing annotation: ", annotation_file)
load(annotation_file)
if (!exists("annot") || !"GOComponent" %in% names(annot)) {
  stop("Annotation_RNAseq.rda must contain 'annot' with a GOComponent column")
}

# Gene symbols occupy the second column in the original annotation object.
membrane_genes <- unique(as.character(
  annot[[2]][grepl("plasma membrane", annot$GOComponent, ignore.case = TRUE)]
))
reference <- celldex::HumanPrimaryCellAtlasData()
dir.create("results", recursive = TRUE, showWarnings = FALSE)

read_samples <- function(tissue) {
  folder <- file.path("data", "scRNAseq", "GSE212966_RAW", tissue)
  sample_dirs <- sort(list.dirs(folder, recursive = FALSE, full.names = TRUE))
  if (length(sample_dirs) != 6) stop("Expected six sample folders in: ", folder)

  samples <- lapply(seq_along(sample_dirs), function(i) {
    path <- sample_dirs[i]
    genes_file <- if (file.exists(file.path(path, "genes.tsv.gz"))) {
      "genes.tsv.gz"
    } else {
      "features.tsv.gz"
    }
    counts <- Seurat::ReadMtx(
      mtx = file.path(path, "matrix.mtx.gz"),
      features = file.path(path, genes_file),
      cells = file.path(path, "barcodes.tsv.gz")
    )
    sample <- Seurat::CreateSeuratObject(counts, project = tissue)
    sample$sample <- paste0(tissue, i)
    sample
  })

  Seurat::merge(samples[[1]], y = samples[-1],
                add.cell.ids = paste0(tissue, seq_along(samples)))
}

process_tissue <- function(tissue) {
  obj <- read_samples(tissue)
  obj[["mitoPercent"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

  obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 &
                  nCount_RNA > 200 & nCount_RNA < 20000 & mitoPercent <= 10)
  obj <- NormalizeData(obj, scale.factor = 10000, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000,
                              verbose = FALSE)
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  obj <- RunPCA(obj, features = VariableFeatures(obj), verbose = FALSE)

  # RPCA provides the integrated space for clustering and visualization.
  obj <- IntegrateLayers(obj, method = RPCAIntegration,
                         orig.reduction = "pca", new.reduction = "integrated.rpca",
                         verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = 1:15,
                       verbose = FALSE)
  obj <- FindClusters(obj, cluster.name = "rpca_clusters", verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:15,
                 reduction.name = "umap.rpca", verbose = FALSE)
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  predictions <- SingleR(test = counts, ref = reference, labels = reference$label.fine)
  obj$singleR.labels <- predictions$labels[match(colnames(obj), rownames(predictions))]

  nk_cells <- colnames(obj)[grepl("NK", obj$singleR.labels)]
  if (!length(nk_cells)) stop("No NK cells annotated in: ", tissue)
  nk_counts <- counts[, nk_cells, drop = FALSE]
  expressed <- rownames(nk_counts)[Matrix::rowSums(nk_counts >= 5) >= 2]
  nk_membrane <- intersect(expressed, membrane_genes)

  sample_qc <- obj@meta.data %>%
    mutate(Tissue = tissue, is_NK = grepl("NK", singleR.labels)) %>%
    group_by(Tissue, sample) %>%
    summarise(Post_QC_cells = n(), NK_cells = sum(is_NK),
              NK_percent = 100 * mean(is_NK), .groups = "drop")

  write.csv(data.frame(Gene = nk_membrane),
            file.path("results", paste0("NK_", tissue, "_membrane_genes.csv")),
            row.names = FALSE)
  write.csv(sample_qc, file.path("results", paste0("NK_", tissue, "_QC.csv")),
            row.names = FALSE)
  saveRDS(obj, file.path("results", paste0("Seurat_", tissue, ".rds")))
  message(tissue, ": ", ncol(obj), " cells; ", length(nk_cells),
          " NK cells; ", length(nk_membrane), " membrane genes")
}

for (tissue in c("PDAC", "ADJ")) process_tissue(tissue)
