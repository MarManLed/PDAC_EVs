# Proteomic preprocessing, exploratory statistics and GSEA
# Run from the repository root: Rscript scripts/01_proteomics.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(openxlsx)
})

input <- "data/proteomics/ID of PDAC cells_JI.xlsx"
out_dir <- "results"
if (!file.exists(input)) stop("Missing proteomic workbook: ", input)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Read the original Excel layout and retain sample measurements and identifiers.
x <- openxlsx::read.xlsx(input, colNames = FALSE, skipEmptyRows = FALSE)
raw <- x[c(6, 9:nrow(x)), c(2:42, 44:46)]
names(raw) <- c(as.character(unlist(raw[1, 1:41])), "Accession", "Name", "Gene")
raw <- raw[-1, , drop = FALSE]
raw <- raw[, !is.na(names(raw)) & !grepl("^NA(\\.|$)", names(raw)), drop = FALSE]

identifiers <- raw[, c("Accession", "Gene", "Name")]
measurements <- raw[, !names(raw) %in% names(identifiers), drop = FALSE]
measurements[] <- lapply(measurements, function(x) as.numeric(gsub(",", ".", as.character(x), fixed = TRUE)))

# Group assignment follows the original workbook column order.
groups <- c(
  "metastatic", "locally_advanced", "metastatic", "metastatic",
  "localized", "localized", rep("metastatic", 4),
  "metastatic", "locally_advanced", rep("metastatic", 2),
  "locally_advanced", "localized", rep("metastatic", 2),
  rep("localized", 4), "metastatic", "localized", "metastatic",
  "localized", "localized", "metastatic", "metastatic",
  "locally_advanced", rep("localized", 4)
)
if (ncol(measurements) != length(groups)) stop("Unexpected number of sample columns")
names(measurements) <- paste0(groups, ".", ave(seq_along(groups), groups, FUN = seq_along))

proteins <- bind_cols(identifiers, measurements) %>%
  pivot_longer(-c(Accession, Gene, Name), names_to = "Sample", values_to = "Abundance") %>%
  mutate(Group = sub("\\.[0-9]+$", "", Sample)) %>%
  filter(Group %in% c("localized", "metastatic"), !is.na(Accession))

# Apply the original missingness limits, then impute within clinical condition.
missingness <- proteins %>%
  group_by(Accession, Group) %>%
  summarise(n_missing = sum(is.na(Abundance)), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = n_missing, values_fill = 0)

retained <- missingness %>%
  filter(localized <= 5, metastatic <= 5, localized + metastatic <= 7) %>%
  pull(Accession)

proteins <- proteins %>%
  filter(Accession %in% retained) %>%
  group_by(Accession, Group) %>%
  mutate(Abundance_imputed = replace(Abundance, is.na(Abundance), mean(Abundance, na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(Log_abundance = log(Abundance_imputed))

if (any(!is.finite(proteins$Log_abundance))) stop("Non-finite log abundance; inspect the input matrix")

# Two-group comparisons; nominal candidates are not FDR-significant discoveries.
protein_stats <- proteins %>%
  group_by(Accession, Gene, Name) %>%
  group_modify(~ {
    localized <- .x$Log_abundance[.x$Group == "localized"]
    metastatic <- .x$Log_abundance[.x$Group == "metastatic"]
    pooled_sd <- sqrt(((length(localized) - 1) * var(localized) +
                       (length(metastatic) - 1) * var(metastatic)) /
                      (length(localized) + length(metastatic) - 2))
    delta <- mean(metastatic) - mean(localized)
    p <- tryCatch(t.test(metastatic, localized)$p.value, error = function(e) NA_real_)
    tibble(mean_localized = mean(localized), mean_metastatic = mean(metastatic),
           difference = delta, effect_size = delta / pooled_sd, p_value = p)
  }) %>%
  ungroup() %>%
  mutate(
    BH_FDR = p.adjust(p_value, method = "BH"),
    exploratory = !is.na(p_value) & p_value < 0.05 & abs(effect_size) > 1,
    direction = if_else(difference > 0, "metastatic", "localized")
  )

if (requireNamespace("qvalue", quietly = TRUE)) {
  protein_stats$q_value <- tryCatch(
    qvalue::qvalue(protein_stats$p_value)$qvalues,
    error = function(e) rep(NA_real_, nrow(protein_stats))
  )
}

openxlsx::write.xlsx(
  list(Protein_statistics = protein_stats, Imputed_abundance = proteins,
       Missingness = missingness),
  file.path(out_dir, "Proteomic_analysis.xlsx"), overwrite = TRUE
)
write.csv(protein_stats, file.path(out_dir, "protein_statistics.csv"), row.names = FALSE)

# Use the strongest absolute effect when multiple accessions share a gene symbol.
gene_directions <- protein_stats %>%
  filter(!is.na(Gene), nzchar(Gene), is.finite(difference)) %>%
  group_by(Gene) %>%
  slice_max(abs(difference), n = 1, with_ties = FALSE) %>%
  ungroup()
write.csv(gene_directions, file.path(out_dir, "secretome_gene_directions.csv"), row.names = FALSE)

# Ranked enrichment: GO and, when supplied, original GMT collections.
ranking <- protein_stats %>%
  filter(!is.na(Gene), nzchar(Gene), is.finite(effect_size)) %>%
  group_by(Gene) %>%
  slice_max(abs(effect_size), n = 1, with_ties = FALSE) %>%
  ungroup()
gene_list <- sort(setNames(ranking$effect_size, ranking$Gene), decreasing = TRUE)

if (requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  go <- lapply(c("BP", "MF", "CC"), function(ontology) {
    tryCatch(
      as.data.frame(clusterProfiler::gseGO(
        geneList = gene_list, ont = ontology, keyType = "SYMBOL",
        OrgDb = org.Hs.eg.db::org.Hs.eg.db, minGSSize = 3, maxGSSize = 100,
        pvalueCutoff = 0.20, pAdjustMethod = "BH", eps = 0,
        by = "fgsea", verbose = FALSE
      )),
      error = function(e) { warning(conditionMessage(e)); data.frame() }
    )
  })
  names(go) <- c("BP", "MF", "CC")
  openxlsx::write.xlsx(go, file.path(out_dir, "GSEA_GO.xlsx"), overwrite = TRUE)
}

gmt_files <- list.files("data/references/gsea", pattern = "\\.gmt$", full.names = TRUE)
if (length(gmt_files) && requireNamespace("clusterProfiler", quietly = TRUE)) {
  for (file in gmt_files) {
    terms <- clusterProfiler::read.gmt(file)
    gsea <- clusterProfiler::GSEA(
      gene_list, TERM2GENE = terms[, c("term", "gene")],
      minGSSize = 3, maxGSSize = 100, pvalueCutoff = 0.20,
      pAdjustMethod = "BH", eps = 0, by = "fgsea", verbose = FALSE
    )
    output <- paste0("GSEA_", tools::file_path_sans_ext(basename(file)), ".csv")
    write.csv(as.data.frame(gsea), file.path(out_dir, output), row.names = FALSE)
  }
}

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Retained accessions: ", length(unique(proteins$Accession)),
        "; exploratory candidates: ", sum(protein_stats$exploratory))
