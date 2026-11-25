# STRING networks, Cytoscape-ready exports and degree-preserving null models
# Run after 01_proteomics.R and 02_single_cell.R.

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
  library(httr)
  library(jsonlite)
  library(igraph)
})

out_dir <- "results/networks"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  "results/secretome_gene_directions.csv",
  "results/NK_PDAC_membrane_genes.csv",
  "results/NK_ADJ_membrane_genes.csv",
  "data/annotations/Exocarta pancreatic cancer cells.xlsx",
  "data/annotations/Annotation_RNAseq.rda"
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("Missing inputs: ", paste(missing_files, collapse = ", "))

# Assemble the three secretome sets and two NK-cell contexts.
directions <- read.csv(required_files[1], check.names = FALSE)
exo <- openxlsx::read.xlsx(required_files[4], startRow = 4, colNames = TRUE)
if (!all(c("Species", "Gene Symbol") %in% names(exo))) {
  stop("Unexpected ExoCarta headers: ", paste(names(exo), collapse = ", "))
}
exo_genes <- unique(trimws(as.character(exo[["Gene Symbol"]][exo$Species == "Homo sapiens"])))

load(required_files[5])
if (!exists("annot") || !"GOComponent" %in% names(annot)) {
  stop("Annotation_RNAseq.rda must contain 'annot' with GOComponent")
}
membrane <- unique(as.character(annot[[2]][grepl("plasma membrane", annot$GOComponent,
                                                  ignore.case = TRUE)]))
secretome <- Reduce(intersect, list(directions$Gene, exo_genes, membrane))
if (!length(secretome)) stop("No genes survived the ExoCarta / membrane intersection")

secretome_sets <- list(
  total = secretome,
  localized = intersect(secretome, directions$Gene[directions$difference < 0]),
  metastatic = intersect(secretome, directions$Gene[directions$difference > 0])
)
nk_sets <- list(
  PDAC = read.csv(required_files[2])$Gene,
  ADJ = read.csv(required_files[3])$Gene
)

# Request high-confidence edges and retain secretome-to-NK associations only.
fetch_string_edges <- function(secretome_genes, nk_genes, cutoff) {
  genes <- unique(c(secretome_genes, nk_genes))
  if (length(genes) < 2) return(data.frame())

  response <- httr::POST(
    "https://string-db.org/api/json/network",
    body = list(identifiers = paste(genes, collapse = "\r"),
                species = 9606, required_score = as.integer(cutoff * 1000)),
    encode = "form"
  )
  httr::stop_for_status(response)
  tab <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))
  if (!nrow(tab)) return(data.frame())

  tab <- tab %>% filter(score >= cutoff, escore >= cutoff | dscore >= cutoff)
  forward <- tab %>%
    filter(preferredName_A %in% secretome_genes, preferredName_B %in% nk_genes) %>%
    transmute(source = preferredName_A, target = preferredName_B,
              combined_score = score, experimental_score = escore, database_score = dscore)
  reverse <- tab %>%
    filter(preferredName_B %in% secretome_genes, preferredName_A %in% nk_genes) %>%
    transmute(source = preferredName_B, target = preferredName_A,
              combined_score = score, experimental_score = escore, database_score = dscore)

  bind_rows(forward, reverse) %>% distinct(source, target, .keep_all = TRUE)
}

# Role-prefixed IDs prevent a shared gene symbol from merging the two compartments.
make_graph <- function(edges) {
  igraph::graph_from_data_frame(
    data.frame(from = paste0("S::", edges$source),
               to = paste0("N::", edges$target)), directed = FALSE
  )
}

# Preserve degree exactly by swapping the NK endpoints of two intergroup edges.
rewire_bipartite <- function(edges, n_swaps) {
  if (nrow(edges) < 2) return(edges)

  existing <- new.env(hash = TRUE, parent = emptyenv())
  edge_key <- function(s, n) paste(s, n, sep = "|")
  for (i in seq_len(nrow(edges))) {
    existing[[edge_key(edges$source[i], edges$target[i])]] <- TRUE
  }

  accepted <- 0L
  attempts <- 0L
  max_attempts <- max(1000L, n_swaps * 100L)

  while (accepted < n_swaps && attempts < max_attempts) {
    attempts <- attempts + 1L
    pair <- sample.int(nrow(edges), 2)
    i <- pair[1]; j <- pair[2]
    if (edges$source[i] == edges$source[j] || edges$target[i] == edges$target[j]) next

    old1 <- edge_key(edges$source[i], edges$target[i])
    old2 <- edge_key(edges$source[j], edges$target[j])
    new1 <- edge_key(edges$source[i], edges$target[j])
    new2 <- edge_key(edges$source[j], edges$target[i])
    if (exists(new1, existing, inherits = FALSE) ||
        exists(new2, existing, inherits = FALSE)) next

    rm(list = c(old1, old2), envir = existing)
    existing[[new1]] <- TRUE
    existing[[new2]] <- TRUE
    tmp <- edges$target[i]
    edges$target[i] <- edges$target[j]
    edges$target[j] <- tmp
    accepted <- accepted + 1L
  }

  if (accepted < n_swaps) warning("Only ", accepted, "/", n_swaps, " swaps completed")
  edges
}

run_null_model <- function(edges, n_random = 1000, seed = 123) {
  set.seed(seed)
  edges <- edges[, c("source", "target"), drop = FALSE]
  observed_graph <- make_graph(edges)
  nodes <- igraph::V(observed_graph)$name
  degrees <- igraph::degree(observed_graph)[nodes]
  observed <- igraph::betweenness(observed_graph, directed = FALSE, normalized = TRUE)[nodes]

  simulated <- replicate(n_random, {
    shuffled <- rewire_bipartite(edges, n_swaps = 20 * nrow(edges))
    graph <- make_graph(shuffled)
    current_degrees <- igraph::degree(graph)[nodes]
    stopifnot(identical(as.integer(current_degrees), as.integer(degrees)))
    as.numeric(igraph::betweenness(graph, directed = FALSE, normalized = TRUE)[nodes])
  })
  if (is.null(dim(simulated))) simulated <- matrix(simulated, nrow = length(nodes))

  empirical_p <- (1 + rowSums(simulated >= as.numeric(observed))) / (n_random + 1)
  data.frame(
    node_id = nodes,
    gene = sub("^[SN]::", "", nodes),
    group = ifelse(startsWith(nodes, "S::"), "Secretome", "NK"),
    degree = as.numeric(degrees),
    observed_betweenness = as.numeric(observed),
    null_mean = rowMeans(simulated),
    null_median = apply(simulated, 1, median),
    null_q025 = apply(simulated, 1, quantile, probs = 0.025),
    null_q975 = apply(simulated, 1, quantile, probs = 0.975),
    empirical_p = empirical_p,
    BH_FDR = p.adjust(empirical_p, method = "BH")
  )
}

# Analyze six context combinations at STRING 0.7, 0.8 and 0.9.
network_summaries <- list()
null_results <- list()
n_random <- as.integer(Sys.getenv("NULL_B", "1000"))
if (is.na(n_random) || n_random < 1) stop("NULL_B must be a positive integer")

for (secretome_name in names(secretome_sets)) {
  for (nk_name in names(nk_sets)) {
    network_name <- paste(secretome_name, nk_name, sep = "_")

    for (cutoff in c(0.7, 0.8, 0.9)) {
      edges <- fetch_string_edges(secretome_sets[[secretome_name]],
                                  nk_sets[[nk_name]], cutoff)
      if (!nrow(edges)) {
        warning("No edges found for ", network_name, " at ", cutoff)
        next
      }

      suffix <- paste0(network_name, "_", as.integer(cutoff * 10))
      nodes <- bind_rows(
        data.frame(id = paste0("S::", unique(edges$source)),
                   gene = unique(edges$source), group = "Secretome"),
        data.frame(id = paste0("N::", unique(edges$target)),
                   gene = unique(edges$target), group = "NK")
      )
      write.csv(edges, file.path(out_dir, paste0(suffix, "_edges.csv")), row.names = FALSE)
      write.csv(nodes, file.path(out_dir, paste0(suffix, "_nodes.csv")), row.names = FALSE)
      network_summaries[[suffix]] <- data.frame(
        network = network_name, score_threshold = cutoff,
        nodes = nrow(nodes), edges = nrow(edges)
      )

      if (cutoff == 0.9) {
        result <- run_null_model(edges, n_random = n_random)
        result$network <- network_name
        null_results[[network_name]] <- result
        write.csv(result, file.path(out_dir, paste0(network_name, "_null.csv")),
                  row.names = FALSE)
      }
    }
  }
}

write.csv(bind_rows(network_summaries), "results/network_threshold_summary.csv", row.names = FALSE)
write.csv(bind_rows(null_results), "results/degree_preserving_null_summary.csv", row.names = FALSE)

# Preserve the original Cytoscape/CytoHubba centrality outputs for figure reproduction.
centrality_files <- list.files("data/references/centrality",
                              pattern = "^Network_VES.*\\.csv$", full.names = TRUE)
dir.create("results/centrality", recursive = TRUE, showWarnings = FALSE)
for (file in centrality_files) {
  file.copy(file, file.path("results/centrality", basename(file)), overwrite = TRUE)
}

message("Network edge/node tables, null-model statistics and historical centrality files exported")
