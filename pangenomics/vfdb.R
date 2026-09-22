# Heatmap para virulencia de Neisseria

library(ComplexHeatmap)
library(grid)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("ComplexHeatmap")

# setwd("/home/user/Escritorio/neisseria/pangenomics/gonorrhoeae_ncbi/complete_genomes/vfdb_results/")

# Read the binary ABRicate matrix
d <- read.delim(
  "vfdb_summary_binary.tsv",
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  comment.char = "",
  quote = "",
  stringsAsFactors = FALSE
)

names(d)[1:2] <- c("isolate", "NUM_FOUND")

# Clean isolate labels
ids <- basename(d$isolate)
ids <- sub(
  "\\.(fsa|fa|fna|fasta)(\\.gz)?$",
  "",
  ids,
  ignore.case = TRUE
)
d
# Optional: remove this suffix from the corrected assembly label
ids <- sub("_corrected$", "", ids)

if (anyDuplicated(ids))
  stop("Duplicate isolate names were produced after cleaning.")

# Convert sample × gene table into gene × sample matrix
m <- as.matrix(d[, -(1:2), drop = FALSE])
storage.mode(m) <- "numeric"

if (anyNA(m) || !all(m %in% c(0, 1)))
  stop("The file is not a valid binary 0/1 matrix.")

rownames(m) <- ids
m <- t(m)

# Gene prevalence
prevalence <- data.frame(
  gene = rownames(m),
  isolates_present = rowSums(m),
  isolates_total = ncol(m),
  prevalence_percent = round(rowMeans(m) * 100, 1),
  classification = ifelse(
    rowSums(m) == ncol(m),
    "Detected in all isolates",
    "Variable"
  )
)

prevalence <- prevalence[
  order(-prevalence$prevalence_percent, prevalence$gene),
]

write.table(
  prevalence,
  "vfdb_gene_prevalence.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Keep only genes that vary among isolates
m_variable <- m[
  rowSums(m) > 0 & rowSums(m) < ncol(m),
  ,
  drop = FALSE
]

message("Genes detected: ", nrow(m))
message("Variable genes: ", nrow(m_variable))
message("Genes detected in all isolates: ",
        sum(rowSums(m) == ncol(m)))

if (nrow(m_variable) == 0)
  stop("No variable genes were detected; a heatmap would not be informative.")

# Jaccard distance for binary profiles
jaccard_distance <- function(x) {
  x <- x > 0
  n <- nrow(x)
  result <- matrix(
    0, n, n,
    dimnames = list(rownames(x), rownames(x))
  )
  
  if (n > 1) {
    for (i in seq_len(n - 1)) {
      for (j in (i + 1):n) {
        union <- sum(x[i, ] | x[j, ])
        distance <- if (union == 0) {
          0
        } else {
          1 - sum(x[i, ] & x[j, ]) / union
        }
        
        result[i, j] <- distance
        result[j, i] <- distance
      }
    }
  }
  
  as.dist(result)
}

row_clustering <- if (nrow(m_variable) > 1) {
  hclust(jaccard_distance(m_variable), method = "average")
} else {
  FALSE
}

column_clustering <- if (ncol(m_variable) > 1) {
  hclust(jaccard_distance(t(m_variable)), method = "average")
} else {
  FALSE
}

make_heatmap <- function() {
  Heatmap(
    m_variable,
    name = "VFDB hit",
    col = c("0" = "#F2F2F2", "1" = "#0072B2"),
    cluster_rows = row_clustering,
    cluster_columns = column_clustering,
    show_row_dend = nrow(m_variable) > 1,
    show_column_dend = ncol(m_variable) > 1,
    rect_gp = gpar(col = "white", lwd = 0.4),
    row_names_gp = gpar(fontsize = 7, fontface = "italic"),
    column_names_gp = gpar(fontsize = 7),
    column_names_rot = 45,
    row_title = "Variable virulence-associated genes",
    column_title = "Neisseria gonorrhoeae isolates",
    heatmap_legend_param = list(
      at = c(0, 1),
      labels = c("Absent", "Present")
    )
  )
}

figure_width  <- min(16, max(8, 4 + 0.22 * ncol(m_variable)))
figure_height <- min(20, max(6, 3 + 0.18 * nrow(m_variable)))
make_heatmap
# Vector version, preferable for publication
pdf(
  "vfdb_variable_heatmap.pdf",
  width = figure_width,
  height = figure_height,
  useDingbats = FALSE
)
draw(make_heatmap())
dev.off()

# High-resolution raster version
tiff(
  "vfdb_variable_heatmap.tiff",
  width = figure_width,
  height = figure_height,
  units = "in",
  res = 600,
  compression = "lzw"
)
draw(make_heatmap())
dev.off()

