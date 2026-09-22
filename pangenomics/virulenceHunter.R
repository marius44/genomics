# ============================================================
# Análisis conjunto de resultados de VirulentHunter
# ============================================================

# Instalar una sola vez si hace falta:
# install.packages(c("tidyverse", "pheatmap", "viridis"))

library(tidyverse)
library(pheatmap)
library(viridis)
library(ggsci)

# ============================================================
# 1. Directorio con los resultados
#
# Estructura esperada:
#
# virulenthunter_results/
# ├── NG01/
# │   └── predict_results.csv
# ├── NG02/
# │   └── predict_results.csv
# ├── NG03/
# │   └── predict_results.csv
# ...
# ============================================================

results_dir <- "/home/user/Escritorio/neisseria/pangenomics/gonorrhoeae_ncbi/complete_genomes/Neisseria/virulenthunter_results/procesadas_mac/"


# ============================================================
# 2. Parámetros
# ============================================================

VF_THRESHOLD <- 0.5
CATEGORY_THRESHOLD <- 0.5


# ============================================================
# 3. Localizar todos los resultados
# ============================================================

files <- list.files(
  results_dir,
  pattern = "predict_results\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

cat("Archivos encontrados:", length(files), "\n")

if (length(files) == 0) {
  stop("No se encontraron archivos predict_results.csv")
}


# ============================================================
# 4. Leer y combinar todos los genomas
# ============================================================

read_virulenthunter <- function(file) {
  
  # El nombre del directorio padre será el nombre de la muestra
  sample_name <- basename(dirname(file))
  
  x <- read_csv(
    file,
    show_col_types = FALSE,
    name_repair = "unique"
  )
  
  # El CSV de pandas incluye una columna de índice innecesaria.
  # Usualmente se llama ...1
  if (grepl("^\\.\\.\\.", names(x)[1])) {
    x <- x[, -1]
  }
  
  x %>%
    mutate(
      sample = sample_name,
      .before = 1
    )
}


all_results <- map_dfr(files, read_virulenthunter)


cat("Genomas:", n_distinct(all_results$sample), "\n")
cat("Proteínas totales:", nrow(all_results), "\n")
str(all_results)

# ============================================================
# 5. Guardar tabla maestra
# ============================================================

write_csv(
  all_results,
  file.path(results_dir, "VirulentHunter_all_results.csv")
)


# ============================================================
# 6. Identificar candidatos a VF
# ============================================================

vf_results <- all_results %>%
  filter(vf_prob >= VF_THRESHOLD)


cat(
  "Proteínas con vf_prob >=",
  VF_THRESHOLD,
  ":",
  nrow(vf_results),
  "\n"
)


write_csv(
  vf_results,
  file.path(results_dir, "VirulentHunter_VF_candidates.csv")
)
vf_results

# ============================================================
# 7. Definir columnas funcionales
# ============================================================

category_columns <- c(
  "Exotoxin",
  "Stress survival",
  "Biofilm",
  "Immune modulation",
  "Invasion",
  "Adherence",
  "Effector delivery system",
  "Nutritional/Metabolic factor",
  "Motility",
  "Antimicrobial activity/Competitive advantage",
  "Others",
  "Post-translational modification",
  "Exoenzyme",
  "Regulation"
)

# Comprobar que estén presentes
missing_categories <- setdiff(category_columns, names(all_results))

if (length(missing_categories) > 0) {
  stop(
    "Faltan estas columnas: ",
    paste(missing_categories, collapse = ", ")
  )
}


# ============================================================
# 8. Categoría principal de cada proteína
#
# Para cada VF candidato, toma la categoría con mayor score.
# ============================================================

vf_primary <- vf_results %>%
  rowwise() %>%
  mutate(
    primary_category = {
      values <- c_across(all_of(category_columns))
      
      if (max(values, na.rm = TRUE) >= CATEGORY_THRESHOLD) {
        category_columns[which.max(values)]
      } else {
        "Unclassified"
      }
    },
    
    primary_score = max(
      c_across(all_of(category_columns)),
      na.rm = TRUE
    )
  ) %>%
  ungroup()


write_csv(
  vf_primary,
  file.path(results_dir, "VirulentHunter_VF_primary_category.csv")
)
vf_primary

# ============================================================
# 9. Número total de VFs por genoma
# ============================================================

vf_per_genome <- all_results %>%
  group_by(sample) %>%
  summarise(
    proteins = n(),
    predicted_VF = sum(vf_prob >= VF_THRESHOLD),
    VF_percent = 100 * predicted_VF / proteins,
    .groups = "drop"
  ) %>%
  arrange(desc(predicted_VF))


write_csv(
  vf_per_genome,
  file.path(results_dir, "VirulentHunter_VF_per_genome.csv")
)

print(vf_per_genome)


# ============================================================
# 10. Gráfico: VFs totales por genoma
# ============================================================

p1 <- ggplot(
  vf_per_genome,
  aes(
    x = reorder(sample, predicted_VF),
    y = predicted_VF
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Genoma",
    y = "Número de factores de virulencia predichos",
    title = "Factores de virulencia predichos por VirulentHunter",
    subtitle = paste0("Umbral vf_prob >= ", VF_THRESHOLD)
  ) +
  theme_bw(base_size = 12)+
  scale_fill_igv()

p1

ggsave(
  file.path(results_dir, "VF_total_por_genoma.pdf"),
  p1,
  width = 8,
  height = 10
)

ggsave(
  file.path(results_dir, "VF_total_por_genoma.png"),
  p1,
  width = 8,
  height = 10,
  dpi = 300
)


# ============================================================
# 11. Transformar categorías a formato largo
# ============================================================

vf_long <- vf_results %>%
  select(
    sample,
    id,
    vf_prob,
    all_of(category_columns)
  ) %>%
  pivot_longer(
    cols = all_of(category_columns),
    names_to = "category",
    values_to = "category_score"
  )
#View(vf_long)

# ============================================================
# 12. Considerar categoría presente si score >= 0.5
#
# IMPORTANTE:
# VirulentHunter es multi-label, por lo que una proteína puede
# pertenecer a más de una categoría.
# ============================================================

vf_category_hits <- vf_long %>%
  filter(category_score >= CATEGORY_THRESHOLD)
#View(vf_category_hits)

write_csv(
  vf_category_hits,
  file.path(results_dir, "VirulentHunter_category_hits.csv")
)


# ============================================================
# 13. Número de VFs de cada categoría por genoma
# ============================================================

category_counts <- vf_category_hits %>%
  count(
    sample,
    category,
    name = "n"
  ) %>%
  complete(
    sample,
    category = category_columns,
    fill = list(n = 0)
  )
#View(category_counts)

write_csv(
  category_counts,
  file.path(results_dir, "VirulentHunter_category_counts.csv")
)


# ============================================================
# 14. Heatmap de número de VFs
# ============================================================

heatmap_matrix <- category_counts %>%
  pivot_wider(
    names_from = category,
    values_from = n,
    values_fill = 0
  ) %>%
  column_to_rownames("sample") %>%
  as.matrix()

#View(heatmap_matrix)
# Mostrar en el panel Plots de RStudio
heatmap_vf <- pheatmap(
  heatmap_matrix,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  fontsize_row = 7,
  fontsize_col = 8,
  angle_col = 45,
  main = "Virulence factor categories"
)

# Guardar directamente en PDF sin abrir un dispositivo gráfico
#pheatmap(
#  heatmap_matrix,
#  scale = "none",
#  cluster_rows = TRUE,
#  cluster_cols = TRUE,
 # border_color = NA,
#  fontsize_row = 7,
#  fontsize_col = 8,
#  angle_col = 45,
#  main = "Virulence factor categories",
#  filename = file.path(
#    results_dir,
#    "VirulentHunter_heatmap_counts.pdf"
#  ),
 # width = 12,
# height = 10
#)



# Guardar PDF
ggsave(
  filename = file.path(
    results_dir,
    "VirulentHunter_heatmap_counts.pdf"
  ),
  plot = heatmap_vf$gtable,
  width = 12,
  height = 10
)

# Guardar también PNG
ggsave(
  filename = file.path(
    results_dir,
    "VirulentHunter_heatmap_counts.pdf"
  ),
  plot = heatmap_vf$gtable,
  width = 12,
  height = 10,
  dpi = 300
)


heatmap_vf

# ============================================================
# 15. Heatmap normalizado
#
# Proporción de cada categoría respecto al total de proteínas
# del proteoma.
# Esto permite comparar genomas con diferentes números de CDS.
# ============================================================

protein_counts <- all_results %>%
  count(
    sample,
    name = "total_proteins"
  )


category_normalized <- category_counts %>%
  left_join(
    protein_counts,
    by = "sample"
  ) %>%
  mutate(
    percentage = 100 * n / total_proteins
  )


norm_matrix <- category_normalized %>%
  select(
    sample,
    category,
    percentage
  ) %>%
  pivot_wider(
    names_from = category,
    values_from = percentage,
    values_fill = 0
  ) %>%
  column_to_rownames("sample") %>%
  as.matrix()

heatmap_vf2 <- pheatmap(
  norm_matrix,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  fontsize_row = 7,
  fontsize_col = 8,
  angle_col = 45,
  main = "Virulence factor categories (% of proteome)"
)


heatmap_vf2


# Guardar PDF
ggsave(
  filename = file.path(
    results_dir,
    "VirulentHunter_heatmap_percentage.pdf"
  ),
  plot = heatmap_vf2$gtable,
  width = 12,
  height = 10
)

# Guardar también PNG
ggsave(
  filename = file.path(
    results_dir,
    "VirulentHunter_heatmap_percentage.png"
  ),
  plot = heatmap_vf2$gtable,
  width = 12,
  height = 10,
  dpi = 300
)

# ============================================================
# 16. Barras apiladas por categoría
# ============================================================

p2 <- ggplot(
  category_counts,
  aes(
    x = sample,
    y = n,
    fill = category
  )
) +
  geom_col() +
  labs(
    x = "Genome",
    y = "Proteins",
    fill = "Category",
    title = "Virulence factors profiles per genome"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )+
  scale_fill_igv()
  
  

p2

ggsave(
  file.path(results_dir, "VF_categories_stacked.pdf"),
  p2,
  width = 14,
  height = 8
)

ggsave(
  file.path(results_dir, "VF_categories_stacked.png"),
  p2,
  width = 14,
  height = 8,
  dpi = 300
)


# ============================================================
# 17. Distribución de vf_prob
# ============================================================

p3 <- ggplot(
  all_results,
  aes(x = vf_prob)
) +
  geom_histogram(
    bins = 50
  ) +
  geom_vline(
    xintercept = VF_THRESHOLD,
    linetype = "dashed"
  ) +
  labs(
    x = "VirulentHunter vf_prob",
    y = "# proteins",
    title = "VirulentHunter score distribution"
  ) +
  theme_bw(base_size = 12)+
  scale_fill_igv()
p3

ggsave(
  file.path(results_dir, "VF_probability_distribution.pdf"),
  p3,
  width = 8,
  height = 6
)


# ============================================================
# 18. Resumen final
# ============================================================

cat("\n========================================\n")
cat("ANÁLISIS TERMINADO\n")
cat("========================================\n")
cat("Genomas:", n_distinct(all_results$sample), "\n")
cat("Proteínas:", nrow(all_results), "\n")
cat("VFs candidatos:", nrow(vf_results), "\n")
cat(
  "Porcentaje global:",
  round(100 * nrow(vf_results) / nrow(all_results), 2),
  "%\n"
)
cat("Resultados guardados en:\n")
cat(results_dir, "\n")

