install.packages("remotes") # if you haven't already
remotes::install_github("AMRverse/AMRgen")


library(AMRgen)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("/home/user/Escritorio/neisseria/pangenomics/gonorrhoeae_ncbi/complete_genomes/")


# ============================================================
# 4. IMPORTAR RESULTADOS DE AMRFINDERPLUS
# ============================================================

geno <- import_amrfp(
  # input_table = "amrfinder_Neisseria_gonorrhoeae.tsv",
  sample_col = "Name"
)


# ============================================================
# 5. REVISAR LOS DATOS IMPORTADOS
# ============================================================

head(geno)

names(geno)

dim(geno)

str(geno)


# Numero de aislamientos que tienen al menos un hit
n_distinct(geno$id)


# Numero de marcadores diferentes
n_distinct(geno$marker)


# Lista de muestras
sort(unique(geno$id))


# ============================================================
# 6. VER TODOS LOS MARCADORES ENCONTRADOS
# ============================================================

geno %>%
  select(
    id,
    marker,
    gene,
    mutation,
    drug_class,
    drug
  ) %>%
  arrange(id, drug_class, marker)


# ============================================================
# 7. FRECUENCIA DE CADA MARCADOR
# ============================================================

marker_frequency <- geno %>%
  distinct(id, marker) %>%
  count(marker, sort = TRUE)


marker_frequency


# ============================================================
# 8. PREVALENCIA DE MARCADORES
# ============================================================
#
# IMPORTANTE:
# AMRFinderPlus no genera filas para muestras sin ningun hit.
#
# Por eso, si sabes el numero TOTAL de genomas analizados,
# colocalo aqui.
#
# Ejemplo:
# total_isolates <- 40
#
# Sustituye 40 por tu numero real.
# ============================================================

total_isolates <- 32


marker_prevalence <- geno %>%
  distinct(id, marker) %>%
  count(marker, sort = TRUE) %>%
  mutate(
    total_isolates = total_isolates,
    prevalence = n / total_isolates,
    prevalence_percent = round(
      100 * prevalence,
      1
    )
  )


marker_prevalence


write.csv(
  marker_prevalence,
  "AMR_marker_prevalence.csv",
  row.names = FALSE
)


# ============================================================
# 9. FRECUENCIA POR CLASE DE ANTIBIOTICO
# ============================================================

drug_class_frequency <- geno %>%
  distinct(id, drug_class) %>%
  count(drug_class, sort = TRUE) %>%
  mutate(
    percent = round(
      100 * n / total_isolates,
      1
    )
  )


drug_class_frequency


write.csv(
  drug_class_frequency,
  "AMR_drug_class_frequency.csv",
  row.names = FALSE
)


# ============================================================
# 10. TABLA DE MARCADORES POR MUESTRA
# ============================================================

markers_per_isolate <- geno %>%
  distinct(id, marker) %>%
  count(
    id,
    name = "n_markers"
  ) %>%
  arrange(desc(n_markers))


markers_per_isolate


write.csv(
  markers_per_isolate,
  "AMR_markers_per_isolate.csv",
  row.names = FALSE
)


# ============================================================
# 11. CREAR MATRIZ MUESTRA x MARCADOR
# ============================================================
#
# 1 = marcador presente
# 0 = marcador ausente
# ============================================================

marker_matrix <- geno %>%
  distinct(id, marker) %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from = marker,
    values_from = present,
    values_fill = 0
  )


marker_matrix


write.csv(
  marker_matrix,
  "AMR_sample_marker_matrix.csv",
  row.names = FALSE
)


# ============================================================
# 12. MATRIZ MUESTRA x CLASE DE ANTIBIOTICO
# ============================================================

class_matrix <- geno %>%
  filter(!is.na(drug_class)) %>%
  distinct(id, drug_class) %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from = drug_class,
    values_from = present,
    values_fill = 0
  )


class_matrix


write.csv(
  class_matrix,
  "AMR_sample_drugclass_matrix.csv",
  row.names = FALSE
)


# ============================================================
# 13. VER DETERMINANTES DE TETRACICLINA
# ============================================================

tet_geno <- geno %>%
  filter(drug_class == "Tetracyclines")


tet_geno %>%
  select(
    id,
    marker,
    gene,
    mutation,
    drug_class
  ) %>%
  arrange(id, marker)


# Frecuencia de marcadores de tetraciclina

tet_frequency <- tet_geno %>%
  distinct(id, marker) %>%
  count(marker, sort = TRUE) %>%
  mutate(
    percent = round(
      100 * n / total_isolates,
      1
    )
  )


tet_frequency


write.csv(
  tet_frequency,
  "Tetracycline_marker_prevalence.csv",
  row.names = FALSE
)


# ============================================================
# 14. BUSCAR rpsJ V57M Y tet(M)
# ============================================================

geno %>%
  filter(
    marker %in% c(
      "rpsJ_V57M",
      "tet(M)"
    )
  ) %>%
  select(
    id,
    marker,
    gene,
    mutation,
    drug_class
  ) %>%
  arrange(marker, id)


# ============================================================
# 15. PREVALENCIA DE rpsJ V57M Y tet(M)
# ============================================================

geno %>%
  filter(
    marker %in% c(
      "rpsJ_V57M",
      "tet(M)"
    )
  ) %>%
  distinct(id, marker) %>%
  count(marker) %>%
  mutate(
    total = total_isolates,
    percent = round(
      100 * n / total_isolates,
      1
    )
  )


# ============================================================
# 16. COMBINACIONES DE MARCADORES DE TETRACICLINA
# ============================================================
#
# Esto permite saber, por ejemplo:
#
# NG01 -> rpsJ_V57M
# NG02 -> rpsJ_V57M + tet(M)
# NG03 -> tet(M)
#
# ============================================================

tet_combinations <- tet_geno %>%
  distinct(id, marker) %>%
  arrange(id, marker) %>%
  group_by(id) %>%
  summarise(
    markers = paste(
      marker,
      collapse = " + "
    ),
    n_markers = n(),
    .groups = "drop"
  )


tet_combinations


# ============================================================
# 17. FRECUENCIA DE CADA COMBINACION DE TETRACICLINA
# ============================================================

tet_combination_frequency <- tet_combinations %>%
  count(
    markers,
    sort = TRUE
  ) %>%
  mutate(
    percent = round(
      100 * n / total_isolates,
      1
    )
  )


tet_combination_frequency


write.csv(
  tet_combination_frequency,
  "Tetracycline_marker_combinations.csv",
  row.names = FALSE
)


# ============================================================
# 18. GRAFICA DE PREVALENCIA DE MARCADORES
# ============================================================

marker_plot_data <- marker_prevalence %>%
  filter(n > 1)


ggplot(
  marker_plot_data,
  aes(
    x = reorder(marker, prevalence_percent),
    y = prevalence_percent
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    x = "AMR marker",
    y = "Prevalence (%)",
    title = "AMR determinants in Neisseria gonorrhoeae"
  ) +
  theme_bw()


# ============================================================
# 19. GRAFICA DE CLASES DE ANTIBIOTICOS
# ============================================================

ggplot(
  drug_class_frequency,
  aes(
    x = reorder(drug_class, percent),
    y = percent
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Antimicrobial class",
    y = "Isolates carrying determinant (%)",
    title = "Genotypic antimicrobial resistance profile"
  ) +
  theme_bw()


# ============================================================
# 20. HEATMAP MUESTRA x MARCADOR
# ============================================================

heatmap_data <- geno %>%
  distinct(id, marker) %>%
  mutate(present = 1)


ggplot(
  heatmap_data,
  aes(
    x = marker,
    y = id,
    fill = factor(present)
  )
) +
  geom_tile() +
  labs(
    x = "AMR marker",
    y = "Isolate",
    fill = "Detected",
    title = "AMR determinant profile"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )


# ============================================================
# 21. HEATMAP SOLO PARA TETRACICLINAS
# ============================================================

tet_heatmap <- tet_geno %>%
  distinct(id, marker) %>%
  mutate(present = 1)


ggplot(
  tet_heatmap,
  aes(
    x = marker,
    y = id,
    fill = factor(present)
  )
) +
  geom_tile() +
  labs(
    x = "Tetracycline resistance marker",
    y = "Isolate",
    fill = "Detected",
    title = "Tetracycline resistance determinants"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )


# ============================================================
# 22. EXPORTAR TABLA COMPLETA PROCESADA POR AMRgen
# ============================================================

write.csv(
  geno,
  "AMRFinderPlus_AMRgen_processed.csv",
  row.names = FALSE
)


# ============================================================
# 23. RESUMEN
# ============================================================

cat("\n")
cat("========================================\n")
cat("RESUMEN GENOTIPICO DE AMR\n")
cat("========================================\n")

cat(
  "Genomas analizados:",
  total_isolates,
  "\n"
)

cat(
  "Genomas con al menos un determinante:",
  n_distinct(geno$id),
  "\n"
)

cat(
  "Marcadores AMR diferentes:",
  n_distinct(geno$marker),
  "\n"
)

cat(
  "Clases de antimicrobianos:",
  n_distinct(
    geno$drug_class,
    na.rm = TRUE
  ),
  "\n"
)

cat("========================================\n")


# ============================================================
# ANALISIS QUE NO DEBEN REALIZARSE SIN DATOS FENOTIPICOS
# ============================================================
#
# NO utilizar:
#
# format_pheno()
# get_binary_matrix()
# solo_ppv()
# amr_ppv()
# amr_logistic()
# concordance()
#
# Tampoco interpretar un determinante genotipico como
# resistencia fenotipica demostrada.
#
# La presencia de un marcador indica un determinante
# genotipico asociado con AMR; no sustituye una prueba
# fenotipica de susceptibilidad.
