# ============================================================
# SANKEY INTEGRADO
# VFDB -> VirulentHunter -> AMRgen / AMRFinderPlus
# Neisseria gonorrhoeae
#
# Diseñado para TUS archivos reales:
#
#   VFDB:
#     vfdb_summary_binary.tsv
#
#   VirulentHunter:
#     VirulentHunter_category_hits.csv
#
#   AMR:
#     amrfinder_Neisseria_gonorrhoeae.tsv
#
# IMPORTANTE:
# Los flujos representan COOCURRENCIA dentro de los mismos
# genomas. NO representan causalidad ni interacción molecular.
# ============================================================



# ============================================================
# 0. PAQUETES
# ============================================================

packages <- c(
  "tidyverse",
  "ggalluvial",
  "AMRgen"
)

missing_packages <- packages[
  !vapply(
    packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  
  stop(
    paste0(
      "Faltan estos paquetes:\n",
      paste(missing_packages, collapse = ", "),
      "\n\nInstálalos antes de continuar."
    )
  )
}

library(tidyverse)
library(ggalluvial)
library(AMRgen)



# ============================================================
# 1. RUTAS
# ============================================================
#
# VFDB está en tu directorio actual.
#
# VirulentHunter está en procesadas_mac.
#
# AMRFinderPlus está en complete_genomes.
#
# Si alguno está en otro sitio, cambia SOLO esta sección.
# ============================================================

setwd("/home/user/Escritorio/neisseria/pangenomics/gonorrhoeae_ncbi/complete_genomes/vfdb_virulenthunter_amrfinder")

VFDB_FILE <- "/home/user/Escritorio/neisseria/pangenomics/gonorrhoeae_ncbi/complete_genomes/vfdb_results/vfdb_summary_binary.tsv"


VH_FILE <- paste0(
  "/home/user/Escritorio/neisseria/pangenomics/",
  "gonorrhoeae_ncbi/complete_genomes/Neisseria/",
  "virulenthunter_results/procesadas_mac/",
  "VirulentHunter_category_hits.csv"
)


AMR_FILE <- paste0(
  "/home/user/Escritorio/neisseria/pangenomics/",
  "gonorrhoeae_ncbi/complete_genomes/",
  "amrfinder_Neisseria_gonorrhoeae.tsv"
)



# ============================================================
# 2. COMPROBAR ARCHIVOS
# ============================================================

cat("\n========================================\n")
cat("COMPROBANDO ARCHIVOS\n")
cat("========================================\n")

cat(
  "VFDB:            ",
  file.exists(VFDB_FILE),
  "  ",
  VFDB_FILE,
  "\n"
)

cat(
  "VirulentHunter:  ",
  file.exists(VH_FILE),
  "  ",
  VH_FILE,
  "\n"
)

cat(
  "AMRFinderPlus:   ",
  file.exists(AMR_FILE),
  "  ",
  AMR_FILE,
  "\n"
)


if (!file.exists(VFDB_FILE)) {
  stop("No se encontró VFDB_FILE.")
}

if (!file.exists(VH_FILE)) {
  stop("No se encontró VH_FILE.")
}

if (!file.exists(AMR_FILE)) {
  stop("No se encontró AMR_FILE.")
}



# ============================================================
# 3. FUNCIÓN PARA LIMPIAR IDs DE LAS MUESTRAS
# ============================================================
#
# Ejemplos:
#
# NG02_final       -> NG02
# NG03.fna         -> NG03
# NG05.fna.faa     -> NG05
# 003-1.fna        -> 003-1
# ATCC.fna         -> ATCC
#
# ============================================================

clean_sample <- function(x) {
  
  x <- basename(as.character(x))
  
  x <- stringr::str_trim(x)
  
  x <- stringr::str_remove(
    x,
    "\\.(csv|tsv|txt|fna|faa|fasta|fa)(\\.(faa|fna))?$"
  )
  
  x <- stringr::str_remove(
    x,
    "_final$"
  )
  
  x
}



# ============================================================
# ============================================================
#
# PARTE A
# VFDB / ABRICATE --summary
#
# ============================================================
# ============================================================



# ============================================================
# 4. IMPORTAR VFDB
# ============================================================

vfdb <- readr::read_tsv(
  VFDB_FILE,
  show_col_types = FALSE,
  name_repair = "minimal"
)


cat("\n========================================\n")
cat("VFDB\n")
cat("========================================\n")

cat(
  "Dimensiones:",
  nrow(vfdb),
  "x",
  ncol(vfdb),
  "\n"
)



# ============================================================
# 5. CREAR / NORMALIZAR COLUMNA sample
# ============================================================

if ("sample" %in% names(vfdb)) {
  
  vfdb <- vfdb %>%
    mutate(
      sample = clean_sample(sample)
    )
  
} else if ("#FILE" %in% names(vfdb)) {
  
  vfdb <- vfdb %>%
    mutate(
      sample = clean_sample(.data[["#FILE"]])
    )
  
} else {
  
  stop(
    "VFDB no contiene ni 'sample' ni '#FILE'."
  )
}


cat(
  "Muestras VFDB:",
  n_distinct(vfdb$sample),
  "\n"
)

print(
  sort(unique(vfdb$sample))
)



# ============================================================
# 6. IDENTIFICAR COLUMNAS DE GENES
# ============================================================

vfdb_gene_cols <- setdiff(
  names(vfdb),
  c(
    "#FILE",
    "NUM_FOUND",
    "sample"
  )
)


cat(
  "\nColumnas de genes VFDB:",
  length(vfdb_gene_cols),
  "\n"
)



# ============================================================
# 7. PASAR VFDB A FORMATO LARGO
# ============================================================

vfdb_hits <- vfdb %>%
  
  pivot_longer(
    cols = all_of(vfdb_gene_cols),
    names_to = "vfdb_gene",
    values_to = "vfdb_count"
  ) %>%
  
  mutate(
    vfdb_count = suppressWarnings(
      as.numeric(vfdb_count)
    )
  ) %>%
  
  filter(
    !is.na(vfdb_count),
    vfdb_count > 0
  ) %>%
  
  select(
    sample,
    vfdb_gene,
    vfdb_count
  ) %>%
  
  distinct()


cat(
  "Asociaciones muestra-gen VFDB:",
  nrow(vfdb_hits),
  "\n"
)


print(
  head(vfdb_hits, 10)
)



# ============================================================
# 8. FRECUENCIA DE GENES VFDB
# ============================================================

vfdb_frequency <- vfdb_hits %>%
  
  distinct(
    sample,
    vfdb_gene
  ) %>%
  
  count(
    vfdb_gene,
    sort = TRUE
  ) %>%
  
  mutate(
    prevalence_percent =
      100 * n /
      n_distinct(vfdb$sample)
  )


cat("\nGenes VFDB más frecuentes:\n")

print(
  head(vfdb_frequency, 20)
)



# ============================================================
# ============================================================
#
# PARTE B
# VIRULENTHUNTER
#
# ============================================================
# ============================================================



# ============================================================
# 9. IMPORTAR VirulentHunter_category_hits.csv
# ============================================================

vh_raw <- readr::read_csv(
  VH_FILE,
  show_col_types = FALSE,
  name_repair = "minimal"
)


cat("\n========================================\n")
cat("VIRULENTHUNTER\n")
cat("========================================\n")

cat(
  "Dimensiones:",
  nrow(vh_raw),
  "x",
  ncol(vh_raw),
  "\n"
)

cat("\nColumnas encontradas:\n")

print(
  names(vh_raw)
)



# ============================================================
# 10. DETECTAR COLUMNA DE MUESTRA
# ============================================================

sample_candidates <- c(
  "sample",
  "Sample",
  "genome",
  "Genome",
  "genome_id",
  "Genome_ID",
  "isolate",
  "Isolate",
  "strain",
  "Strain",
  "file",
  "File",
  "filename",
  "Filename"
)


vh_sample_col <- intersect(
  sample_candidates,
  names(vh_raw)
)


# Si no coincide exactamente, buscar por patrón

if (length(vh_sample_col) == 0) {
  
  vh_sample_col <- grep(
    "sample|genome|isolate|strain|file",
    names(vh_raw),
    ignore.case = TRUE,
    value = TRUE
  )
}


if (length(vh_sample_col) == 0) {
  
  stop(
    paste0(
      "No pude identificar la columna de muestra en ",
      "VirulentHunter_category_hits.csv.\n\n",
      "Columnas:\n",
      paste(names(vh_raw), collapse = ", ")
    )
  )
}


vh_sample_col <- vh_sample_col[1]


cat(
  "\nColumna de muestra detectada:",
  vh_sample_col,
  "\n"
)



# ============================================================
# 11. DETECTAR SI CATEGORY_HITS ESTÁ EN FORMATO LARGO
# ============================================================

category_candidates <- c(
  "category",
  "Category",
  "categories",
  "Categories",
  "vf_category",
  "VF_category",
  "VF_Category",
  "primary_category",
  "Primary_category",
  "Primary_Category"
)


vh_category_col <- intersect(
  category_candidates,
  names(vh_raw)
)


if (length(vh_category_col) == 0) {
  
  vh_category_col <- grep(
    "categor",
    names(vh_raw),
    ignore.case = TRUE,
    value = TRUE
  )
}



# ============================================================
# 12A. CASO 1:
# category_hits ya está en formato largo
#
# sample | category
# NG02   | Adherence
# NG02   | Biofilm
# ...
# ============================================================
# ============================================================
# 12. CREAR TABLA sample -> vf_category
# ============================================================

if (length(vh_category_col) > 0) {
  
  # ----------------------------------------------------------
  # CASO A:
  # VirulentHunter_category_hits.csv está en formato largo
  #
  # Ejemplo:
  # sample   category
  # NG02     Adherence
  # NG02     Biofilm
  # ----------------------------------------------------------
  
  vh_category_col <- vh_category_col[1]
  
  cat(
    "Formato VirulentHunter detectado: LARGO\n"
  )
  
  cat(
    "Columna de categoría:",
    vh_category_col,
    "\n"
  )
  
  vh_categories <- vh_raw %>%
    
    transmute(
      
      sample = clean_sample(
        .data[[vh_sample_col]]
      ),
      
      vf_category = as.character(
        .data[[vh_category_col]]
      )
    ) %>%
    
    # Eliminar corchetes y comillas si las categorías
    # fueron guardadas como texto tipo lista
    
    mutate(
      vf_category = stringr::str_replace_all(
        vf_category,
        "[\\[\\]\"']",
        ""
      )
    ) %>%
    
    # Si una celda contiene varias categorías,
    # separarlas en filas distintas
    
    tidyr::separate_rows(
      vf_category,
      sep = "\\s*[;,|]+\\s*"
    ) %>%
    
    mutate(
      vf_category = stringr::str_trim(
        vf_category
      )
    ) %>%
    
    filter(
      !is.na(sample),
      sample != "",
      !is.na(vf_category),
      vf_category != ""
    ) %>%
    
    distinct(
      sample,
      vf_category
    )
  
  
} else {
  
  # ----------------------------------------------------------
  # CASO B:
  # VirulentHunter_category_hits.csv está en formato ancho
  #
  # Ejemplo:
  #
  # sample   Adherence   Biofilm   Exotoxin
  # NG02         1          1          0
  # NG03         1          0          1
  # ----------------------------------------------------------
  
  cat(
    "Formato VirulentHunter detectado: ANCHO\n"
  )
  
  
  # Columnas que NO son categorías
  
  exclude_cols <- c(
    vh_sample_col,
    "id",
    "ID",
    "protein",
    "Protein",
    "vf_prob",
    "probability",
    "Probability",
    "count",
    "Count",
    "total",
    "Total",
    "n"
  )
  
  
  possible_category_cols <- setdiff(
    names(vh_raw),
    exclude_cols
  )
  
  
  # Identificar columnas numéricas/lógicas
  
  category_is_usable <- vapply(
    
    possible_category_cols,
    
    function(nm) {
      
      x <- vh_raw[[nm]]
      
      is.numeric(x) ||
        is.integer(x) ||
        is.logical(x)
    },
    
    FUN.VALUE = logical(1)
  )
  
  
  numeric_category_cols <- possible_category_cols[
    category_is_usable
  ]
  
  
  if (length(numeric_category_cols) == 0) {
    
    stop(
      paste0(
        "No pude identificar automáticamente las categorías ",
        "de VirulentHunter.\n\n",
        "Columnas disponibles:\n",
        paste(
          names(vh_raw),
          collapse = ", "
        )
      )
    )
  }
  
  
  cat(
    "\nColumnas interpretadas como categorías:\n"
  )
  
  print(
    numeric_category_cols
  )
  
  
  vh_categories <- vh_raw %>%
    
    mutate(
      sample = clean_sample(
        .data[[vh_sample_col]]
      )
    ) %>%
    
    mutate(
      across(
        all_of(numeric_category_cols),
        ~ suppressWarnings(
          as.numeric(.x)
        )
      )
    ) %>%
    
    select(
      sample,
      all_of(numeric_category_cols)
    ) %>%
    
    pivot_longer(
      
      cols = all_of(
        numeric_category_cols
      ),
      
      names_to = "vf_category",
      
      values_to = "category_hit"
    ) %>%
    
    filter(
      !is.na(category_hit),
      category_hit > 0
    ) %>%
    
    distinct(
      sample,
      vf_category
    )
}


# ============================================================
# 13. COMPROBAR RESULTADO
# ============================================================

cat(
  "\nNúmero de muestras VirulentHunter:",
  n_distinct(vh_categories$sample),
  "\n"
)

cat(
  "Número de categorías:",
  n_distinct(vh_categories$vf_category),
  "\n"
)


cat(
  "\nCategorías detectadas:\n"
)

print(
  vh_categories %>%
    count(
      vf_category,
      sort = TRUE
    )
)


cat(
  "\nPrimeras 20 asociaciones muestra-categoría:\n"
)

print(
  head(
    vh_categories,
    20
  )
)

# ============================================================
# ============================================================
#
# PARTE C
# AMRFINDERPLUS / AMRgen
#
# ============================================================
# ============================================================



# ============================================================
# 14. IMPORTAR AMRFINDERPLUS
# ============================================================

amr <- AMRgen::import_amrfp(
  input_table = AMR_FILE,
  sample_col = "Name"
)


cat("\n========================================\n")
cat("AMRFINDERPLUS / AMRgen\n")
cat("========================================\n")

cat(
  "Filas importadas:",
  nrow(amr),
  "\n"
)


# ============================================================
# 15. NORMALIZAR IDs
# ============================================================

amr <- amr %>%
  mutate(
    sample =
      clean_sample(id)
  )


cat(
  "Muestras AMRgen:",
  n_distinct(amr$sample),
  "\n"
)



# ============================================================
# 16. EXTRAER CLASES DE RESISTENCIA
# ============================================================

amr_classes <- amr %>%
  
  transmute(
    sample,
    amr_class =
      as.character(
        drug_class
      )
  ) %>%
  
  mutate(
    amr_class =
      stringr::str_replace_all(
        amr_class,
        "[\\[\\]\"']",
        ""
      )
  ) %>%
  
  tidyr::separate_rows(
    amr_class,
    sep = "\\s*[;,|]+\\s*"
  ) %>%
  
  mutate(
    amr_class =
      stringr::str_trim(
        amr_class
      )
  ) %>%
  
  filter(
    !is.na(sample),
    sample != "",
    !is.na(amr_class),
    amr_class != ""
  ) %>%
  
  distinct(
    sample,
    amr_class
  )


cat("\nClases AMR detectadas:\n")

print(
  amr_classes %>%
    count(
      amr_class,
      sort = TRUE
    )
)



# ============================================================
# ============================================================
#
# PARTE D
# COMPROBAR QUE LOS IDs COINCIDEN
#
# ============================================================
# ============================================================



# ============================================================
# 17. LISTAS DE MUESTRAS
# ============================================================

samples_vfdb <- sort(
  unique(
    vfdb_hits$sample
  )
)


samples_vh <- sort(
  unique(
    vh_categories$sample
  )
)


samples_amr <- sort(
  unique(
    amr$sample
  )
)



# ============================================================
# 18. INTERSECCIÓN
# ============================================================

common_samples <- Reduce(
  intersect,
  list(
    samples_vfdb,
    samples_vh,
    samples_amr
  )
)


cat("\n========================================\n")
cat("COINCIDENCIA DE MUESTRAS\n")
cat("========================================\n")

cat(
  "VFDB:",
  length(samples_vfdb),
  "\n"
)

cat(
  "VirulentHunter:",
  length(samples_vh),
  "\n"
)

cat(
  "AMRgen:",
  length(samples_amr),
  "\n"
)

cat(
  "\nMuestras comunes:",
  length(common_samples),
  "\n"
)


print(
  common_samples
)



# ============================================================
# 19. MOSTRAR IDS QUE NO COINCIDEN
# ============================================================

cat(
  "\nVFDB que NO están en VirulentHunter:\n"
)

print(
  setdiff(
    samples_vfdb,
    samples_vh
  )
)


cat(
  "\nVirulentHunter que NO están en VFDB:\n"
)

print(
  setdiff(
    samples_vh,
    samples_vfdb
  )
)


cat(
  "\nVFDB que NO están en AMRgen:\n"
)

print(
  setdiff(
    samples_vfdb,
    samples_amr
  )
)


cat(
  "\nAMRgen que NO están en VFDB:\n"
)

print(
  setdiff(
    samples_amr,
    samples_vfdb
  )
)



# ============================================================
# 20. CONTROL
# ============================================================

if (length(common_samples) < 2) {
  
  stop(
    paste0(
      "Hay menos de dos muestras comunes entre los tres ",
      "análisis. Revisa los nombres antes de continuar."
    )
  )
}



# ============================================================
# ============================================================
#
# PARTE E
# SELECCIONAR GENES VFDB INFORMATIVOS
#
# ============================================================
# ============================================================



# ============================================================
# 21. PREVALENCIA EN LAS MUESTRAS COMUNES
# ============================================================

vfdb_common_frequency <- vfdb_hits %>%
  
  filter(
    sample %in%
      common_samples
  ) %>%
  
  distinct(
    sample,
    vfdb_gene
  ) %>%
  
  count(
    vfdb_gene,
    sort = TRUE
  ) %>%
  
  mutate(
    
    prevalence =
      n /
      length(common_samples),
    
    prevalence_percent =
      prevalence * 100,
    
    # Máximo cuando prevalencia = 0.5
    # Sirve para seleccionar genes variables/informativos.
    
    variability_score =
      prevalence *
      (1 - prevalence)
  )



# ============================================================
# 22. QUITAR GENES UBIQUOS
# ============================================================

vfdb_variable <- vfdb_common_frequency %>%
  
  filter(
    n > 0,
    n < length(common_samples)
  ) %>%
  
  arrange(
    desc(variability_score),
    desc(n)
  )


cat(
  "\nGenes VFDB variables:",
  nrow(vfdb_variable),
  "\n"
)


print(
  vfdb_variable
)



# ============================================================
# 23. ELEGIR LOS 10 GENES MÁS INFORMATIVOS
# ============================================================
#
# Puedes aumentar TOP_VFDB después si la figura queda limpia.
# ============================================================

TOP_VFDB <- 10


if (nrow(vfdb_variable) > 0) {
  
  selected_vfdb_genes <- vfdb_variable %>%
    
    slice_head(
      n = TOP_VFDB
    ) %>%
    
    pull(
      vfdb_gene
    )
  
} else {
  
  warning(
    paste(
      "Todos los genes VFDB son ubicuos.",
      "Se usarán los 10 genes más frecuentes."
    )
  )
  
  
  selected_vfdb_genes <- vfdb_common_frequency %>%
    
    slice_head(
      n = TOP_VFDB
    ) %>%
    
    pull(
      vfdb_gene
    )
}


cat(
  "\nGenes VFDB utilizados en el Sankey:\n"
)

print(
  selected_vfdb_genes
)



# ============================================================
# ============================================================
#
# PARTE F
# CREAR RUTAS DEL SANKEY
#
# ============================================================
# ============================================================



# ============================================================
# 24. VFDB
# ============================================================

vfdb_sankey <- vfdb_hits %>%
  
  filter(
    sample %in%
      common_samples,
    
    vfdb_gene %in%
      selected_vfdb_genes
  ) %>%
  
  distinct(
    sample,
    vfdb_gene
  )



# ============================================================
# 25. VIRULENTHUNTER
# ============================================================

vh_sankey <- vh_categories %>%
  
  filter(
    sample %in%
      common_samples
  ) %>%
  
  distinct(
    sample,
    vf_category
  )



# ============================================================
# 26. AMR
# ============================================================

amr_sankey <- amr_classes %>%
  
  filter(
    sample %in%
      common_samples
  ) %>%
  
  distinct(
    sample,
    amr_class
  )



# ============================================================
# 27. CREAR TODAS LAS COOCURRENCIAS DENTRO DEL MISMO GENOMA
# ============================================================
#
# Ejemplo:
#
# NG02:
#
# pilE
# Adherence
# Tetracyclines
#
# genera una ruta:
#
# pilE -> Adherence -> Tetracyclines
#
# ============================================================

paths <- vfdb_sankey %>%
  
  inner_join(
    vh_sankey,
    by = "sample",
    relationship = "many-to-many"
  ) %>%
  
  inner_join(
    amr_sankey,
    by = "sample",
    relationship = "many-to-many"
  )


cat(
  "\nRutas crudas:",
  nrow(paths),
  "\n"
)



# ============================================================
# 28. NORMALIZAR EL PESO POR GENOMA
# ============================================================
#
# MUY IMPORTANTE:
#
# Un genoma con:
#
# 10 genes VFDB
# 8 categorías VH
# 5 clases AMR
#
# produciría muchísimas combinaciones.
#
# Por eso cada genoma debe contribuir en TOTAL con peso = 1.
#
# ============================================================

paths <- paths %>%
  
  group_by(
    sample
  ) %>%
  
  mutate(
    weight =
      1 / n()
  ) %>%
  
  ungroup()



# ============================================================
# 29. COMPROBAR NORMALIZACIÓN
# ============================================================

weight_check <- paths %>%
  
  group_by(
    sample
  ) %>%
  
  summarise(
    total_weight =
      sum(weight),
    .groups = "drop"
  )


cat(
  "\nComprobación de pesos:\n"
)

print(
  weight_check
)



# ============================================================
# 30. AGREGAR RUTAS IGUALES
# ============================================================

sankey_data <- paths %>%
  
  group_by(
    vfdb_gene,
    vf_category,
    amr_class
  ) %>%
  
  summarise(
    
    weight =
      sum(weight),
    
    n_isolates =
      n_distinct(sample),
    
    .groups = "drop"
  ) %>%
  
  arrange(
    desc(n_isolates),
    desc(weight)
  )



# ============================================================
# 31. GUARDAR TABLA COMPLETA
# ============================================================

write.csv(
  sankey_data,
  "Sankey_all_cooccurrences.csv",
  row.names = FALSE
)



# ============================================================
# 32. FILTRAR RUTAS MUY RARAS PARA LA FIGURA
# ============================================================
#
# Con 24 genomas empezamos requiriendo presencia en >=2.
#
# Si queda demasiado vacía:
#
# MIN_ISOLATES <- 1
#
# Si queda demasiado saturada:
#
# MIN_ISOLATES <- 3
#
# ============================================================

MIN_ISOLATES <- 2


sankey_plot_data <- sankey_data %>%
  
  filter(
    n_isolates >=
      MIN_ISOLATES
  )


cat(
  "\nRutas totales:",
  nrow(sankey_data),
  "\n"
)

cat(
  "Rutas mostradas en la figura:",
  nrow(sankey_plot_data),
  "\n"
)


if (nrow(sankey_plot_data) == 0) {
  
  stop(
    paste0(
      "Ninguna ruta cumple MIN_ISOLATES = ",
      MIN_ISOLATES,
      ". Cambia MIN_ISOLATES a 1."
    )
  )
}



# ============================================================
# ============================================================
#
# PARTE G
# SANKEY / ALLUVIAL
#
# ============================================================
# ============================================================



# ============================================================
# 33. CONSTRUIR FIGURA
# ============================================================

p_sankey <- ggplot(
  
  sankey_plot_data,
  
  aes(
    
    y =
      weight,
    
    axis1 =
      vfdb_gene,
    
    axis2 =
      vf_category,
    
    axis3 =
      amr_class
  )
  
) +
  
  
  geom_alluvium(
    
    aes(
      fill =
        vf_category
    ),
    
    width =
      0.10,
    
    alpha =
      0.75
  ) +
  
  
  geom_stratum(
    
    width =
      0.14,
    
    colour =
      "grey30",
    
    linewidth =
      0.35
  ) +
  
  
  geom_text(
    
    stat =
      "stratum",
    
    aes(
      label =
        after_stat(stratum)
    ),
    
    size =
      3
  ) +
  
  
  scale_x_discrete(
    
    limits = c(
      "VFDB genes",
      "VirulentHunter",
      "AMR classes"
    ),
    
    expand =
      c(
        0.08,
        0.08
      )
  ) +
  
  
  labs(
    
    title =
      "Integrated virulence and antimicrobial resistance profiles",
    
    subtitle =
      "Neisseria gonorrhoeae",
    
    x =
      NULL,
    
    y =
      "Normalized co-occurrence weight",
    
    fill =
      "Virulence category",
    
    caption =
      paste0(
        "Flows represent genomic co-occurrence within the same isolate. ",
        "Each isolate contributes a total normalized weight of 1. ",
        "Connections do not imply causal or direct molecular associations."
      )
  ) +
  
  
  theme_minimal(
    base_size =
      12
  ) +
  
  
  theme(
    
    panel.grid =
      element_blank(),
    
    axis.text.x =
      element_text(
        face = "bold",
        size = 12
      ),
    
    axis.text.y =
      element_blank(),
    
    axis.ticks.y =
      element_blank(),
    
    legend.position =
      "right",
    
    plot.title =
      element_text(
        face = "bold"
      )
  )

p_sankey

# ============================================================
# 34. MOSTRAR
# ============================================================

print(
  p_sankey
)



# ============================================================
# 35. EXPORTAR
# ============================================================

ggsave(
  
  filename =
    "Sankey_VFDB_VirulentHunter_AMR.pdf",
  
  plot =
    p_sankey,
  
  width =
    16,
  
  height =
    10
)


ggsave(
  
  filename =
    "Sankey_VFDB_VirulentHunter_AMR.png",
  
  plot =
    p_sankey,
  
  width =
    16,
  
  height =
    10,
  
  dpi =
    600
)



# ============================================================
# 36. EXPORTAR LOS DATOS UTILIZADOS EN LA FIGURA
# ============================================================

write.csv(
  
  sankey_plot_data,
  
  "Sankey_plotted_cooccurrences.csv",
  
  row.names =
    FALSE
)



# ============================================================
# 37. RESUMEN FINAL
# ============================================================

cat("\n")
cat("========================================\n")
cat("RESUMEN FINAL\n")
cat("========================================\n")

cat(
  "Muestras VFDB:",
  length(samples_vfdb),
  "\n"
)

cat(
  "Muestras VirulentHunter:",
  length(samples_vh),
  "\n"
)

cat(
  "Muestras AMRgen:",
  length(samples_amr),
  "\n"
)

cat(
  "Muestras comunes:",
  length(common_samples),
  "\n"
)

cat(
  "Genes VFDB seleccionados:",
  length(selected_vfdb_genes),
  "\n"
)

cat(
  "Categorías VirulentHunter:",
  n_distinct(vh_sankey$vf_category),
  "\n"
)

cat(
  "Clases AMR:",
  n_distinct(amr_sankey$amr_class),
  "\n"
)

cat(
  "Rutas totales:",
  nrow(sankey_data),
  "\n"
)

cat(
  "Rutas mostradas:",
  nrow(sankey_plot_data),
  "\n"
)

cat("\nArchivos generados:\n")

cat(
  "  Sankey_VFDB_VirulentHunter_AMR.pdf\n"
)

cat(
  "  Sankey_VFDB_VirulentHunter_AMR.png\n"
)

cat(
  "  Sankey_all_cooccurrences.csv\n"
)

cat(
  "  Sankey_plotted_cooccurrences.csv\n"
)

cat("========================================\n")

######################################################
# Segunda parte
# ============================================================
# ANÁLISIS DE ASOCIACIÓN VIRULENCIA <-> AMR
# Neisseria gonorrhoeae
#
# Requiere que ya existan:
#
#   vfdb_hits
#   vh_categories
#   amr_classes
#   common_samples
#
# Produce:
#
# 1) Fisher exact test
# 2) Odds ratio
# 3) Coeficiente phi
# 4) Corrección Benjamini-Hochberg (FDR)
#
# 5) Sankey:
#       VFDB genes -> AMR classes
#
# 6) Sankey:
#       VirulentHunter categories -> AMR classes
#
# IMPORTANTE:
# Una asociación indica coocurrencia estadística entre
# características genómicas. No demuestra causalidad.
# ============================================================


# ============================================================
# 0. PAQUETES
# ============================================================

library(tidyverse)
library(ggalluvial)


# Opcional, para juntar figuras al final
#if (!requireNamespace("patchwork", quietly = TRUE)) {
#  install.packages("patchwork")
#}

library(patchwork)


# ============================================================
# 1. COMPROBAR OBJETOS NECESARIOS
# ============================================================

required_objects <- c(
  "vfdb_hits",
  "vh_categories",
  "amr_classes",
  "common_samples"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1)
  )
]

if (length(missing_objects) > 0) {
  
  stop(
    paste0(
      "Faltan estos objetos:\n",
      paste(missing_objects, collapse = ", "),
      "\n\nEjecuta primero el script de integración."
    )
  )
}


cat(
  "\nNúmero de aislamientos comunes:",
  length(common_samples),
  "\n"
)


# ============================================================
# 2. TABLA BASE DE MUESTRAS
# ============================================================

base_samples <- tibble(
  sample = common_samples
)


# ============================================================
# ============================================================
#
# PARTE A
# MATRIZ BINARIA VFDB
#
# ============================================================
# ============================================================


# ============================================================
# 3. VFDB: AISLADO x GEN
# ============================================================

vfdb_binary <- vfdb_hits %>%
  
  filter(
    sample %in% common_samples
  ) %>%
  
  distinct(
    sample,
    vfdb_gene
  ) %>%
  
  mutate(
    present = 1L
  ) %>%
  
  pivot_wider(
    names_from = vfdb_gene,
    values_from = present,
    values_fill = 0
  )


# Añadir cualquier muestra que pudiera no tener hits

vfdb_binary <- base_samples %>%
  
  left_join(
    vfdb_binary,
    by = "sample"
  ) %>%
  
  mutate(
    across(
      -sample,
      ~ replace_na(.x, 0L)
    )
  )


cat(
  "\nMatriz VFDB:",
  nrow(vfdb_binary),
  "aislamientos x",
  ncol(vfdb_binary) - 1,
  "genes\n"
)


# ============================================================
# ============================================================
#
# PARTE B
# MATRIZ BINARIA VIRULENTHUNTER
#
# ============================================================
# ============================================================


# ============================================================
# 4. VIRULENTHUNTER: AISLADO x CATEGORÍA
# ============================================================

vh_binary <- vh_categories %>%
  
  filter(
    sample %in% common_samples
  ) %>%
  
  distinct(
    sample,
    vf_category
  ) %>%
  
  mutate(
    present = 1L
  ) %>%
  
  pivot_wider(
    names_from = vf_category,
    values_from = present,
    values_fill = 0
  )


vh_binary <- base_samples %>%
  
  left_join(
    vh_binary,
    by = "sample"
  ) %>%
  
  mutate(
    across(
      -sample,
      ~ replace_na(.x, 0L)
    )
  )


cat(
  "Matriz VirulentHunter:",
  nrow(vh_binary),
  "aislamientos x",
  ncol(vh_binary) - 1,
  "categorías\n"
)


# ============================================================
# ============================================================
#
# PARTE C
# MATRIZ BINARIA AMR
#
# ============================================================
# ============================================================


# ============================================================
# 5. AMR: AISLADO x CLASE
# ============================================================

amr_binary <- amr_classes %>%
  
  filter(
    sample %in% common_samples
  ) %>%
  
  distinct(
    sample,
    amr_class
  ) %>%
  
  mutate(
    present = 1L
  ) %>%
  
  pivot_wider(
    names_from = amr_class,
    values_from = present,
    values_fill = 0
  )


amr_binary <- base_samples %>%
  
  left_join(
    amr_binary,
    by = "sample"
  ) %>%
  
  mutate(
    across(
      -sample,
      ~ replace_na(.x, 0L)
    )
  )


cat(
  "Matriz AMR:",
  nrow(amr_binary),
  "aislamientos x",
  ncol(amr_binary) - 1,
  "clases\n"
)


# ============================================================
# 6. GUARDAR MATRICES
# ============================================================

write.csv(
  vfdb_binary,
  "VFDB_binary_for_association.csv",
  row.names = FALSE
)

write.csv(
  vh_binary,
  "VirulentHunter_binary_for_association.csv",
  row.names = FALSE
)

write.csv(
  amr_binary,
  "AMR_binary_for_association.csv",
  row.names = FALSE
)


# ============================================================
# ============================================================
#
# PARTE D
# FUNCIÓN ESTADÍSTICA
#
# ============================================================
# ============================================================


# ============================================================
# 7. FUNCIÓN PARA FISHER + OR + PHI
# ============================================================

test_binary_association <- function(x, y) {
  
  x <- as.integer(x)
  y <- as.integer(y)
  
  tab <- table(
    factor(x, levels = c(0, 1)),
    factor(y, levels = c(0, 1))
  )
  
  
  # ----------------------------------------------------------
  # Tabla:
  #
  #             AMR-
  #             AMR+
  #
  # VF-    d     c
  # VF+    b     a
  #
  # ----------------------------------------------------------
  
  d <- unname(tab[1, 1])
  c <- unname(tab[1, 2])
  b <- unname(tab[2, 1])
  a <- unname(tab[2, 2])
  
  
  fisher_result <- fisher.test(tab)
  
  
  denominator <- sqrt(
    (a + b) *
      (c + d) *
      (a + c) *
      (b + d)
  )
  
  
  phi <- if (
    denominator == 0
  ) {
    NA_real_
  } else {
    (a * d - b * c) / denominator
  }
  
  
  tibble(
    
    n_both_present = a,
    
    n_feature_only = b,
    
    n_amr_only = c,
    
    n_neither = d,
    
    odds_ratio =
      unname(
        fisher_result$estimate
      ),
    
    p_value =
      fisher_result$p.value,
    
    phi =
      phi
  )
}


# ============================================================
# ============================================================
#
# PARTE E
# VFDB <-> AMR
#
# ============================================================
# ============================================================


# ============================================================
# 8. TODAS LAS COMBINACIONES VFDB x AMR
# ============================================================

vfdb_names <- setdiff(
  names(vfdb_binary),
  "sample"
)


amr_names <- setdiff(
  names(amr_binary),
  "sample"
)


vfdb_amr_results <- tidyr::crossing(
  
  vfdb_gene =
    vfdb_names,
  
  amr_class =
    amr_names
  
) %>%
  
  rowwise() %>%
  
  mutate(
    
    result = list(
      
      test_binary_association(
        
        vfdb_binary[[vfdb_gene]],
        
        amr_binary[[amr_class]]
      )
    )
  ) %>%
  
  ungroup() %>%
  
  unnest(
    result
  ) %>%
  
  mutate(
    
    p_adjusted =
      p.adjust(
        p_value,
        method = "BH"
      ),
    
    abs_phi =
      abs(phi)
    
  ) %>%
  
  arrange(
    p_adjusted,
    desc(abs_phi)
  )


# ============================================================
# 9. GUARDAR TODOS LOS RESULTADOS
# ============================================================

write.csv(
  vfdb_amr_results,
  "VFDB_AMR_all_associations.csv",
  row.names = FALSE
)


cat(
  "\nTop asociaciones VFDB <-> AMR:\n"
)

print(
  head(
    vfdb_amr_results,
    20
  )
)


# ============================================================
# ============================================================
#
# PARTE F
# VIRULENTHUNTER <-> AMR
#
# ============================================================
# ============================================================


# ============================================================
# 10. TODAS LAS COMBINACIONES VH x AMR
# ============================================================

vh_names <- setdiff(
  names(vh_binary),
  "sample"
)


vh_amr_results <- tidyr::crossing(
  
  vf_category =
    vh_names,
  
  amr_class =
    amr_names
  
) %>%
  
  rowwise() %>%
  
  mutate(
    
    result = list(
      
      test_binary_association(
        
        vh_binary[[vf_category]],
        
        amr_binary[[amr_class]]
      )
    )
  ) %>%
  
  ungroup() %>%
  
  unnest(
    result
  ) %>%
  
  mutate(
    
    p_adjusted =
      p.adjust(
        p_value,
        method = "BH"
      ),
    
    abs_phi =
      abs(phi)
    
  ) %>%
  
  arrange(
    p_adjusted,
    desc(abs_phi)
  )


# ============================================================
# 11. GUARDAR
# ============================================================

write.csv(
  vh_amr_results,
  "VirulentHunter_AMR_all_associations.csv",
  row.names = FALSE
)


cat(
  "\nTop asociaciones VirulentHunter <-> AMR:\n"
)

print(
  head(
    vh_amr_results,
    20
  )
)


# ============================================================
# ============================================================
#
# PARTE G
# ASOCIACIONES ESTADÍSTICAMENTE SIGNIFICATIVAS
#
# ============================================================
# ============================================================


# ============================================================
# 12. CRITERIOS
# ============================================================
#
# Con solamente ~24 genomas, Fisher + FDR puede ser
# bastante conservador.
#
# Para considerar una asociación confirmatoria usamos:
#
# FDR < 0.05
# phi > 0
# al menos 3 aislamientos con ambas características
#
# ============================================================

FDR_LIMIT <- 0.05
MIN_BOTH <- 3


vfdb_amr_significant <- vfdb_amr_results %>%
  
  filter(
    !is.na(phi),
    phi > 0,
    n_both_present >= MIN_BOTH,
    p_adjusted < FDR_LIMIT
  )


vh_amr_significant <- vh_amr_results %>%
  
  filter(
    !is.na(phi),
    phi > 0,
    n_both_present >= MIN_BOTH,
    p_adjusted < FDR_LIMIT
  )


cat(
  "\n========================================\n"
)

cat(
  "ASOCIACIONES SIGNIFICATIVAS\n"
)

cat(
  "========================================\n"
)

cat(
  "VFDB <-> AMR:",
  nrow(vfdb_amr_significant),
  "\n"
)

cat(
  "VirulentHunter <-> AMR:",
  nrow(vh_amr_significant),
  "\n"
)


# ============================================================
# 13. GUARDAR SIGNIFICATIVAS
# ============================================================

write.csv(
  vfdb_amr_significant,
  "VFDB_AMR_significant_FDR005.csv",
  row.names = FALSE
)

write.csv(
  vh_amr_significant,
  "VirulentHunter_AMR_significant_FDR005.csv",
  row.names = FALSE
)


# ============================================================
# ============================================================
#
# PARTE H
# PREPARAR FIGURAS
#
# ============================================================
# ============================================================


# ============================================================
# 14. DECIDIR QUÉ MOSTRAR
# ============================================================
#
# Si hay asociaciones FDR < 0.05:
# -> mostramos esas asociaciones.
#
# Si NO hay suficientes, mostramos asociaciones exploratorias
# con mayor phi.
#
# Esto NO cambia el análisis estadístico.
# Sólo permite visualizar patrones candidatos.
#
# ============================================================

TOP_EXPLORATORY <- 15


if (nrow(vfdb_amr_significant) >= 2) {
  
  vfdb_plot_data <- vfdb_amr_significant %>%
    
    mutate(
      evidence =
        "FDR < 0.05"
    )
  
  VFDB_PLOT_TYPE <-
    "FDR-significant associations"
  
} else {
  
  warning(
    paste(
      "Hay pocas o ninguna asociación VFDB-AMR",
      "significativa después de FDR.",
      "La figura VFDB será EXPLORATORIA."
    )
  )
  
  vfdb_plot_data <- vfdb_amr_results %>%
    
    filter(
      !is.na(phi),
      phi > 0,
      n_both_present >= MIN_BOTH
    ) %>%
    
    arrange(
      desc(phi),
      p_value
    ) %>%
    
    slice_head(
      n = TOP_EXPLORATORY
    ) %>%
    
    mutate(
      evidence =
        "Exploratory"
    )
  
  VFDB_PLOT_TYPE <-
    "Top exploratory associations"
}



if (nrow(vh_amr_significant) >= 2) {
  
  vh_plot_data <- vh_amr_significant %>%
    
    mutate(
      evidence =
        "FDR < 0.05"
    )
  
  VH_PLOT_TYPE <-
    "FDR-significant associations"
  
} else {
  
  warning(
    paste(
      "Hay pocas o ninguna asociación VirulentHunter-AMR",
      "significativa después de FDR.",
      "La figura VirulentHunter será EXPLORATORIA."
    )
  )
  
  vh_plot_data <- vh_amr_results %>%
    
    filter(
      !is.na(phi),
      phi > 0,
      n_both_present >= MIN_BOTH
    ) %>%
    
    arrange(
      desc(phi),
      p_value
    ) %>%
    
    slice_head(
      n = TOP_EXPLORATORY
    ) %>%
    
    mutate(
      evidence =
        "Exploratory"
    )
  
  VH_PLOT_TYPE <-
    "Top exploratory associations"
}


# ============================================================
# 15. PESO PARA EL SANKEY
# ============================================================
#
# Utilizamos phi como magnitud del flujo.
#
# phi:
#
# 0       = sin asociación
# ~0.1    = débil
# ~0.3    = moderada
# ~0.5+   = fuerte
#
# El ancho NO representa número de genomas.
# Representa fuerza de asociación.
#
# ============================================================

vfdb_plot_data <- vfdb_plot_data %>%
  
  mutate(
    flow_weight =
      phi
  )


vh_plot_data <- vh_plot_data %>%
  
  mutate(
    flow_weight =
      phi
  )


# ============================================================
# ============================================================
#
# PARTE I
# SANKEY VFDB -> AMR
#
# ============================================================
# ============================================================


# ============================================================
# 16. FIGURA VFDB
# ============================================================

p_vfdb_amr <- ggplot(
  
  vfdb_plot_data,
  
  aes(
    y =
      flow_weight,
    
    axis1 =
      vfdb_gene,
    
    axis2 =
      amr_class
  )
  
) +
  
  geom_alluvium(
    
    aes(
      fill =
        amr_class
    ),
    
    width =
      0.14,
    
    alpha =
      0.75
  ) +
  
  geom_stratum(
    
    width =
      0.18,
    
    colour =
      "grey30",
    
    linewidth =
      0.35
  ) +
  
  geom_text(
    
    stat =
      "stratum",
    
    aes(
      label =
        after_stat(stratum)
    ),
    
    size =
      3.2
  ) +
  
  scale_x_discrete(
    
    limits =
      c(
        "VFDB genes",
        "AMR classes"
      ),
    
    expand =
      c(
        0.12,
        0.12
      )
  ) +
  
  labs(
    
    title =
      "VFDB–AMR associations",
    
    subtitle =
      VFDB_PLOT_TYPE,
    
    x =
      NULL,
    
    y =
      "Association strength (phi)",
    
    fill =
      "AMR class",
    
    caption =
      paste0(
        "Links represent positive genomic associations. ",
        "Fisher's exact test with Benjamini-Hochberg correction. ",
        "Link width represents phi coefficient."
      )
  ) +
  
  theme_minimal(
    base_size = 11
  ) +
  
  theme(
    
    panel.grid =
      element_blank(),
    
    axis.text.x =
      element_text(
        face = "bold"
      ),
    
    axis.text.y =
      element_blank(),
    
    axis.ticks.y =
      element_blank(),
    
    legend.position =
      "right"
  )
p_vfdb_amr

# ============================================================
# ============================================================
#
# PARTE J
# SANKEY VIRULENTHUNTER -> AMR
#
# ============================================================
# ============================================================


# ============================================================
# 17. FIGURA VIRULENTHUNTER
# ============================================================

p_vh_amr <- ggplot(
  
  vh_plot_data,
  
  aes(
    y =
      flow_weight,
    
    axis1 =
      vf_category,
    
    axis2 =
      amr_class
  )
  
) +
  
  geom_alluvium(
    
    aes(
      fill =
        amr_class
    ),
    
    width =
      0.14,
    
    alpha =
      0.75
  ) +
  
  geom_stratum(
    
    width =
      0.18,
    
    colour =
      "grey30",
    
    linewidth =
      0.35
  ) +
  
  geom_text(
    
    stat =
      "stratum",
    
    aes(
      label =
        after_stat(stratum)
    ),
    
    size =
      3.1
  ) +
  
  scale_x_discrete(
    
    limits =
      c(
        "VirulentHunter",
        "AMR classes"
      ),
    
    expand =
      c(
        0.12,
        0.12
      )
  ) +
  
  labs(
    
    title =
      "VirulentHunter–AMR associations",
    
    subtitle =
      VH_PLOT_TYPE,
    
    x =
      NULL,
    
    y =
      "Association strength (phi)",
    
    fill =
      "AMR class",
    
    caption =
      paste0(
        "Links represent positive genomic associations. ",
        "Fisher's exact test with Benjamini-Hochberg correction. ",
        "Link width represents phi coefficient."
      )
  ) +
  
  theme_minimal(
    base_size = 11
  ) +
  
  theme(
    
    panel.grid =
      element_blank(),
    
    axis.text.x =
      element_text(
        face = "bold"
      ),
    
    axis.text.y =
      element_blank(),
    
    axis.ticks.y =
      element_blank(),
    
    legend.position =
      "right"
  )
p_vh_amr

# ============================================================
# 18. MOSTRAR POR SEPARADO
# ============================================================

print(
  p_vfdb_amr
)

print(
  p_vh_amr
)


# ============================================================
# ============================================================
#
# PARTE K
# FIGURA COMBINADA
#
# ============================================================
# ============================================================


# ============================================================
# 19. COMBINAR PANELES
# ============================================================

p_combined <- (
  
  p_vfdb_amr +
    
    p_vh_amr +
    
    plot_layout(
      ncol = 2,
      guides = "collect"
    )
  
) &
  
  theme(
    legend.position = "bottom"
  )


p_combined


# ============================================================
# 20. EXPORTAR
# ============================================================

ggsave(
  
  "Sankey_associations_VFDB_AMR.pdf",
  
  p_vfdb_amr,
  
  width = 9,
  
  height = 8
)


ggsave(
  
  "Sankey_associations_VFDB_AMR.png",
  
  p_vfdb_amr,
  
  width = 9,
  
  height = 8,
  
  dpi = 600
)


ggsave(
  
  "Sankey_associations_VirulentHunter_AMR.pdf",
  
  p_vh_amr,
  
  width = 10,
  
  height = 8
)


ggsave(
  
  "Sankey_associations_VirulentHunter_AMR.png",
  
  p_vh_amr,
  
  width = 10,
  
  height = 8,
  
  dpi = 600
)


ggsave(
  
  "Sankey_associations_combined.pdf",
  
  p_combined,
  
  width = 18,
  
  height = 8
)


ggsave(
  
  "Sankey_associations_combined.png",
  
  p_combined,
  
  width = 18,
  
  height = 8,
  
  dpi = 600
)


# ============================================================
# ============================================================
#
# PARTE L
# RESUMEN ESTADÍSTICO
#
# ============================================================
# ============================================================


cat("\n")
cat("========================================\n")
cat("RESUMEN DE ASOCIACIONES\n")
cat("========================================\n")

cat(
  "Aislamientos:",
  length(common_samples),
  "\n\n"
)

cat(
  "Pruebas VFDB-AMR:",
  nrow(vfdb_amr_results),
  "\n"
)

cat(
  "VFDB-AMR FDR < 0.05:",
  nrow(vfdb_amr_significant),
  "\n\n"
)

cat(
  "Pruebas VirulentHunter-AMR:",
  nrow(vh_amr_results),
  "\n"
)

cat(
  "VirulentHunter-AMR FDR < 0.05:",
  nrow(vh_amr_significant),
  "\n\n"
)

cat(
  "Tipo de figura VFDB:",
  VFDB_PLOT_TYPE,
  "\n"
)

cat(
  "Tipo de figura VirulentHunter:",
  VH_PLOT_TYPE,
  "\n"
)

cat("========================================\n")


# ============================================================
# 21. MOSTRAR LAS ASOCIACIONES MÁS FUERTES
# ============================================================

cat(
  "\nTOP VFDB-AMR POR PHI:\n"
)

print(
  
  vfdb_amr_results %>%
    
    filter(
      !is.na(phi),
      phi > 0
    ) %>%
    
    arrange(
      desc(phi)
    ) %>%
    
    select(
      vfdb_gene,
      amr_class,
      n_both_present,
      odds_ratio,
      phi,
      p_value,
      p_adjusted
    ) %>%
    
    head(20)
)


cat(
  "\nTOP VIRULENTHUNTER-AMR POR PHI:\n"
)

print(
  
  vh_amr_results %>%
    
    filter(
      !is.na(phi),
      phi > 0
    ) %>%
    
    arrange(
      desc(phi)
    ) %>%
    
    select(
      vf_category,
      amr_class,
      n_both_present,
      odds_ratio,
      phi,
      p_value,
      p_adjusted
    ) %>%
    
    head(20)
)


###############################33
####################################3
#################################

# ============================================================
# ¿CUÁNTO VARÍAN REALMENTE LAS CATEGORÍAS VIRULENTHUNTER?
# ============================================================

vh_prevalence <- vh_categories %>%
  filter(sample %in% common_samples) %>%
  distinct(sample, vf_category) %>%
  count(vf_category, name = "n_present") %>%
  mutate(
    total = length(common_samples),
    prevalence = n_present / total,
    prevalence_percent = 100 * prevalence
  ) %>%
  arrange(desc(prevalence_percent))

vh_prevalence


# ============================================================
# VARIABILIDAD DE LAS CLASES AMR
# ============================================================

amr_prevalence <- amr_classes %>%
  filter(sample %in% common_samples) %>%
  distinct(sample, amr_class) %>%
  count(amr_class, name = "n_present") %>%
  mutate(
    total = length(common_samples),
    prevalence = n_present / total,
    prevalence_percent = 100 * prevalence
  ) %>%
  arrange(desc(prevalence_percent))

amr_prevalence
#################################################################
############################################3
##############################################3
##############################################33

# ============================================================
# VIRULENTHUNTER CUANTITATIVO vs AMRgen / AMRFinderPlus
# Neisseria gonorrhoeae
#
# OBJETIVO:
#
# Sustituir el análisis binario:
#
#   "¿Está presente Biofilm? sí/no"
#
# por:
#
#   "¿Cuánto representa Biofilm en cada aislamiento?"
#
# y relacionarlo con:
#
#   "¿Cuántos determinantes AMR de cada clase tiene?"
#
#
# ANÁLISIS:
#
# 1. Conteos VirulentHunter por categoría
# 2. Proporción relativa de cada categoría
# 3. Número de marcadores AMR por clase
# 4. Spearman:
#
#      categoría VirulentHunter
#               vs
#      carga de marcadores AMR
#
# 5. Wilcoxon exploratorio para clases AMR variables
# 6. Heatmap de rho
# 7. Sankey reducido con asociaciones más fuertes
#
# ============================================================



# ============================================================
# 0. PAQUETES
# ============================================================

library(tidyverse)
library(ggalluvial)
library(AMRgen)



# ============================================================
# 1. ARCHIVOS
# ============================================================

VH_COUNTS_FILE <- paste0(
  "/home/user/Escritorio/neisseria/pangenomics/",
  "gonorrhoeae_ncbi/complete_genomes/Neisseria/",
  "virulenthunter_results/procesadas_mac/",
  "VirulentHunter_category_counts.csv"
)


AMR_FILE <- paste0(
  "/home/user/Escritorio/neisseria/pangenomics/",
  "gonorrhoeae_ncbi/complete_genomes/",
  "amrfinder_Neisseria_gonorrhoeae.tsv"
)



# ============================================================
# 2. COMPROBAR ARCHIVOS
# ============================================================

cat(
  "VirulentHunter:",
  file.exists(VH_COUNTS_FILE),
  "\n"
)

cat(
  "AMRFinderPlus:",
  file.exists(AMR_FILE),
  "\n"
)


if (!file.exists(VH_COUNTS_FILE)) {
  stop("No se encontró VirulentHunter_category_counts.csv")
}

if (!file.exists(AMR_FILE)) {
  stop("No se encontró el TSV combinado de AMRFinderPlus")
}



# ============================================================
# 3. FUNCIÓN PARA LIMPIAR IDs
# ============================================================

clean_sample <- function(x) {
  
  x <- basename(
    as.character(x)
  )
  
  x <- stringr::str_trim(x)
  
  x <- stringr::str_remove(
    x,
    "\\.(csv|tsv|txt|fna|faa|fasta|fa)(\\.(faa|fna))?$"
  )
  
  x <- stringr::str_remove(
    x,
    "_final$"
  )
  
  x
}



# ============================================================
# ============================================================
#
# PARTE A
# VIRULENTHUNTER CATEGORY COUNTS
#
# ============================================================
# ============================================================



# ============================================================
# 4. IMPORTAR
# ============================================================

vh_raw <- readr::read_csv(
  VH_COUNTS_FILE,
  show_col_types = FALSE,
  name_repair = "minimal"
)


cat("\n======================================\n")
cat("VIRULENTHUNTER CATEGORY COUNTS\n")
cat("======================================\n")

cat(
  "Dimensiones:",
  nrow(vh_raw),
  "x",
  ncol(vh_raw),
  "\n\n"
)

print(
  names(vh_raw)
)



# ============================================================
# 5. DETECTAR COLUMNA DE MUESTRA
# ============================================================

sample_candidates <- c(
  "sample",
  "Sample",
  "genome",
  "Genome",
  "genome_id",
  "Genome_ID",
  "isolate",
  "Isolate",
  "strain",
  "Strain",
  "file",
  "File",
  "filename",
  "Filename"
)


vh_sample_col <- intersect(
  sample_candidates,
  names(vh_raw)
)


if (length(vh_sample_col) == 0) {
  
  vh_sample_col <- grep(
    "sample|genome|isolate|strain|file",
    names(vh_raw),
    ignore.case = TRUE,
    value = TRUE
  )
}


# Si aún no encuentra nada, usamos la primera columna

if (length(vh_sample_col) == 0) {
  
  warning(
    "No identifiqué automáticamente la columna de muestra. Se usará la primera columna."
  )
  
  vh_sample_col <- names(vh_raw)[1]
  
} else {
  
  vh_sample_col <- vh_sample_col[1]
}


cat(
  "\nColumna de muestra:",
  vh_sample_col,
  "\n"
)



# ============================================================
# 6. DETERMINAR SI LA TABLA ESTÁ EN FORMATO LARGO
# ============================================================

category_candidates <- c(
  "category",
  "Category",
  "vf_category",
  "VF_category",
  "VF_Category"
)


count_candidates <- c(
  "count",
  "Count",
  "counts",
  "Counts",
  "n",
  "N",
  "hits",
  "Hits",
  "n_hits",
  "VF_count"
)


category_col <- intersect(
  category_candidates,
  names(vh_raw)
)


count_col <- intersect(
  count_candidates,
  names(vh_raw)
)



# ============================================================
# 7A. SI YA ESTÁ EN FORMATO LARGO
# ============================================================

if (
  length(category_col) > 0 &&
  length(count_col) > 0
) {
  
  category_col <- category_col[1]
  count_col <- count_col[1]
  
  
  cat(
    "Formato detectado: LARGO\n"
  )
  
  
  vh_counts <- vh_raw %>%
    
    transmute(
      
      sample =
        clean_sample(
          .data[[vh_sample_col]]
        ),
      
      vf_category =
        as.character(
          .data[[category_col]]
        ),
      
      category_count =
        as.numeric(
          .data[[count_col]]
        )
    )
  
  
  # ============================================================
  # 7B. SI ESTÁ EN FORMATO ANCHO
  # ============================================================
  
} else {
  
  
  cat(
    "Formato detectado: ANCHO\n"
  )
  
  
  exclude_cols <- c(
    vh_sample_col,
    "total",
    "Total",
    "TOTAL",
    "VF_total",
    "n_total"
  )
  
  
  possible_category_cols <- setdiff(
    names(vh_raw),
    exclude_cols
  )
  
  
  numeric_category_cols <- possible_category_cols[
    vapply(
      possible_category_cols,
      function(nm) {
        
        is.numeric(
          vh_raw[[nm]]
        )
        
      },
      FUN.VALUE = logical(1)
    )
  ]
  
  
  if (length(numeric_category_cols) == 0) {
    
    stop(
      paste0(
        "No pude identificar las columnas de categorías.\n\n",
        paste(
          names(vh_raw),
          collapse = ", "
        )
      )
    )
  }
  
  
  cat(
    "\nCategorías numéricas detectadas:\n"
  )
  
  print(
    numeric_category_cols
  )
  
  
  vh_counts <- vh_raw %>%
    
    transmute(
      
      sample =
        clean_sample(
          .data[[vh_sample_col]]
        ),
      
      across(
        all_of(
          numeric_category_cols
        )
      )
    ) %>%
    
    pivot_longer(
      
      cols =
        -sample,
      
      names_to =
        "vf_category",
      
      values_to =
        "category_count"
    )
}



# ============================================================
# 8. LIMPIAR
# ============================================================

vh_counts <- vh_counts %>%
  
  mutate(
    
    vf_category =
      stringr::str_trim(
        vf_category
      ),
    
    category_count =
      as.numeric(
        category_count
      )
    
  ) %>%
  
  filter(
    !is.na(sample),
    sample != "",
    !is.na(vf_category),
    vf_category != "",
    !is.na(category_count)
  ) %>%
  
  group_by(
    sample,
    vf_category
  ) %>%
  
  summarise(
    category_count =
      sum(category_count),
    .groups = "drop"
  )



# ============================================================
# 9. CALCULAR COMPOSICIÓN RELATIVA
# ============================================================
#
# VirulentHunter es multietiqueta.
#
# Por eso esto NO significa "% de proteínas del genoma".
#
# Es:
#
# porcentaje relativo de las asignaciones de categorías
# dentro de cada aislamiento.
#
# ============================================================

vh_quant <- vh_counts %>%
  
  group_by(
    sample
  ) %>%
  
  mutate(
    
    total_category_assignments =
      sum(category_count),
    
    category_percent =
      ifelse(
        total_category_assignments > 0,
        100 *
          category_count /
          total_category_assignments,
        NA_real_
      )
  ) %>%
  
  ungroup()



cat(
  "\nMuestras VirulentHunter:",
  n_distinct(vh_quant$sample),
  "\n"
)

cat(
  "Categorías:",
  n_distinct(vh_quant$vf_category),
  "\n"
)



# ============================================================
# ============================================================
#
# PARTE B
# AMRFINDERPLUS / AMRgen
#
# ============================================================
# ============================================================



# ============================================================
# 10. IMPORTAR
# ============================================================

amr <- AMRgen::import_amrfp(
  input_table = AMR_FILE,
  sample_col = "Name"
)



# ============================================================
# 11. LIMPIAR IDs
# ============================================================

amr <- amr %>%
  
  mutate(
    sample =
      clean_sample(id)
  )



# ============================================================
# 12. SEPARAR CLASES AMR
# ============================================================

amr_markers <- amr %>%
  
  transmute(
    
    sample,
    
    marker =
      as.character(marker),
    
    amr_class =
      as.character(drug_class)
  ) %>%
  
  mutate(
    
    amr_class =
      stringr::str_replace_all(
        amr_class,
        "[\\[\\]\"']",
        ""
      )
  ) %>%
  
  tidyr::separate_rows(
    amr_class,
    sep = "\\s*[;,|]+\\s*"
  ) %>%
  
  mutate(
    amr_class =
      stringr::str_trim(
        amr_class
      )
  ) %>%
  
  filter(
    !is.na(sample),
    sample != "",
    !is.na(marker),
    marker != "",
    !is.na(amr_class),
    amr_class != ""
  ) %>%
  
  distinct(
    sample,
    amr_class,
    marker
  )



# ============================================================
# 13. MUESTRAS COMUNES
# ============================================================

common_samples <- intersect(
  unique(vh_quant$sample),
  unique(amr_markers$sample)
)


common_samples <- sort(
  common_samples
)


cat("\n======================================\n")
cat("MUESTRAS COMUNES\n")
cat("======================================\n")

cat(
  "VirulentHunter:",
  n_distinct(vh_quant$sample),
  "\n"
)

cat(
  "AMRgen:",
  n_distinct(amr_markers$sample),
  "\n"
)

cat(
  "Comunes:",
  length(common_samples),
  "\n\n"
)


print(
  common_samples
)



# ============================================================
# 14. LIMITAR A MUESTRAS COMUNES
# ============================================================

vh_quant <- vh_quant %>%
  
  filter(
    sample %in%
      common_samples
  )


amr_markers <- amr_markers %>%
  
  filter(
    sample %in%
      common_samples
  )



# ============================================================
# ============================================================
#
# PARTE C
# CARGA AMR CUANTITATIVA
#
# ============================================================
# ============================================================



# ============================================================
# 15. MARCADORES POR CLASE Y AISLADO
# ============================================================

amr_class_burden <- amr_markers %>%
  
  count(
    sample,
    amr_class,
    name =
      "n_amr_markers"
  )



# ============================================================
# 16. AÑADIR CEROS
# ============================================================

all_amr_classes <- sort(
  unique(
    amr_markers$amr_class
  )
)


amr_class_burden <- tidyr::expand_grid(
  
  sample =
    common_samples,
  
  amr_class =
    all_amr_classes
  
) %>%
  
  left_join(
    amr_class_burden,
    by = c(
      "sample",
      "amr_class"
    )
  ) %>%
  
  mutate(
    
    n_amr_markers =
      replace_na(
        n_amr_markers,
        0L
      ),
    
    amr_present =
      n_amr_markers > 0
  )



# ============================================================
# 17. CARGA AMR TOTAL POR AISLADO
# ============================================================

amr_total_burden <- amr_markers %>%
  
  distinct(
    sample,
    marker
  ) %>%
  
  count(
    sample,
    name =
      "n_amr_total"
  )


amr_total_burden <- tibble(
  sample =
    common_samples
) %>%
  
  left_join(
    amr_total_burden,
    by =
      "sample"
  ) %>%
  
  mutate(
    n_amr_total =
      replace_na(
        n_amr_total,
        0L
      )
  )



# ============================================================
# 18. REVISAR VARIABILIDAD AMR
# ============================================================

amr_class_variability <- amr_class_burden %>%
  
  group_by(
    amr_class
  ) %>%
  
  summarise(
    
    n_present =
      sum(
        amr_present
      ),
    
    n_absent =
      sum(
        !amr_present
      ),
    
    min_markers =
      min(
        n_amr_markers
      ),
    
    max_markers =
      max(
        n_amr_markers
      ),
    
    distinct_marker_counts =
      n_distinct(
        n_amr_markers
      ),
    
    .groups =
      "drop"
  )


cat(
  "\nVariabilidad de clases AMR:\n"
)

print(
  amr_class_variability
)



# ============================================================
# ============================================================
#
# PARTE D
# SPEARMAN:
# VIRULENTHUNTER vs CARGA AMR POR CLASE
#
# ============================================================
# ============================================================



# ============================================================
# 19. FUNCIÓN SEGURA PARA SPEARMAN
# ============================================================

safe_spearman <- function(x, y) {
  
  ok <- complete.cases(
    x,
    y
  )
  
  x <- x[ok]
  y <- y[ok]
  
  n <- length(x)
  
  
  if (
    n < 5 ||
    dplyr::n_distinct(x) < 2 ||
    dplyr::n_distinct(y) < 2
  ) {
    
    return(
      tibble(
        n = n,
        rho = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  
  result <- suppressWarnings(
    
    cor.test(
      x,
      y,
      method =
        "spearman",
      exact =
        FALSE
    )
  )
  
  
  tibble(
    
    n =
      n,
    
    rho =
      unname(
        result$estimate
      ),
    
    p_value =
      result$p.value
  )
}



# ============================================================
# 20. UNIR DATOS CUANTITATIVOS
# ============================================================

spearman_input <- vh_quant %>%
  
  select(
    sample,
    vf_category,
    category_count,
    category_percent
  ) %>%
  
  inner_join(
    amr_class_burden,
    by =
      "sample",
    relationship =
      "many-to-many"
  )



# ============================================================
# 21. CORRELACIÓN UTILIZANDO PROPORCIÓN DE VIRULENCIA
# ============================================================

vh_amr_spearman <- spearman_input %>%
  
  group_by(
    vf_category,
    amr_class
  ) %>%
  
  group_modify(
    
    ~ safe_spearman(
      .x$category_percent,
      .x$n_amr_markers
    )
    
  ) %>%
  
  ungroup() %>%
  
  mutate(
    
    p_adjusted =
      p.adjust(
        p_value,
        method =
          "BH"
      ),
    
    abs_rho =
      abs(
        rho
      )
  ) %>%
  
  arrange(
    p_adjusted,
    desc(abs_rho)
  )



# ============================================================
# 22. GUARDAR
# ============================================================

write.csv(
  vh_amr_spearman,
  "VirulentHunter_AMR_Spearman_quantitative.csv",
  row.names =
    FALSE
)



# ============================================================
# 23. MOSTRAR TOP
# ============================================================

cat(
  "\n======================================\n"
)

cat(
  "TOP CORRELACIONES VH-AMR\n"
)

cat(
  "======================================\n"
)


print(
  
  vh_amr_spearman %>%
    
    filter(
      !is.na(rho)
    ) %>%
    
    arrange(
      desc(abs_rho)
    ) %>%
    
    head(25)
  
)



# ============================================================
# ============================================================
#
# PARTE E
# VIRULENCIA vs CARGA AMR TOTAL
#
# ============================================================
# ============================================================



# ============================================================
# 24. CORRELACIÓN CON TOTAL DE MARCADORES AMR
# ============================================================

vh_total_amr_input <- vh_quant %>%
  
  select(
    sample,
    vf_category,
    category_percent
  ) %>%
  
  left_join(
    amr_total_burden,
    by =
      "sample"
  )



vh_total_amr <- vh_total_amr_input %>%
  
  group_by(
    vf_category
  ) %>%
  
  group_modify(
    
    ~ safe_spearman(
      .x$category_percent,
      .x$n_amr_total
    )
    
  ) %>%
  
  ungroup() %>%
  
  mutate(
    
    p_adjusted =
      p.adjust(
        p_value,
        method =
          "BH"
      ),
    
    abs_rho =
      abs(
        rho
      )
  ) %>%
  
  arrange(
    p_adjusted,
    desc(abs_rho)
  )


cat(
  "\nVirulentHunter vs carga AMR total:\n"
)


print(
  vh_total_amr
)


write.csv(
  vh_total_amr,
  "VirulentHunter_total_AMR_burden_Spearman.csv",
  row.names =
    FALSE
)



# ============================================================
# ============================================================
#
# PARTE F
# WILCOXON PARA CLASES AMR VARIABLES
#
# ============================================================
# ============================================================



# ============================================================
# 25. REQUERIR AL MENOS 3 + Y 3 -
# ============================================================
#
# Con tus datos:
#
# Rifamicinas probablemente será la más apropiada.
#
# Quinolones:
# 11 positivos / 2 negativos
#
# Sulfonamidas:
# 12 positivos / 1 negativo
#
# Estas últimas tienen poco poder como comparación binaria.
#
# ============================================================

MIN_GROUP_SIZE <- 3


eligible_amr_classes <- amr_class_variability %>%
  
  filter(
    n_present >=
      MIN_GROUP_SIZE,
    n_absent >=
      MIN_GROUP_SIZE
  )


cat(
  "\nClases elegibles para Wilcoxon:\n"
)


print(
  eligible_amr_classes
)



# ============================================================
# 26. FUNCIÓN WILCOXON + RANK-BISERIAL
# ============================================================

safe_wilcox <- function(x, present) {
  
  x1 <- x[
    present
  ]
  
  x0 <- x[
    !present
  ]
  
  
  n1 <- length(x1)
  n0 <- length(x0)
  
  
  if (
    n1 < MIN_GROUP_SIZE ||
    n0 < MIN_GROUP_SIZE
  ) {
    
    return(
      tibble(
        n_positive = n1,
        n_negative = n0,
        median_positive = NA_real_,
        median_negative = NA_real_,
        median_difference = NA_real_,
        rank_biserial = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  
  test <- suppressWarnings(
    
    wilcox.test(
      x1,
      x0,
      exact =
        FALSE
    )
  )
  
  
  all_values <- c(
    x1,
    x0
  )
  
  
  ranks <- rank(
    all_values,
    ties.method =
      "average"
  )
  
  
  U <- sum(
    ranks[
      seq_along(
        x1
      )
    ]
  ) -
    n1 *
    (n1 + 1) /
    2
  
  
  rank_biserial <-
    2 * U /
    (n1 * n0) -
    1
  
  
  tibble(
    
    n_positive =
      n1,
    
    n_negative =
      n0,
    
    median_positive =
      median(
        x1
      ),
    
    median_negative =
      median(
        x0
      ),
    
    median_difference =
      median(x1) -
      median(x0),
    
    rank_biserial =
      rank_biserial,
    
    p_value =
      test$p.value
  )
}



# ============================================================
# 27. PREPARAR DATOS
# ============================================================

wilcox_input <- vh_quant %>%
  
  select(
    sample,
    vf_category,
    category_percent
  ) %>%
  
  inner_join(
    amr_class_burden,
    by =
      "sample",
    relationship =
      "many-to-many"
  ) %>%
  
  filter(
    amr_class %in%
      eligible_amr_classes$amr_class
  )



# ============================================================
# 28. WILCOXON
# ============================================================

if (
  nrow(
    eligible_amr_classes
  ) > 0
) {
  
  
  vh_amr_wilcox <- wilcox_input %>%
    
    group_by(
      vf_category,
      amr_class
    ) %>%
    
    group_modify(
      
      ~ safe_wilcox(
        .x$category_percent,
        .x$amr_present
      )
      
    ) %>%
    
    ungroup() %>%
    
    mutate(
      
      p_adjusted =
        p.adjust(
          p_value,
          method =
            "BH"
        ),
      
      abs_effect =
        abs(
          rank_biserial
        )
    ) %>%
    
    arrange(
      p_adjusted,
      desc(abs_effect)
    )
  
  
  print(
    vh_amr_wilcox
  )
  
  
  write.csv(
    vh_amr_wilcox,
    "VirulentHunter_AMR_Wilcoxon.csv",
    row.names =
      FALSE
  )
  
  
} else {
  
  
  cat(
    "\nNinguna clase AMR tiene suficientes + y - para Wilcoxon.\n"
  )
  
  
  vh_amr_wilcox <- tibble()
  
}



# ============================================================
# ============================================================
#
# PARTE G
# HEATMAP DE SPEARMAN
#
# ============================================================
# ============================================================



# ============================================================
# 29. HEATMAP
# ============================================================

heatmap_data <- vh_amr_spearman %>%
  
  filter(
    !is.na(rho)
  )



p_heatmap <- ggplot(
  
  heatmap_data,
  
  aes(
    x =
      amr_class,
    y =
      vf_category,
    fill =
      rho
  )
  
) +
  
  geom_tile(
    colour =
      "white",
    linewidth =
      0.3
  ) +
  
  geom_text(
    
    aes(
      label =
        sprintf(
          "%.2f",
          rho
        )
    ),
    
    size =
      3
  ) +
  
  scale_fill_gradient2(
    midpoint =
      0,
    limits =
      c(-1, 1)
  ) +
  
  labs(
    
    title =
      "VirulentHunter categories vs AMR determinant burden",
    
    subtitle =
      paste0(
        "Spearman correlation; n = ",
        length(common_samples),
        " isolates"
      ),
    
    x =
      "AMR class",
    
    y =
      "VirulentHunter category",
    
    fill =
      "Spearman rho",
    
    caption =
      paste0(
        "VirulentHunter values represent relative category composition. ",
        "AMR values represent number of unique AMRFinderPlus markers per class."
      )
  ) +
  
  theme_minimal(
    base_size =
      11
  ) +
  
  theme(
    
    panel.grid =
      element_blank(),
    
    axis.text.x =
      element_text(
        angle =
          45,
        hjust =
          1
      )
  )


print(
  p_heatmap
)


ggsave(
  "VirulentHunter_AMR_Spearman_heatmap.pdf",
  p_heatmap,
  width =
    11,
  height =
    8
)


ggsave(
  "VirulentHunter_AMR_Spearman_heatmap.png",
  p_heatmap,
  width =
    11,
  height =
    8,
  dpi =
    600
)



# ============================================================
# ============================================================
#
# PARTE H
# SANKEY DE ASOCIACIONES CUANTITATIVAS
#
# ============================================================
# ============================================================



# ============================================================
# 30. ELEGIR RELACIONES PARA MOSTRAR
# ============================================================
#
# Primero intentamos:
#
# FDR < 0.05
#
# Si no existen, se muestran las 15 correlaciones de mayor
# magnitud como figura EXPLORATORIA.
#
# ============================================================

significant_quant <- vh_amr_spearman %>%
  
  filter(
    !is.na(rho),
    p_adjusted <
      0.05
  )


TOP_EXPLORATORY <- 15


if (
  nrow(
    significant_quant
  ) >= 2
) {
  
  
  sankey_quant <- significant_quant
  
  
  SANKEY_SUBTITLE <-
    "FDR-significant quantitative associations"
  
  
} else {
  
  
  warning(
    paste(
      "No hay suficientes asociaciones que sobrevivan FDR.",
      "El Sankey cuantitativo será EXPLORATORIO."
    )
  )
  
  
  sankey_quant <- vh_amr_spearman %>%
    
    filter(
      !is.na(rho)
    ) %>%
    
    arrange(
      desc(
        abs(rho)
      )
    ) %>%
    
    slice_head(
      n =
        TOP_EXPLORATORY
    )
  
  
  SANKEY_SUBTITLE <-
    "Top exploratory quantitative associations"
  
}



# ============================================================
# 31. PREPARAR
# ============================================================

sankey_quant <- sankey_quant %>%
  
  mutate(
    
    flow_weight =
      abs(
        rho
      ),
    
    direction =
      case_when(
        
        rho > 0 ~
          "Positive",
        
        rho < 0 ~
          "Negative",
        
        TRUE ~
          "No association"
      )
  )



# ============================================================
# 32. SANKEY
# ============================================================

p_sankey_quant <- ggplot(
  
  sankey_quant,
  
  aes(
    
    y =
      flow_weight,
    
    axis1 =
      vf_category,
    
    axis2 =
      amr_class
  )
  
) +
  
  geom_alluvium(
    
    aes(
      fill =
        direction
    ),
    
    width =
      0.14,
    
    alpha =
      0.75
  ) +
  
  geom_stratum(
    
    width =
      0.18,
    
    colour =
      "grey30",
    
    linewidth =
      0.35
  ) +
  
  geom_text(
    
    stat =
      "stratum",
    
    aes(
      label =
        after_stat(
          stratum
        )
    ),
    
    size =
      3
  ) +
  
  scale_x_discrete(
    
    limits =
      c(
        "VirulentHunter",
        "AMR classes"
      ),
    
    expand =
      c(
        0.12,
        0.12
      )
  ) +
  
  labs(
    
    title =
      "Virulence–AMR quantitative associations",
    
    subtitle =
      SANKEY_SUBTITLE,
    
    x =
      NULL,
    
    y =
      "|Spearman rho|",
    
    fill =
      "Direction",
    
    caption =
      paste0(
        "Link width represents the absolute Spearman correlation. ",
        "Exploratory links should not be interpreted as statistically ",
        "significant unless they survive FDR correction."
      )
  ) +
  
  theme_minimal(
    base_size =
      11
  ) +
  
  theme(
    
    panel.grid =
      element_blank(),
    
    axis.text.x =
      element_text(
        face =
          "bold"
      ),
    
    axis.text.y =
      element_blank(),
    
    axis.ticks.y =
      element_blank()
  )


print(
  p_sankey_quant
)


ggsave(
  "Sankey_VirulentHunter_AMR_quantitative.pdf",
  p_sankey_quant,
  width =
    11,
  height =
    8
)


ggsave(
  "Sankey_VirulentHunter_AMR_quantitative.png",
  p_sankey_quant,
  width =
    11,
  height =
    8,
  dpi =
    600
)



# ============================================================
# 33. EXPORTAR DATOS DEL SANKEY
# ============================================================

write.csv(
  sankey_quant,
  "Sankey_VirulentHunter_AMR_quantitative_data.csv",
  row.names =
    FALSE
)



# ============================================================
# ============================================================
#
# PARTE I
# RESUMEN FINAL
#
# ============================================================
# ============================================================


cat("\n")
cat("======================================\n")
cat("RESUMEN\n")
cat("======================================\n")

cat(
  "Aislamientos comunes:",
  length(common_samples),
  "\n"
)

cat(
  "Categorías VirulentHunter:",
  n_distinct(
    vh_quant$vf_category
  ),
  "\n"
)

cat(
  "Clases AMR:",
  length(
    all_amr_classes
  ),
  "\n"
)

cat(
  "Correlaciones evaluadas:",
  nrow(
    vh_amr_spearman
  ),
  "\n"
)

cat(
  "Correlaciones FDR < 0.05:",
  sum(
    vh_amr_spearman$p_adjusted < 0.05,
    na.rm =
      TRUE
  ),
  "\n"
)

cat(
  "Clases elegibles para Wilcoxon:",
  nrow(
    eligible_amr_classes
  ),
  "\n"
)

cat("\nArchivos principales:\n")

cat(
  "  VirulentHunter_AMR_Spearman_quantitative.csv\n"
)

cat(
  "  VirulentHunter_total_AMR_burden_Spearman.csv\n"
)

cat(
  "  VirulentHunter_AMR_Spearman_heatmap.pdf\n"
)

cat(
  "  Sankey_VirulentHunter_AMR_quantitative.pdf\n"
)

cat("======================================\n")

####################################################3
#################################################
########################################
# ============================================================
# HEATMAP MODIFICADO:
# rho + significancia FDR
# ============================================================

# Requiere que ya exista:
# vh_amr_spearman
# con columnas:
#   vf_category
#   amr_class
#   rho
#   p_value
#   p_adjusted


# ============================================================
# 1. PREPARAR DATOS
# ============================================================

heatmap_data <- vh_amr_spearman %>%
  filter(!is.na(rho)) %>%
  mutate(
    
    # estrellas por FDR
    sig_label = case_when(
      !is.na(p_adjusted) & p_adjusted < 0.001 ~ "***",
      !is.na(p_adjusted) & p_adjusted < 0.01  ~ "**",
      !is.na(p_adjusted) & p_adjusted < 0.05  ~ "*",
      TRUE ~ ""
    ),
    
    significant = !is.na(p_adjusted) & p_adjusted < 0.05,
    
    # etiqueta dentro de la celda
    plot_label = ifelse(
      sig_label == "",
      sprintf("%.2f", rho),
      paste0(sprintf("%.2f", rho), "\n", sig_label)
    )
  )


# ============================================================
# 2. ORDEN DE FACTORES
# ============================================================
# Puedes cambiar estos vectores si quieres otro orden

vf_order <- c(
  "Stress survival",
  "Regulation",
  "Post-translational modification",
  "Others",
  "Nutritional/Metabolic factor",
  "Motility",
  "Invasion",
  "Immune modulation",
  "Exotoxin",
  "Exoenzyme",
  "Effector delivery system",
  "Biofilm",
  "Antimicrobial activity/Competitive advantage",
  "Adherence"
)

amr_order <- c(
  "Beta-lactams",
  "Cephalosporins (3rd gen.)",
  "Macrólidos",
  "Quinolones",
  "Rifamicinas",
  "Sulfonamidas",
  "Tetracyclines"
)

heatmap_data <- heatmap_data %>%
  mutate(
    vf_category = factor(vf_category, levels = vf_order),
    amr_class   = factor(amr_class, levels = amr_order)
  )


# ============================================================
# 3. CONTEO DE CELDAS SIGNIFICATIVAS
# ============================================================

n_sig <- heatmap_data %>%
  filter(significant) %>%
  nrow()

cat(
  "\nNúmero de asociaciones con FDR < 0.05:",
  n_sig,
  "\n"
)


# ============================================================
# 4. HEATMAP PRINCIPAL
# ============================================================

p_heatmap_sig <- ggplot(
  heatmap_data,
  aes(
    x = amr_class,
    y = vf_category,
    fill = rho
  )
) +
  
  # mosaico base
  geom_tile(
    colour = "white",
    linewidth = 0.4
  ) +
  
  # remarcar celdas significativas
  geom_tile(
    data = heatmap_data %>% filter(significant),
    fill = NA,
    colour = "black",
    linewidth = 1.0
  ) +
  
  # texto dentro de celdas
  geom_text(
    aes(label = plot_label),
    size = 3.2,
    lineheight = 0.9
  ) +
  
  scale_fill_gradient2(
    low = "#b3584a",
    mid = "white",
    high = "#5e4fa2",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman rho"
  ) +
  
  labs(
    title = "VirulentHunter categories vs AMR determinant burden",
    subtitle = paste0(
      "Spearman correlation; n = ",
      length(common_samples),
      " isolates"
    ),
    x = "AMR class",
    y = "VirulentHunter category",
    caption = paste0(
      "Cell values are Spearman rho. ",
      "Black borders indicate associations with Benjamini-Hochberg FDR < 0.05. ",
      "Asterisks: * FDR < 0.05, ** FDR < 0.01, *** FDR < 0.001. ",
      "VirulentHunter values represent relative category composition. ",
      "AMR values represent number of unique AMRFinderPlus markers per class."
    )
  ) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "right"
  )


print(p_heatmap_sig)


# ============================================================
# 5. EXPORTAR
# ============================================================

ggsave(
  "VirulentHunter_AMR_Spearman_heatmap_FDR.pdf",
  p_heatmap_sig,
  width = 12,
  height = 8.5
)

ggsave(
  "VirulentHunter_AMR_Spearman_heatmap_FDR.png",
  p_heatmap_sig,
  width = 12,
  height = 8.5,
  dpi = 600
)


# ============================================================
# 6. TABLA DE APOYO SOLO CON LAS SIGNIFICATIVAS
# ============================================================

heatmap_significant_table <- heatmap_data %>%
  filter(significant) %>%
  arrange(p_adjusted, desc(abs(rho))) %>%
  select(
    vf_category,
    amr_class,
    rho,
    p_value,
    p_adjusted,
    sig_label
  )

print(heatmap_significant_table, n = Inf)

write.csv(
  heatmap_significant_table,
  "VirulentHunter_AMR_Spearman_significant_only.csv",
  row.names = FALSE
)

####################################################################
####################################################################
##############################################################3

# ============================================================
# CHORD DIAGRAM
# VirulentHunter categories vs AMR determinant burden
#
# Requiere que ya exista:
#
#   vh_amr_spearman
#
# con las columnas:
#
#   vf_category
#   amr_class
#   rho
#   p_value
#   p_adjusted
#
# ============================================================


# ============================================================
# 0. PAQUETES
# ============================================================

if (!requireNamespace("circlize", quietly = TRUE)) {
  install.packages("circlize")
}

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

library(circlize)
library(dplyr)


# ============================================================
# 1. PARÁMETROS
# ============================================================

USE_ONLY_FDR_SIGNIFICANT <- TRUE

FDR_THRESHOLD <- 0.05

MIN_ABS_RHO <- 0.40

TOP_EXPLORATORY <- 20


# Colores de enlaces

COL_POS <- "#2CBCC3"
COL_NEG <- "#F28E85"


# Colores de sectores

COL_VF <- "#D9D9D9"
COL_AMR <- "#BDBDBD"


# ============================================================
# 2. PREPARAR TABLA
# ============================================================

chord_df <- vh_amr_spearman %>%
  
  filter(
    !is.na(rho)
  ) %>%
  
  mutate(
    
    abs_rho =
      abs(rho),
    
    direction =
      if_else(
        rho >= 0,
        "Positive",
        "Negative"
      )
  )


# ============================================================
# 3. SELECCIONAR ASOCIACIONES
# ============================================================

if (USE_ONLY_FDR_SIGNIFICANT) {
  
  chord_use <- chord_df %>%
    
    filter(
      !is.na(p_adjusted),
      p_adjusted < FDR_THRESHOLD,
      abs_rho >= MIN_ABS_RHO
    )
  
  
  if (nrow(chord_use) >= 2) {
    
    chord_subtitle <- paste0(
      "FDR-significant associations (BH-adjusted p < ",
      FDR_THRESHOLD,
      ")"
    )
    
  } else {
    
    warning(
      paste(
        "Hay menos de dos asociaciones FDR-significativas.",
        "Se mostrarán asociaciones exploratorias."
      )
    )
    
    
    chord_use <- chord_df %>%
      
      filter(
        abs_rho >= MIN_ABS_RHO
      ) %>%
      
      arrange(
        desc(abs_rho)
      ) %>%
      
      slice_head(
        n = TOP_EXPLORATORY
      )
    
    
    chord_subtitle <-
      "Top exploratory quantitative associations"
  }
  
} else {
  
  chord_use <- chord_df %>%
    
    filter(
      abs_rho >= MIN_ABS_RHO
    ) %>%
    
    arrange(
      desc(abs_rho)
    ) %>%
    
    slice_head(
      n = TOP_EXPLORATORY
    )
  
  
  chord_subtitle <-
    "Top exploratory quantitative associations"
}


if (nrow(chord_use) == 0) {
  
  stop(
    "No existen asociaciones que cumplan los criterios para generar el chord."
  )
}


# ============================================================
# 4. ORDEN DE LAS CATEGORÍAS
# ============================================================

vf_order <- c(
  
  "Adherence",
  
  "Antimicrobial activity/Competitive advantage",
  
  "Biofilm",
  
  "Effector delivery system",
  
  "Exoenzyme",
  
  "Exotoxin",
  
  "Immune modulation",
  
  "Invasion",
  
  "Motility",
  
  "Nutritional/Metabolic factor",
  
  "Others",
  
  "Post-translational modification",
  
  "Regulation",
  
  "Stress survival"
)


amr_order <- c(
  
  "Beta-lactams",
  
  "Cephalosporins (3rd gen.)",
  
  "Macrólidos",
  
  "Quinolones",
  
  "Rifamicinas",
  
  "Sulfonamidas",
  
  "Tetracyclines"
)


# ============================================================
# 5. CONSERVAR SÓLO LAS CATEGORÍAS PRESENTES
# ============================================================

vf_present <- vf_order[
  vf_order %in%
    unique(chord_use$vf_category)
]


amr_present <- amr_order[
  amr_order %in%
    unique(chord_use$amr_class)
]


sector_order <- c(
  vf_present,
  amr_present
)


# ============================================================
# 6. TABLA PARA circlize
# ============================================================

plot_df <- chord_use %>%
  
  transmute(
    
    from =
      vf_category,
    
    to =
      amr_class,
    
    value =
      abs_rho,
    
    rho =
      rho,
    
    direction =
      direction,
    
    p_adjusted =
      p_adjusted
  )


# ============================================================
# 7. COLORES DE ENLACES
# ============================================================

link_cols <- ifelse(
  
  plot_df$direction == "Positive",
  
  adjustcolor(
    COL_POS,
    alpha.f = 0.80
  ),
  
  adjustcolor(
    COL_NEG,
    alpha.f = 0.80
  )
)


# ============================================================
# 8. COLORES DE SECTORES
# ============================================================

grid_cols <- c(
  
  setNames(
    rep(
      COL_VF,
      length(vf_present)
    ),
    vf_present
  ),
  
  setNames(
    rep(
      COL_AMR,
      length(amr_present)
    ),
    amr_present
  )
)


# ============================================================
# 9. ESPACIOS ENTRE SECTORES
# ============================================================
#
# Se deja un espacio grande entre:
#
# VirulentHunter | AMR
#
# ============================================================

gap_vector <- c(
  
  if (length(vf_present) > 1) {
    rep(
      2,
      length(vf_present) - 1
    )
  },
  
  12,
  
  if (length(amr_present) > 1) {
    rep(
      2,
      length(amr_present) - 1
    )
  },
  
  12
)


# ============================================================
# 10. FUNCIÓN PARA DIBUJAR EL CHORD
# ============================================================

draw_chord <- function() {
  
  # Limpiar diagramas circulares anteriores
  
  circos.clear()
  
  
  # Márgenes del dispositivo gráfico
  
  par(
    mar = c(
      4,
      4,
      5,
      4
    )
  )
  
  
  # Configuración de circlize
  
  circos.par(
    
    start.degree = 90,
    
    gap.degree = gap_vector,
    
    track.margin = c(
      0.01,
      0.01
    ),
    
    canvas.xlim = c(
      -1.35,
      1.35
    ),
    
    canvas.ylim = c(
      -1.35,
      1.35
    )
  )
  
  
  # ==========================================================
  # CHORD
  # ==========================================================
  
  chordDiagram(
    
    x = plot_df[
      ,
      c(
        "from",
        "to",
        "value"
      )
    ],
    
    order =
      sector_order,
    
    grid.col =
      grid_cols,
    
    col =
      link_cols,
    
    transparency =
      0.10,
    
    annotationTrack =
      "grid",
    
    preAllocateTracks =
      list(
        track.height = 0.14
      ),
    
    directional =
      0,
    
    link.sort =
      TRUE,
    
    link.decreasing =
      FALSE
  )
  
  
  # ==========================================================
  # ETIQUETAS EXTERNAS
  # ==========================================================
  
  circos.trackPlotRegion(
    
    track.index = 1,
    
    panel.fun = function(x, y) {
      
      sector_name <- get.cell.meta.data(
        "sector.index"
      )
      
      
      xcenter <- get.cell.meta.data(
        "xcenter"
      )
      
      
      ylim <- get.cell.meta.data(
        "ylim"
      )
      
      
      circos.text(
        
        x =
          xcenter,
        
        y =
          ylim[1] + 0.1,
        
        labels =
          sector_name,
        
        facing =
          "clockwise",
        
        niceFacing =
          TRUE,
        
        adj =
          c(
            0,
            0.5
          ),
        
        cex =
          0.70
      )
    },
    
    bg.border =
      NA
  )
  
  
  # ==========================================================
  # TÍTULO
  # ==========================================================
  
  title(
    
    main =
      "VirulentHunter–AMR associations",
    
    sub =
      chord_subtitle,
    
    cex.main =
      1.35,
    
    cex.sub =
      0.95,
    
    line =
      1
  )
  
  
  # ==========================================================
  # LEYENDA
  # ==========================================================
  
  legend(
    
    "topleft",
    
    legend = c(
      "Positive association",
      "Negative association"
    ),
    
    fill = c(
      COL_POS,
      COL_NEG
    ),
    
    border =
      NA,
    
    bty =
      "n",
    
    cex =
      0.80
  )
  
  
  # ==========================================================
  # TEXTO INFERIOR
  # ==========================================================
  
  mtext(
    
    paste0(
      "Link width represents |Spearman rho|. ",
      "Only associations passing the selected FDR and effect-size ",
      "criteria are shown."
    ),
    
    side =
      1,
    
    line =
      2,
    
    cex =
      0.70
  )
}


# ============================================================
# 11. MOSTRAR AHORA EN RSTUDIO
# ============================================================
#
# IMPORTANTE:
#
# Esta llamada NO está dentro de pdf() ni png().
# Por eso debe aparecer directamente en la pestaña Plots.
#
# ============================================================

draw_chord()


# ============================================================
# 12. GUARDAR PDF
# ============================================================

pdf(
  
  file =
    "Chord_VirulentHunter_AMR_significant.pdf",
  
  width =
    12,
  
  height =
    12
)


draw_chord()


dev.off()


# ============================================================
# 13. GUARDAR PNG
# ============================================================

png(
  
  filename =
    "Chord_VirulentHunter_AMR_significant.png",
  
  width =
    3600,
  
  height =
    3600,
  
  res =
    300
)


draw_chord()


dev.off()


# ============================================================
# 14. VOLVER A MOSTRAR EN RSTUDIO
# ============================================================
#
# dev.off() cerró el dispositivo PDF/PNG.
# Esta llamada vuelve a dibujar el chord en Plots.
#
# ============================================================

draw_chord()


# ============================================================
# 15. EXPORTAR DATOS UTILIZADOS
# ============================================================

write.csv(
  
  plot_df,
  
  "Chord_VirulentHunter_AMR_significant_data.csv",
  
  row.names =
    FALSE
)


# ============================================================
# 16. MOSTRAR ASOCIACIONES INCLUIDAS
# ============================================================

cat("\n")
cat("========================================\n")
cat("ASOCIACIONES INCLUIDAS EN EL CHORD\n")
cat("========================================\n")


print(
  plot_df,
  n = Inf
)


# ============================================================
# 17. CONFIRMAR ARCHIVOS
# ============================================================

cat("\n")
cat("========================================\n")
cat("ARCHIVOS GENERADOS\n")
cat("========================================\n")


cat(
  
  "PDF: ",
  
  file.exists(
    "Chord_VirulentHunter_AMR_significant.pdf"
  ),
  
  "\n"
)


cat(
  
  "PNG: ",
  
  file.exists(
    "Chord_VirulentHunter_AMR_significant.png"
  ),
  
  "\n"
)


cat(
  
  "CSV: ",
  
  file.exists(
    "Chord_VirulentHunter_AMR_significant_data.csv"
  ),
  
  "\n"
)


cat(
  "\nDirectorio:\n",
  getwd(),
  "\n"
)

