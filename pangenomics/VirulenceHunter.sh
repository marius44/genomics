#!/usr/bin/env bash
# Se corre en una Mac o máquina con GPUs potentes
set -u

# ============================================================
# VirulentHunter batch runner
# Procesa secuencialmente todos los proteomas .faa
# ============================================================

# Directorio del repositorio VirulentHunter
VH_DIR="/Users/UAML/VirulentHunter"

# Directorio donde están los proteomas
INPUT_DIR="/Users/UAML/VirulentHunter"

# Directorio donde se guardarán los resultados
RESULTS_DIR="${INPUT_DIR}/virulenthunter_results"

# Patrón de archivos a procesar
PATTERN="*.faa"

# Modelo ESM2
ESM_MODEL="facebook/esm2_t30_150M_UR50D"

# Entorno conda
CONDA_ENV="virulenthunter-mps"


# ------------------------------------------------------------
# Preparación
# ------------------------------------------------------------

mkdir -p "${RESULTS_DIR}"

# Permite fallback a CPU si alguna operación no está soportada
# por MPS.
export PYTORCH_ENABLE_MPS_FALLBACK=1

# Cargar conda
source /opt/miniconda3/etc/profile.d/conda.sh
conda activate "${CONDA_ENV}"

cd "${VH_DIR}" || {
    echo "ERROR: No pude entrar a ${VH_DIR}"
    exit 1
}


echo "============================================================"
echo " VirulentHunter batch"
echo "============================================================"
echo "Input:      ${INPUT_DIR}"
echo "Resultados: ${RESULTS_DIR}"
echo "Modelo:     ${ESM_MODEL}"
echo "Inicio:     $(date)"
echo "============================================================"
echo


# ------------------------------------------------------------
# Procesamiento
# ------------------------------------------------------------

shopt -s nullglob

FILES=("${INPUT_DIR}"/${PATTERN})

if [ "${#FILES[@]}" -eq 0 ]; then
    echo "ERROR: No encontré archivos ${PATTERN} en:"
    echo "${INPUT_DIR}"
    exit 1
fi

TOTAL=${#FILES[@]}
COUNT=0

for FAA in "${FILES[@]}"; do

    COUNT=$((COUNT + 1))

    FILE=$(basename "${FAA}")

    # Quitar extensiones habituales:
    # NG41.fna.faa -> NG41
    SAMPLE="${FILE%.faa}"
    SAMPLE="${SAMPLE%.fna}"

    OUTDIR="${RESULTS_DIR}/${SAMPLE}"
    LOGFILE="${OUTDIR}/${SAMPLE}_virulenthunter.log"
    RESULTFILE="${OUTDIR}/predict_results.csv"

    echo
    echo "============================================================"
    echo "[${COUNT}/${TOTAL}] ${SAMPLE}"
    echo "Archivo: ${FAA}"
    echo "Inicio:  $(date)"
    echo "============================================================"

    mkdir -p "${OUTDIR}"

    # Saltar muestras terminadas
    if [ -s "${RESULTFILE}" ]; then
        echo "Resultado existente:"
        echo "${RESULTFILE}"
        echo "SKIP ${SAMPLE}"
        continue
    fi

    START=$(date +%s)

    python predict.py \
        -i "${FAA}" \
        -o "${OUTDIR}/" \
        --base_esm_path "${ESM_MODEL}" \
        2>&1 | tee "${LOGFILE}"

    EXITCODE=${PIPESTATUS[0]}

    END=$(date +%s)
    ELAPSED=$((END - START))

    HOURS=$((ELAPSED / 3600))
    MINUTES=$(((ELAPSED % 3600) / 60))
    SECONDS=$((ELAPSED % 60))

    if [ "${EXITCODE}" -eq 0 ] && [ -s "${RESULTFILE}" ]; then
        echo
        echo "OK: ${SAMPLE}"
        echo "Tiempo: ${HOURS}h ${MINUTES}m ${SECONDS}s"
        echo "Resultado: ${RESULTFILE}"
    else
        echo
        echo "ERROR procesando ${SAMPLE}"
        echo "Código de salida: ${EXITCODE}"
        echo "Log: ${LOGFILE}"
        echo
        echo "El batch continuará con la siguiente muestra."
    fi

done


echo
echo "============================================================"
echo " Batch terminado"
echo " Fecha: $(date)"
echo " Resultados: ${RESULTS_DIR}"
echo "============================================================"
