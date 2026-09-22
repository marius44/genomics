#!/usr/bin/env bash

set -u
set -o pipefail

source /opt/conda/etc/profile.d/conda.sh
conda activate amrfinder

DB="/data/databases/amrfinder/2026-08-07.1"
OUTDIR="amrfinder_results"

mkdir -p "$OUTDIR"

shopt -s nullglob
GENOMES=(*.fna)

if (( ${#GENOMES[@]} == 0 )); then
    echo "ERROR: no se encontraron archivos .fna en:"
    pwd
    conda deactivate
    exit 1
fi

echo "Genomas encontrados: ${#GENOMES[@]}"
echo

PROCESSED=0
FAILED=0

for IN in "${GENOMES[@]}"; do
    SAMPLE="${IN%.fna}"

    echo "=========================================="
    echo "Procesando: $SAMPLE"
    echo "Archivo: $IN"
    echo "=========================================="

    if amrfinder \
        --threads 40 \
        --print_node \
        -O Neisseria_gonorrhoeae \
        --name "$SAMPLE" \
        -n "$IN" \
        --database "$DB" \
        -o "$OUTDIR/${SAMPLE}_amrfp.tsv" \
        2> "$OUTDIR/${SAMPLE}_amrfp.log"
    then
        echo "Terminado: $SAMPLE"
        ((PROCESSED+=1))
    else
        echo "ERROR al procesar: $SAMPLE"
        ((FAILED+=1))
    fi

    echo
done

echo "=========================================="
echo "Ejecución terminada"
echo "Genomas encontrados: ${#GENOMES[@]}"
echo "Procesados correctamente: $PROCESSED"
echo "Fallidos: $FAILED"
echo "Resultados: $OUTDIR"
echo "=========================================="

# ============================================================
# CONSOLIDAR RESULTADOS PARA AMRgen
# ============================================================

COMBINED="$OUTDIR/amrfinder_Neisseria_gonorrhoeae.tsv"
SAMPLE_LIST="$OUTDIR/amrfinder_samples.txt"

# Lista completa de genomas analizados, incluidos los que no tuvieron hits
printf '%s\n' "${GENOMES[@]%.fna}" > "$SAMPLE_LIST"

# Combinar resultados conservando una sola cabecera
awk '
    FNR == 1 && NR != 1 {next}
    {print}
' "$OUTDIR"/*_amrfp.tsv > "$COMBINED"

echo "Tabla consolidada: $COMBINED"
echo "Lista de muestras: $SAMPLE_LIST"
echo "Muestras incluidas: $(wc -l < "$SAMPLE_LIST")"
echo "Hits consolidados: $(( $(wc -l < "$COMBINED") - 1 ))"

conda deactivate

if (( FAILED > 0 )); then
    exit 1
fi
