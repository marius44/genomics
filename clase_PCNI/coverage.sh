#!/bin/bash

#coverage_ont_batch.sh
#
# Calcula automáticamente coverage para genomas ensamblados
# usando lecturas Oxford Nanopore (ONT).
#
# Requisitos:
#   minimap2
#   samtools
#
# Entrada:
#   FASTA y FASTQ en el directorio actual.
#
# Ejemplos compatibles:
#
#   005A-004-JEJ.fasta
#   005A-004-JEJ.fastq.gz
#
#   005A-007-AJL.fasta
#   005A-007-AJL_001.fastq.gz
#   005A-007-AJL_002.fastq.gz
#
# Salida principal:
#   coverage_results/coverage_NCBI.csv
#
# ============================================================

set -uo pipefail

# ------------------------------------------------------------
# CONFIGURACION
# ------------------------------------------------------------

DIR="."

OUTDIR="coverage_results"
BAMDIR="${OUTDIR}/bam"
LOGDIR="${OUTDIR}/logs"
TMPDIR="${OUTDIR}/tmp"

CSV="${OUTDIR}/coverage_NCBI.csv"

# Porcentaje de MemAvailable que se permitira usar para sort.
# No usar 100%; Linux, minimap2 y samtools necesitan margen.
RAM_PERCENT=85

# Preset de minimap2.
# Para ONT convencional:
PRESET="${PRESET:-map-ont}"

# Para ONT Q20+/R10.4.1 muy preciso puedes ejecutar:
# PRESET=lr:hq bash coverage_ont_batch.sh


# ------------------------------------------------------------
# COMPROBAR PROGRAMAS
# ------------------------------------------------------------

for PROGRAM in minimap2 samtools awk find; do
    if ! command -v "$PROGRAM" >/dev/null 2>&1; then
        echo "ERROR: No se encontro $PROGRAM"
        exit 1
    fi
done


# ------------------------------------------------------------
# CREAR DIRECTORIOS
# ------------------------------------------------------------

mkdir -p "$OUTDIR" "$BAMDIR" "$LOGDIR" "$TMPDIR"


# ------------------------------------------------------------
# DETECTAR CPU
# ------------------------------------------------------------

THREADS=$(nproc)

if [[ "$THREADS" -lt 1 ]]; then
    THREADS=1
fi


# ------------------------------------------------------------
# DETECTAR RAM DISPONIBLE
#
# /proc/meminfo reporta MemAvailable en kB.
#
# samtools sort -m especifica RAM POR THREAD.
# Por eso:
#
# RAM por thread =
#   (MemAvailable * 0.85) / numero de threads
#
# ------------------------------------------------------------

MEM_AVAILABLE_KB=$(awk '/MemAvailable:/ {print $2}' /proc/meminfo)

if [[ -z "${MEM_AVAILABLE_KB}" ]]; then
    echo "ERROR: No se pudo determinar MemAvailable."
    exit 1
fi

TOTAL_AVAILABLE_MB=$(( MEM_AVAILABLE_KB / 1024 ))

USABLE_MB=$(( TOTAL_AVAILABLE_MB * RAM_PERCENT / 100 ))

SORT_MEM_MB=$(( USABLE_MB / THREADS ))

# Evitar valores demasiado pequenos.
if [[ "$SORT_MEM_MB" -lt 256 ]]; then
    SORT_MEM_MB=256
fi


# ------------------------------------------------------------
# INFORMACION DEL SERVIDOR
# ------------------------------------------------------------

echo
echo "============================================================"
echo " COBERTURA DE GENOMAS ONT"
echo "============================================================"
echo
echo "CPU disponibles:          ${THREADS}"
echo "RAM disponible:           ${TOTAL_AVAILABLE_MB} MB"
echo "RAM asignable a sort:     ${USABLE_MB} MB (${RAM_PERCENT}%)"
echo "RAM por thread de sort:   ${SORT_MEM_MB} MB"
echo "Preset minimap2:          ${PRESET}"
echo "Directorio de resultados: ${OUTDIR}"
echo
echo "============================================================"
echo


# ------------------------------------------------------------
# CREAR CSV
# ------------------------------------------------------------

echo \
"sample,fasta,fastq_files,contigs,assembly_bp,total_reads,primary_mapped_reads,mapped_percent,total_read_bases,raw_coverage_x,mean_depth_x,covered_bases_percent,ncbi_coverage,status" \
> "$CSV"


# ------------------------------------------------------------
# CONTADOR DE MUESTRAS
# ------------------------------------------------------------

N=0
SUCCESS=0
FAILED=0


# ------------------------------------------------------------
# RECORRER TODOS LOS FASTA
# ------------------------------------------------------------

while IFS= read -r -d '' FASTA
do

    N=$((N + 1))

    FASTA_NAME=$(basename "$FASTA")

    # Quitar extension .fasta / .fa / .fna
    SAMPLE="${FASTA_NAME%.*}"

    # Eliminar espacios accidentales al principio/final.
    # Esto corrige casos como:
    # "005A-TAG-037 .fasta"
    SAMPLE=$(printf '%s' "$SAMPLE" | \
        sed 's/^[[:space:]]*//;s/[[:space:]]*$//')

    # Nombre seguro para BAM/logs
    SAFE_SAMPLE=$(printf '%s' "$SAMPLE" | \
        sed 's/[^A-Za-z0-9._-]/_/g')

    BAM="${BAMDIR}/${SAFE_SAMPLE}.sorted.bam"
    LOG="${LOGDIR}/${SAFE_SAMPLE}.log"
    TMPPREFIX="${TMPDIR}/${SAFE_SAMPLE}"

    echo
    echo "------------------------------------------------------------"
    echo "Muestra: $SAMPLE"
    echo "FASTA:   $FASTA_NAME"
    echo "------------------------------------------------------------"


    # --------------------------------------------------------
    # BUSCAR FASTQ CORRESPONDIENTES
    #
    # Busca archivos que COMIENCEN con el nombre de la muestra:
    #
    # sample.fastq.gz
    # sample_001.fastq.gz
    # sample_run1.fastq.gz
    # etc.
    # --------------------------------------------------------

    FASTQS=()

    while IFS= read -r -d '' FQ
    do
        FASTQS+=("$FQ")
    done < <(
        find -L "$DIR" -maxdepth 1 -type f \
        \( \
            -iname "${SAMPLE}*.fastq.gz" -o \
            -iname "${SAMPLE}*.fq.gz"    -o \
            -iname "${SAMPLE}*.fastq"    -o \
            -iname "${SAMPLE}*.fq" \
        \) \
        -print0 | sort -z
    )


    # --------------------------------------------------------
    # INFORMACION DEL ENSAMBLAJE
    # --------------------------------------------------------

    CONTIGS=$(awk '
        /^>/ {n++}
        END {print n+0}
    ' "$FASTA")

    ASSEMBLY_BP=$(awk '
        /^>/ {next}
        {
            gsub(/[[:space:]]/, "", $0)
            n += length($0)
        }
        END {print n+0}
    ' "$FASTA")


    echo "Contigs:           $CONTIGS"
    echo "Tamano ensamblaje: $ASSEMBLY_BP bp"


    # --------------------------------------------------------
    # SI NO EXISTEN FASTQ
    # --------------------------------------------------------

    if [[ "${#FASTQS[@]}" -eq 0 ]]; then

        echo "ADVERTENCIA: No se encontro FASTQ para $SAMPLE"

        echo \
"\"$SAMPLE\",\"$FASTA_NAME\",NA,$CONTIGS,$ASSEMBLY_BP,NA,NA,NA,NA,NA,NA,NA,NA,NO_FASTQ" \
        >> "$CSV"

        FAILED=$((FAILED + 1))

        continue
    fi


    echo "FASTQ encontrados: ${#FASTQS[@]}"

    for FQ in "${FASTQS[@]}"; do
        echo "  -> $(basename "$FQ")"
    done


    # Crear una cadena de nombres para el CSV
    FASTQ_NAMES=""

    for FQ in "${FASTQS[@]}"; do

        NAME=$(basename "$FQ")

        if [[ -z "$FASTQ_NAMES" ]]; then
            FASTQ_NAMES="$NAME"
        else
            FASTQ_NAMES="${FASTQ_NAMES};${NAME}"
        fi

    done


    # --------------------------------------------------------
    # MAPEO ONT -> ENSAMBLAJE
    #
    # minimap2 usa todos los CPU.
    #
    # samtools sort:
    #   -@ numero de threads
    #   -m memoria maxima aproximada POR THREAD
    #
    # --------------------------------------------------------

    echo
    echo "Mapeando lecturas ONT..."

    {
        echo "================================================="
        echo "Sample: $SAMPLE"
        echo "FASTA:  $FASTA"
        echo "Preset: $PRESET"
        echo "Threads: $THREADS"
        echo "RAM/thread samtools sort: ${SORT_MEM_MB}M"
        echo "================================================="
    } > "$LOG"


    if ! minimap2 \
        -ax "$PRESET" \
        -t "$THREADS" \
        "$FASTA" \
        "${FASTQS[@]}" \
        2>> "$LOG" \
        | samtools sort \
            -@ "$THREADS" \
            -m "${SORT_MEM_MB}M" \
            -T "$TMPPREFIX" \
            -o "$BAM" \
            - \
            2>> "$LOG"
    then

        echo "ERROR durante minimap2/samtools para $SAMPLE"
        echo "Revisa: $LOG"

        echo \
"\"$SAMPLE\",\"$FASTA_NAME\",\"$FASTQ_NAMES\",$CONTIGS,$ASSEMBLY_BP,NA,NA,NA,NA,NA,NA,NA,NA,MAPPING_ERROR" \
        >> "$CSV"

        FAILED=$((FAILED + 1))

        rm -f "$BAM"

        continue
    fi


    # --------------------------------------------------------
    # VERIFICAR BAM
    # --------------------------------------------------------

    if ! samtools quickcheck "$BAM"; then

        echo "ERROR: BAM invalido para $SAMPLE"

        echo \
"\"$SAMPLE\",\"$FASTA_NAME\",\"$FASTQ_NAMES\",$CONTIGS,$ASSEMBLY_BP,NA,NA,NA,NA,NA,NA,NA,NA,BAM_ERROR" \
        >> "$CSV"

        FAILED=$((FAILED + 1))

        continue
    fi


    # --------------------------------------------------------
    # INDEXAR BAM
    # --------------------------------------------------------

    echo "Indexando BAM..."

    samtools index -@ "$THREADS" "$BAM"


    # --------------------------------------------------------
    # TOTAL DE LECTURAS
    #
    # -F 0x900:
    # excluye SECONDARY y SUPPLEMENTARY.
    #
    # De esta forma contamos cada read una sola vez.
    # --------------------------------------------------------

    TOTAL_READS=$(samtools view \
        -@ "$THREADS" \
        -c \
        -F 0x900 \
        "$BAM")


    # --------------------------------------------------------
    # LECTURAS PRIMARIAS MAPEADAS
    #
    # 0x904 =
    # 0x900 secondary/supplementary
    # +
    # 0x004 unmapped
    # --------------------------------------------------------

    MAPPED_READS=$(samtools view \
        -@ "$THREADS" \
        -c \
        -F 0x904 \
        "$BAM")


    # --------------------------------------------------------
    # PORCENTAJE DE MAPPING
    # --------------------------------------------------------

    MAP_PERCENT=$(awk \
        -v m="$MAPPED_READS" \
        -v t="$TOTAL_READS" \
        'BEGIN {
            if (t > 0)
                printf "%.2f", (m/t)*100
            else
                printf "0.00"
        }')


    # --------------------------------------------------------
    # TOTAL DE BASES EN LAS LECTURAS
    #
    # Se utilizan solamente registros primarios para evitar
    # contar una misma lectura varias veces.
    # --------------------------------------------------------

    TOTAL_READ_BASES=$(samtools view \
        -@ "$THREADS" \
        -F 0x900 \
        "$BAM" \
        | awk '
            $10 != "*" {
                total += length($10)
            }
            END {
                printf "%.0f\n", total+0
            }
        ')


    # --------------------------------------------------------
    # COBERTURA TEORICA
    #
    # total de bases FASTQ / tamano ensamblaje
    # --------------------------------------------------------

    RAW_COVERAGE=$(awk \
        -v b="$TOTAL_READ_BASES" \
        -v g="$ASSEMBLY_BP" \
        'BEGIN {
            if (g > 0)
                printf "%.2f", b/g
            else
                printf "0.00"
        }')


    # --------------------------------------------------------
    # COBERTURA REAL POR MAPEO
    #
    # -aa incluye TODAS las posiciones del ensamblaje,
    # tambien las posiciones con profundidad cero.
    #
    # Calculamos simultaneamente:
    #
    # 1) mean depth
    # 2) porcentaje de bases cubiertas >=1x
    #
    # --------------------------------------------------------

    DEPTH_RESULT=$(samtools depth -aa "$BAM" | \
        awk '
        {
            sum += $3
            n++

            if ($3 > 0)
                covered++
        }
        END {
            if (n > 0) {
                printf "%.2f %.2f\n",
                    sum/n,
                    (covered/n)*100
            }
            else {
                printf "0.00 0.00\n"
            }
        }')


    MEAN_DEPTH=$(echo "$DEPTH_RESULT" | awk '{print $1}')
    COVERED_PERCENT=$(echo "$DEPTH_RESULT" | awk '{print $2}')


    # --------------------------------------------------------
    # FORMATO PARA NCBI
    # --------------------------------------------------------

    NCBI_COVERAGE="${MEAN_DEPTH}x"


    # --------------------------------------------------------
    # MOSTRAR RESULTADOS
    # --------------------------------------------------------

    echo
    echo "Resultados:"
    echo
    echo "  Contigs:                  $CONTIGS"
    echo "  Tamano del ensamblaje:    $ASSEMBLY_BP bp"
    echo "  Lecturas ONT:             $TOTAL_READS"
    echo "  Lecturas mapeadas:        $MAPPED_READS"
    echo "  Mapping:                  ${MAP_PERCENT}%"
    echo "  Bases en lecturas:        $TOTAL_READ_BASES"
    echo "  Cobertura teorica:        ${RAW_COVERAGE}x"
    echo "  Profundidad media:        ${MEAN_DEPTH}x"
    echo "  Bases cubiertas >=1x:     ${COVERED_PERCENT}%"
    echo
    echo "  >>> NCBI coverage:        $NCBI_COVERAGE"
    echo


    # --------------------------------------------------------
    # ESCRIBIR CSV
    # --------------------------------------------------------

    echo \
"\"$SAMPLE\",\"$FASTA_NAME\",\"$FASTQ_NAMES\",$CONTIGS,$ASSEMBLY_BP,$TOTAL_READS,$MAPPED_READS,$MAP_PERCENT,$TOTAL_READ_BASES,$RAW_COVERAGE,$MEAN_DEPTH,$COVERED_PERCENT,\"$NCBI_COVERAGE\",OK" \
    >> "$CSV"


    SUCCESS=$((SUCCESS + 1))

done < <(
    find -L "$DIR" -maxdepth 1 -type f \
    \( \
        -iname "*.fasta" -o \
        -iname "*.fa"    -o \
        -iname "*.fna" \
    \) \
    -print0 | sort -z
)


# ------------------------------------------------------------
# FINAL
# ------------------------------------------------------------

echo
echo "============================================================"
echo " PROCESAMIENTO TERMINADO"
echo "============================================================"
echo
echo "Genomas detectados:   $N"
echo "Procesados OK:         $SUCCESS"
echo "Con error/sin FASTQ:   $FAILED"
echo
echo "CSV:"
echo "  $CSV"
echo
echo "BAM:"
echo "  $BAMDIR"
echo
echo "Logs:"
echo "  $LOGDIR"
echo
echo "============================================================"

# Activar ambiente
#conda activate coverage

# 1. Mapear las lecturas ONT contra el ensamblaje
#minimap2 -ax map-ont -t 16 \
#    005-AAC-051.fasta \
#    005-AAC-051.fastq.gz | \
#samtools sort -@ 8 -o 005A-004-JEJ.sorted.bam

#conda deactivate

# 2. Indexar el BAM
#samtools index 005A-004-JEJ.sorted.bam

# 3. Revisar porcentaje de lecturas mapeadas
#samtools flagstat 005A-004-JEJ.sorted.bam

# 4. Revisar estadísticas de cobertura por contig
#samtools coverage 005A-004-JEJ.sorted.bam

# 5. Calcular la cobertura promedio de TODO el ensamblaje
#    -aa incluye también posiciones con cobertura cero
#samtools depth -aa 005A-004-JEJ.sorted.bam | \
#awk '{sum += $3; n++} END {if (n > 0) printf "%.2fx\n", sum/n}'
