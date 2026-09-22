#!/bin/bash
# Primero debemos activar conda
# Para hacerlo sin un error CondaError: Run 'conda init' before 'conda activate'
# debemos declarar desde donde se llama a conda
source /opt/conda/etc/profile.d/conda.sh
conda activate abricate

find "$(pwd)" -maxdepth 1 -type f -name "*.fna" | sort > genomes_to_screen.txt

echo "Creando carpeta de resultados"
mkdir -p vfdb_results

echo "Buscando factores de virulencia "
abricate --version > vfdb_results/abricate_version.txt 2>&1
abricate --list > vfdb_results/database_versions.tsv

printf 'database\tvfdb\nmin_identity\t80\nmin_coverage\t80\nthreads\t36\n' \
> vfdb_results/run_parameters.tsv

abricate \
  --db vfdb \
  --minid 80 \
  --mincov 80 \
  --threads 36 \
  --nopath \
  --fofn genomes_to_screen.txt \
  > vfdb_results/vfdb_all_hits.tsv \
  2> vfdb_results/vfdb_run.log


echo "Haciendo resumen"
abricate --summary vfdb_results/vfdb_all_hits.tsv > vfdb_results/vfdb_summary_coverage.tsv

echo "Haciendo matriz"

awk '
BEGIN {
    FS=OFS="\t"
}
NR == 1 {
    print
    next
}
{
    for (i=3; i<=NF; i++) {
        $i = ($i == "." ? 0 : 1)
    }
    print
}
' vfdb_results/vfdb_summary_coverage.tsv \
> vfdb_results/vfdb_summary_binary.tsv
rm genomes_to_screen.txt

echo "Hecho"
