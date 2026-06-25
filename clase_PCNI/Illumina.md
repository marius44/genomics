# Pipeline para ensamblar un genoma con lecturas cortas (illumina)
## Primero haremos un Analisis de calidad con fastqc
fastqc -t 40 -o output_fastqc

## Se hace el recorte y control de calidad con fastp para generar los archivos clean con los que alimentaremos al ensamblador  
fastp     -i SU_MUESTRA_R1_001.fastq     -I SU_MUESTRA_R2_001.fastq     -o SU_MUESTRA_R1.clean.fastq.gz     -O SU_MUESTRA_R2.clean.fastq.gz     --detect_adapter_for_pe     --correction     --trim_poly_g     --html fastp.html     --json fastp.json     -w 16

## Ensamble con unicycler/spades
unicycler -1 SU_MUESTRA_R1.clean.fastq.gz -2 SU_MUESTRA_R2.clean.fastq.gz -o out_SU_MUESTRA -t 15 --keep 3 --verbosity 2 

## Quast
## Activar el ambiente conda
conda activate quast

quast.py *.fasta -o quast -t 10 --circos 

## CheckM2
## Activar el ambiente conda
conda activate checkm2
## Correr checkm2
checkm2 predict --threads 10 -x fasta --force --input . --output-directory checkm2/

## Anotacion
### Activar el entorno conda de prokka
conda activate prokka
prokka SU_MUESTRA.fasta --outdir anotacion_SU_MUESTRA --prefix SU_MUESTRA_ --force --cpus 10

# Hypro anotacion extra
 nextflow run hoelzer-lab/hypro -r 0.0.4 -profile local,conda --fasta ./NG19.fna --database uniprotkb --output ./hypro_results --customdb ~/Escritorio/Fagos_prokka/nextflow-autodownload-databases/uniprotkb/uniprotkb.fasta

