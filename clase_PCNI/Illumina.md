# Pipeline para ensamblar un genoma con lecturas cortas (illumina)
# Primero haremos un Analisis de calidad con AfterQC 
pypy /home/user/AfterQC/after.py -1 R1.fastq -2 R2.fastq (-d dir ) -f -1 -t -1
### -d permite colocar la carpeta con todos los fastq que querramos analizar sin tener que hacer una lista o hacerlo uno por uno

## Ensamble con unicycler/spades
unicycler -1 R1.fastq -2 R2.fq -o out_folder --spades_path /home/user/SPAdes-3.15.4-Linux/bin/spades.py -t 15 --keep 3 --verbosity 2 

## Calidad de los ensambles con Quast
python /home/user/quast-5.2.0/quast.py archivo.fasta -o /directorio_salida/ -t 20 

## Anotacion prokka 
docker run --rm --cpus 10 -v /ruta/donde/estan/los/fasta:/data staphb/prokka:latest prokka TU_ENSAMBLE.fasta --outdir anotacion --prefix PREFIJO_PARA_DIFERENCIAR_MUESTRAS --force --cpus 10

## Anotacion bakta
bash /home/user/bakta-podman.sh --db /home/user/Escritorio/databases/bakta/ --output anotacio_bakta TU_ENSAMBLE.fasta


# Pipeline para ensamblar un genoma con lecturas largas (ONT)
## Concatenar todos los fastq que se generaron por corrida
## Eg.
zcat barcode01/*.gz >> Nombre_muestra.fastq
mkdir fastqc

# Revisión de calidad
fastqc -t 44 -o fastqc *.gz


# Hybracter (FiltLong/porechop/Ensamble)
## Se debe hacer un archivo input_SU_MUESTRA.csv
## Debe tener dos columnas separadas por comas. Serán el nombre de la muestra y el nombre del archivo fastq
## E.g.
s_aureus_sample1,sample1_long_read.fastq.gz

# Crear una carpeta con el nombre del genoma que será ensamblado
mkdir SU_MUESTRA_hybracter

# Activar el entorno conda de hybracter
conda activate hybracter

# Se coloca el comando para ensamblar con hybracter
hybracter long --input SU_MUESTRA.csv --databases /data/databases/hybracter --output SU_MUESTRA_hybracter -t 30 --auto

# Quast
quast.py *.fasta -o quast -t 5 --circos 

# CheckM2
conda activate checkm2
checkm2 predict --threads 4 -x fasta --force --input . --output-directory checkm2/

# Anotacion
## Activar el entorno conda de prokka
conda activate prokka
prokka NG40_final.fasta --outdir anotacion_NG18 --prefix NG29_ --force --cpus 4

# Hypro anotacion extra
 nextflow run hoelzer-lab/hypro -r 0.0.4 -profile local,conda --fasta ./NG19.fna --database uniprotkb --output ./hypro_results --customdb ~/Escritorio/Fagos_prokka/nextflow-autodownload-databases/uniprotkb/uniprotkb.fasta
