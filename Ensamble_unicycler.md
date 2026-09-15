# Ensamble de genomas de secuencias pareadas (Illumina)
## Primero haremos un Analisis de calidad con Fastp 
fastp -1 R1.fastq -2 R2.fastq (-d dir ) -f -1 -t -1
#-d permite colocar la carpeta con todos los fastq que querramos analizar sin tener que hacer una lista o hacerlo uno por uno

## Ensamble con unicycler/spades
unicycler -1 R1.fastq -2 R2.fq -o out_folder -t 15 --keep 3 --verbosity 2 

## Calidad de los ensambles con Quast
conda activate quast
quast.py archivo.fasta -o /directorio_salida/ -t 20 

## Anotacion prokka 
prokka TU_ENSAMBLE.fasta --outdir anotacion --prefix PREFIJO_PARA_DIFERENCIAR_MUESTRAS --force --cpus 10

## Anotacion bakta
bash /home/user/bakta-podman.sh --db /home/user/Escritorio/databases/bakta/ --output anotacio_bakta TU_ENSAMBLE.fasta
