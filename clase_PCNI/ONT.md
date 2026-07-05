# Revisión de calidad
fastqc -t 44 -o fastqc *.gz

# Hybracter (FiltLong/porechop/Ensamble)
## Se debe hacer un archivo input_SU_MUESTRA.csv (Use el comando nano)

nano input_SU_MUESTRA.csv

## Debe tener dos columnas separadas por comas. Serán el nombre de la muestra y el nombre del archivo fastq
## Ejemplo
s_aureus_sample1,sample1_long_read.fastq.gz

## Recuerde que debe guardar los cambios con la combinación de teclas "Ctrl + o" y cerrar con "Ctrl + x"

# Crear una carpeta con el nombre del genoma que será ensamblado

mkdir SU_MUESTRA_hybracter

# Activar el entorno conda de hybracter

conda activate hybracter

# Se coloca el comando para ensamblar con hybracter

hybracter long --input SU_MUESTRA.csv --databases /data/databases/hybracter --output SU_MUESTRA_hybracter -t 30 --auto

conda deactivate

# Quast

conda activate quast

quast.py *.fasta -o quast -t 5 --circos 

conda deactivate 

# CheckM2
conda activate checkm2

checkm2 predict --threads 4 -x fasta --force --input . --output-directory checkm2/

conda deactivate 

# Anotacion
## Activar el entorno conda de prokka

conda activate prokka

prokka NG40_final.fasta --outdir anotacion_NG18 --prefix NG29_ --force --cpus 4

conda deactivate 

# Hypro anotacion extra
 nextflow run hoelzer-lab/hypro -r 0.0.4 -profile local,conda --fasta ./NG19.fna --database uniprotkb --output ./hypro_results --customdb ~/Escritorio/Fagos_prokka/nextflow-autodownload-databases/uniprotkb/uniprotkb.fasta
