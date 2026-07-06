# Hybracter (FiltLong/porechop/Ensamble)

## Para este ensamble usaremos los archivos con secuencias cortas y ya limpias de Illumina (clean) y el archivo con secuencias largas. Copielos a su carpeta

## Se debe hacer un archivo input_SU_MUESTRA.csv (Use el comando nano)

nano input_SU_MUESTRA.csv

## Debe tener cuatro columnas separadas por comas. Serán el nombre de la muestra, el nombre del archivo fastq de secuencias largas, la secuencia R1.clean y finalmente la secuencia R2.clean
## Ejemplo
Hjan13,Hjan13.fastq,Hjan13_R1.clean.fastq.gz,Hjan13_R2.clean.fastq.gz

## Recuerde que debe guardar los cambios con la combinación de teclas "Ctrl + o" y cerrar con "Ctrl + x"

# Crear una carpeta con el nombre del genoma que será ensamblado

mkdir SU_MUESTRA_hybracter

# Activar el entorno conda de hybracter

conda activate hybracter

# Se coloca el comando para ensamblar con hybracter

hybracter hybrid --input hybrid.csv --databases /data/databases/hybracter --output hybrid -t 30 --auto

conda deactivate

# Quast

conda activate quast

cd FINAL_OUTPUT/incomplete/

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
