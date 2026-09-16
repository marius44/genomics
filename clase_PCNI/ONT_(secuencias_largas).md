## Revisión de calidad con fastqc
Se hace una carpeta para colocar los resultados y se corre el programa
```bash
mkdir fastqc
```
Se corre fastqc
```bash
fastqc -t 24 -o fastqc *.fastq
```

Al finalizar, entrar a la carpeta fastqc y revisar manualmente cada métrica de calidad. Complementariamente se puede condensar toda la información con multiqc. Para hacerlo en la carpeta fastqc correr

```bash
multiqc .
```

## Ensamble con Hybracter.
Integra remoción de lecturas cortas (FiltLong), corte de adaptadores (porechop) y ensamble (flye). Se debe hacer un archivo input_SU_MUESTRA.csv (Use el comando nano)
```bash
nano input_SU_MUESTRA.csv
```
Debe tener dos columnas separadas por comas. Serán el nombre de la muestra y el nombre del archivo fastq.

Ejemplo
```
s_aureus_sample1,sample1_long_read.fastq
```
**Recuerde que debe guardar los cambios con la combinación de teclas "Ctrl + o" y cerrar con "Ctrl + x"**

Crear una carpeta con el nombre del genoma que será ensamblado
```
mkdir SU_MUESTRA_hybracter
```
Activar el entorno conda de hybracter
```
conda activate hybracter
```
Ensamble con hybracter
```
hybracter long --input SU_MUESTRA.csv --databases /data/databases/hybracter --output SU_MUESTRA_hybracter -t 30 --auto

conda deactivate
```
# Análisis del ensamble generado
Se realiza con Quast
```bash
conda activate quast
```
Se entra a la carpeta de cada ensamble generado
```bash
cd SU_MUESTRA/FINAL_OUTPUT/incomplete/
```
```bash
quast.py *.fasta -o quast -t 5 --circos 
```
Al finalizar el análisis de todos los ensambles se desactiva el entorno

```bash
conda deactivate 
```
# Analisis de completitud del ensamble
Se hace con CheckM2
```bash
conda activate checkm2
```

Entrar a la carpeta correspondiente de cada ensamble y usar el ensamble en formato fasta

```bash
checkm2 predict --threads 4 -x fasta --force --input . --output-directory checkm2/
```

Al finalizar el análisis de todos los ensambles se desactiva el entorno

```bash
conda deactivate 
```
# Anotacion
Se debe activar el entorno conda de prokka
```bash
conda activate prokka
```
Se entra a la carpeta en donde esta cada ensamble en formato fasta y se corre (en el grupo tenemos un script automático para la anotación solicitarlo si son muchos genomas):
```bash
prokka NG40_final.fasta --outdir anotacion_NG18 --prefix NG29_ --force --cpus 4

Al finalizar el análisis de todos los ensambles se desactiva el entorno

```bash
conda deactivate 
```
