#!/bin/bash
# Ensamble automatico de ONT

# Lo primero es concatenar todos los fastq que se generaron por corrida
# Eg.
# cat barcode01/*.gz > ./barcode01.fastq.gz
# Revisión de calidad
fastqc -t 40 -o fastqc *.gz

ls *.fastq.gz > adaptadores.txt
cat adaptadores.txt | while read line; do
    porechop_abi --ab_initio -i $line.fastq.gz -o $line.fastq
done    
# Revisión de calidad
fastqc -t 40 -o fastqc_porechop_abi *.fastq

#Creación de lista_genomas
ls *.fastq > genomas.txt

cat genomas.txt | while read line;do #GENOMAS.TXT ES UNA LISTA CON LAS CARPETAS DONDE ESTAN LAS LECTURAS
    echo "Entrando a " $line
    cd $line
    echo "Puliendo " $line
    seqkit seq -g -m 1000 *.fastq -o $line.filtered.fastq # QUITA READS MENORES A 1000 bp
    echo "Ensamblando " $line
    canu -p $line -d assembly_canu_$line genomeSize=2.3m corThreads=20 -nanopore-raw $line.filtered.fastq # ENSAMBLA CON LAS CORRECCIONES
    cd ..
done
