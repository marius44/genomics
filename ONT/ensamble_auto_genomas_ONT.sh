#!/bin/bash
# Ensamble automatico de ONT

cat genomas.txt | while read line;do #GENOMAS.TXT ES UNA LISTA CON LAS CARPETAS DONDE ESTAN LAS LECTURAS
    echo "Entrando a " $line
    cd $line
    echo "Puliendo " $line
    seqkit seq -g -m 1000 *.fastq.gz -o $line.filtered.fastq # QUITA READS MENORES A 1000 bp
    echo "Ensamblando " $line
    canu -p $line -d assembly_canu_$line genomeSize=2.3m corThreads=20 -nanopore-raw $line.filtered.fastq # ENSAMBLA CON LAS CORRECCIONES
    cd ..
done
