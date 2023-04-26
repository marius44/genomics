#!/bin/bash
## Este script va a generar Bases de datos de Anvi'o. Además estimará los metabolismos de varios genomas.
# Las bases de datos generadas se pueden acomplar a otros analisis posteriores como la estimación de
# pangenomas.

#La version de anvio cargada actualmente (junio 2022) es la 7.1
# La línea esta comentada porque de momento no es posible activarlo desde el script
#conda activate anvio-7.1

#Crea lista_genomas.txt solo con los fasta en la carpeta
echo "Creando lista de genomas a analizar: "
ls *.fasta >> lista_genomas.txt

# De la lista de genomas, los lee y les da formato para anvio 7
echo "Haciendo los archivos fasta compatibles con Anvi'o: "
cat lista_genomas.txt | while read line;do
   anvi-script-reformat-fasta --seq-type NT "$line" -o $line.nts.fasta -l 0 --simplify-names
done

# crea una lista con los fasta con formato correcto para leerse en sl siguiente paso
echo "Creando lista de genomas con formato correcto: "
ls *.fasta.nts.fasta >> lista.nts.txt

# Crea las db necesarias para anvio
echo "Creando bases de datos de Anvi'o: "
cat lista.nts.txt | while read line;do
   anvi-gen-contigs-database -f "$line" -o $line.contigs.db
done

# crea una lista con los contigs.db con formato correcto para leerse en sl siguiente paso
echo "Creando lista de bases de datos: "
ls *.contigs.db >> lista.db.txt

# Anotación
echo "Calculando hidden Markov models (hmms). Para anotar genes encontrados: "
cat lista.db.txt | while read line;do
    anvi-run-hmms -c "$line" -T 37
done

# KEGG 
echo "Comparando contra KEGG: "
cat lista.db.txt | while read line;do
    anvi-run-kegg-kofams -c "$line" -T 37 --just-do-it --skip-bitscore-heuristic
done

# Metabolismo
echo "Estimando Metabolismos: "
cat lista.db.txt | while read line;do
    anvi-estimate-metabolism -c "$line" -O "$line"
done

echo "Borrando listas: "
rm lista.db.txt
rm lista_genomas.txt
rm lista.nts.txt
