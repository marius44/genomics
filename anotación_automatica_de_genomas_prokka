#!/bin/bash
# Esto hace una lista con elementos especificos de una extensión

# Pedimos al usuario que ingrese el género de los organismos que va a anotar
echo "Introduce el genero del organismo analizado: "
read genero
echo "El genero es: " $genero
fasta_DIR=.

# Esto hace una lista con elementos especificos de una extensión
# Se alimenta de archivos en formato fasta
echo "Haciendo lista de genomas "
ls *.fasta >> lista_genomas.txt
# Si tienes *.nts.fasta de corrección de formato generado por anvio.
# Corrección disponible en https://anvio.org/help/main/programs/anvi-script-reformat-fasta/
# anvi-script-reformat-fasta Halomonas_salina_VN10_sedA.fasta -o Halomonas_salina_VN10_sedA.nts.fasta --min-len 1000 --seq-type NT --simplify-names --report-file Hsalina-report.txt
# Esto es útil cuando Prokka tiene errores por nombres en los contigs
#ls *nts.fasta >> lista_genomas.txt

# Usamos un bucle while. Toma los las lineas de la lista generada en y da la orden.
# Entre comillas la orden solo toma el nombre y se lo pega a algo, funciona bien para nombrar archivos generados.
# Sin comillas queda como input directamente para el comando.
cat lista_genomas.txt | while read line;do
    echo "procesando " $line
    #docker run --rm -v $fasta_DIR:/data staphb/prokka:latest prokka /data/"$line" --outdir /data/$genero --prefix "$line" --force --cpus 40
    prokka "$line" --outdir $genero --prefix "$line" --force --cpus 40
    echo "Listo, ya quedó anotado " $line
done

echo "Se anotaron todos sus genomas."

# Borrar la lista de genomas
rm lista_genomas.txt

# Pedimos al usuario que ingrese el archivo para hacer metricas
chmod -R 777 $genero
cd $genero

#### HAY QUE CAMBIAR PERMISOS DE LA CARPETA QUE GENERA; LO HACE POR DOCKER
echo "Se calculan las metricas de: " *.gbk

# localizar ribosomales
echo "Numero de ribosomales" >> metricas.txt
grep -E "S ribosomal RNA" *.gbk | grep -v -E -c "note=\"aligned only" >> metricas.txt

echo "5S rRNA" >> metricas.txt
grep -c "product=\"5S ribosomal RNA" *.gbk >> metricas.txt

echo "16s rRNA" >> metricas.txt
grep -c "product=\"16S ribosomal RNA" *.gbk >> metricas.txt

echo "23s rRNA" >> metricas.txt
grep -c "product=\"23S ribosomal RNA" *.gbk >> metricas.txt
echo " " >> metricas.txt

# CDS
echo "Numero de CDS" >> metricas.txt
grep -c "CDS             " *.gbk >> metricas.txt
echo " " >> metricas.txt

# Contigs
echo "Numero de contigs" >> metricas.txt
grep -c "LOCUS" *.gbk >> metricas.txt
echo " " >> metricas.txt

# Metricas de ensamble
#echo "Calidad de ensamble" >> metricas.txt
#more ../quast/report.txt >> metricas.txt

echo "Hecho."
