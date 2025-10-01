#!/bin/bash
# Ensamble automatico de ONT

# Concatenar todos los fastq que se generaron por corrida
# Eg.
# cat barcode01/*.gz > ./barcode01.fastq.gz
mkdir fastqc

# Revisión de calidad
echo "Vamos a echale un ojo a la calidad"
fastqc -t 44 -o fastqc *.gz

# Omito el uso de porechop_abi porque desde MinKNOW ya se recortaron los adaptadores

#Trimming usando ProwlerTrimmer 
# ONT da reads con calidades variables
# Las variables en este programa van entre comillas.

mkdir prowler_trim
python ~/ProwlerTrimmer/TrimmerLarge.py   -f "NG36_gel.fastq"   -i "./"   -o "prowler_trim/"   -w 300   -l 500   -c "T"   -g "F0"   -m "D"   -q 4    -r ".fastq"

# -f, 	--file,		filename:	The name of the file you want to trim, wihtout the folderpath"
#-i, 	--infolder, 	inFolder:	The folderpath where your file to be trimmed is located (default = cwd)
#-o, 	--outfolder,	outFolder:	The folderpath where your want to save the trimmed file (default = cwd)
#-w, 	--windowsize,	windowSize:	Change the size of the trimming window (default= 100bp)
#-l, 	--minlen,	minLen:		Change the minimum acceptable numer of bases in a read (default=100)
#-m, 	--trimmode,	mode:		Select trimming algorithm: S for static  or D for dynamic (default=S)
#-q, 	--qscore,	Qcutoff:	Select the phred quality score trimming threshold (default=7)
#-d, 	--datamax,	maxDataMB:	Select a maximum data subsample in MB (default=0, entire file)
#-r, 	--outformat,	outMode:	Select output format of trimmed file (fastq or fasta) (default=.fastq)
#-c, 	--clip,		clipping:	Select L to clip leading Ns, T to trim trialing Ns and LT to trim both (default=LT)
#-g, 	--fragments,	fragments:	Select fragmentation mode (default=U0)


# Revisión de calidad
cd prowler_trim/
fastqc -t 40 -o fastqc *.fastq

# Creación de lista_genomas
ls *porechoped.fastq > porechoped_genomas.txt

cat porechoped_genomas.txt | while read line;do #GENOMAS.TXT ES UNA LISTA CON LAS CARPETAS DONDE ESTAN LAS LECTURAS
    echo "Puliendo " $line
    echo "Ensamblando " $line
    canu -p $line -d assembly_canu_$line genomeSize=2.3m corThreads=20 -nanopore-raw $line.filtered.fastq.gz # ENSAMBLA CON LAS CORRECCIONES
done

## Calidad de los ensambles con Quast
python /home/user/quast-5.2.0/quast.py archivo.fasta -o /directorio_salida/ -t 20 

## Checar contaminacion con CheckM2
micromamba activate checkm2
checkm2
checkm2 predict --threads 40 -x fasta --force --input . --output-directory checkm2/

## Anotacion prokka 
docker run --rm --cpus 10 -v /ruta/donde/estan/los/fasta:/data staphb/prokka:latest prokka TU_ENSAMBLE.fasta --outdir anotacion --prefix PREFIJO_PARA_DIFERENCIAR_MUESTRAS --force --cpus 10

