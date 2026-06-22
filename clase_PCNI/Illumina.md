# Pipeline para ensamblar un genoma con lecturas cortas (illumina)
# Primero haremos un Analisis de calidad con AfterQC 
fastp     -i Hjan13SSA2_S295_L001_R1_001.fastq     -I Hjan13SSA2_S295_L001_R2_001.fastq     -o Hjan13_R1.clean.fastq.gz     -O Hjan13_R2.clean.fastq.gz     --detect_adapter_for_pe     --correction     --trim_poly_g     --html fastp.html     --json fastp.json     -w 16

## Ensamble con unicycler/spades
unicycler -1 R1.fastq -2 R2.fq -o out_folder --spades_path /home/user/SPAdes-3.15.4-Linux/bin/spades.py -t 15 --keep 3 --verbosity 2 

## Quast
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

