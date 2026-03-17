#!/bin/bash
#La version de anvio cargada actualmente (marzo 2026) es la 9
#conda activate anvio-9
#Crea lista_genomas.txt solo con los fasta en la carpeta
#ls CP*.fasta >> lista_genomas.txt

# De la lista de genomas, los lee y les da formato para anvio 7
#cat lista_genomas.txt | while read line;do
#    anvi-script-reformat-fasta --seq-type NT "$line" -o $line.nts.fasta -l 0 --simplify-names
#done

# crea una lista con los fasta con formato correcto para leerse en sl siguiente paso
#ls *.fasta.nts.fasta >> lista.nts.txt

# Crea las db necesarias para anvio
#cat lista.nts.txt | while read line;do
#    anvi-gen-contigs-database -f "$line" -o $line.contigs.db
#done

# crea una lista con los contigs.db con formato correcto para leerse en sl siguiente paso
#ls *.contigs.db >> lista.db.txt

#HMMS
#cat lista.db.txt | while read line;do
#    anvi-run-hmms -c "$line" -T 40
#done

# import taxonomy
#cat lista.db.txt | while read line;do
#    anvi-export-functions -c "$line" -o "$line".txt
#done


#kegg
#cat lista.db.txt | while read line;do
#    anvi-run-kegg-kofams -c "$line" -T 35 --just-do-it --skip-bitscore-heuristic
#done

# Metabolismo
#cat lista.db.txt | while read line;do
#    anvi-estimate-metabolism -c "$line" -O "$line"
#done

#rm lista.db.txt
#rm lista_genomas.txt

# Para concatenar genes
#anvi-get-sequences-for-hmm-hits --external-genomes external_genomes.txt \
#                                -o concatenated-proteins.fa \
#                                --hmm-source Bacteria_71 \
#                                --gene-names Ribosomal_L1,Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6 \
#                                --return-best-hit \
#                                --get-aa-sequences \
#                                --concatenate

# Arbol
anvi-gen-phylogenomic-tree -f concatenated-proteins.fa \
                           -o phylogenomic-tree.txt --just-do-it

anvi-gen-phylogenomic-profiles \      
  -i GENOMES.db \
  -t phylogenomic-tree.txt \
  -o PHYLO_PROFILE



