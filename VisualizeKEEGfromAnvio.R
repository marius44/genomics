library("pathview")
library("ggplot2")
library("ggsci")
library("gridExtra")
library("dplyr")

setwd("/home/user/Escritorio/mario/PGAP_NCBI/Hven4/")

kegg_table_Hven4 <- read.table("/home/user/Escritorio/mario/PGAP_NCBI/Hven4/Hven4_kegg-metabolism_modules.txt", sep = "\t", header = T)
kegg_table_Hven4[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Hven4 <- mutate(kegg_table_Hven4, genome_name="Hven4")
str(kegg_table_Hven4)
kegg_table_Hven7 <- read.table("/home/user/Escritorio/mario/PGAP_NCBI/Hven7/kegg-metabolism_modules.txt", sep = "\t", header = T)
kegg_table_Hven7[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Hven7 <- mutate(kegg_table_Hven7, genome_name="Hven7")
str(kegg_table_Hven7)
kegg_table_Hven9 <- read.table("/home/user/Escritorio/mario/PGAP_NCBI/Hven9/kegg-metabolism_modules.txt", sep = "\t", header = T)
kegg_table_Hven9[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Hven9 <- mutate(kegg_table_Hven9, genome_name="Hven9")
str(kegg_table_Hven9)
kegg_table_Hjan13 <- read.table("/home/user/Escritorio/mario/PGAP_NCBI/Hven4/kegg-metabolism_modules_Hjan13.txt", sep = "\t", header = T)
kegg_table_Hjan13[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Hjan13 <- mutate(kegg_table_Hjan13, genome_name="Hjan13")
str(kegg_table_Hjan13)
kegg_table_Hjan14 <- read.table("/home/user/Escritorio/mario/PGAP_NCBI/Hven4/kegg-metabolism_modules.txt", sep = "\t", header = T)
kegg_table_Hjan14[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Hjan14 <- mutate(kegg_table_Hjan14, genome_name="Hjan14")
str(kegg_table_Hjan14)
kegg_table_bins_metagenome <- read.table("/home/user/Escritorio/mario/metagenomes_WMS_isabel/Analisis_2021/MA_MB/MA/MA_anvio/kegg_contigs2_modules.txt", sep = "\t", header = T)
kegg_table_bins_metagenome[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_bins_metagenome <- mutate(kegg_table_bins_metagenome, genome_name="Metagenome")
str(kegg_table_bins_metagenome)

kegg_table_Halomonas_maxbin_001 <- read.table("~/Escritorio/mario/PGAP_NCBI/MAGs/Halomonas_maxbin.001/kegg-metabolism_modules.txt", sep = "\t", header = T)
kegg_table_Halomonas_maxbin_001[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Halomonas_maxbin_001 <- mutate(kegg_table_Halomonas_maxbin_001, genome_name="Halomonas_maxbin_001")
str(kegg_table_Halomonas_maxbin_001)

kegg_table_Halomonas_metabat2_2 <- read.table("~/Escritorio/mario/PGAP_NCBI/MAGs/Halomonas_metabat2.2/kegg_Halomonas_metabat_modules.txt", sep = "\t", header = T)
kegg_table_Halomonas_metabat2_2[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Halomonas_metabat2_2 <- mutate(kegg_table_Halomonas_metabat2_2, genome_name="Halomonas_metabat2.2")
str(kegg_table_Halomonas_metabat2_2)

kegg_table_nitriluruptor_maxbin_009 <- read.table("~/Escritorio/mario/PGAP_NCBI/MAGs/Nitriliruptor_maxbin.009/Nitriliruptor_maxbin.009_kegg_modules.txt", sep = "\t", header = T)
kegg_table_nitriluruptor_maxbin_009[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_nitriluruptor_maxbin_009 <- mutate(kegg_table_nitriluruptor_maxbin_009, genome_name="Nitriliruptor_maxbin.009")
str(kegg_table_nitriluruptor_maxbin_009)

kegg_table_Lacimicrobium_maxbin.004 <- read.table("~/Escritorio/mario/PGAP_NCBI/MAGs/Lacimicrobium_maxbin.004/Lacimicrobium_maxbin.004_kegg_modules.txt", sep = "\t", header = T)
kegg_table_Lacimicrobium_maxbin.004[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Lacimicrobium_maxbin.004 <- mutate(kegg_table_Lacimicrobium_maxbin.004, genome_name="Lacimicrobium_maxbin.004")
str(kegg_table_Lacimicrobium_maxbin.004)

kegg_table_Lacimicrobium_metabat2.30 <- read.table("~/Escritorio/mario/PGAP_NCBI/MAGs/Lacimicrobium_metabat2.30/Lacimicrobium_metabat2.30_kegg__modules.txt", sep = "\t", header = T)
kegg_table_Lacimicrobium_metabat2.30[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Lacimicrobium_metabat2.30 <- mutate(kegg_table_Lacimicrobium_metabat2.30, genome_name="Lacimicrobium_metabat2.30")
str(kegg_table_Lacimicrobium_metabat2.30)

kegg_table_Synechococcus_metabat2.14 <- read.table("~/Escritorio/mario/PGAP_NCBI/MAGs/Synechococcus_metabat2.14/Synechococcus_metabat2.14_kegg_modules.txt", sep = "\t", header = T)
kegg_table_Synechococcus_metabat2.14[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Synechococcus_metabat2.14 <- mutate(kegg_table_Synechococcus_metabat2.14, genome_name="Synechococcus_metabat2.14")
str(kegg_table_Synechococcus_metabat2.14)



#################################################################"""""""""""2
kegg_Hven4 <- ggplot(kegg_table_Hven4, 
                     aes(x=module_category, 
                         y=module_subcategory,
                         #size=module_completeness,
                         color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("KEGG metabolism modules", subtitle = "Hven4") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

kegg_Hven4_chido <- kegg_Hven4 + theme(legend.position="bottom", legend.box = "horizontal")
kegg_Hven4_chido
ggsave("kegg_Hven4.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)

###
kegg_table7 <- ggplot(kegg_table_Hven7, 
                     aes(x=module_category, 
                         y=module_subcategory,
                         #size=module_completeness,
                         color=module_completeness)) +
  geom_point(size=3) +
  ggtitle("KEGG metabolism modules", subtitle = "Hven7") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

kegg_Hven7_chido <- kegg_table7 + theme(legend.position="bottom", legend.box = "horizontal")
kegg_Hven7_chido
ggsave("kegg_Hven7.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)


###
###

kegg_table9 <- ggplot(kegg_table_Hven9, 
                      aes(x=module_category, 
                          y=module_subcategory,
                          #size=module_completeness,
                          color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("KEGG metabolism modules", subtitle = "Hven9") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

kegg_Hven9_chido <- kegg_table9 + theme(legend.position="bottom", legend.box = "horizontal")
kegg_Hven9_chido
ggsave("kegg_Hven9.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)

#############
kegg_table_max1 <- read.table("/home/user/Escritorio/mario/PGAP_NCBI/MAGs/Halomonas_maxbin.001/kegg-metabolism_modules.txt", sep = "\t", header = T)
str(kegg_table_max1)
kegg_max1 <- ggplot(kegg_table_max1, 
                      aes(x=module_category, 
                          y=module_subcategory,
                          #size=module_completeness,
                          color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("KEGG metabolism modules", subtitle = "maxbin1") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

kegg_max1_chido <- kegg_max1 + theme(legend.position="bottom", legend.box = "horizontal")
kegg_max1_chido
ggsave("kegg_max1.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)



kegg_Hjan13 <- ggplot(kegg_table_Hjan13, 
                     aes(x=module_category, 
                         y=module_subcategory,
                         #size=module_completeness,
                         color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("KEGG metabolism modules", subtitle = "Hjan13") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

kegg_Hjan13_chido <- kegg_Hjan13 + theme(legend.position="bottom", legend.box = "horizontal")
kegg_Hjan13_chido
ggsave("kegg_Hjan13.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)


kegg_Hjan14 <- ggplot(kegg_table_Hjan14, 
                      aes(x=module_category, 
                          y=module_subcategory,
                          #size=module_completeness,
                          color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("KEGG metabolism modules", subtitle = "Hjan14") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

kegg_Hjan14_chido <- kegg_Hjan14 + theme(legend.position="bottom", legend.box = "horizontal")
kegg_Hjan14_chido
ggsave("kegg_Hjan14.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)




###

grid.arrange(kegg_Hven9_chido, kegg_Hven7_chido,kegg_Hven4_chido)

#
#kegg_kofam_table_Hven4 <- read.table("kegg-metabolism_kofam_hits.txt", sep = "\t", header = T)
#str(kegg_kofam_table_Hven4)
#kofam_table <- ggplot(kegg_kofam_table_Hven4, aes(x=ko)) +
#  geom_bar() + 
#  theme(axis.text.x = element_text(angle = 90))+
#  scale_fill_igv()
#kofam_table
#ggsave("kegg_table.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)



#################################3
#Vamos a filtrar con dplyr

Lipids_Hven7 <- filter(kegg_table_Hven7, module_category == "Lipid metabolism")
xeno_Hven7 <- filter(kegg_table_Hven7, module_category == "Xenobiotics biodegradation")
xeno_Hven4 <- filter(kegg_table_Hven4, module_category == "Xenobiotics biodegradation")
xeno_Hven9 <- filter(kegg_table_Hven9, module_category == "Xenobiotics biodegradation")
xeno_Hven9 <- filter(kegg_table_Hven9, module_category == "Xenobiotics biodegradation")

xeno_maxbin <- filter(kegg_table_max1, module_category == "Xenobiotics biodegradation")


#
polyketides_Hven9 <- filter(kegg_table_Hven9, module_category == "Biosynthesis of other secondary metabolites")
polyketides_Hven7 <- filter(kegg_table_Hven7, module_category == "Biosynthesis of other secondary metabolites")
polyketides_Hven7_plot <- ggplot(polyketides_Hven7, 
                      aes(x=module_category, 
                          y=module_subcategory,
                          #size=module_completeness,
                          color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("poliketydes", subtitle = "Hven7") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()
polyketides_Hven7_plot
ggsave("polyke_hven7.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)

xeno_Hven7_plot <- polyketides_Hven7_plot + theme(legend.position="bottom", legend.box = "horizontal")
xeno_Hven7_plot
#


xeno_Hven7_plot <- ggplot(xeno_Hven7,
                          aes(x=module_category, 
                          y=module_subcategory,
                          #size=module_completeness,
                          color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("Xenobiotics", subtitle = "Hven7") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

xeno_Hven7_plot <- xeno_Hven7_plot + theme(legend.position="bottom", legend.box = "horizontal")
xeno_Hven7_plot
ggsave("xeno_hven7.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)
##############33
##################
xeno_maxbin1_plot <- ggplot(xeno_maxbin,
                          aes(x=module_name, 
                              y=module_completeness,
                              fill=as.factor(module_completeness))) +
  geom_col() + 
  ggtitle("Xenobiotics", subtitle = "Maxbin1") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()+ 
  legend

xeno_maxbin1_plot <- xeno_maxbin1_plot + theme(legend.position="bottom", legend.box = "horizontal")
xeno_maxbin1_plot
ggsave("xeno_maxbin1.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)




####################################33


modulos_completos_maxbin1 <- xeno_maxbin <- filter(kegg_table_max1, module_is_complete == "True")

modulos_completos_maxbin1_a <- ggplot(modulos_completos_maxbin1, 
                      aes(x=module_category, 
                          y=module_subcategory,
                          #size=module_completeness,
                          color=module_completeness)) +
  geom_bar() + 
  ggtitle("KEGG metabolism modules", subtitle = "MAG Maxbin1") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

modulos_completos_maxbin1_a <- modulos_completos_maxbin1_a + theme(legend.position="bottom", legend.box = "horizontal")
modulos_completos_maxbin1_a

xeno_Hven7


####################
#PARA VER TODOS EN UNA SOLA GRAFICA
#UNIR DATA FRAMES
#Con tidyverse
#install.packages("tidyverse")
library("tidyverse")
tablas_unidas_Hjan13_y_Hjan13 <- tablas_unidas_Hjan13_y_Hjan13 + theme(legend.position="bottom", legend.box = "horizontal")

str(kegg_table_Hjan13)
str(kegg_table_Hjan13)
tablas_unidas_Hjan13_y_Hjan13 <- ggplot(tablas_unidas, 
                                        aes(x=genome_name, 
                                            y=module_subcategory,
                                            #size=module_completeness,
                                            color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("KEGG metabolism modules", subtitle = "UNidas") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

tablas_unidas_Hjan13_y_Hjan13


tablas_unidas_Hven4_Hven7_Hven9 <- bind_rows(kegg_table_Hven4, kegg_table_Hven7, kegg_table_Hven9)
str(tablas_unidas_Hven4_Hven7_Hven9)
#View(tablas_unidas_Hven4_Hven7_Hven9)
tablas_unidas_Hven4_Hven7_Hven9 <- ggplot(tablas_unidas_Hven4_Hven7_Hven9, 
                                          aes(x=genome_name, 
                                              y=module_subcategory,
                                              #size=module_completeness,
                                              color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("KEGG metabolism modules", subtitle = "UNidas") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()


tablas_unidas_Hven4_Hven7_Hven9

all_tables <- bind_rows(kegg_table_Hven4, kegg_table_Hven7, kegg_table_Hven9, kegg_table_Hjan13, kegg_table_Hjan14)
str(all_tables)
#View(all_tables)
all_tables_kegg <- ggplot(all_tables,
                          aes(x=genome_name, 
                          y=module_subcategory,
                          #size=module_completeness,
                          color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("KEGG metabolism modules", subtitle = "All samples") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()
all_tables_kegg

#Otros graficos
#MEJOR PORQUE ES UN HEATMAP
#TAMBIÉN PUSE ELA LEYENDA CORRIDA
all_tables_kegg_tiles <- ggplot(all_tables,
                          aes(x=genome_name, 
                              y=module_subcategory,
                              fill=module_completeness)) +
  geom_tile() + 
  ggtitle("KEGG metabolism modules", subtitle = "All samples") +
  guides(fill=guide_legend(title="Completeness (fraction)"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  xlab("Strain") +
  ylab("Module subcategory") +
  scale_fill_gradient( low = "#132B43",
                       high = "#56B1F7",
                       space = "Lab",
                       guide = FALSE,
                       aesthetics = "fill") #Esto es para que el color sea un gradiente dependiente de los valores de fill
all_tables_kegg_tiles 


library("pathview")
library("ggplot2")
library("ggsci")
library("gridExtra")
library("dplyr")

setwd("/home/user/Escritorio/mario/PGAP_NCBI/Hven4/")

kegg_table_Hven4 <- read.table("/home/user/Escritorio/mario/PGAP_NCBI/Hven4/Hven4_kegg-metabolism_modules.txt", sep = "\t", header = T)
str(kegg_table_Hven4)
kegg_Hven4 <- ggplot(kegg_table_Hven4, 
                     aes(x=module_category, 
                         y=module_subcategory,
                         #size=module_completeness,
                         color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("KEGG metabolism modules", subtitle = "Hven4") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

kegg_Hven4_chido <- kegg_Hven4 + theme(legend.position="bottom", legend.box = "horizontal")
kegg_Hven4_chido
ggsave("kegg_Hven4.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)

###
kegg_table_Hven7 <- read.table("/home/user/Escritorio/mario/PGAP_NCBI/Hven7/kegg-metabolism_modules.txt", sep = "\t", header = T)
str(kegg_table_Hven7)
kegg_table7 <- ggplot(kegg_table_Hven7, 
                      aes(x=module_category, 
                          y=module_subcategory,
                          #size=module_completeness,
                          color=module_completeness)) +
  geom_point(size=3) +
  ggtitle("KEGG metabolism modules", subtitle = "Hven7") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

kegg_Hven7_chido <- kegg_table7 + theme(legend.position="bottom", legend.box = "horizontal")
kegg_Hven7_chido
ggsave("kegg_Hven7.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)


###
###
kegg_table_Hven9 <- read.table("/home/user/Escritorio/mario/PGAP_NCBI/Hven9/kegg-metabolism_modules.txt", sep = "\t", header = T)
str(kegg_table_Hven9)
kegg_table9 <- ggplot(kegg_table_Hven9, 
                      aes(x=module_category, 
                          y=module_subcategory,
                          #size=module_completeness,
                          color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("KEGG metabolism modules", subtitle = "Hven9") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

kegg_Hven9_chido <- kegg_table9 + theme(legend.position="bottom", legend.box = "horizontal")
kegg_Hven9_chido
ggsave("kegg_Hven9.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)

#############
kegg_table_max1 <- read.table("/home/user/Escritorio/mario/PGAP_NCBI/MAGs/Halomonas_maxbin.001/kegg-metabolism_modules.txt", sep = "\t", header = T)
str(kegg_table_max1)
kegg_max1 <- ggplot(kegg_table_max1, 
                    aes(x=module_category, 
                        y=module_subcategory,
                        #size=module_completeness,
                        color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("KEGG metabolism modules", subtitle = "maxbin1") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

kegg_max1_chido <- kegg_max1 + theme(legend.position="bottom", legend.box = "horizontal")
kegg_max1_chido
ggsave("kegg_max1.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)


kegg_table_Hjan13 <- read.table("/home/user/Escritorio/mario/PGAP_NCBI/Hven4/kegg-metabolism_modules_Hjan13.txt", sep = "\t", header = T)
str(kegg_table_Hjan13)
kegg_Hjan13 <- ggplot(kegg_table_Hjan13, 
                      aes(x=module_category, 
                          y=module_subcategory,
                          #size=module_completeness,
                          color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("KEGG metabolism modules", subtitle = "Hjan13") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

kegg_Hjan13_chido <- kegg_Hjan13 + theme(legend.position="bottom", legend.box = "horizontal")
kegg_Hjan13_chido
ggsave("kegg_Hjan13.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)

kegg_table_Hjan14 <- read.table("/home/user/Escritorio/mario/PGAP_NCBI/Hven4/kegg-metabolism_modules.txt", sep = "\t", header = T)
str(kegg_table_Hjan14)
kegg_Hjan14 <- ggplot(kegg_table_Hjan14, 
                      aes(x=module_category, 
                          y=module_subcategory,
                          #size=module_completeness,
                          color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("KEGG metabolism modules", subtitle = "Hjan14") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

kegg_Hjan14_chido <- kegg_Hjan14 + theme(legend.position="bottom", legend.box = "horizontal")
kegg_Hjan14_chido
ggsave("kegg_Hjan14.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)




###

grid.arrange(kegg_Hven9_chido, kegg_Hven7_chido,kegg_Hven4_chido)

#
#kegg_kofam_table_Hven4 <- read.table("kegg-metabolism_kofam_hits.txt", sep = "\t", header = T)
#str(kegg_kofam_table_Hven4)
#kofam_table <- ggplot(kegg_kofam_table_Hven4, aes(x=ko)) +
#  geom_bar() + 
#  theme(axis.text.x = element_text(angle = 90))+
#  scale_fill_igv()
#kofam_table
#ggsave("kegg_table.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)



#################################3
#Vamos a filtrar con dplyr

Lipids_Hven7 <- filter(kegg_table_Hven7, module_category == "Lipid metabolism")
xeno_Hven7 <- filter(kegg_table_Hven7, module_category == "Xenobiotics biodegradation")
xeno_Hven4 <- filter(kegg_table_Hven4, module_category == "Xenobiotics biodegradation")
xeno_Hven9 <- filter(kegg_table_Hven9, module_category == "Xenobiotics biodegradation")
xeno_Hven9 <- filter(kegg_table_Hven9, module_category == "Xenobiotics biodegradation")

xeno_maxbin <- filter(kegg_table_max1, module_category == "Xenobiotics biodegradation")


#
polyketides_Hven9 <- filter(kegg_table_Hven9, module_category == "Biosynthesis of other secondary metabolites")
polyketides_Hven7 <- filter(kegg_table_Hven7, module_category == "Biosynthesis of other secondary metabolites")
polyketides_Hven7_plot <- ggplot(polyketides_Hven7, 
                                 aes(x=module_category, 
                                     y=module_subcategory,
                                     #size=module_completeness,
                                     color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("poliketydes", subtitle = "Hven7") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()
polyketides_Hven7_plot
ggsave("polyke_hven7.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)

xeno_Hven7_plot <- polyketides_Hven7_plot + theme(legend.position="bottom", legend.box = "horizontal")
xeno_Hven7_plot
#


xeno_Hven7_plot <- ggplot(xeno_Hven7,
                          aes(x=module_category, 
                              y=module_subcategory,
                              #size=module_completeness,
                              color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("Xenobiotics", subtitle = "Hven7") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

xeno_Hven7_plot <- xeno_Hven7_plot + theme(legend.position="bottom", legend.box = "horizontal")
xeno_Hven7_plot
ggsave("xeno_hven7.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)
##############33
##################
xeno_maxbin1_plot <- ggplot(xeno_maxbin,
                            aes(x=module_name, 
                                y=module_completeness,
                                fill=as.factor(module_completeness))) +
  geom_col() + 
  ggtitle("Xenobiotics", subtitle = "Maxbin1") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()+ 
  legend

xeno_maxbin1_plot <- xeno_maxbin1_plot + theme(legend.position="bottom", legend.box = "horizontal")
xeno_maxbin1_plot
ggsave("xeno_maxbin1.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)




####################################33


modulos_completos_maxbin1 <- xeno_maxbin <- filter(kegg_table_max1, module_is_complete == "True")

modulos_completos_maxbin1_a <- ggplot(modulos_completos_maxbin1, 
                                      aes(x=module_category, 
                                          y=module_subcategory,
                                          #size=module_completeness,
                                          color=module_completeness)) +
  geom_bar() + 
  ggtitle("KEGG metabolism modules", subtitle = "MAG Maxbin1") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

modulos_completos_maxbin1_a <- modulos_completos_maxbin1_a + theme(legend.position="bottom", legend.box = "horizontal")
modulos_completos_maxbin1_a

xeno_Hven7


####################
#PARA VER TODOS EN UNA SOLA GRAFICA
#UNIR DATA FRAMES
#Con tidyverse
#install.packages("tidyverse")
library("tidyverse")
tablas_unidas_Hjan13_y_Hjan13 <- tablas_unidas_Hjan13_y_Hjan13 + theme(legend.position="bottom", legend.box = "horizontal")

str(kegg_table_Hjan13)
str(kegg_table_Hjan13)
tablas_unidas_Hjan13_y_Hjan13 <- ggplot(tablas_unidas, 
                                        aes(x=genome_name, 
                                            y=module_subcategory,
                                            #size=module_completeness,
                                            color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("KEGG metabolism modules", subtitle = "UNidas") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

tablas_unidas_Hjan13_y_Hjan13


tablas_unidas_Hven4_Hven7_Hven9 <- bind_rows(kegg_table_Hven4, kegg_table_Hven7, kegg_table_Hven9)
str(tablas_unidas_Hven4_Hven7_Hven9)
#View(tablas_unidas_Hven4_Hven7_Hven9)
tablas_unidas_Hven4_Hven7_Hven9 <- ggplot(tablas_unidas_Hven4_Hven7_Hven9, 
                                          aes(x=genome_name, 
                                              y=module_subcategory,
                                              #size=module_completeness,
                                              color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("KEGG metabolism modules", subtitle = "UNidas") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()


tablas_unidas_Hven4_Hven7_Hven9

all_tables <- bind_rows(kegg_table_Hven4, kegg_table_Hven7, kegg_table_Hven9, kegg_table_Hjan13, kegg_table_Hjan14)
str(all_tables)
#View(all_tables)
all_tables_kegg <- ggplot(all_tables,
                          aes(x=genome_name, 
                              y=module_subcategory,
                              #size=module_completeness,
                              color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("KEGG metabolism modules", subtitle = "All samples") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()
all_tables_kegg

#Otros graficos
#MEJOR PORQUE ES UN HEATMAP
#TAMBIÉN PUSE ELA LEYENDA CORRIDA
all_tables_kegg_tiles <- ggplot(all_tables,
                                aes(x=genome_name, 
                                    y=module_subcategory,
                                    fill=module_completeness)) +
  geom_tile() + 
  ggtitle("KEGG metabolism modules", subtitle = "All samples") +
  guides(fill=guide_legend(title="Completeness (fraction)"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  xlab("Strain") +
  ylab("Module subcategory") +
  scale_fill_gradient( low = "#132B43",
                       high = "#56B1F7",
                       space = "Lab",
                       guide = FALSE,
                       aesthetics = "fill") #Esto es para que el color sea un gradiente dependiente de los valores de fill
all_tables_kegg_tiles 

library("pathview")
library("ggplot2")
library("ggsci")
library("gridExtra")
library("dplyr")

setwd("/home/user/Escritorio/mario/PGAP_NCBI/Hven4/")

kegg_table_Hven4 <- read.table("/home/user/Escritorio/mario/PGAP_NCBI/Hven4/Hven4_kegg-metabolism_modules.txt", sep = "\t", header = T)
str(kegg_table_Hven4)
kegg_Hven4 <- ggplot(kegg_table_Hven4, 
                     aes(x=module_category, 
                         y=module_subcategory,
                         #size=module_completeness,
                         color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("KEGG metabolism modules", subtitle = "Hven4") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

kegg_Hven4_chido <- kegg_Hven4 + theme(legend.position="bottom", legend.box = "horizontal")
kegg_Hven4_chido
ggsave("kegg_Hven4.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)

###
kegg_table_Hven7 <- read.table("/home/user/Escritorio/mario/PGAP_NCBI/Hven7/kegg-metabolism_modules.txt", sep = "\t", header = T)
str(kegg_table_Hven7)
kegg_table7 <- ggplot(kegg_table_Hven7, 
                      aes(x=module_category, 
                          y=module_subcategory,
                          #size=module_completeness,
                          color=module_completeness)) +
  geom_point(size=3) +
  ggtitle("KEGG metabolism modules", subtitle = "Hven7") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

kegg_Hven7_chido <- kegg_table7 + theme(legend.position="bottom", legend.box = "horizontal")
kegg_Hven7_chido
ggsave("kegg_Hven7.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)


###
###
kegg_table_Hven9 <- read.table("/home/user/Escritorio/mario/PGAP_NCBI/Hven9/kegg-metabolism_modules.txt", sep = "\t", header = T)
str(kegg_table_Hven9)
kegg_table9 <- ggplot(kegg_table_Hven9, 
                      aes(x=module_category, 
                          y=module_subcategory,
                          #size=module_completeness,
                          color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("KEGG metabolism modules", subtitle = "Hven9") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

kegg_Hven9_chido <- kegg_table9 + theme(legend.position="bottom", legend.box = "horizontal")
kegg_Hven9_chido
ggsave("kegg_Hven9.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)

#############
kegg_table_max1 <- read.table("/home/user/Escritorio/mario/PGAP_NCBI/MAGs/Halomonas_maxbin.001/kegg-metabolism_modules.txt", sep = "\t", header = T)
str(kegg_table_max1)
kegg_max1 <- ggplot(kegg_table_max1, 
                    aes(x=module_category, 
                        y=module_subcategory,
                        #size=module_completeness,
                        color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("KEGG metabolism modules", subtitle = "maxbin1") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

kegg_max1_chido <- kegg_max1 + theme(legend.position="bottom", legend.box = "horizontal")
kegg_max1_chido
ggsave("kegg_max1.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)


kegg_table_Hjan13 <- read.table("/home/user/Escritorio/mario/PGAP_NCBI/Hven4/kegg-metabolism_modules_Hjan13.txt", sep = "\t", header = T)
str(kegg_table_Hjan13)
kegg_Hjan13 <- ggplot(kegg_table_Hjan13, 
                      aes(x=module_category, 
                          y=module_subcategory,
                          #size=module_completeness,
                          color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("KEGG metabolism modules", subtitle = "Hjan13") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

kegg_Hjan13_chido <- kegg_Hjan13 + theme(legend.position="bottom", legend.box = "horizontal")
kegg_Hjan13_chido
ggsave("kegg_Hjan13.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)

kegg_table_Hjan14 <- read.table("/home/user/Escritorio/mario/PGAP_NCBI/Hven4/kegg-metabolism_modules.txt", sep = "\t", header = T)
str(kegg_table_Hjan14)
kegg_Hjan14 <- ggplot(kegg_table_Hjan14, 
                      aes(x=module_category, 
                          y=module_subcategory,
                          #size=module_completeness,
                          color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("KEGG metabolism modules", subtitle = "Hjan14") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

kegg_Hjan14_chido <- kegg_Hjan14 + theme(legend.position="bottom", legend.box = "horizontal")
kegg_Hjan14_chido
ggsave("kegg_Hjan14.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)




###

grid.arrange(kegg_Hven9_chido, kegg_Hven7_chido,kegg_Hven4_chido)

#
#kegg_kofam_table_Hven4 <- read.table("kegg-metabolism_kofam_hits.txt", sep = "\t", header = T)
#str(kegg_kofam_table_Hven4)
#kofam_table <- ggplot(kegg_kofam_table_Hven4, aes(x=ko)) +
#  geom_bar() + 
#  theme(axis.text.x = element_text(angle = 90))+
#  scale_fill_igv()
#kofam_table
#ggsave("kegg_table.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)



#################################3
#Vamos a filtrar con dplyr

Lipids_Hven7 <- filter(kegg_table_Hven7, module_category == "Lipid metabolism")
xeno_Hven7 <- filter(kegg_table_Hven7, module_category == "Xenobiotics biodegradation")
xeno_Hven4 <- filter(kegg_table_Hven4, module_category == "Xenobiotics biodegradation")
xeno_Hven9 <- filter(kegg_table_Hven9, module_category == "Xenobiotics biodegradation")
xeno_Hven9 <- filter(kegg_table_Hven9, module_category == "Xenobiotics biodegradation")

xeno_maxbin <- filter(kegg_table_max1, module_category == "Xenobiotics biodegradation")


#
polyketides_Hven9 <- filter(kegg_table_Hven9, module_category == "Biosynthesis of other secondary metabolites")
polyketides_Hven7 <- filter(kegg_table_Hven7, module_category == "Biosynthesis of other secondary metabolites")
polyketides_Hven7_plot <- ggplot(polyketides_Hven7, 
                                 aes(x=module_category, 
                                     y=module_subcategory,
                                     #size=module_completeness,
                                     color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("poliketydes", subtitle = "Hven7") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()
polyketides_Hven7_plot
ggsave("polyke_hven7.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)

xeno_Hven7_plot <- polyketides_Hven7_plot + theme(legend.position="bottom", legend.box = "horizontal")
xeno_Hven7_plot
#


xeno_Hven7_plot <- ggplot(xeno_Hven7,
                          aes(x=module_category, 
                              y=module_subcategory,
                              #size=module_completeness,
                              color=module_completeness)) +
  geom_point(size=3) + 
  ggtitle("Xenobiotics", subtitle = "Hven7") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

xeno_Hven7_plot <- xeno_Hven7_plot + theme(legend.position="bottom", legend.box = "horizontal")
xeno_Hven7_plot
ggsave("xeno_hven7.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)
##############33
##################
xeno_maxbin1_plot <- ggplot(xeno_maxbin,
                            aes(x=module_name, 
                                y=module_completeness,
                                fill=as.factor(module_completeness))) +
  geom_col() + 
  ggtitle("Xenobiotics", subtitle = "Maxbin1") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()+ 
  legend

xeno_maxbin1_plot <- xeno_maxbin1_plot + theme(legend.position="bottom", legend.box = "horizontal")
xeno_maxbin1_plot
ggsave("xeno_maxbin1.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)




####################################33


modulos_completos_maxbin1 <- xeno_maxbin <- filter(kegg_table_max1, module_is_complete == "True")

modulos_completos_maxbin1_a <- ggplot(modulos_completos_maxbin1, 
                                      aes(x=module_category, 
                                          y=module_subcategory,
                                          #size=module_completeness,
                                          color=module_completeness)) +
  geom_bar() + 
  ggtitle("KEGG metabolism modules", subtitle = "MAG Maxbin1") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

modulos_completos_maxbin1_a <- modulos_completos_maxbin1_a + theme(legend.position="bottom", legend.box = "horizontal")
modulos_completos_maxbin1_a

xeno_Hven7


####################
#PARA VER TODOS EN UNA SOLA GRAFICA
#UNIR DATA FRAMES
#Con tidyverse
#install.packages("tidyverse")
#library("tidyverse")
#tablas_unidas_Hjan13_y_Hjan13 <- tablas_unidas_Hjan13_y_Hjan13 + theme(legend.position="bottom", legend.box = "horizontal")

#str(kegg_table_Hjan13)
#str(kegg_table_Hjan13)
#tablas_unidas_Hjan13_y_Hjan13 <- ggplot(tablas_unidas, 
 #                                       aes(x=genome_name, 
  #                                          y=module_subcategory,
   #                                         #size=module_completeness,
    #                                        color=module_completeness)) +
#  geom_point(size=3) + 
 # ggtitle("KEGG metabolism modules", subtitle = "UNidas") +
  #guides(fill=guide_legend(title=" "))+
#  theme_bw() +
 # theme(axis.text.x = element_text(angle = 90))+
  #labs(ylab="module", xlab="category") +
  #scale_fill_igv()

#tablas_unidas_Hjan13_y_Hjan13


#tablas_unidas_Hven4_Hven7_Hven9 <- bind_rows(kegg_table_Hven4, kegg_table_Hven7, kegg_table_Hven9)
#str(tablas_unidas_Hven4_Hven7_Hven9)
#View(tablas_unidas_Hven4_Hven7_Hven9)
#tablas_unidas_Hven4_Hven7_Hven9 <- ggplot(tablas_unidas_Hven4_Hven7_Hven9, 
 #                                         aes(x=genome_name, 
  #                                            y=module_subcategory,
   #                                           #size=module_completeness,
    #                                          color=module_completeness)) +
  #geom_point(size=3) + 
  #ggtitle("KEGG metabolism modules", subtitle = "UNidas") +
  #guides(fill=guide_legend(title=" "))+
  #theme_bw() +
  #theme(axis.text.x = element_text(angle = 90))+
  #labs(ylab="module", xlab="category") +
  #scale_fill_igv()


#tablas_unidas_Hven4_Hven7_Hven9
###################################################################
# UNI TODAS LAS TABLAS (GENOMAS, MAGs y METAGENOMA)
all_tables <- bind_rows(kegg_table_Hven4, kegg_table_Hven7, kegg_table_Hven9, 
                        kegg_table_Hjan13, kegg_table_Hjan14, kegg_table_bins_metagenome, 
                        kegg_table_Halomonas_maxbin_001, kegg_table_Halomonas_metabat2_2,
                        kegg_table_Lacimicrobium_maxbin.004, kegg_table_Lacimicrobium_metabat2.30,
                        kegg_table_nitriluruptor_maxbin_009, kegg_table_Synechococcus_metabat2.14)
str(all_tables)
#solo metagenomas
kegg_table_bins_metagenome_tiles <- ggplot(kegg_table_bins_metagenome, 
                                           aes(x=genome_name, 
                                               y=module_subcategory,
                                               fill=module_completeness)) +
  geom_tile(color="white",
            lwd = 0.3,
            linetype = 1) + 
  ggtitle("KEGG metabolism modules", subtitle = "All samples") +
  guides(fill=guide_legend(title="Completeness (fraction)"), guide_colourbar(barwidth = 0.5,
                                                                             barheight = 20))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  xlab("Strain") +
  ylab("Module subcategory") +
  scale_fill_gradient( low = "#132B43",
                       high = "#56B1F7",
                       space = "Lab",
                       guide = FALSE,
                       aesthetics = "fill") #Esto es para que el color sea un gradiente dependiente de los valores de fill

kegg_table_bins_metagenome_tiles


#$%#%"
# PLOT TODAS LAS TABLAS
#Primero cambiaremos el orden de los valores de los nivels para que los ordene com yo quiero y no en orden alfabetico

all_tables$genome_name <- factor(all_tables$genome_name, 
                                 levels = c("Hven4", "Hven7", "Hven9", "Hjan13", "Hjan14", "Metagenome",
                                            "Halomonas_maxbin_001", "Halomonas_metabat2.2",
                                            "Nitriliruptor_maxbin.009", "Lacimicrobium_maxbin.004", 
                                            "Lacimicrobium_metabat2.30", "Synechococcus_metabat2.14"))

all_tables_kegg_tiles <- ggplot(all_tables,
                                aes(x=genome_name, 
                                    y=module_subcategory,
                                    fill=module_completeness)) +
  geom_tile(color="white", #esto hace el heatmap y pone divisiones blancas
            lwd = 0.3,
            linetype = 1) + 
  ggtitle("KEGG metabolism modules", subtitle = "All samples") +
  guides(fill=guide_legend(title="Metabolic module completeness (fraction)"), guide_colourbar(barwidth = 0.5,
                                                                             barheight = 20))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))+
  xlab(" ") +
  ylab("Module subcategory") +
  scale_fill_gradient( low = "#132B43", # #FFCC00 Amarillo http://derekogle.com/NCGraphing/resources/colors
                       high = "#56B1F7", # #FF3300 Rojo
                       space = "Lab",
                       guide = FALSE,
                       aesthetics = "fill") #Esto es para que el color sea un gradiente dependiente de los valores de fill

all_tables_kegg_tiles + coord_fixed()






ggsave("GENOMES_VS_METAGENOME_VS_MAGS2.pdf", path = ".", width = 10, height = 10, device = "pdf", dpi = 800)


