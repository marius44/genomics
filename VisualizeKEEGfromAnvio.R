library("ggplot2")
library("ggsci")
library("gridExtra")
library("dplyr")
library("extrafont")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("pathview")

setwd("/home/user/Escritorio/mario/genomes_Isabel_WGS/pange_halomonas/anvio_metabolism/")
font_import(paths = "/home/user/miniconda3/pkgs/cctyper-1.2.3-py39_0/info/recipe/data/",prompt = F)

kegg_table_Halomonas_venusta_4743 <- read.table("Halomonas_venusta_DSM_4743_modules.txt", sep = "\t", header = T)
kegg_table_Halomonas_venusta_4743[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Halomonas_venusta_4743 <- mutate(kegg_table_Halomonas_venusta_4743, genome_name="H. venusta DSM 4743")
#View(kegg_table_Halomonas_venusta_4743)

kegg_table_Hven4 <- read.table("Halomonas_venusta_Hven4_modules.txt", sep = "\t", header = T)
kegg_table_Hven4[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Hven4 <- mutate(kegg_table_Hven4, genome_name="→ H. venusta Hven4")
str(kegg_table_Hven4)

kegg_table_Hven7 <- read.table("Halomonas_venusta_Hven7_modules.txt", sep = "\t", header = T)
kegg_table_Hven7[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Hven7 <- mutate(kegg_table_Hven7, genome_name="→ H. venusta Hven7")
str(kegg_table_Hven7)

kegg_table_Hven9 <- read.table("Halomonas_venusta_Hven9_modules.txt", sep = "\t", header = T)
kegg_table_Hven9[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Hven9 <- mutate(kegg_table_Hven9, genome_name="→ H. venusta Hven9")
str(kegg_table_Hven9)

kegg_table_Hjan13 <- read.table("Halomonas_janggokensis_Hjan13_modules.txt", sep = "\t", header = T)
kegg_table_Hjan13[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Hjan13 <- mutate(kegg_table_Hjan13, genome_name="→ H. janggokensis Hjan13")
str(kegg_table_Hjan13)

kegg_table_Hjan14 <- read.table("Halomonas_janggokensis_Hjan14_modules.txt", sep = "\t", header = T)
kegg_table_Hjan14[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Hjan14 <- mutate(kegg_table_Hjan14, genome_name="→ H. janggokensis Hjan14")
str(kegg_table_Hjan14)


kegg_table_Halomonas_aestuarii <- read.table("Halomonas_aestuarii_modules.txt", sep = "\t", header = T)
kegg_table_Halomonas_aestuarii[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Halomonas_aestuarii <- mutate(kegg_table_Halomonas_aestuarii, genome_name="H. aestuarii Hb3")
select(kegg_table_Halomonas_aestuarii, -starts_with("warnings"))

str(kegg_table_Halomonas_aestuarii)


kegg_table_Halomonas_alkaliphila <- read.table("Halomonas_alkaliphila_modules.txt", sep = "\t", header = T)
kegg_table_Halomonas_alkaliphila[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Halomonas_alkaliphila <- mutate(kegg_table_Halomonas_alkaliphila, genome_name="H. alkaliphila X3")
str(kegg_table_Halomonas_alkaliphila)

kegg_table_Halomonas_axialensis <- read.table("Halomonas_axialensis_modules.txt", sep = "\t", header = T)
kegg_table_Halomonas_axialensis[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Halomonas_axialensis <- mutate(kegg_table_Halomonas_axialensis, genome_name="H. axialensis Althf1")
str(kegg_table_Halomonas_axialensis)

kegg_table_Halomonas_beimenensis2 <- read.table("Halomonas_beimenensis_modules.txt", sep = "\t", header = T)
kegg_table_Halomonas_beimenensis2[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Halomonas_beimenensis2 <- mutate(kegg_table_Halomonas_beimenensis2, genome_name="H. beimenensis NTU-11")
str(kegg_table_Halomonas_beimenensis2)

kegg_table_Halomonas_beimenensis <- read.table("Halomonas_beimenensis_NTU-111_modules.txt", sep = "\t", header = T)
kegg_table_Halomonas_beimenensis[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Halomonas_beimenensis <- mutate(kegg_table_Halomonas_beimenensis, genome_name="H. beimenensis")
str(kegg_table_Halomonas_beimenensis)

kegg_table_Halomonas_campaniensis <- read.table("Halomonas_campaniensis_modules.txt", sep = "\t", header = T)
kegg_table_Halomonas_campaniensis[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Halomonas_campaniensis <- mutate(kegg_table_Halomonas_campaniensis, genome_name="H. campaniensis LS21")
str(kegg_table_Halomonas_campaniensis)


kegg_table_Halomonas_campisalis <- read.table("Halomonas_campisalis_modules.txt", sep = "\t", header = T)
kegg_table_Halomonas_campisalis[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Halomonas_campisalis <- mutate(kegg_table_Halomonas_campisalis, genome_name="H. campisalis SS10-MC5")
str(kegg_table_Halomonas_campisalis)

kegg_table_Halomonas_chromatireducens <- read.table("Halomonas_chromatireducens_modules.txt", sep = "\t", header = T)
kegg_table_Halomonas_chromatireducens[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Halomonas_chromatireducens <- mutate(kegg_table_Halomonas_chromatireducens, genome_name="H. chromatireducens AGD 8-3")
str(kegg_table_Halomonas_chromatireducens)

kegg_table_Halomonas_elongata <- read.table("Halomonas_elongata_modules.txt", sep = "\t", header = T)
kegg_table_Halomonas_elongata[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Halomonas_elongata <- mutate(kegg_table_Halomonas_elongata, genome_name="H. elongata DSM 2581")
str(kegg_table_Halomonas_elongata)

kegg_table_Halomonas_huangheensis <- read.table("Halomonas_huangheensis_modules.txt", sep = "\t", header = T)
kegg_table_Halomonas_huangheensis[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Halomonas_huangheensis <- mutate(kegg_table_Halomonas_huangheensis, genome_name="H. huangheensis BJGMM-B45")
str(kegg_table_Halomonas_huangheensis)

kegg_table_Halomonas_hydrothermalis <- read.table("Halomonas_hydrothermalis_modules.txt", sep = "\t", header = T)
kegg_table_Halomonas_hydrothermalis[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Halomonas_hydrothermalis <- mutate(kegg_table_Halomonas_hydrothermalis, genome_name="H. hydrothermalis Y2")
str(kegg_table_Halomonas_hydrothermalis)

kegg_table_Halomonas_meridiana <- read.table("Halomonas_meridiana_modules.txt", sep = "\t", header = T)
kegg_table_Halomonas_meridiana[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Halomonas_meridiana <- mutate(kegg_table_Halomonas_meridiana, genome_name="H. meridiana SCSIO 43005")
str(kegg_table_Halomonas_meridiana)

kegg_table_Halomonas_olivaria <- read.table("Halomonas_olivaria_modules.txt", sep = "\t", header = T)
kegg_table_Halomonas_olivaria[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Halomonas_olivaria <- mutate(kegg_table_Halomonas_olivaria, genome_name="H. olivaria TYRC17")
str(kegg_table_Halomonas_olivaria)

kegg_table_Halomonas_piezotolerans <- read.table("Halomonas_piezotolerans_modules.txt", sep = "\t", header = T)
kegg_table_Halomonas_piezotolerans[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Halomonas_piezotolerans <- mutate(kegg_table_Halomonas_piezotolerans, genome_name="H. piezotolerans NBT06E8")
str(kegg_table_Halomonas_piezotolerans)

kegg_table_Halomonas_sulfidaeris <- read.table("Halomonas_sulfidaeris_modules.txt", sep = "\t", header = T)
kegg_table_Halomonas_sulfidaeris[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Halomonas_sulfidaeris <- mutate(kegg_table_Halomonas_sulfidaeris, genome_name="H. sulfidaeris ATCC BAA-803")
str(kegg_table_Halomonas_sulfidaeris)

kegg_table_Halomonas_sulfidivorans <- read.table("Halomonas_sulfidivorans_modules.txt", sep = "\t", header = T)
kegg_table_Halomonas_sulfidivorans[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Halomonas_sulfidivorans <- mutate(kegg_table_Halomonas_sulfidivorans, genome_name="H. sulfidivorans MCCC 1A13718")
str(kegg_table_Halomonas_sulfidivorans)

kegg_table_Halomonas_sulfidoxydans <- read.table("Halomonas_sulfidoxydans3_modules.txt", sep = "\t", header = T)
kegg_table_Halomonas_sulfidoxydans[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Halomonas_sulfidoxydans <- mutate(kegg_table_Halomonas_sulfidoxydans, genome_name="H. sulfidoxydans MCCC 1A11059")
str(kegg_table_Halomonas_sulfidoxydans)

kegg_table_Halomonas_titanicae <- read.table("Halomonas_titanicae_modules.txt", sep = "\t", header = T)
kegg_table_Halomonas_titanicae[2] <- NULL #ELIMINA LA COLUMNA 2 que no tiene nombre de metagenomas y tiene un punto
kegg_table_Halomonas_titanicae <- mutate(kegg_table_Halomonas_titanicae, genome_name="H. titanicae SOB56")
str(kegg_table_Halomonas_titanicae)


# UNI TODAS LAS TABLAS (GENOMAS, MAGs y METAGENOMA)
all_tables <- bind_rows(kegg_table_Halomonas_venusta_4743, kegg_table_Hven4, kegg_table_Hven7, kegg_table_Hven9, 
                        kegg_table_Hjan13, kegg_table_Hjan14, kegg_table_Halomonas_aestuarii,
                        kegg_table_Halomonas_alkaliphila, kegg_table_Halomonas_axialensis,
                        kegg_table_Halomonas_beimenensis2, kegg_table_Halomonas_beimenensis,
                        kegg_table_Halomonas_campaniensis, kegg_table_Halomonas_campisalis,
                        kegg_table_Halomonas_chromatireducens, kegg_table_Halomonas_elongata,
                        kegg_table_Halomonas_huangheensis, kegg_table_Halomonas_hydrothermalis,
                        kegg_table_Halomonas_meridiana, kegg_table_Halomonas_olivaria,
                        kegg_table_Halomonas_piezotolerans, kegg_table_Halomonas_sulfidaeris,
                        kegg_table_Halomonas_sulfidivorans, kegg_table_Halomonas_sulfidoxydans,
                        kegg_table_Halomonas_titanicae)



str(all_tables)
#View(all_tables)

#######


#$%#%"
# PLOT TODAS LAS TABLAS
#Primero cambiaremos el orden de los valores de los nivels para que los ordene com yo quiero y no en orden alfabetico

all_tables$genome_name <- factor(all_tables$genome_name, 
                                 levels = c("H. venusta DSM 4743", "→ H. venusta Hven4", "→ H. venusta Hven7", 
                                            "→ H. venusta Hven9", "→ H. janggokensis Hjan13", 
                                            "→ H. janggokensis Hjan14", "H. aestuarii Hb3",
                                            "H. alkaliphila X3", "H. axialensis Althf1",
                                            "H. beimenensis NTU-11", "H. beimenensis", 
                                            "H. campaniensis LS21", "H. campisalis SS10-MC5",
                                            "H. chromatireducens AGD 8-3", "H. elongata DSM 2581",
                                            "H. huangheensis BJGMM-B45", "H. hydrothermalis Y2",
                                            "H. meridiana SCSIO 43005", "H. olivaria TYRC17",
                                            "H. piezotolerans NBT06E8", "H. sulfidaeris ATCC BAA-803",
                                            "H. sulfidivorans MCCC 1A13718", "H. sulfidoxydans MCCC 1A11059",
                                            "H. titanicae SOB56"))

# Unir todas las muestras
all_tables_kegg_tiles <- ggplot(all_tables,
                                aes(x=genome_name, 
                                    y=module_subcategory,
                                    fill=module_completeness)) +
  geom_tile(color="white", #esto hace el heatmap y pone divisiones blancas
            lwd = 0.3,
            linetype = 1) + 
  ggtitle("KEGG metabolism modules", subtitle = "All strains") +
  guides(fill=guide_legend(title="Metabolic \nmodule \ncompleteness \n(fraction)"), 
         guide_colourbar(barwidth = 0.5,barheight = 20))+
  theme_classic() +
  theme(axis.text.x = element_text(face="plain", angle = 45, hjust = 1, size = 8),axis.text.y = element_text(size = 11),)+

   #     axis.text.y = element_text(size = 9), 
    #    axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
     #   legend.position = "none") +
  xlab(" ") +
  ylab("Module subcategory") +
  scale_fill_gradient( low = "#33FFFF", # #FFCC00 Amarillo http://derekogle.com/NCGraphing/resources/colors
                       high = "#6600FF", # #FF3300 Rojo
                       space = "Lab",
                       guide = FALSE,
                       aesthetics = "fill") #Esto es para que el color sea un gradiente dependiente de los valores de fill

all_tables_kegg_tiles + coord_fixed()





ggsave("GENOMES_all4.png", path = ".", width = 10, height = 10, device = "png", dpi = 1000)
ggsave("GENOMES_all4.pdf", path = ".", width = 10, height = 10, device = "pdf", dpi = 800)


######
#quiero ver solo los de capacidad metabolica
metabolic_capacity <- all_tables %>% 
  filter(module_subcategory  == "Metabolic capacity")
#View(metabolic_capacity)

metabolic_capacity_plot <- ggplot(metabolic_capacity,
                                  aes(x=genome_name, 
                                      y=module_name,
                                      fill=module_completeness)) +
  geom_tile(color="white", #esto hace el heatmap y pone divisiones blancas
            lwd = 0.3,
            linetype = 1) + 
  ggtitle("Metabolic capacities", subtitle = "All samples") +
  guides(fill=guide_legend(title="Metabolic \nmodule \ncompleteness \n(fraction)"), guide_colourbar(barwidth = 0.5,
                                                                                              barheight = 20))+
  theme_classic() +
  theme(axis.text.x = element_text(face="plain", angle = 45, hjust = 1, size = 8),axis.text.y = element_text(size = 11),)+
  
  #     axis.text.y = element_text(size = 9), 
  #    axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
  #   legend.position = "none") +
  xlab(" ") +
  ylab("Module subcategory") +
  scale_fill_gradient( low = "#33FFFF", # #FFCC00 Amarillo http://derekogle.com/NCGraphing/resources/colors
                       high = "#6600FF", # #FF3300 Rojo
                       space = "Lab",
                       guide = FALSE,
                       aesthetics = "fill") #Esto es para que el color sea un gradiente dependiente de los valores de fill


metabolic_capacity_plot + coord_fixed()

ggsave("Cap_metabolicas_all4.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)
ggsave("Cap_metabolicas_all4.pdf", path = ".", width = 10, height = 10, device = "pdf", dpi = 800)

#####
#fatty acids

fatty_acids_all <- all_tables %>% 
  filter(module_subcategory  == "Fatty acid metabolism")
#View(fatty_acids_all)

fatty_acids_all_plot <- ggplot(fatty_acids_all,
                               aes(x=genome_name, 
                                   y=module_name,
                                   fill=module_completeness)) +
  geom_tile(color="white", #esto hace el heatmap y pone divisiones blancas
            lwd = 0.3,
            linetype = 1) + 
  ggtitle("Fatty acid metabolism", subtitle = "All samples") +
  guides(fill=guide_legend(title="Metabolic \nmodule \ncompleteness \n(fraction)"), guide_colourbar(barwidth = 0.5,
                                                                                              barheight = 20))+
  theme_classic() +
  theme(axis.text.x = element_text(face="plain", angle = 45, hjust = 1, size = 8),axis.text.y = element_text(size = 11),)+
  
  #     axis.text.y = element_text(size = 9), 
  #    axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
  #   legend.position = "none") +
  xlab(" ") +
  ylab("Module subcategory") +
  scale_fill_gradient( low = "#33FFFF", # #FFCC00 Amarillo http://derekogle.com/NCGraphing/resources/colors
                       high = "#6600FF", # #FF3300 Rojo
                       space = "Lab",
                       guide = FALSE,
                       aesthetics = "fill") #Esto es para que el color sea un gradiente dependiente de los valores de fill

fatty_acids_all_plot + coord_fixed()

ggsave("fatty_acids_all4.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)
ggsave("fatty_acids_all4.pdf", path = ".", width = 10, height = 10, device = "pdf", dpi = 800)
#####
#aromaticos
aromatics_all <- all_tables %>% 
  filter(module_subcategory  == "Aromatics degradation")
#View(aromatics_all)

aromatics_all_plot <- ggplot(aromatics_all,
                             aes(x=genome_name, 
                                 y=module_name,
                                 fill=module_completeness)) +
  geom_tile(color="white", #esto hace el heatmap y pone divisiones blancas
            lwd = 0.3,
            linetype = 1) + 
  ggtitle("Aromatics metabolism", subtitle = "All samples") +
  guides(fill=guide_legend(title="Metabolic \nmodule \ncompleteness \n(fraction)"), guide_colourbar(barwidth = 0.5,
                                                                                              barheight = 20))+
  theme_classic() +
  theme(axis.text.x = element_text(face="plain", angle = 45, hjust = 1, size = 8),axis.text.y = element_text(size = 11),)+
  
  #     axis.text.y = element_text(size = 9), 
  #    axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
  #   legend.position = "none") +
  xlab(" ") +
  ylab("Module subcategory") +
  scale_fill_gradient( low = "#33FFFF", # #FFCC00 Amarillo http://derekogle.com/NCGraphing/resources/colors
                       high = "#6600FF", # #FF3300 Rojo
                       space = "Lab",
                       guide = FALSE,
                       aesthetics = "fill") #Esto es para que el color sea un gradiente dependiente de los valores de fill

aromatics_all_plot + coord_fixed()

ggsave("aromatics_all_all4.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)
ggsave("aromatics_all_all4.pdf", path = ".", width = 10, height = 10, device = "pdf", dpi = 800)

### 
#central_metabolism


central_metabolism_all <- all_tables %>% 
  filter(module_subcategory  == "Central carbohydrate metabolism")
#View(central_metabolism_all)

central_metabolism_all_plot <- ggplot(central_metabolism_all,
                                      aes(x=genome_name, 
                                          y=module_name,
                                          fill=module_completeness)) +
  geom_tile(color="white", #esto hace el heatmap y pone divisiones blancas
            lwd = 0.3,
            linetype = 1) + 
  ggtitle("Central carbohydrate metabolism", subtitle = "All samples") +
  guides(fill=guide_legend(title="Metabolic \nmodule \ncompleteness \n(fraction)"), guide_colourbar(barwidth = 0.5,
                                                                                              barheight = 20))+
  theme_classic() +
  theme(axis.text.x = element_text(face="plain", angle = 45, hjust = 1, size = 8),axis.text.y = element_text(size = 11),)+
  
  #     axis.text.y = element_text(size = 9), 
  #    axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
  #   legend.position = "none") +
  xlab(" ") +
  ylab("Module subcategory") +
  scale_fill_gradient( low = "#33FFFF", # #FFCC00 Amarillo http://derekogle.com/NCGraphing/resources/colors
                       high = "#6600FF", # #FF3300 Rojo
                       space = "Lab",
                       guide = FALSE,
                       aesthetics = "fill") #Esto es para que el color sea un gradiente dependiente de los valores de fill

central_metabolism_all_plot + coord_fixed()

ggsave("central_metabolism_all4.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)
ggsave("central_metabolism_all4.pdf", path = ".", width = 10, height = 10, device = "pdf", dpi = 800)

#########
#Lipid

Lipid_all <- all_tables %>% 
  filter(module_subcategory  == "Lipid metabolism")
#View(Lipid_all)

Lipid_all_plot <- ggplot(Lipid_all,
                         aes(x=genome_name, 
                             y=module_name,
                             fill=module_completeness)) +
  geom_tile(color="white", #esto hace el heatmap y pone divisiones blancas
            lwd = 0.3,
            linetype = 1) + 
  ggtitle("Lipid metabolism", subtitle = "All samples") +
  guides(fill=guide_legend(title="Metabolic \nmodule \ncompleteness \n(fraction)"), guide_colourbar(barwidth = 0.5,
                                                                                              barheight = 20))+
  theme_classic() +
  theme(axis.text.x = element_text(face="plain", angle = 45, hjust = 1, size = 8),axis.text.y = element_text(size = 11),)+
  
  #     axis.text.y = element_text(size = 9), 
  #    axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
  #   legend.position = "none") +
  xlab(" ") +
  ylab("Module subcategory") +
  scale_fill_gradient( low = "#33FFFF", # #FFCC00 Amarillo http://derekogle.com/NCGraphing/resources/colors
                       high = "#6600FF", # #FF3300 Rojo
                       space = "Lab",
                       guide = FALSE,
                       aesthetics = "fill") #Esto es para que el color sea un gradiente dependiente de los valores de fill


Lipid_all_plot + coord_fixed()

ggsave("Lipid metabolism_all4.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)
ggsave("Lipid metabolism_all4.pdf", path = ".", width = 10, height = 10, device = "pdf", dpi = 800)

#########
#nitrogen

Nitrogen_all <- all_tables %>% 
  filter(module_subcategory  == "Nitrogen metabolism")
#View(Nitrogen_all)

Nitrogen_all_plot <- ggplot(Nitrogen_all,
                            aes(x=genome_name, 
                                y=module_name,
                                fill=module_completeness)) +
  geom_tile(color="white", #esto hace el heatmap y pone divisiones blancas
            lwd = 0.3,
            linetype = 1) + 
  ggtitle("Nitrogen metabolism", subtitle = "All samples") +
  guides(fill=guide_legend(title="Metabolic \nmodule \ncompleteness \n(fraction)"), guide_colourbar(barwidth = 0.5,
                                                                                              barheight = 20))+
  theme_classic() +
  theme(axis.text.x = element_text(face="plain", angle = 45, hjust = 1, size = 8),axis.text.y = element_text(size = 11),)+
  
  #     axis.text.y = element_text(size = 9), 
  #    axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
  #   legend.position = "none") +
  xlab(" ") +
  ylab("Module subcategory") +
  scale_fill_gradient( low = "#33FFFF", # #FFCC00 Amarillo http://derekogle.com/NCGraphing/resources/colors
                       high = "#6600FF", # #FF3300 Rojo
                       space = "Lab",
                       guide = FALSE,
                       aesthetics = "fill") #Esto es para que el color sea un gradiente dependiente de los valores de fill


Nitrogen_all_plot + coord_fixed()

ggsave("Nitrogen metabolism_all4.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)
ggsave("Nitrogen metabolism_all4.pdf", path = ".", width = 10, height = 10, device = "pdf", dpi = 800)

#########
#Lipopolysaccharide metabolism

Lipopolysaccharide_metabolism_all <- all_tables %>% 
  filter(module_subcategory  == "Lipopolysaccharide metabolism")
#View(Nitrogen_all)

Lipopolysaccharide_metabolism_all_plot <- ggplot(Lipopolysaccharide_metabolism_all,
                                                 aes(x=genome_name, 
                                                     y=module_name,
                                                     fill=module_completeness)) +
  geom_tile(color="white", #esto hace el heatmap y pone divisiones blancas
            lwd = 0.3,
            linetype = 1) + 
  ggtitle("Lipopolysaccharide metabolism", subtitle = "All samples") +
  guides(fill=guide_legend(title="Metabolic \nmodule \ncompleteness \n(fraction)"), guide_colourbar(barwidth = 0.5,
                                                                                              barheight = 20))+
  theme_classic() +
  theme(axis.text.x = element_text(face="plain", angle = 45, hjust = 1, size = 8),axis.text.y = element_text(size = 11),)+
  
  #     axis.text.y = element_text(size = 9), 
  #    axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
  #   legend.position = "none") +
  xlab(" ") +
  ylab("Module subcategory") +
  scale_fill_gradient( low = "#33FFFF", # #FFCC00 Amarillo http://derekogle.com/NCGraphing/resources/colors
                       high = "#6600FF", # #FF3300 Rojo
                       space = "Lab",
                       guide = FALSE,
                       aesthetics = "fill") #Esto es para que el color sea un gradiente dependiente de los valores de fill


Lipopolysaccharide_metabolism_all_plot + coord_fixed()

ggsave("Lipopolysaccharide_metabolism_all4.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)
ggsave("Lipopolysaccharide_metabolism_all4.pdf", path = ".", width = 10, height = 10, device = "pdf", dpi = 800)


##########
#polyketides

Polyketides_metabolism_all <- all_tables %>% 
  filter(module_category  == "Biosynthesis of terpenoids and polyketides")
View(Polyketides_metabolism_all)

Polyketides_metabolism_all_plot <- ggplot(Polyketides_metabolism_all,
                                                 aes(x=genome_name, 
                                                     y=module_name,
                                                     fill=module_completeness)) +
  geom_tile(color="white", #esto hace el heatmap y pone divisiones blancas
            lwd = 0.3,
            linetype = 1) + 
  ggtitle("Biosynthesis of terpenoids and polyketides", subtitle = "All samples") +
  guides(fill=guide_legend(title="Metabolic \nmodule \ncompleteness \n(fraction)"), guide_colourbar(barwidth = 0.5,
                                                                                                    barheight = 20))+
  theme_classic() +
  theme(axis.text.x = element_text(face="plain", angle = 45, hjust = 1, size = 8),axis.text.y = element_text(size = 11),)+
  
  #     axis.text.y = element_text(size = 9), 
  #    axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
  #   legend.position = "none") +
  xlab(" ") +
  ylab("Module subcategory") +
  scale_fill_gradient( low = "#33FFFF", # #FFCC00 Amarillo http://derekogle.com/NCGraphing/resources/colors
                       high = "#6600FF", # #FF3300 Rojo
                       space = "Lab",
                       guide = FALSE,
                       aesthetics = "fill") #Esto es para que el color sea un gradiente dependiente de los valores de fill


Polyketides_metabolism_all_plot + coord_fixed()

ggsave("Polyketides_metabolism_all4.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)
ggsave("Polyketides_metabolism_all4.pdf", path = ".", width = 10, height = 10, device = "pdf", dpi = 800)

##########
#

Terpenoid_metabolism_all <- all_tables %>% 
  filter(module_category  == "Biosynthesis of terpenoids and polyketides")
View(Terpenoid_metabolism_all)

Terpenoid_metabolism_all_plot <- ggplot(Terpenoid_metabolism_all,
                                aes(x=genome_name, 
                                    y=module_name,
                                    fill=module_completeness)) +
  geom_tile(color="white", #esto hace el heatmap y pone divisiones blancas
            lwd = 0.3,
            linetype = 1) + 
  ggtitle("Biosynthesis of terpenoids and polyketides", subtitle = "All samples") +
  guides(fill=guide_legend(title="Metabolic \nmodule \ncompleteness \n(fraction)"), guide_colourbar(barwidth = 0.5,
                                                                                                    barheight = 20))+
  theme_classic() +
  theme(axis.text.x = element_text(face="plain", angle = 45, hjust = 1, size = 8),axis.text.y = element_text(size = 11),)+
  
  #     axis.text.y = element_text(size = 9), 
  #    axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
  #   legend.position = "none") +
  xlab(" ") +
  ylab("Module subcategory") +
  scale_fill_gradient( low = "#33FFFF", # #FFCC00 Amarillo http://derekogle.com/NCGraphing/resources/colors
                       high = "#6600FF", # #FF3300 Rojo
                       space = "Lab",
                       guide = FALSE,
                       aesthetics = "fill") #Esto es para que el color sea un gradiente dependiente de los valores de fill


Terpenoid_metabolism_all_plot + coord_fixed()

ggsave("X_metabolism_all4.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)
ggsave("X_metabolism_all4.pdf", path = ".", width = 10, height = 10, device = "pdf", dpi = 800)

############
#polyamines

Polyamines_metabolism_all <- all_tables %>% 
  filter(module_subcategory  == "Polyamine biosynthesis")
#View(Polyketides_metabolism_all)

Polyamines_metabolism_all_plot <- ggplot(Polyamines_metabolism_all,
                                aes(x=genome_name, 
                                    y=module_name,
                                    fill=module_completeness)) +
  geom_tile(color="white", #esto hace el heatmap y pone divisiones blancas
            lwd = 0.3,
            linetype = 1) + 
  ggtitle("Polyamine biosynthesis", subtitle = "All samples") +
  guides(fill=guide_legend(title="Metabolic \nmodule \ncompleteness \n(fraction)"), guide_colourbar(barwidth = 0.5,
                                                                                                    barheight = 20))+
  theme_classic() +
  theme(axis.text.x = element_text(face="plain", angle = 45, hjust = 1, size = 8),axis.text.y = element_text(size = 11),)+
  
  #     axis.text.y = element_text(size = 9), 
  #    axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
  #   legend.position = "none") +
  xlab(" ") +
  ylab("Module subcategory") +
  scale_fill_gradient( low = "#33FFFF", # #FFCC00 Amarillo http://derekogle.com/NCGraphing/resources/colors
                       high = "#6600FF", # #FF3300 Rojo
                       space = "Lab",
                       guide = FALSE,
                       aesthetics = "fill") #Esto es para que el color sea un gradiente dependiente de los valores de fill


Polyamines_metabolism_all_plot + coord_fixed()

ggsave("Polyamines_metabolism_all4.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)
ggsave("Polyamines_metabolism_all4.pdf", path = ".", width = 10, height = 10, device = "pdf", dpi = 800)

##############
#serine

Serine_metabolism_all <- all_tables %>% 
  filter(module_subcategory  == "Serine and threonine metabolism")
#View(Polyketides_metabolism_all)

Serine_metabolism_all_plot <- ggplot(Serine_metabolism_all,
                                         aes(x=genome_name, 
                                             y=module_name,
                                             fill=module_completeness)) +
  geom_tile(color="white", #esto hace el heatmap y pone divisiones blancas
            lwd = 0.3,
            linetype = 1) + 
  ggtitle("Serine and threonine biosynthesis", subtitle = "All samples") +
  guides(fill=guide_legend(title="Metabolic \nmodule \ncompleteness \n(fraction)"), guide_colourbar(barwidth = 0.5,
                                                                                                    barheight = 20))+
  theme_classic() +
  theme(axis.text.x = element_text(face="plain", angle = 45, hjust = 1, size = 8),axis.text.y = element_text(size = 11),)+
  
  #     axis.text.y = element_text(size = 9), 
  #    axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
  #   legend.position = "none") +
  xlab(" ") +
  ylab("Module subcategory") +
  scale_fill_gradient( low = "#33FFFF", # #FFCC00 Amarillo http://derekogle.com/NCGraphing/resources/colors
                       high = "#6600FF", # #FF3300 Rojo
                       space = "Lab",
                       guide = FALSE,
                       aesthetics = "fill") #Esto es para que el color sea un gradiente dependiente de los valores de fill


Serine_metabolism_all_plot + coord_fixed()

ggsave("Serine_metabolism_all4.png", path = ".", width = 10, height = 10, device = "png", dpi = 800)
ggsave("Serine_metabolism_all4.pdf", path = ".", width = 10, height = 10, device = "pdf", dpi = 800)

