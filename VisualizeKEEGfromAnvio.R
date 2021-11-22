library("pathview")
library("ggplot2")
library("ggsci")
library("gridExtra")
library("dplyr")

setwd("/home/user/Escritorio/mario/PGAP_NCBI/Hven4/")

kegg_table_Hven4 <- read.table("/home/user/Escritorio/mario/PGAP_NCBI/Hven4/Hven4_kegg-metabolism_modules.txt", sep = "\t", header = T) #en mi primer intento hubo  un error
str(kegg_table_Hven4)
kegg_Hven4 <- ggplot(kegg_table_Hven4, 
                     aes(x=module_category, 
                         y=module_subcategory,
                         size=module_completeness,
                         color=module_completeness)) +
  geom_point() + 
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
                         size=module_completeness,
                         color=module_completeness)) +
  geom_point() + 
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
                          size=module_completeness,
                          color=module_completeness)) +
  geom_point() + 
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
                          size=module_completeness,
                          color=module_completeness)) +
  geom_point() + 
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
                         size=module_completeness,
                         color=module_completeness)) +
  geom_point() + 
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
                          size=module_completeness,
                          color=module_completeness)) +
  geom_point() + 
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

#
polyketides_Hven9 <- filter(kegg_table_Hven9, module_category == "Biosynthesis of other secondary metabolites")
polyketides_Hven7 <- filter(kegg_table_Hven7, module_category == "Biosynthesis of other secondary metabolites")
polyketides_Hven7_plot <- ggplot(polyketides_Hven7, 
                      aes(x=module_category, 
                          y=module_subcategory,
                          #size=module_completeness,
                          color=module_completeness)) +
  geom_point() + 
  ggtitle("poliketydes", subtitle = "Hven7") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

xeno_Hven7_plot <- polyketides_Hven7_plot + theme(legend.position="bottom", legend.box = "horizontal")
polyketides_Hven7_plot
#


xeno_Hven7_plot <- ggplot(xeno_Hven7,
                          aes(x=module_category, 
                          y=module_subcategory,
                          #size=module_completeness,
                          color=module_completeness)) +
  geom_point() + 
  ggtitle("Xenobiotics", subtitle = "Hven7") +
  guides(fill=guide_legend(title=" "))+
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 90))+
  labs(ylab="module", xlab="category") +
  scale_fill_igv()

xeno_Hven7_plot <- xeno_Hven7_plot + theme(legend.position="bottom", legend.box = "horizontal")
xeno_Hven7_plot


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



