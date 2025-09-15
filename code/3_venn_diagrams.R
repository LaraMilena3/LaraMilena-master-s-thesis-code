##############################################
# Standard Tasks at the beginning
##############################################

rm(list=ls()) #remove ALL objects 
Sys.setenv(LANG = "en") #Let's keep stuff in English
Sys.setlocale("LC_ALL","English")
options(scipen=999) #avoid scientific annotation

##############################################
# Version Number
##############################################
vers<-1

##############################################
# Set constants
##############################################

debbuging = TRUE

############
# Set the current workdirectory, i.e. where all files are kept
############

#The folder were the data comes from
inwd <- "r_data"
#The folder where tables go to
outwd <- "r_output"
#The folder where the pdfs go
outwd2 = "graphs/"

#############################
#packages
####################

library(ggvenn)
library("eulerr")
require(VennDiagram)

#######################################
#main body
###############################

df = read.csv(paste0(inwd, "/data_leth_subl.csv"))

#receive all chemicals that have been tested at least once
chemicals <- data.frame(CAS = df$CAS, pesticide = df$pesticide, pharmaceutical = df$pharmaceutical, industrial = df$industrial)
chemicals <- unique(chemicals)

#make a venn diagram about how many chemicals have been tested in all the groups and their overlap
plot1 <- chemicals %>%
  ggplot() +
  geom_venn(aes(A = pesticide, B = pharmaceutical, C = industrial#, D = reach
  ))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        axis.ticks=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


#print plot1 as a pdf into the outwd2 folder if debbuging == FALSE
if (!debbuging){
  ggsave(file= paste0(outwd2, "venn_diagram_chemicals.pdf"),
         width = 7,
         plot = plot1,
         height = 7)
}

#make a venn diagram about how many tests have been conducted in all the groups and their overlap
plot2 <- df %>%
  ggplot() +
  geom_venn(aes(A = pesticide, B = pharmaceutical, C = industrial
  ),text_size = 4.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        axis.ticks=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

#print plot2 as a pdf into the outwd2 folder if debbuging == FALSE
if (!debbuging){
  ggsave(paste0(outwd2, "venn_diagram_tests.pdf"),
         plot = plot2,
         width = 7,
         height = 7)
}


