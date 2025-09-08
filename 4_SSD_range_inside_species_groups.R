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

#if debbuging = FALSE, the pdf will be saved to outwd2
debbuging = TRUE

#the colours for the boxplots
colours <- c("#E2F6B6", "#90C819", 
             "#F5CCD7", "#A82445", 
             "#ADEBFF", "#0099CC", 
             "#FF9933", "#C16200")

#the colours for the legend
legend_colours1 <- c("#E2F6B6", "#F5CCD7","#ADEBFF", "#FF9933")
legend_colours2 <- c("#90C819", "#A82445", "#0099CC", "#C16200")


############
# Set the current workdirectory, i.e. where all files are kept
############

#the folder were the data comes from
inwd <- "r_data"
#the folder where tables go to
outwd <- "r_output"
#the folder where the pdfs go
outwd2 = "graphs/"

####################
# packages, imports, functions
################


library(tidyverse)
library(dplyr)



# Function to add legend with the colours
draw_legend <- function(x, y, colours) {
  for (i in seq_along(colours)) {
    rect(x + (i - 1) * 0.35,
         y - 0.3 / 2,
         x + (i - 1) * 0.35 + 0.3,
         y + 0.3 / 2,
         col = colours[i], 
         border = NA)
  }
}

################
# Imports
################

df <- read.csv(paste0(outwd, "/effect_of_chemicals_number.csv")) #check that it is a data frame

#########################
#main body
#################

#Only keep the groups that we defined as acute and chronic
groups <- c("algae.EC10.leth", "crustaceans.EC10.subl.long", "fish.EC10.subl.long", "molluscs.EC10.subl.long", "algae.EC50.leth", "crustaceans.EC50.leth", "fish.EC50.leth", "molluscs.EC50.leth")
df <- df[df$group %in% groups,]
#Split the group name into it's components
df[ , 9:12] <- str_split_fixed(df$group, "\\.", 4)

#Calculate the range between HC95 and HC05 (and the log10 of it)
df$range <- df$HC95/df$HC05
df$range_log <- log10(df$range)
#define the endpoint as a factor
df$V10 <- factor(df$V10, levels = c("EC50", "EC10"))

#When printing the pdf save it to outwd2
if (!debbuging){
  pdf(file= paste0(outwd2, "inside_group_difference.pdf"),
      width=9,
      height=6
  )
}
#Define the parameters for the plotting
par(mar = c(1.5, 5, 0.1, 0.1), 
    mgp = c(2.5, 0.5, 0),
    las = 1,
    cex.lab = 1.3,
    cex.main = 1.8,
    cex = 1.5,
    tcl = -0.3)

boxplot(df$range_log~df$V10 + df$V9,
        at = c(1.1,1.9, 3.6,4.4, 6.1,6.9, 8.6,9.4),
        xaxt = 'n',
        yaxt = "n",
        col = colours,

        ylab = expression(frac(HC95,HC05)),
        xlab = "")

#Add axes to the bottom and the left
axis(1, at = c(1.5, 4, 6.5, 9), labels = c("algae", "crustacean", "fish", "molluscs"))
axis(2, at = axTicks(2), label = paste0("1e+", axTicks(2)))

#Define the starting point of the legend
legend_x <- 7.4
legend_y <- 8

# Draw the legend for acute data
draw_legend(legend_x, legend_y, legend_colours1)
text(legend_x + 1.6, legend_y+0.05, "acute", adj = 0)

# Draw the legend for chronic data
draw_legend(legend_x, legend_y - 0.5, legend_colours2)
text(legend_x + 1.6, legend_y - 0.45, "chronic", adj = 0)

#When printing the pdf, the saving device has to be closed
if (!debbuging){
  dev.off()
}


