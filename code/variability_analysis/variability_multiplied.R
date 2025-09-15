##############################################
# Standard Tasks at the beginning
##############################################

rm(list=ls()) #remove ALL objects 
Sys.setenv(LANG = "en") #Let's keep stuff in English
Sys.setlocale("LC_ALL","English")
options(scipen=999) #avoid scientific annotation


##############################################
# Set constants
##############################################
#if debbuging = FALSE the pdf will be printed into outwd2
debbuging = TRUE 

#to more easily change what is plotted here are several constants that are defined which will be important later
#If any of the following constants are misspelled the last possible option of those in the # is chosen

#the weighting of sampling that will be conducted
weight = 10 #numeric(for log) / sqrt / none / occurence
#what endpoint wants to be analyzed
effect = "EC50" #EC10 / EC50
#if innerspecies variation should be defined as IQR or max/min
range = "IQR" #IQR / maxmin
#which species should be analyzed
species = "algae" #fish / algae / molluscs / crustaceans
#how many iterations are called
sample_size = 10000 #whatever you want
#to create an density plot a band width can be chosen
band_width = 0.15 #numeric / SJ / nrd / nrd0

parts = FALSE #TRUE/FALSE

############
# Set the current workdirectory, i.e. where all files are kept
############


#The folder were the data comes from
inwd <- "r_data"
#The folder where tables go to
outwd <- "r_output"
#The folder where the pdfs go
outwd2 = "graphs/"

####################
# packages, imports, functions
################

library(RColorBrewer)
library(tidyverse)
library(dplyr)
library(ssdtools)


################
# Imports
################

df_list_raw <- readRDS(paste0(inwd, "/data_list_leth_subl.RData")) #check that it is a data frame
effect_data_raw <- read.csv(paste0(outwd, "/effect_of_chemicals.csv"))





#choose the species group and endpoint
#it first checks to see which species group is wanted
#then it checks if EC10 or EC50 is wanted
#If the species group or endpoint is misspelled it uses EC50 and crustaceans as defaults
if(species == "fish"){
  if(effect == "EC10"){
    effect_data <- effect_data_raw[effect_data_raw$group == "fish.EC10.subl.long",]
    df <- df_list_raw[["fish.EC10.subl.long"]]
  } else{
    effect_data <- effect_data_raw[effect_data_raw$group == "fish.EC50.leth",]
    df <- df_list_raw[["fish.EC50.leth"]]
    effect = "EC50"
  }
} else if(species == "algae"){
  if(effect == "EC10"){
    effect_data <- effect_data_raw[effect_data_raw$group == "algae.EC10.leth",]
    df <- df_list_raw[["algae.EC10.leth"]]
  } else{
    effect_data <- effect_data_raw[effect_data_raw$group == "algae.EC50.leth",]
    df <- df_list_raw[["algae.EC50.leth"]]
    effect = "EC50"
  }
} else if(species == "molluscs"){
  if(effect == "EC10"){
    effect_data <- effect_data_raw[effect_data_raw$group == "molluscs.EC10.subl.long",]
    df <- df_list_raw[["molluscs.EC10.subl.long"]]
  } else{
    effect_data <- effect_data_raw[effect_data_raw$group == "molluscs.EC50.leth",]
    df <- df_list_raw[["molluscs.EC50.leth"]]
    effect = "EC50"
  }
} else{
  species = "crustaceans"
  if(effect == "EC10"){
    effect_data <- effect_data_raw[effect_data_raw$group == "crustaceans.EC10.subl.long",]
    df <- df_list_raw[["crustaceans.EC10.subl.long"]]
  } else{
    effect_data <- effect_data_raw[effect_data_raw$group == "crustaceans.EC50.leth",]
    df <- df_list_raw[["crustaceans.EC50.leth"]]
    effect = "EC50"
  }
} 

#to define the between species variation the range of the SSD is calculated
#to calculate the range the HC95 is divided by the HC05
effect_data$range <- effect_data$HC95/effect_data$HC05



#organisms that have been tested a lot with a specific chemical should be more represented in the evaluation
#for every chemical organism pair the sum of all tests is calculated 
weight_df <- aggregate(count~organism+CAS, data=df, sum)


#remove organism CAS pairs that only have 3 or less tests
#pairs with 4 or more tests have a more stable range between HC95 and HC05
#pairs with 2 or 3 tests have in the median lower ranges than tests with more than 3
weight_df = weight_df[weight_df$count>3,]

#to not only sample the most common pair a weighting function can we used
#the df weight_df with the weight every pair gets is changed by the weight factor "weight"
if(is.numeric(weight)){
  weight_df$count = log(weight_df$count+1, base = weight)
  weight = paste0("log", weight)
} else if(weight == "none"){
  weight_df$count = 1
}else if (weight == "sqrt"){
  weight_df$count = sqrt(weight_df$count)
}else{
  weight = "occurence"
}

#creates a new column with organism name and Cas number
df$merge <- paste(df$organism, df$CAS, sep = " ")  # Adds a space between the two columns

#the variation inside species can either be determined by the Inter quartile range or by max/min
#this calculates the innerspecies variation in either way depending which is has been asked for
if (range == "IQR"){
  temp <- data.frame(organism = tapply(X = df$organism,
                                       INDEX = df$merge,
                                       FUN = function(x) unique(x)),
                     CAS = tapply(X = df$CAS,
                                  INDEX = df$merge,
                                  FUN = function(x) unique(x)),
                     range = tapply(X = df$mgperL,
                                    INDEX = df$merge,
                                    FUN = function(x) quantile(x, prob=c(.75))/quantile(x, prob=c(.25))))
} else{
  temp <- data.frame(organism = tapply(X = df$organism,
                                       INDEX = df$merge,
                                       FUN = function(x) unique(x)),
                     CAS = tapply(X = df$CAS,
                                  INDEX = df$merge,
                                  FUN = function(x) unique(x)),
                     range = tapply(X = df$mgperL,
                                    INDEX = df$merge,
                                    FUN = function(x) max(x)/min(x)))
  #as before if it got misspelled max/min is used
  range = "max/min"
}

#add the weighting to the innerspecies variation
df <- merge(temp, weight_df, by = c("organism", "CAS"))

#chose a random organism chemical pair weighted by the weight calculated before
#how many are chosen depends on the constant sample_size
results <- df[sample(1:nrow(df), size = sample_size, replace = TRUE, prob = df$count), ]
#multiply the innerspecies variation with the between species variation
results$range <- results$range * c(sample(effect_data$range, size = sample_size, replace = TRUE))
#the values are logged 
results$range_log <- log10(results$range)


#to visualize the two variations combined a density graph is created
#to make a density estimate a band width has to be chosen
#this code calculates the density estimate with the provided bandwidth of the constants
if(is.numeric(band_width) | band_width == "SJ" | band_width == "nrd"){
  results_dens = density(results$range_log,
                         bw = band_width)
} else{
  results_dens = density(results$range_log,
                         bw = "nrd0")
}


#this saves the pdf that is going to be created into the folder defined as outwd2
if (!debbuging){
    pdf(file= paste0(outwd2, weight, "_", species, "_", effect,"_", range, "_two_ranges_multilplied.pdf"),
        width = 9,
        height = 6)
}

#setting the graphical parameters, such as margins and font sizes
par(mar = c(3.1, 3, 0.2, 0.1), 
    mgp = c(2, 0.5, 0),
    las = 1,
    mfcol = c(1,1),
    cex.lab = 1.3,
    cex.main = 2,
    cex = 1.5,
    bty = "n",
    tcl = -0.3)


# make an empty plot
plot(results_dens,
     main = "",
     xaxt = "n",
     xlab = "",
     xlim = c(-0.2, 8.5))

#give the plot a base colour
polygon(results_dens, col = rgb(0.78, 0.89, 1, alpha = 0.6))

#add an x axis
axis(side = 1, line = -0.5, at = c(axTicks(1),1), labels = paste0("1e+", c(axTicks(1), 1)))

#############################
#here we define thresholds, calculate their x and y positions and the percentage of iterations left of it
#there are two kinds of thresholds
# at a value of 10
# the value at which 50% and 95% of the iterations have a lower value
##################

#threshold at 10
#calculate how many iterations have a value below 10
ten <- length(results$range_log[results$range_log <= 1])
#calculate the y position
#this will be used to draw it on the graph
position_1 <- approx(results_dens$x, results_dens$y, xout = 1)$y
#calculate the percentage of iterations that have a value below ten
percent_1 <- paste0(signif((ten/sample_size) * 100,digits = 3), "%")

#threshold where 50% have a lower value
#calculate the x position
median_x = median(results$range_log)
#calculate the y position
median_y <- approx(results_dens, xout = median_x)$y

#threshold where 95% have a lower value
#calculate the x position
p95_x <- quantile(results$range_log, 0.95)
#calculate the y position
p95_y <- approx(results_dens, xout = p95_x)$y

#the next three functions colour in the parts of the graph right of the threshold
#with increasingly darkening colours

#the threshold at 10
polygon(c(results_dens$x[results_dens$x >= 1 ], 1),
        c(results_dens$y[results_dens$x >= 1 ], 0),
        col = rgb(1, 0.9, 0.9, alpha = 1),
        border = 1)

#the threshold for 50%
polygon(c(results_dens$x[results_dens$x >= median_x ], median_x),
        c(results_dens$y[results_dens$x >= median_x ], 0),
        col = rgb(1, 0.8, 0.8, alpha = 1),
        border = 1)

#the threshold for 95%
polygon(c(results_dens$x[results_dens$x >= p95_x ], p95_x),
        c(results_dens$y[results_dens$x >= p95_x ], 0),
        col = rgb(1, 0.4, 0.4, alpha = 1),
        border = 1)

#add a label of the percentage of iterations that are below that threshold
text(x = 1, y = position_1, labels = percent_1, pos = 3, xpd = TRUE)
text(x = median_x, y = median_y, labels = "50%", pos = 3, xpd = TRUE)
text(x = p95_x, y = p95_y, labels = "95%", pos = 3, xpd = TRUE)

#this adds the value at which the 50% and 95% thresholds are
axis(side = 1, line = 1, at = c(median_x, p95_x), 
     labels = paste0("1e+", c(signif(median_x, digits = 2), signif(p95_x, digits = 2))), 
     tick = FALSE)
axis(side = 1, line = -0.5, at = c(median_x, p95_x), 
     labels = FALSE,
     tcl = -2)


if (!debbuging){
  dev.off()
}