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
vers<-3


##############################################
# Set constants
##############################################

#if debbuging = FALSE the pdf will be printed into outwd2
debbuging = TRUE 

#If any of the following constants are misspelled the last possible option of those in the # is chosen

#Defining the group that is plotted
effect = "EC50" #EC10 / EC50
species = "fish" #fish / algae / /molluscs / crustaceans

#The parameters of the density plot
weight = 10 #lnumeric (for log) / sqrt / none / occurence
sample_size = 100 #numeric (letters are not allowed)
band_width = 0.15 #numeric / SJ / nrd / nrd0


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


library(tidyverse)
library(dplyr)
library(ssdtools)


################
# Imports
################

df_list_raw <- readRDS(paste0(inwd, "/data_list_leth_subl.RData")) #check that it is a data frame


############################################
#main body
####################################

#choose the species group and endpoint that was defined in the consonants
if(species == "fish"){
  if(effect == "EC10"){
    df <- df_list_raw[["fish.EC10.subl.long"]]
  } else{
    df <- df_list_raw[["fish.EC50.leth"]]
    effect = "EC50"
  }
} else if(species == "algae"){
  if(effect == "EC10"){
    df <- df_list_raw[["algae.EC10.leth"]]
  } else{
    df <- df_list_raw[["algae.EC50.leth"]]
    effect = "EC50"
  }
} else if(species == "molluscs"){
  if(effect == "EC10"){
    df <- df_list_raw[["molluscs.EC10.subl.long"]]
  } else{
    df <- df_list_raw[["molluscs.EC50.leth"]]
    effect = "EC50"
  }
} else{
  species = "crustaceans"
  if(effect == "EC10"){
    df <- df_list_raw[["crustaceans.EC10.subl.long"]]
  } else{
    df <- df_list_raw[["crustaceans.EC50.leth"]]
    effect = "EC50"
  }
} 

#get the sum of the tests of one species with one chemical
weight_df <- aggregate(count~organism+CAS, data=df, sum)
weight_df$needed <- 1

#get rid of the chemicals with less than 6 species tested
temp <- aggregate(needed~CAS, data =weight_df, sum)
temp <- temp[temp$needed>=6,]
weight_df <- weight_df[weight_df$CAS %in% temp$CAS,]
weight_df$needed <- NULL

#multiply the amount of tests with one species and one chemicals 
#with the tests from the other species with the same chemical
weight_df <- aggregate(count~CAS, data = weight_df, prod)

#multiply the weight by the weighting factor defined in the constants
if(is.numeric(weight)){
  weight_df$count = log(weight_df$count+1, base = weight)
  weight = paste0("log", weight)
} else if(weight == "none"){
  weight_df$count = 1
} else if (weight == "sqrt"){
  weight_df$count = sqrt(weight_df$count)
} else{
  weight = "occurence"
}

###############################################
#make the SSDs with the random species selection
#####################################

#make an empty data frame where the results will be saved to
results <-  data.frame(matrix(ncol = 5, nrow = 0))
colnames(results) <- c("group", "chemical", "HC05", "HC95", "fit")

#make the SSD as often as defined in constants with random tests for every species
for(i in sample(weight_df$CAS, size = sample_size, replace = TRUE, prob = weight_df$count)){
  df2 <- df[df$CAS == i,]
  df2 = data.frame("mgperL" = tapply(df2$mgperL, df2$organism, function (x) x[sample(length(x), 1)]))
  
  SSD_fit_all <- ssd_fit_dists(data = df2, 
                               left =  "mgperL",
                               dists = c("burrIII3",
                                         "gamma",
                                         "lgumbel",
                                         "llogis", 
                                         "lnorm",
                                         "weibull"))
  
  
  #Examine the goodness of fit for the selected models, and put in in a new dataframe
  Fit_compare <- data.frame(ssd_gof(x = SSD_fit_all, 
                                    pvalue = TRUE))
  
  Best_fit_model <- Fit_compare[Fit_compare$delta == 0, "dist"]
  
  SSD_fit_best <- ssd_fit_dists(data = df2, 
                                left =  "mgperL",
                                dists = Best_fit_model)
  
  #Calcuate the HC05 from the curve, and out it in a new dataframe
  # using the best model 
  
  EC10_HC05_best <- data.frame(ssd_hc(x = SSD_fit_best,
                                      percent = c(5,95)))
  
  results[nrow(results) + 1,] <- c(paste0(species, "_", effect),
                                   i,
                                   EC10_HC05_best[1,3],
                                   EC10_HC05_best[2,3],
                                   names(SSD_fit_best))
}
results$HC95 <- as.numeric(results$HC95)
results$HC05 <- as.numeric(results$HC05)

results$range <- results$HC95/results$HC05
results$range_log <- log10(results$range)


if(is.numeric(band_width) | band_width == "SJ" | band_width == "nrd"){
  results_dens = density(results$range_log,
                         bw = band_width)
} else{
  results_dens = density(results$range_log,
                         bw = "nrd0")
}


if (!debbuging){
  pdf(file= paste0(outwd2, weight, "_", species, "_", effect,"_random_SSD.pdf"),
      width = 9,
      height = 6)
}

#define the graphical parameters
par(mar = c(3.1, 3, 0.2, 0.1), 
    mgp = c(2, 0.5, 0),
    las = 1,
    mfcol = c(1,1),
    cex.lab = 1.3,
    cex.main = 2,
    cex = 1.5,
    bty = "n",
    tcl = -0.3)

#make an plot with the outline of the density plot
plot(results_dens,
     main = "",
     xaxt = "n",
     xlab = "",
     xlim = c(-0.2, 8.5))

#give the density plot a nice base colour
polygon(results_dens, col = rgb(0.78, 0.89, 1, alpha = 0.6))

#add an x axis
axis(side = 1, line = -0.5, at = c(axTicks(1),1), labels = paste0("1e+", c(axTicks(1), 1)))


#calculate the amount of SSDs that have a range below 10
ten <- length(results$range_log[results$range_log <= 1])
percent_1 <- paste0(signif((ten/sample_size) * 100,digits = 3), "%")

#approximate the y position
position_1 <- approx(results_dens$x, results_dens$y, xout = 1)$y

#calculate the median range between the HC05 and HC95
median_x = median(results$range_log)
median_y <- approx(results_dens, xout = median_x)$y

#calculate the 95th percentile of the range between the HC05 and HC95
p95_x <- quantile(results$range_log, 0.95)
p95_y <- approx(results_dens, xout = p95_x)$y

#colour the three values that were calculated in
polygon(c(results_dens$x[results_dens$x >= 1 ], 1),
        c(results_dens$y[results_dens$x >= 1 ], 0),
        col = rgb(1, 0.9, 0.9, alpha = 1),
        border = 1)

polygon(c(results_dens$x[results_dens$x >= median_x ], median_x),
        c(results_dens$y[results_dens$x >= median_x ], 0),
        col = rgb(1, 0.8, 0.8, alpha = 1),
        border = 1)

polygon(c(results_dens$x[results_dens$x >= p95_x ], p95_x),
        c(results_dens$y[results_dens$x >= p95_x ], 0),
        col = rgb(1, 0.4, 0.4, alpha = 1),
        border = 1)

#add the labels
text(x = 1, y = position_1, labels = percent_1, pos = 3, xpd = TRUE)
text(x = median_x, y = median_y, labels = "50%", pos = 3, xpd = TRUE)
text(x = p95_x, y = p95_y, labels = "95%", pos = 3, xpd = TRUE)


#add the median range as a label at the bottom
axis(side = 1, line = 1, at = c(median_x, p95_x), 
     labels = paste0("1e+", c(signif(median_x, digits = 2), signif(p95_x, digits = 2))), 
     tick = FALSE)

#add the 95th percentile range as a label at the bottom
axis(side = 1, line = -0.5, at = c(median_x, p95_x), 
     labels = FALSE,
     tcl = -2)


if (!debbuging){
  dev.off()
}





