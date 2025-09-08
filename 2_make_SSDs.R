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


############
# Set the current workdirectory, i.e. where all files are kept
############

inwd <- "r_data"
outwd <- "r_output"


####################
# packages, imports, functions
################


require(tidyverse)
library(ssdtools)


################
# Imports
################

df_list_raw <- readRDS(paste0(inwd, "/data_list_leth_subl.RData")) #check that it is a data frame
df_list = list()
for (i in 1:length(df_list_raw)) {
  try(df_list[[i]] <- aggregate(mgperL ~ organism+CAS, data = df_list_raw[[i]], function (x) exp(mean(log(x)))))
}
names(df_list) <- names(df_list_raw)

names_df_list <- names(df_list)


#########################################################
# get a list of all chemicals that have 6 tests
###################################################

for (i in 1:length(names_df_list)){
  
  CAS_list = df_list[[i]]
  
  CAS_list$count <- 1
  
  CAS_list <- aggregate(count ~ CAS, data = CAS_list, "sum")
  
  CAS_list_6 <- CAS_list[CAS_list$count >= 6,]
  assign(paste("CAS_list",i, sep=""), CAS_list_6$CAS)
}#make the chemical lists for all the chemicals



###############################################
#make the SSDs
#####################################

effect_all <- data.frame(matrix(ncol = 5, nrow = 0))

colnames(effect_all) <- c("group", "chemical", "HC05", "HC95", "fit")


for (i in 1:20){
  start.time <- Sys.time()
  for (n in get(paste0("CAS_list", i))){
    df2 <- df_list[[i]][df_list[[i]]$CAS == n,]
    if (length(unique(df2$mgperL)) > 1){
      
      SSD_fit_all <- ssd_fit_dists(data = df2, 
                                   left =  "mgperL",
                                   dists = c("burrIII3",
                                             "gamma",
                                             "lgumbel",
                                             "llogis", 
                                             "lnorm",
                                             "weibull"))
      
      #Although we can generally see by eye which one os the best fit, lets choose the model thats
      # the best fit for our data based on the goodness of fit statistics:
      
      #Examine the goodness of fit for the selected models, and put in in a new dataframe
      Fit_compare <- data.frame(ssd_gof(x = SSD_fit_all, 
                                        pvalue = TRUE))
      
      Best_fit_model <- Fit_compare[Fit_compare$delta == 0, "dist"]
      
      SSD_fit_best <- ssd_fit_dists(data = df2, 
                                    left =  "mgperL",
                                    dists = Best_fit_model)
      
      #Calcuate the HC05 from the curve, and out it in a new dataframe
      #   1.  using the best model 
      
      EC10_HC05_best <- data.frame(ssd_hc(x = SSD_fit_best,
                                          percent = c(5,95)))
      
      effect_all[nrow(effect_all) + 1,] <- c(names(df_list[i]),
                                             n,
                                             EC10_HC05_best[1,3],
                                             EC10_HC05_best[2,3],
                                             SSD_fit_best)
    }
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
}




effect_all <- apply(effect_all,2,as.character)
write.csv(x = effect_all, file = paste0(outwd, "/effect_of_chemicals.csv"))
