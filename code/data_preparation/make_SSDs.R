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
vers<-2


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

df_list_raw <- readRDS(paste0(inwd, "/data_list_leth_subl.RData")) 
df_list = list()
#we are only interested in one concentration per chemcial with one species
#therefore all rows that had the same chemical and species are taken together with the geometric mean
for (i in 1:length(df_list_raw)) {
  try(df_list[[i]] <- aggregate(mgperL ~ organism+CAS, data = df_list_raw[[i]], function (x) exp(mean(log(x)))))
}
names(df_list) <- names(df_list_raw)

names_df_list <- names(df_list)


#########################################################
# get a list of all chemicals that have 6 tests
#the SSD tool used in this requires 6 different species to function
#therefore all other chemicals have to be removed
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

#make a data frame where the results will be saved in
effect_all <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(effect_all) <- c("group", "chemical", "HC05", "HC95", "fit")

#run it 20 times because we have 20 different data frames (one for every group)
for (i in 1:20){
  #this timer with the time code at the end is, so it can be seen how long it took
  start.time <- Sys.time()
  for (n in get(paste0("CAS_list", i))){
    df2 <- df_list[[i]][df_list[[i]]$CAS == n,]
    #SSDtools fails if all species have the same concentration
    #therefore it has to be checked that there is at last 2 unique values
    if (length(unique(df2$mgperL)) > 1){

      #using several models to fit the SSD so we can afterwards check which one fits the best
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
      
      #Calcuate the HC05  and the HC95 from the curve
      EC10_HC05_best <- data.frame(ssd_hc(x = SSD_fit_best,
                                          percent = c(5,95)))

      #saving it to the data frame created earlier (attaching another row)
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



#saving the results
#saving everything as characters helped to resolve a problem with saving
effect_all <- apply(effect_all,2,as.character)
write.csv(x = effect_all, file = paste0(outwd, "/effect_of_chemicals.csv"))
