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
vers<-5

##############################################
# Set constants
##############################################

#should this be saved? if yes = TRUE
FALSE -> printing


############
# Set the current workdirectory, i.e. where all files are kept
############

#the folder where the data comes from
inwd <- "r_data"
#the folder the data gets saved too
#since this file creates a smaller standardised csv file, it gets saved back into r_data
outwd <- "r_data"


################
# packages, imports, functions
################

library("tidyverse")
library("dplyr")
library("data.table")
library("webchem")
library("stringr")

#this function will be used to later to check if a test duration was durring a specific timeframe
#it returns either 1 or 0 which will be saved inside a column
time.interest = function(data, x, y){
  ifelse(data >= x & data <= y, 1, 0)
}

################
# Imports
################

data <- read.csv(paste0(inwd, "/raw_data.csv")) 
source_chemical <- read.csv(paste0(inwd, "/source_with_cas.csv")) 

####################################################################################################3
#main body
#####################################################################################
#################################
#removing duplicates and giving an index
##########################
data <- data[!duplicated(data), ]
setDT(data, keep.rownames = TRUE)[]
data$rn <- as.integer(data$rn)
data$milena_index <- data$rn + 1000000
data$rn <- NULL
data <- as.data.frame(data)

################################
#renaming  and redefining stuff
####################################

names(data)[1] <- "organism"
data$CAS <- as.character(data$CAS)


#there are differnt kind of field tests. 
#for this project that difference is not important, therefore all field tests get the same name
#if it's not declared as a field tests it's going to be a lab test, therefore all other cases are defined as LAB
data$test_location <- ifelse(grepl("FIELD", data$test_location), "FIELD", "LAB")

########################################
#remove rows that are not needed for the analysis
################################

#only tests with a clearly specified concentration can be used for further analysis
data = data[data$conc_sign == "=",]

#tests that did not specify a concentration or that reported the concentration as 0 can not be used
#a concentration of 0 can not be used because that doesn't make sense
data <- data[!is.na(data$mgperL),]
data <- data[data$mgperL != 0, ]

#####################################
#we can only need rows that have one tested chemical
#but some rows have the chemical information in the source column
#this part of the code checkes all unique words in the source column if they are a chemical

#if this section of the code has been run before it can get skipped, since it saves a file of the results
##############################

#find all the chemicals with a cas number in the source column
#and save it in the CAS column
no_CAS <- data[is.na(data$CAS),]
no_CAS <- unique(no_CAS$source)

#this checks if the words in the source column are a chemical
source_temp = cir_query(no_CAS,
                match = "first",
                representation = "cas"
)

#this adds the CAS number to one data frame
source_temp <- gsub("-", "", source_temp$CAS)
source_chemical <- source_temp[!is.na(source_temp),]

#save it as a file so this step can be skipped next time
write.csv(x = source_chemical, file = paste0(inwd, "/source_with_cas.csv"))

rm(no_CAS)
rm(source_temp)

###################################
#adding a cas number where a singular chemical is in the source field
######################

#to increase the speed all rows that have substance or reaction in the source column and no cas number or chemical name are removed
data <- data[!is.na(data$CAS) | !is.na(data$chemical_name) | !grepl(pattern = "ubstance|eaction", x = data$source), ]

#give the CAS column in the source_chemical dataframe a different name as in the main dataframe
names(source_chemical) <- c("query", "CAS2")
#replace NA values of Data with the CAS values from the source column
data <- data %>%
  left_join(source_chemical, by = c("source" = "query")) %>%
  mutate(CAS = ifelse(is.na(CAS), CAS2, CAS)) %>%
  select(-CAS2)

########################################
#remove additional rows that are not needed for the analysis
################################

#remove rows where the chemical is not specified
data <- data[!is.na(data$CAS),]

#only include the species groups we are interested in
data <- data[grepl("algae|crustaceans|fish|molluscs", data$species_group),]#only keep the species groups, we are intersted in

#remove unusable species
data$organism <- str_trim(data$organism, side = "right")
data <- data[grepl(" ", data$organism), ]
data <- data[!grepl("sp\\.", data$organism), ]
data <- data[!grepl("species", data$organism),] #remove rows where several species where tested
data <- data[!grepl("algae", data$organism),] #remove rows where the species is only specified as algae

#only keep the effects we are interested in
data <- data[grepl("DVP|GRO|ITX|MOR|MPH|POP|REP", data$effect),]


#only use EC10 and EC50 values and comparable

#EC10
data$endpoint[data$endpoint == "NOEC"|
                data$endpoint == "EC05"|
                data$endpoint == "NEC (NO EFFECT CONCENTRATION)"|
                data$endpoint == "EC07"|
                data$endpoint == "HIGHEST CONCENTRATION WITH 0% MORTALITY"|
                data$endpoint == "NOAEL"|
                data$endpoint == "HIGHEST CONCENTRATION WITHOUT OBSERVED EFFECT"|
                data$endpoint == "THE HIGHEST CONCENTRATION WITOUT OBSERVED EFFECTS"|
                data$endpoint == "NOELR"|
                data$endpoint == "IC07"|
                data$endpoint == "EC08"|
                data$endpoint == "MAX. CONC. CAUSING NO MORTALIT"|
                data$endpoint == "EL05"] <- "EC10"

#EC50
data$endpoint[data$endpoint == "TL50"|
                data$endpoint == "EBC50"| #reduction in biomass growth(EbC50)
                data$endpoint == "IGC 50"| #growth inhibitory concentration
                data$endpoint == "TLM (50% SURVIVAL POINT)"|
                data$endpoint == "TLM CONCENTRATION OF OIL CAUSING 50% MORTALITY"] <- "EC50"                                                                                                                 

data = data[data$endpoint == "EC50" | data$endpoint == "EC10",]

#remove columns that are not needed (to decrease size)
colums_to_remove <- c("rn",
                      "REACH_reliability",
                      "Conc._based_on",
                      "Type_of_Study",
                      "title",
                      "author",
                      "publication_year",
                      "common_name", 
                      "REACH_dossier",
                      "EC_Number",
                      "conc_sign",
                      "REACH_adequacy_of_study",
                      "REACH_reference_type",
                      "REACH_source_detailed_species_grouping",
                      "reference_number",
                      "organism_habitat",
                      "result_id",
                      "ncbi_taxid",
                      "database")
for (i in colums_to_remove){
  data[[i]] <- NULL
}
rm(colums_to_remove)

########################################
#rename species (because they have several names)
###############################
data$organism[data$organism == "pseudokirchneriella subcapitata"] <- "raphidocelis subcapitata"


#######################################################
#add columns
#######################################################

#if it's a lethal or sublethal effect
data$lethality <- ifelse(data$species_group == "algae" |
                           (data$species_group == "crustaceans" & data$effect == "ITX") |
                           data$effect == "MOR" |
                           data$effect == "POP",
                         "lethal", "sublethal")

#have the test duration in hours
data$Duration_Value = data$duration_days * 24

#add a boolean if a row contains a species we defined as a test species
data$test_species <-data$organism %in% c("daphnia magna", "ceriodaphnia dubia", "pimephales promelas", "oncorhynchus mykiss", "lepomis macrochirus", "danio rerio", "cyprinus carpio", "raphidocelis subcapitata", "desmodesmus subspicatus")

#######################################
#adding the information what group of chemical it is
#pesticide, pharmaceutical or industrial chemical
###############################

#defining pesticides
pesticide <- read.csv(paste0(inwd, "/pesticide_information_data.csv"))

#defining pharmaceuticals
pharmaceutical <- read.csv(paste0(inwd, "/Drugbank_library_2022-07-28.csv"))
names(pharmaceutical)[1] = "CAS"

pharmaceutical$CAS <- gsub("-", "", pharmaceutical$CAS)

#defining industrial chemicals
industrial <- read.csv(paste0(inwd, "/ECHA_industrial_chemicals.csv"), skip = 2) 
names(industrial)[3] = "CAS"

industrial$CAS <- gsub("-", "", industrial$CAS)

#adding it to the data frame
chemicals <- data.frame(CAS = unique(data$CAS))

chemicals$pesticide <- chemicals$CAS %in% pesticide$CAS
chemicals$pharmaceutical <- chemicals$CAS %in% pharmaceutical$CAS
chemicals$industrial <- chemicals$CAS %in% industrial$CAS

data <- merge(data, chemicals, by = "CAS", all.x = TRUE)

####################################
#add mode of actions to pesticides and pharmaceuticals
############################

#add the mode of actions of pesticides (fungicide, herbicide etc.)
df_pesticide <- read.csv(paste0(inwd, "/2024-03-26_MOA_data_with_ToC.csv"))

df_pesticide <- data.frame("CAS" = df_pesticide$CAS,"PPDB_MoA" = df_pesticide$PPDB_MoA)

df_pesticide$CAS <- as.character(df_pesticide$CAS)
data <- merge(data, df_pesticide, by = "CAS", all.x = TRUE)

data$PPDB_MoA <- tolower(data$PPDB_MoA)

#add the mode of actions of pharmaceuticals (atc code)
df_pharma <- read.csv(paste0(inwd, "/Drugbank_library_2022-07-28.csv"))
df_pharma = df_pharma[df_pharma$ATC_codes_drugbank != "",]
df_pharma <- df_pharma[!duplicated(df_pharma), ]
df_pharma <- df_pharma[df_pharma$Name_drugbank != "dalteparin",]
df_pharma = df_pharma[df_pharma$ï..CAS_drugbank != "Not Available", ]

df_pharma <- data.frame("CAS"= df_pharma$ï..CAS_drugbank,
                        ATC_codes_drugbank = df_pharma$ATC_codes_drugbank,
                        Solubility_drugbank = df_pharma$Solubility_drugbank,
                        MOA_drugbank = df_pharma$MOA_drugbank)

df_pharma$CAS <- gsub("-", "", df_pharma$CAS)

data <- merge(data, df_pharma, by = "CAS", all.x = TRUE)


df_pest <- read.csv(paste0(inwd, "/pesticide_information_data.csv"))


data <- merge(data, df_pest, by = "CAS", all.x = TRUE)


#############################################
#split the data frame into groups by endpoint and species group
#######################################

#split by endpoint and species groups
df_list <- split(data, list(data$endpoint, data$species_group))


names_df_list <- names(df_list)

######################################
#seperate the chemicals into time frames
#there are two timeframes a short one and a long one
#algae only has a short timeframe
#everything outside the timeframes gets removed
#except field tests that are longer than the longest time frame (they are still added in that case)
###############################3
y_values <- c(96,96,168,96,168,96,168,96)

#define the short timeframe and set the longer timeframe as 0 for now
for (i in 1:length(y_values)){
  df_list[[i]]$timeframe1 = sapply(df_list[[i]]$Duration_Value, x=24, y = y_values[i], FUN = time.interest)
  df_list[[i]]$timeframe2 = 0
}

#algae only have one timeframe therefore all field tests are included in the short timeframe
for(i in 1){
  df_list[[i]]$timeframe1 <- ifelse(grepl("FIELD", df_list[[i]]$test_location), 1, df_list[[i]]$timeframe1)
}

#add the longer timeframe to the groups that have a long timeframe
df_list[[3]]$timeframe2 <- ifelse(df_list[[3]]$Duration_Value >= 336 & grepl("FIELD", df_list[[3]]$test_location) | df_list[[3]]$Duration_Value >= 336 & df_list[[3]]$Duration_Value <= 744, 1,0)
df_list[[5]]$timeframe2 <- ifelse(df_list[[5]]$Duration_Value >= 336 & grepl("FIELD", df_list[[5]]$test_location) | df_list[[5]]$Duration_Value >= 336 & df_list[[5]]$Duration_Value <= 744, 1,0)
df_list[[7]]$timeframe2 <- ifelse(df_list[[7]]$Duration_Value >= 336 & grepl("FIELD", df_list[[7]]$test_location) | df_list[[7]]$Duration_Value >= 336 & df_list[[7]]$Duration_Value <= 744, 1,0)

#make a new list of data tables with the time frame information
df_list_timeframe <- list()
for (i in names_df_list){
  df_list_timeframe = append(df_list_timeframe, list(df_list[[i]][df_list[[i]]$timeframe1 == 1,]))
  if (sum(df_list[[i]]$timeframe2 != 0)){
    df_list_timeframe = append(df_list_timeframe, list(df_list[[i]][df_list[[i]]$timeframe2 == 1,]))
  }
}

#set the names of data frames
names(df_list_timeframe) <- c("EC10.algae","EC50.algae",
                              "EC10.crustaceans.short","EC10.crustaceans.long","EC50.crustaceans",
                              "EC10.fish.short","EC10.fish.long","EC50.fish", 
                              "EC10.molluscs.short","EC10.molluscs.long","EC50.molluscs")
names_timeframe_list = names(df_list_timeframe)
##########################################################
#separate the efects into lethal and sublethal
###############################################

#receive all the rows with lethal effects
lethal_df_list <- lapply(df_list_timeframe, function(df) {
  df[df$lethality == "lethal", ]
})

for (i in 1:length(names_timeframe_list)){
  names(lethal_df_list)[i] <- paste0(names_timeframe_list[i] , ".leth")
}

#receive all the rows with sublethal effects
sublethal_df_list <- lapply(df_list_timeframe, function(df) {
  df[df$lethality == "sublethal", ]
})

for (i in 1:length(names_timeframe_list)){
  names(sublethal_df_list)[i] <- paste0(names_timeframe_list[i] , ".subl")
}

#remove data frames that don't have any values in it
sublethal_df_list <- sublethal_df_list[sapply(sublethal_df_list, nrow) > 0]

#combine the lists
df_list_leth_subl <- c(lethal_df_list, sublethal_df_list)


#manually change the name so the order is more useful
names(df_list_leth_subl) <- c(
  "algae.EC10.leth",             "algae.EC50.leth",
  "crustaceans.EC10.leth.short", "crustaceans.EC10.leth.long", 
  "crustaceans.EC50.leth",       "fish.EC10.leth.short", 
  "fish.EC10.leth.long",         "fish.EC50.leth",      
  "molluscs.EC10.leth.short",    "molluscs.EC10.leth.long",    
  "molluscs.EC50.leth",          "crustaceans.EC10.subl.short",
  "crustaceans.EC10.subl.long",  "crustaceans.EC50.subl",      
  "fish.EC10.subl.short",        "fish.EC10.subl.long",        
  "fish.EC50.subl",              "molluscs.EC10.subl.short",   
  "molluscs.EC10.subl.long",     "molluscs.EC50.subl")

df_list_leth_subl <- df_list_leth_subl[sort(names(df_list_leth_subl))]

#add the group name into the data frame as a column
for (i in names(df_list_leth_subl)){
  df_list_leth_subl[[i]]$group <- i
}

#turn it into a data frame
df_leth_subl <- do.call(rbind, df_list_leth_subl)



###########################################################
#saving the data
####################################################


#TRUE -> printing


if (printing == TRUE){
  #save the data frame as a .csv file:
  write.csv(x = df_leth_subl, file = paste0(outwd, "/data_leth_subl.csv"))
  saveRDS(df_list_leth_subl, file=paste0(outwd, "/data_list_leth_subl.RData"))
}

