##############################################
# Standard Tasks at the beginning
##############################################

rm(list=ls()) #remove ALL objects 
#cat("\014") # clear console window prior to new run (can be runed somtimes, but i don't like to loss my output every time)
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

#should this be printed? if yes = TRUE
FALSE -> drucken


############
# Set the current workdirectory, i.e. where all files are kept
############

#This defines a string which contains the filepath for the working directory.
# These folders must already have been created beforehand.
# You dont have to call them "data" or "r_output", 
# if you want to use a different naming convention its up to you...

inwd <- "r_data"
inwd2 <- "r_output"
outwd <- "r_data"
outwd2 <- "r_output"


################
# packages, imports, functions
################

library("tidyverse")
library("dplyr")
library("data.table")
library("webchem")
library("stringr")

zeit.interesse = function(data, x, y){
  ifelse(data >= x & data <= y, 1, 0)
}

################
# Imports
################

#First thing, import the data.
# the filetype of the data is an Excel Workbook with an ".xlsx" suffix at the end of the filename.
# To import the data, we are using "read_xlsx".
#
# For other filetypes (e.g. ".csv"), there are other options such as read.csv() or read_csv()


data <- read.csv(paste0(inwd, "/raw_data.csv")) #check that it is a data frame
source_chemical <- read.csv(paste0(inwd, "/source_with_cas.csv")) #check that it is a data frame

data$test_location <- ifelse(grepl("FIELD", data$test_location), "FIELD", "LAB")

test <- data[grepl("FIELD", data$test_location),]
test <- test[test$species_group == "algae",]
#first <- as.data.frame(sort(table(test$endpoint), decreasing = TRUE))
#data <- test
#I use this constantly for data exploration 
data$count <- 1 #can be useful, for example for aggregate

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
#renaming stuff
####################################
names(data)[1] <- "organism"
data$CAS <- as.character(data$CAS)


#######################################################
#the steps in mikaels paper
#################################################

#step 2: klimisch score
#data <- data[data$REACH_reliability == "2 (reliable with restrictions)" | data$REACH_reliability == "1 (reliable without restriction)",]
#data <- data[!grepl('other|4|3', data$REACH_reliability),]

#step 3
#step 3.1: Numeric qualifiers other than "="
data = data[data$conc_sign == "=",]

#step3.3: not needed
#step 3.4: no mgperL added
data <- data[!is.na(data$mgperL),]

data <- data[data$mgperL != 0, ]

#####################################################
#add the CAS numbers where the source is a name of a chemical
#and then remove the values that still don't have a CAS number
#############################################
data <- data[!is.na(data$CAS) | !is.na(data$chemical_name) | !grepl(pattern = "ubstance|eaction", x = data$source), ]

names(source_chemical) <- c("query", "CAS2")
#replace NA values of Data with the CAS values from source
data <- data %>%
  left_join(source_chemical, by = c("source" = "query")) %>%
  mutate(CAS = ifelse(is.na(CAS), CAS2, CAS)) %>%
  select(-CAS2)


#step 1: remove values where the chemical is not specified

data <- data[!is.na(data$CAS),]

#I needed this for somewhere else and don't want it to break
#hydro <- data[!is.na(data$chemical_name)]


#####################################
#find all the source chemicals that have a cas number
##############################
#tested <- data[is.na(data$CAS),]
#tested <- unique(tested$source)

#abc = cir_query(tested,
#                match = "first",
#                representation = "cas"
#)

#abc <- gsub("-", "", abc$CAS)

#source_chemical <- abc[!is.na(abc),]
#rm(tested)
#rm(abc)

############################################
#do this now, so I have to remove less endpoints
####################################
data <- data[grepl("algae|crustaceans|fish|molluscs", data$species_group),]#only keep the species groups, we are intersted in

#step 3.2: remove unusable species rows
#data[c('genus', 'species')] <- str_split_fixed(data$organism, ' ', 2)
data$organism <- str_trim(data$organism, side = "right")
data <- data[grepl(" ", data$organism), ]
data <- data[!grepl("sp\\.", data$organism), ]
data <- data[!grepl("species", data$organism),]#remove rows where several species where tested
data <- data[!grepl("algae", data$organism),]#remove rows where the species is only specified as algae

data <- data[grepl("DVP|GRO|ITX|MOR|MPH|POP|REP", data$effect),]#only keep the effects we are interested in

#step 3.5: only use EC10 and EC50 values and comparable
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
                data$endpoint == "EBC50"| #eduction in biomass growth(EbC50)
                data$endpoint == "IGC 50"| #growth inhibitory concentration
                data$endpoint == "TLM (50% SURVIVAL POINT)"|
                data$endpoint == "TLM CONCENTRATION OF OIL CAUSING 50% MORTALITY"] <- "EC50"                                                                                                                 

data = data[data$endpoint == "EC50" | data$endpoint == "EC10",]


########################################
#rename species (for example because they have several names)
###############################
data$organism[data$organism == "pseudokirchneriella subcapitata"] <- "raphidocelis subcapitata"

##########################################################
#remove rows that are not needed
#remove test with several species
#should I also remove columns to make the life of my computer easier?
###########################################################



definitely_remove <- c("rn",
                       "REACH_reliability",
                       "Conc._based_on")

colums_to_remove <- c(definitely_remove,
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
table(data$test_location)


df <- read.csv(paste0(inwd, "/2024-03-26_MOA_data_with_ToC.csv")) #check that it is a data frame

df <- data.frame("CAS" = df$CAS,"PPDB_MoA" = df$PPDB_MoA)

df$CAS <- as.character(df$CAS)
data <- merge(data, df, by = "CAS", all.x = TRUE)

data$PPDB_MoA <- tolower(data$PPDB_MoA)

df <- read.csv(paste0(inwd, "/Drugbank_library_2022-07-28.csv"))
df = df[df$ATC_codes_drugbank != "",]
df <- df[!duplicated(df), ]
df <- df[df$Name_drugbank != "dalteparin",]
df = df[df$ï..CAS_drugbank != "Not Available", ]

df <- data.frame("CAS"= df$ï..CAS_drugbank,
                   ATC_codes_drugbank = df$ATC_codes_drugbank,
                   Solubility_drugbank = df$Solubility_drugbank,
                   MOA_drugbank = df$MOA_drugbank)

df$CAS <- gsub("-", "", df$CAS)

data <- merge(data, df, by = "CAS", all.x = TRUE)


pesticide <- read.csv(paste0(inwd, "/pesticide_information_data.csv"))


data <- merge(data, pesticide, by = "CAS", all.x = TRUE)

#write.csv(x = pesticides, file = paste0(outwd, "/pesticides_with_MOA.csv"))

#######################################################
#add useful things
#######################################################


data$lethality <- ifelse(data$species_group == "algae" |
                           (data$species_group == "crustaceans" & data$effect == "ITX") |
                           data$effect == "MOR" |
                           data$effect == "POP",
                         "lethal", "sublethal")

data$Duration_Value = data$duration_days * 24
data$test_species <-data$organism %in% c("daphnia magna", "ceriodaphnia dubia", "pimephales promelas", "oncorhynchus mykiss", "lepomis macrochirus", "danio rerio", "cyprinus carpio", "raphidocelis subcapitata", "desmodesmus subspicatus")

#######################################
#adding the information what gorup of chemical it is
###############################
chemicals <- data.frame(CAS = unique(data$CAS))


#pesticide data
#pesticide <- read.csv(paste0(inwd, "/2024-03-26_MOA_data_with_ToC.csv")) #check that it is a data frame

#pharmaceuticals
pharmaceutical <- read.csv(paste0(inwd, "/Drugbank_library_2022-07-28.csv"))
names(pharmaceutical)[1] = "CAS"

pharmaceutical$CAS <- gsub("-", "", pharmaceutical$CAS)

#industrial chemicals
industrial <- read.csv(paste0(inwd, "/ECHA_industrial_chemicals.csv"), skip = 2) #downloaded on 2025.01.16
names(industrial)[3] = "CAS"

industrial$CAS <- gsub("-", "", industrial$CAS)

##############################
#adding it to the list
#########################

chemicals$pesticide <- chemicals$CAS %in% pesticide$CAS
chemicals$pharmaceutical <- chemicals$CAS %in% pharmaceutical$CAS
chemicals$industrial <- chemicals$CAS %in% industrial$CAS

data <- merge(data, chemicals, by = "CAS", all.x = TRUE)
data$test_location <- ifelse(grepl("FIELD", data$test_location), "FIELD", data$test_location)

#############################################
#make the objects that are suposed to be saved
#######################################

df_list <- split(data, list(data$endpoint, data$species_group))


names_df_list <- names(df_list)


data_EC50 <- data[data$endpoint == "EC50",]
data_EC10 <- data[data$endpoint == "EC10",]
words <- unlist(strsplit(data$PPDB_MoA, " "))
word_count2 <- data.frame(sort(table(words), decreasing = TRUE))

#########################################################
# seperate the chemicals into time groups
###################################################


y_values <- c(96,96,168,96,168,96,168,96)


for (i in 1:length(y_values)){
  df_list[[i]]$timeframe1 = sapply(df_list[[i]]$Duration_Value, x=24, y = y_values[i], FUN = zeit.interesse)
  df_list[[i]]$timeframe2 = 0
}

for(i in 1){
  df_list[[i]]$timeframe1 <- ifelse(grepl("FIELD", df_list[[i]]$test_location), 1, df_list[[i]]$timeframe1)
}

df_list[[3]]$timeframe2 <- ifelse(df_list[[3]]$Duration_Value >= 336 & grepl("FIELD", df_list[[3]]$test_location) | df_list[[3]]$Duration_Value >= 336 & df_list[[3]]$Duration_Value <= 774, 1,0)
df_list[[5]]$timeframe2 <- ifelse(df_list[[5]]$Duration_Value >= 336 & grepl("FIELD", df_list[[5]]$test_location) | df_list[[5]]$Duration_Value >= 336 & df_list[[5]]$Duration_Value <= 774, 1,0)
df_list[[7]]$timeframe2 <- ifelse(df_list[[7]]$Duration_Value >= 336 & grepl("FIELD", df_list[[7]]$test_location) | df_list[[7]]$Duration_Value >= 336 & df_list[[7]]$Duration_Value <= 774, 1,0)

df_list_timeframe <- list()
for (i in names_df_list){
  df_list_timeframe = append(df_list_timeframe, list(df_list[[i]][df_list[[i]]$timeframe1 == 1,]))
  if (sum(df_list[[i]]$timeframe2 != 0)){
    df_list_timeframe = append(df_list_timeframe, list(df_list[[i]][df_list[[i]]$timeframe2 == 1,]))
  }
}
names(df_list_timeframe) <- c("EC10.algae","EC50.algae",
                              "EC10.crustaceans.short","EC10.crustaceans.long","EC50.crustaceans",
                              "EC10.fish.short","EC10.fish.long","EC50.fish", 
                              "EC10.molluscs.short","EC10.molluscs.long","EC50.molluscs")
names_timeframe_list = names(df_list_timeframe)
##########################################################
#separate the efects into lethal and sublethal
###############################################

lethal_df_list <- lapply(df_list_timeframe, function(df) {
  df[df$lethality == "lethal", ]
})

for (i in 1:length(names_timeframe_list)){
  names(lethal_df_list)[i] <- paste0(names_timeframe_list[i] , ".leth")
}

sublethal_df_list <- lapply(df_list_timeframe, function(df) {
  df[df$lethality == "sublethal", ]
})

#remove lists that don't have any values in it (not needed anymore)
#sublethal_df_list <- sublethal_df_list[sapply(sublethal_df_list, nrow) > 0]


for (i in 1:length(names_timeframe_list)){
  names(sublethal_df_list)[i] <- paste0(names_timeframe_list[i] , ".subl")
}

df_list_leth_subl <- c(lethal_df_list, sublethal_df_list)

for (i in names(df_list_leth_subl)){
  if (nrow(df_list_leth_subl[[i]]) == 0){
    df_list_leth_subl[[i]] <- NULL
  }
}


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
for (i in names(df_list_leth_subl)){
  df_list_leth_subl[[i]]$group <- i
}

df_leth_subl <- do.call(rbind, df_list_leth_subl)

sum(is.na(df_leth_subl$ATC_codes_drugbank))

sum(is.na(df_leth_subl$PPDB_MoA))


#######################################
#adding the information if reach sufficient
###############################

chemicals <- data.frame(CAS = unique(df_leth_subl$CAS))


df_leth_subl <- merge(df_leth_subl, chemicals, by = "CAS", all.x = TRUE)


###########################################################
#saving the data
####################################################


#TRUE -> drucken


if (drucken == TRUE){
  
  #save the data frame as a .csv file:
  write.csv(x = data, file = paste0(outwd, "/data_all.csv"))
  write.csv(x = data_EC50, file = paste0(outwd, "/data_EC50.csv"))
  write.csv(x = data_EC10, file = paste0(outwd, "/data_EC10.csv"))
  write.csv(x = df_leth_subl, file = paste0(outwd, "/data_leth_subl.csv"))
  
  
  saveRDS(df_list, file=paste0(outwd, "/data_list.RData"))
  saveRDS(df_list_timeframe, file=paste0(outwd, "/data_list_timeframe.RData"))
  saveRDS(df_list_leth_subl, file=paste0(outwd, "/data_list_leth_subl.RData"))
}



#df <- do.call(rbind, df_list_timeframe)
#hello = as.data.frame(sort(table(df$PPDB_MoA), decreasing = TRUE))
#write.csv(x = hello, file = paste0(outwd, "/mode_of_actions.csv"))

#as a .xlsx file (only needed if asked for):
#require(writexl)
#write_xlsx(x = ECOTOX_data, path = paste0(outwd, "/tidy_data_diclofenac_v", vers, ".xlsx"))



  

