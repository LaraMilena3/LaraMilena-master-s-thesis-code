##############################################
# Standard Tasks at the beginning
##############################################

rm(list=ls()) #remove ALL objects 
Sys.setenv(LANG = "en") #to keep the language to english
Sys.setlocale("LC_ALL","English")
options(scipen=999) #avoid scientific annotation

##############################################
# Version Number
##############################################
vers<-2

##############################################
# Set constants
##############################################
debbuging = TRUE

par(tck = -0.025, las = 1, bty = 'o')

############
# Set the current workdirectory, i.e. where all files are kept
############
inwd <- "r_data"
outwd <- "r_output"
outwd2 <- "graphs/"

  ####################
# packages, imports, functions
################
library(tidyverse)
library(dplyr)


################
# Imports
################
df_list <- readRDS(paste0(inwd, "/data_list_leth_subl.RData")) #check that it is a data frame
for (i in 1:length(df_list)) {
  df_list[[i]]$test_location <- ifelse(grepl("FIELD", df_list[[i]]$test_location), "FIELD", "LAB")
  df_list[[i]] <- aggregate(mgperL ~ CAS+test_location+organism+species_group, data = df_list[[i]], function (x) exp(mean(log(x))))
  df_list[[i]]$group <- names(df_list[i])
}

df <- do.call(rbind, df_list)
test <- df[grepl("FIELD", df$test_location),]
test <- test[test$species_group == "algae",]
table(df$test_location)

df$count <- 1

results <- data.frame(CAS = NULL, ratio = NULL, group = NULL)
i = 3
for (i in 1:20){
  df2 <- df[df$group == names(df_list)[i],]
  
  #the way to only compare by cas not also by species
  
  #common_chemicals <- df2 %>%
  #  group_by(CAS) %>%
  #  filter(n_distinct(test_location) == 2) %>%  # Ensure the `chemical` appears in both groups
  #  pull(CAS) %>%
  #  unique()
  
  #Find common `chemical` values that exist in both groups
  common_chemicals <- df2 %>%
    group_by(CAS,organism) %>%
    filter(n_distinct(test_location) == 2) %>%  # Ensure the `chemical` appears in both groups
    pull(CAS) %>%
    unique()
  
  # 3. Filter the data to only keep rows with common `chemical` values
  df2 <- df2 %>%
    filter(CAS %in% common_chemicals)
  
  df2 <- df2[duplicated(df2[, c("CAS", "organism")]) | duplicated(df2[, c("CAS", "organism")], fromLast = TRUE), ]

  df2 <- df2[order(df2$test_location, decreasing = TRUE), ]
  
  df2 <- df2 %>%
    group_by(CAS, organism) %>%
    mutate(lab_tox = mgperL[1]) %>%
    ungroup()
  
  if (sum(grepl("FIELD", df2$test_location)) >= 1){
    test <- aggregate(test_location ~ CAS+organism+species_group+lab_tox, data = df2, function(x) x[1])
    df2 <- aggregate(mgperL ~ CAS+organism+species_group+lab_tox, data = df2, function(x) x[2]/x[1])
    df2$group <- paste0(names(df_list)[i], " (", nrow(df2), ")")
    
    results <- rbind(results, df2) 
  }
}
plot(mgperL~lab_tox,
     data = results,
     log = "xy")

lab_correlation = ggplot(results, aes(y = mgperL, x = lab_tox, colour = species_group)) +
  geom_point(size = 3) +
  labs(y = "lab to field ratio", x = "Lab efect concentration(mgperL)") +
  scale_color_manual(values = c("#90C819","#A82445", "#0099CC", "#C16200")) +
  geom_smooth(method='lm', aes(group = 1), colour = "black") +
  scale_x_continuous(trans='log10') +  
  scale_y_continuous(trans='log10') +
  theme_minimal(base_size = 20) 

if (!debbuging){
  ggsave(paste0(outwd2, "lab_field_corr.pdf"),
         plot = lab_correlation,
         width = 10,
         height = 7)
}
cor.test(log10(results$mgperL), log10(results$lab_tox), method = 'pearson')
cor.test(results$mgperL, results$lab_tox, method = 'spearman')



results$mgperL_log <- log10(results$mgperL)
results$group <- as.factor(results$group)
results$species_group <- as.factor(results$species_group)

bp <- boxplot(mgperL_log ~ species_group,
              data = results)


if (!debbuging){
  pdf(file= paste0(outwd2, "lab_field_box.pdf"),
      width = 9,
      height = 5)
}

par(mar = c(3.5, 5.5, 0.2, 10), las = 1, mfcol = c(1,1), mgp = c(4,1,0))

plot(1,
     type = "n", 
     xlab = "",
     ylab = "ratio",
     yaxt = "n",
     xaxt = "n",
     #main = "  common chemicals (environment/lab)",
     xlim = c(0.5, 4.5),
     ylim = c(-2, 4),
     cex.lab = 1.5,
     cex.main = 1.5)

rect(0,0,5,5.5, col = "#F7FCEE", border = NA)
rect(0,-2.5,5,0, col = "#FFECEB", border = NA)
segments(-1, 0, 5, 0, lwd =2, lty = 2)

#abline(h = axTicks(2), col = "lightgrey", lty = 2, lwd = 1.5)
#abline(h = c(0,1,2,3), col = "darkgrey", lty = 2, lwd = 1.5)


Q1 <- bp$stats[2,]  # First quartile (Q1)
Q3 <- bp$stats[4,]
tapply(10**results$ratio_log, results$group, FUN = function(x) median(x))

10**bp$stats[3,]
rect(0.6,Q1[1],1.4,Q3[1], col = "#90C819", border = NA)
rect(1.6,Q1[2],2.4,Q3[2], col = "#A82445", border = NA)
rect(2.6,Q1[3],3.4,Q3[3], col = "#0099CC", border = NA)
rect(3.6,Q1[4],4.4,Q3[4], col = "#C16200", border = NA)

#abline(abline(h = 1, col = "orange3", lwd = 2))

boxplot(mgperL_log ~ species_group,
        data = results,
        add = TRUE,
        col = NA,
        yaxt = "n",
        lwd = 2,
        xaxt = "n"
        #log = "y"
)

#points(mgperL_log ~ species_group,
#       data = results)

axis(2, at = seq(from = -2, to = 5, by = 1),
     labels = c(paste0("1e", seq(from = -2, to = -1, by = 1)),
                paste0("1e+", seq(from = 0, to = 5, by = 1))),
     cex.axis = 1.5)

axis(1, at = c(1,2,3,4),
     labels = c("algae\n(2)", "crustaceans\n(26)", "fish\n(25)", "molluscs\n(5)"),
     cex.axis = 1.5,
     mgp = c(3, 2.3, 0))

axis(side = 4,
     at = c(-1, 3),
     labels = c("Field more toxic", "Field less toxic"),
     cex.axis = 1.5)


if (!debbuging){
  dev.off()
}

sum(results$mgperL < 0.2)
sum(grepl("long", results$group))



################################
#make the densitiy plots
#########################3
par(bty = "n", mgp = c(3, 1, 0), mar = c(5,5,2.5,0.2), las = 1, mfcol = c(1,1))

results_dens = density(results$mgperL_log,
                       bw = 0.15)

plot(results_dens,
     main = "Lab to Field Variation",
     xaxt = "n",
     xlim = c(-2,7))

axis(side = 1, at = axTicks(1), labels = ifelse(axTicks(1)>=0, 
                                                paste0("1e+", axTicks(1)), 
                                                paste0("1e", axTicks(1))))

#give it a nice base colour

polygon(results_dens, col = rgb(0.78, 0.89, 1, alpha = 0.6))

median_x = median(results$mgperL_log)

median_y <- approx(results_dens, xout = median_x)$y



segments(median_x, 0, median_x, median_y, col = "darkblue", lwd = 2, lty = 2)


#add a label
#text(x = 1, y = position_1, labels = percent_1, pos = 3, xpd = TRUE)
text(x = median_x, y = median_y, labels = signif(10**median_x, digits = 2), pos = 3, xpd = TRUE)


#axis(side = 1, at = c(1,2,3,4), labels = c(1,2,3,4))

