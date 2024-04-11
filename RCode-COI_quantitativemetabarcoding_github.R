# 
# R script for:
# Quantitative assessment of reef foraminifera community from metabarcoding data
# Elsa B. Girard a,b, Emilie A. Didaskalou c, Andi M. A. Pratama d, Carolina Rattner a, Raphaël Morard e, Willem Renema a,b

# a Naturalis Biodiversity Center, Darwinweg 2, 2333 CR Leiden, the Netherlands
# b IBED, University of Amsterdam, Sciencepark 904, 1098 XH Amsterdam, the Netherlands
# c CML, University of Leiden, Einsteinweg 2, 2333 CC Leiden, the Netherlands
# d Marine Science Department, Faculty of Marine Science and Fisheries, Hasanuddin University, Jl. Perintis Kemerdekaan Km. 10 Tamalenrea, Makassar 90245, Indonesia
# e MARUM, University of Bremen, Leobener Str. 8, 28359 Bremen, Germany




options(rstudio.help.showDataPreview = FALSE)

rm(list=ls())
setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
setSessionTimeLimit(cpu = Inf, elapsed = Inf)

#package needed
library(stringr)
library(pheatmap)
library(vegan)
library(ggplot2)
library(reshape2)
library(dplyr)
library(car)
library(patchwork)
library(ggvenn)
library(ggbreak)
library(scales)
library(cowplot)
library(ggpubr)
library(gghalves)
library(tidyverse)
library(rstatix)
library(devtools)
library(propagate)
library(ggpol)

"%not%" <- Negate("%in%") #to create a "not in" sign, the opposite of %in%
"is.not.na" <- Negate("is.na")

######### ddPCR and CT-scanning ############################

#--------load and arrange dataset----

df <- read.csv("~/Downloads/Dataset/Data_area_volume_genecopies_allspecimens_corrected.csv")
df$area_mm2 <- df$area_um2/1000000
df$volume_chambers_mm3 <- df$volume_chambers_um3/1000000000
df$volume_shell_mm3 <- df$volume_shell_um3/1000000000
df$volume_residuals_mm3 <- df$volume_residuals_um3/1000000000
df$total_gene_copies_in_specimen_million <- df$total_gene_copies_in_specimen/1000000

#make new total shell volume depending on the residual category
df$final_volume_shell_mm3 <- NA
for (i in 1:nrow(df)) {
  if (is.na(df$residual_category[i])) next
  if (df$residual_category[i] == "shell") {
    df$final_volume_shell_mm3[i] <- df$volume_shell_mm3[i] + df$volume_residuals_mm3[i]
  } else {
    df$final_volume_shell_mm3[i] <- df$volume_shell_mm3[i]
  }
}

DNA <- df[df$purpose == "DNA",]
DNA$area_mm2_log <- log(DNA$area_mm2)
DNA$volume_mm3_log <- log(DNA$volume_chambers_mm3)
DNA$total_gene_copies_in_specimen_million_log <- log(DNA$total_gene_copies_in_specimen_million)
DNA$density <- DNA$total_gene_copies_in_specimen_million/DNA$area_mm2

CTscan <- df[df$purpose == "ct-scan",]
CTscan$surf_vol_ratio <- CTscan$area_mm2 / CTscan$volume_shell_mm3

DNA1 <- DNA %>% group_by(species) %>% slice_sample(n=30)
DNA1 <- DNA1[,colSums(is.na(DNA1))==0]
DNA1$calculated_volume_chambers_mm3 <- NA

depth <- read.csv("~/Downloads/Dataset/ddpcr_sample_depth.csv")

DNA5 <- merge(DNA1, depth, by = "sample_id", all=TRUE)
DNA5 <- DNA5[DNA5$species %not% NA,]
DNA1 <- DNA5

test <- DNA1 %>% group_by(species) %>% summarize(mean_copies = mean(total_gene_copies_in_specimen_million), mean_density = mean(density))

#--------fig4 - Show graph with a 30 specimen subset two models-----------------------

p1 <- ggplot(DNA, aes(x=area_mm2, y=total_gene_copies_in_specimen_million)) +
  geom_point(aes(col = species), alpha = 0.3, size = 1) +
  stat_smooth(aes(col = species), method = "lm", formula = y~x, se = FALSE) + 
  stat_cor(aes(col = species)) +
  scale_color_manual(values =c("Amphisorus SpL" = "#00B6EB","Amphistegina lessonii" = "#FB61D7","Baculogypsinoides spinosus" = "#F8766D",
                               "Calcarina spengleri" = "#E6194B","Heterostegina depressa" = "#000075",
                               "Neorotalia gaimardi" = "#53B400","Operculina ammonoides" = "#A58AFF")) +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  theme_classic2() + 
  ylab("Number of gene copies (millions)") + 
  xlab("Surface area (mm2)") +
  labs(tag = "A") +
  theme(legend.position = "none",strip.background = element_blank(), strip.text.x = element_blank())

p2 <- ggplot(DNA, aes(y=density, x=species)) +
  geom_boxplot(aes(col = species), size = 0.5) +
  scale_color_manual(values =c("Amphisorus SpL" = "#00B6EB","Amphistegina lessonii" = "#FB61D7","Baculogypsinoides spinosus" = "#F8766D",
                               "Calcarina spengleri" = "#E6194B","Heterostegina depressa" = "#000075",
                               "Neorotalia gaimardi" = "#53B400","Operculina ammonoides" = "#A58AFF")) +
  theme_classic2()+
  scale_y_continuous(trans='log10') +
  ylab("Copies (millions)/surface area (mm2)") + 
  labs(tag = "B") +
  theme(legend.position = "none",legend.title = element_blank(),strip.background = element_blank(), strip.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.title.x = element_blank()) 


layout <- "
A
A
A
B
"
p1 + p2 + plot_layout(design = layout)

#to check for significance between boxplot
ggplot(DNA, aes(y=density, x=species)) +
  geom_boxplot(aes(col = species), size = 0.5) +
  scale_color_manual(values =c("Amphisorus SpL" = "#00B6EB","Amphistegina lessonii" = "#FB61D7","Baculogypsinoides spinosus" = "#F8766D",
                               "Calcarina spengleri" = "#E6194B","Heterostegina depressa" = "#000075",
                               "Neorotalia gaimardi" = "#53B400","Operculina ammonoides" = "#A58AFF")) +
  theme_classic2()+
  scale_y_continuous(trans='log10') +
  geom_signif(comparisons = list(c("Amphisorus SpL", "Calcarina spengleri")), 
              map_signif_level=TRUE, 
              y_position = 0.4, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Amphisorus SpL", "Heterostegina depressa")), 
              map_signif_level=TRUE,
              y_position = 0.6, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Amphistegina lessonii", "Calcarina spengleri")), 
              map_signif_level=TRUE,
              y_position = 0.8, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Amphistegina lessonii", "Heterostegina depressa")), 
              map_signif_level=TRUE,
              y_position = 1, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Amphistegina lessonii", "Neorotalia gaimardi")), 
              map_signif_level=TRUE,
              y_position = 1.2, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Baculogypsinoides spinosus", "Calcarina spengleri")), 
              map_signif_level=TRUE,
              y_position = 1.4, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Baculogypsinoides spinosus", "Heterostegina depressa")), 
              map_signif_level=TRUE,
              y_position = 1.6, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Baculogypsinoides spinosus", "Neorotalia gaimardi")), 
              map_signif_level=TRUE,
              y_position = 1.8, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Calcarina spengleri", "Heterostegina depressa")), 
              map_signif_level=TRUE,
              y_position = 2.0, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Calcarina spengleri", "Neorotalia gaimardi")), 
              map_signif_level=TRUE,
              y_position = 2.2, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Calcarina spengleri", "Operculina ammonoides")), 
              map_signif_level=TRUE,
              y_position = 2.4, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Heterostegina depressa", "Neorotalia gaimardi")), 
              map_signif_level=TRUE,
              y_position = 2.6, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Heterostegina depressa", "Operculina ammonoides")), 
              map_signif_level=TRUE,
              y_position = 2.8, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Neorotalia gaimardi", "Operculina ammonoides")), 
              map_signif_level=TRUE,
              y_position = 3, tip_length = 0.01, vjust = .1) +
  ylab("Copies (millions)/surface area (mm2)") + 
  labs(tag = "B") +
  theme(legend.position = "none",legend.title = element_blank(),strip.background = element_blank(), strip.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.title.x = element_blank()) 



#--------Correction factors: Get powerlaw + linear regression coefficient for all species surface vs copies with 999 permutation------------------

#surface vs gene copies: extract a and b parameter from y ~ a*x^b (powerlaw regression) and parameter a from y ~ a*x + 0 (linear regression)

species <- c(1:6993)
a <- c(1:6993) #coef a powerlaw
b <- c(1:6993) #coef b powerlaw
se_a <- c(1:6993) #error a powerlaw
se_b <- c(1:6993) #error b powerlaw
pr_a <- c(1:6993) #prob a powerlaw
pr_b <- c(1:6993) #prob b powerlaw
log_rsquared <- c(1:6993) #r squared powerlaw
log_rsquared_adj <- c(1:6993) #r squared adjusted powerlaw
log_pvalue <- c(1:6993) #p-value powerlaw
lma <- c(1:6993) #coef a linear
se_lma <- c(1:6993) #error a linear
pr_lma <- c(1:6993) #prob a linear
lm_rsquared <- c(1:6993) #r squared linear
lm_rsquared_adj <- c(1:6993) #r squared adjusted linear
lm_pvalue <- c(1:6993) #p-value linear
regression_result <- data.frame(species, b, a, se_b, se_a, pr_b, pr_a,
                                log_rsquared, log_rsquared_adj,log_pvalue,
                                lma, se_lma, pr_lma,lm_rsquared, lm_rsquared_adj, lm_pvalue)
regression_result$species[1:6993] <- unique(DNA1$species)

DNA2 <- DNA[, c("species", "total_gene_copies_in_specimen_million", "area_mm2", "total_gene_copies_in_specimen_million_log", "area_mm2_log")]

w <- 0

#define function to extract overall p-value of model
#extract overall p-value of model
#overall_p(model)
overall_p <- function(my_model) {
  f <- summary(my_model)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}



for (v in 1:999) {
  
  DNA3 <- DNA2 %>% group_by(species) %>% slice_sample(n=30)
  
  for (i in 1:7) {
    
    w <- w + 1
    
    model <- lm(total_gene_copies_in_specimen_million_log ~ area_mm2_log,               #log(y) ~ a * log(x) + b
                data = DNA3[DNA3$species == regression_result[w,1],])
    
    regression_result[w, 2] <- coef(summary(model))[[1,1]]
    regression_result[w, 3] <- coef(summary(model))[[2,1]]
    regression_result[w, 4] <- coef(summary(model))[[1,2]]
    regression_result[w, 5] <- coef(summary(model))[[2,2]]
    regression_result[w, 6] <- coef(summary(model))[[1,4]]
    regression_result[w, 7] <- coef(summary(model))[[2,4]]
    regression_result[w, 8] <- summary(model)$r.squared
    regression_result[w, 9] <- summary(model)$adj.r.squared
    regression_result[w, 10] <- overall_p(model)
    
    model <- lm(total_gene_copies_in_specimen_million ~ area_mm2 + 0,                   #y ~ lma * x
                data = DNA3[DNA3$species == regression_result[[w,1]],])
    
    regression_result[w, 11] <- coef(summary(model))[[1,1]]
    regression_result[w, 12] <- coef(summary(model))[[1,2]]
    regression_result[w, 13] <- coef(summary(model))[[1,4]]
    regression_result[w, 14] <- summary(model)$r.squared
    regression_result[w, 15] <- summary(model)$adj.r.squared
    regression_result[w, 16] <- overall_p(model)
    
  }
  
  
}

reg_coef_mean <- regression_result %>% 
  group_by(species) %>% 
  summarize(mean_a = mean(a, na.rm=T), sd_a = sd(a, na.rm = T), mean_se_a = mean(se_a, na.rm=T), mean_pr_a = mean(pr_a, na.rm=T),
            mean_b = mean(b, na.rm=T), sd_b = sd(b, na.rm = T), mean_se_b = mean(se_b, na.rm=T),  mean_pr_b = mean(pr_b, na.rm=T), 
            mean_rsq_log = mean(log_rsquared, na.rm=T), mean_rsq_adj_log = mean(log_rsquared_adj, na.rm=T), mean_pvalue_log = mean(log_pvalue, na.rm=T),
            mean_lma = mean(lma, na.rm=T), sd_lma = sd(lma, na.rm = T), mean_se_lma = mean(se_lma, na.rm=T), mean_pr_lma = mean(pr_lma, na.rm=T), 
            mean_rsq_lm = mean(lm_rsquared, na.rm=T), mean_rsq_adj_lm = mean(lm_rsquared_adj, na.rm=T), mean_pvalue_lm = mean(lm_pvalue, na.rm=T))

write.csv(reg_coef_mean, "C:/Temp/01-PhD-project/Chapter3-Quantitativemetabarcoding/_Manuscript-Quantitativemetabarcoding+ddPCR/Size-vs-genecopy/ddpcr_powerlaw_linear_coefficients_permutations999.csv", row.names=FALSE)




#--------fig3 - Morphology weight, cell volume and shell volume vs surface area ratio estimation--------

ctscanh <- CTscan[CTscan$species == "Heterostegina depressa",]


p1 <- ggplot(CTscan, aes(x=area_mm2, y=volume_chambers_mm3)) +
  geom_point(aes(col = species), alpha = 0.3, size = 1) +
  stat_smooth(aes(col = species), method = "lm", formula = y~x, 
              method.args = list(start= c(a=1, b=1)),
              se=FALSE) +
  stat_cor(aes(col = species))+
  scale_color_manual(values =c("Amphisorus SpL" = "#00B6EB","Amphistegina lessonii" = "#FB61D7","Baculogypsinoides spinosus" = "#F8766D",
                               "Calcarina spengleri" = "#E6194B","Heterostegina depressa" = "#000075",
                               "Neorotalia gaimardi" = "#53B400","Operculina ammonoides" = "#A58AFF")) +
  theme_classic2()+ 
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  ylab("Cell volume (mm3)") + 
  xlab("Surface area (mm2)") +
  labs(tag = "A") +
  theme(legend.position = "none",legend.title = element_blank(),strip.background = element_blank(), strip.text.x = element_blank())

p2 <- ggplot(CTscan, aes(x=area_mm2, y=final_volume_shell_mm3)) +
  geom_point(aes(col = species), alpha = 0.3, size = 1) +
  stat_smooth(aes(col = species), method = "lm", formula = y~x, 
              method.args = list(start= c(a=1, b=1)),
              se=FALSE) +
  stat_cor(aes(col = species))+
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  scale_color_manual(values =c("Amphisorus SpL" = "#00B6EB","Amphistegina lessonii" = "#FB61D7","Baculogypsinoides spinosus" = "#F8766D",
                               "Calcarina spengleri" = "#E6194B","Heterostegina depressa" = "#000075",
                               "Neorotalia gaimardi" = "#53B400","Operculina ammonoides" = "#A58AFF")) +
  theme_classic2()+ 
  ylab("Shell volume (mm3)") + 
  xlab("Surface area (mm2)") +
  labs(tag = "B") +
  theme(legend.position = "none",legend.title = element_blank(),strip.background = element_blank(), strip.text.x = element_blank())

p3 <- ggplot(CTscan, aes(y=surf_vol_ratio, x=species)) +
  geom_boxplot(aes(col = species), size = 0.5) +
  scale_color_manual(values =c("Amphisorus SpL" = "#00B6EB","Amphistegina lessonii" = "#FB61D7","Baculogypsinoides spinosus" = "#F8766D",
                               "Calcarina spengleri" = "#E6194B","Heterostegina depressa" = "#000075",
                               "Neorotalia gaimardi" = "#53B400","Operculina ammonoides" = "#A58AFF")) +
  theme_classic2()+
  #scale_y_continuous(trans='log10') +
  ylab("Surface area to shell volume ratio") + 
  labs(tag = "C") +
  theme(legend.position = "none",legend.title = element_blank(),strip.background = element_blank(), strip.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.title.x = element_blank()) 

p4 <- ggplot(ctscanh[ctscanh$species == "Heterostegina depressa",], aes(x=final_volume_shell_mm3, y=weight_mg)) +
  geom_point(col = "#000075", alpha = 0.3, size = 2) +
  stat_smooth(col = "#000075", method = "lm", formula = y~x+0, 
              se=FALSE, fullrange=TRUE) +
  stat_cor(col = "#000075")+
  theme_classic2()+ 
  ylab("Weight (mg)") + 
  xlab("Shell volume (mm3)") +
  labs(tag = "D") +
  theme(legend.position = "none")

p5 <- ggplot(ctscanh[ctscanh$species == "Heterostegina depressa",], aes(x=area_mm2, y=weight_mg)) +
  geom_point(col = "#000075", alpha = 0.3, size = 2) +
  stat_smooth(col = "#000075", method = "lm", formula = y~x+0, 
              se=FALSE, fullrange=TRUE) +
  stat_cor(col = "#000075")+
  theme_classic2()+ 
  ylab("Weight (mg)") + 
  xlab("Surface area (mm2)") +
  labs(tag = "E") +
  theme(legend.position = "none")


layout <- "
AB
AB
AB
CC
DE
"
p1 + p2 + p3 + p4 + p5 + plot_layout(design = layout)

#to calculate differences
ggplot(CTscan, aes(y=surf_vol_ratio, x=species)) +
  geom_boxplot(aes(col = species), size = 0.5) +
  scale_color_manual(values =c("Amphisorus SpL" = "#00B6EB","Amphistegina lessonii" = "#FB61D7","Baculogypsinoides spinosus" = "#F8766D",
                               "Calcarina spengleri" = "#E6194B","Heterostegina depressa" = "#000075",
                               "Neorotalia gaimardi" = "#53B400","Operculina ammonoides" = "#A58AFF")) +
  geom_signif(comparisons = list(c("Amphisorus SpL", "Amphistegina lessonii")), 
              map_signif_level=TRUE, 
              y_position = 1, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Amphisorus SpL", "Baculogypsinoides spinosus")), 
              map_signif_level=TRUE, 
              y_position = 2, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Amphisorus SpL", "Calcarina spengleri")), 
              map_signif_level=TRUE, 
              y_position = 3, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Amphisorus SpL", "Heterostegina depressa")), 
              map_signif_level=TRUE,
              y_position = 4, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Amphisorus SpL", "Neorotalia gaimardi")), 
              map_signif_level=TRUE, 
              y_position = 5, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Amphisorus SpL", "Operculina ammonoides")), 
              map_signif_level=TRUE, 
              y_position = 6, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Amphistegina lessonii", "Baculogypsinoides spinosus")), 
              map_signif_level=TRUE,
              y_position = 7, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Amphistegina lessonii", "Calcarina spengleri")), 
              map_signif_level=TRUE,
              y_position = 8, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Amphistegina lessonii", "Heterostegina depressa")), 
              map_signif_level=TRUE,
              y_position = 9, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Amphistegina lessonii", "Neorotalia gaimardi")), 
              map_signif_level=TRUE,
              y_position = 10, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Amphistegina lessonii", "Operculina ammonoides")), 
              map_signif_level=TRUE,
              y_position = 11, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Baculogypsinoides spinosus", "Calcarina spengleri")), 
              map_signif_level=TRUE,
              y_position = 12, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Baculogypsinoides spinosus", "Heterostegina depressa")), 
              map_signif_level=TRUE,
              y_position = 13, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Baculogypsinoides spinosus", "Neorotalia gaimardi")), 
              map_signif_level=TRUE,
              y_position = 14, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Baculogypsinoides spinosus", "Operculina ammonoides")), 
              map_signif_level=TRUE,
              y_position = 15, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Calcarina spengleri", "Heterostegina depressa")), 
              map_signif_level=TRUE,
              y_position = 16, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Calcarina spengleri", "Neorotalia gaimardi")), 
              map_signif_level=TRUE,
              y_position = 17, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Calcarina spengleri", "Operculina ammonoides")), 
              map_signif_level=TRUE,
              y_position = 18, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Heterostegina depressa", "Neorotalia gaimardi")), 
              map_signif_level=TRUE,
              y_position = 19, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Heterostegina depressa", "Operculina ammonoides")), 
              map_signif_level=TRUE,
              y_position = 20, tip_length = 0.01, vjust = .1) +
  geom_signif(comparisons = list(c("Neorotalia gaimardi", "Operculina ammonoides")), 
              map_signif_level=TRUE,
              y_position = 21, tip_length = 0.01, vjust = .1) +
  theme_classic2()+
  #scale_y_continuous(trans='log10') +
  ylab("Surface area to shell volume ratio") + 
  labs(tag = "C") +
  theme(legend.position = "none",legend.title = element_blank(),strip.background = element_blank(), strip.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.title.x = element_blank()) 
















######### METABARCODING PART ############################

#--------NEW - prep data - Filter dataset 0.1%-----------------------------------------------
df <- read.csv("~/Downloads/Dataset/Forams_run_apscale_ESV_table_filtered_target.csv", row.names = 1)

#df[row,column]
df1 <- df[,-c(209)]

#transformation switching rows and columns
df1 <- t(df1)

df2 <- as.data.frame(df1)

#remove samples with less than 1000 reads
df2 <- df2[rowSums(df2) > 1000,]
otudf <- df2

sum(otudf)

#transform raw values into percentages 
df2 <- as.data.frame(t(apply(df2, 1, function(x) x/sum(x))))
#make sure row sum equal 1
rowSums(df2)

#remove asv < 0.1% of reads of the total reads but keeping reads number
for (i in 1:ncol(otudf)) {
  
  for (x in 1:nrow(otudf)) {
    
    if (df2[x,i] < 0.001) { #to account for cross contamination and barcode switching during sequencing
      
      otudf[x,i] <- 0
    }
  }
}

rowSums(otudf)
colSums(otudf)

#keep ESVs that still have reads
otudf <- as.data.frame(otudf[,colSums(otudf) > 0]) #from 2959 asv to 532 asv

sum(otudf)

21145944/21522300*100

#transform back data (switch rows and columns)
datafilt <- as.data.frame(t(otudf))

#create a column with ESVs
df$ESVs <- row.names(df) #original data
datafilt$ESVs <- row.names(datafilt) #filtered data

#match ESVs from filtered data with original metadata
datafilt1 <- plyr::match_df(df[,c(209,210)], datafilt, on = "ESVs")

#combine datasets
datafilt1 <- cbind(datafilt1[,c(1,2)], datafilt)

#rearrange column order
datafilt2 <- datafilt1[, c(2,1,3:209)]

colnames(datafilt2)

write.csv(datafilt2, "~/Downloads/Dataset/foram_data_apscale_qualitycheck.csv")

#--------NEW - prep data - create dataset with metadata + blast---------------------------------------------------------------

#load data
df <- read.csv("~/Downloads/Dataset/foram_data_apscale_qualitycheck.csv", row.names = 1)
metadata <- read.csv("~/Downloads/Dataset/METADATA_ForamSamples-bulkvspicked.csv")
taxo <- read.csv("~/Downloads/Dataset/blast_forams_galaxy.csv")

#transform the dataset without ESVs info
otudf <- as.data.frame(t(df[,-c(1:2)]))
otudf$baseclear_ID_real <- row.names(otudf) #create column with baseclear ID real names
otudf <- otudf[,c(ncol(otudf), 1:(ncol(otudf)-1))] #rearrange the columns 
row.names(otudf) <- NULL #remove row names

#join metadata with asv data
df1 <- merge(metadata[,c(3,7,10,16,17)], otudf, by = "baseclear_ID_real")

#join esv metadata (blast)
df_m <- melt(df1, id.vars = c("baseclear_ID_real","field.nmbr.","island", "type","new.name"), 
             value.name = "reads", variable.name = "ESVs")
df2 <- merge(df_m, df[,c(1:2)], by = "ESVs", all = TRUE)
df3 <- merge(df2, taxo[,c(1:6)], by = "ESVs", all = TRUE)

df4 <- df3[df3$reads > 0,]

write.csv(df4, "~/Downloads/Dataset/workingdataset-quantitativemetabarcoding_new.csv", row.names=FALSE)

#--------NEW - prep data - create 99.4% ID datasets-----------------------------------------------------------

df <- read.csv("~/Downloads/Dataset/workingdataset-quantitativemetabarcoding_new.csv")

df <- df[!is.na(df$ESVs),]

df1 <- df
df1$only_species_level <- NA

for (i in 1:nrow(df1)) {
      
      if (is.na(df1[i,10])) next #in case we have an NA value, stop loop and go to next row

      if (df1[i,10] < 99.4) {
      
           df1[i,14] <- "Z_other_foraminifera"
      
            } else { 
      
      df1[i,14] <- df1[i,9] 
      
            }}

df1$species[is.na(df1$identity_percentage)] <- "Z_not_foraminifera"
df1$only_species_level[is.na(df1$identity_percentage)] <- "Z_not_foraminifera"

write.csv(df1, "~/Downloads/Dataset/workingdataset-quantitativemetabarcoding-994_new.csv", row.names=FALSE)

#__________________stats on dataset___________________
df <- read.csv("~/Downloads/Dataset/workingdataset-quantitativemetabarcoding-994_new.csv")

#______all samples together________________________

df1 <- df[,c(1,14)]
df1$counts <- 1
df2 <- unique(df1)

#total number of ESVs
sum(df2$counts)

#number of ESVs assigned to species level
sum(df2$counts[df2$only_species_level %not% c("Z_other_foraminifera","Z_not_foraminifera")])
#percentage of ESVs assigned to species level
sum(df2$counts[df2$only_species_level %not% c("Z_other_foraminifera","Z_not_foraminifera")])/sum(df2$counts) *100
#percentage of reads ESVs assigned to species level
sum(df$reads[df$only_species_level %not% c("Z_other_foraminifera","Z_not_foraminifera")])/sum(df$reads) *100

#percentage of reads of ESVs assigned to foraminifera
sum(df$reads[df$only_species_level != "Z_not_foraminifera"])/sum(df$reads) *100
#percentage of ESVs assigned to foraminifera
sum(df2$counts[df2$only_species_level != "Z_not_foraminifera"])/sum(df2$counts) *100
#number of ESVs assigned to foraminifera
sum(df2$counts[df2$only_species_level != "Z_not_foraminifera"])

#percentage of reads of ESVs not assigned to species level
sum(df$reads[df$only_species_level == "Z_other_foraminifera"])/sum(df$reads) *100

#to calculate percentage of ESVs assigned to none-LBF

lbf <- c("Sorites_sp2_Spermonde","Sorites_sp1_Spermonde","Peneroplis_sp2_&_Peneroplis_pertusus_5117_&_Dendritina_ambigua_Spermonde",
         "Peneroplis_sp1_Spermonde","Parasorites_sp","Operculina_complanata","Operculina_ammonoides","Nummulites_venosus",                                                                
         "Neorotalia_gaimardi_&_Baculogypsina_sphaerulata_Spermonde","Neorotalia_calcar_Spermonde","Heterostegina_depressa_sp2_Spermonde",                                              
         "Heterostegina_depressa_sp1_Spermonde","Elphidium_Rik_9327_9337_9350_Spermonde","Calcarina_spengleri_Rik_9352_consensus0_1732_Spermonde",                            
         "Calcarina_sp._Spermonde","Calcarina_sp._5164_Spermonde","Calcarina_hispida_Ambon_&_Calcarina_spengleri_&_Calcarina_sp._5247_Spermonde",      
         "Calcarina_hispida_&_Calcarina_sp_5163_Spermonde","Baculogypsinoides_spinosus_8078_8079_8080",                                         
         "Amphistegina_radiata_Spermonde","Amphistegina_papillosa_Spermonde_3741","Amphistegina_lobifera","Amphistegina_lessonii",                                                             
         "Amphistegina_bicirculata_Rik_9379_9380_9381_9382_Eilat","Amphisorus_&_Amphisorus_SpL_&_Amphisorus_SpS_Spermonde")

lbf1 <- c(lbf, "Z_other_foraminifera")
lbf2 <- c(lbf, "Z_other_foraminifera","Z_not_foraminifera")

df3 <- df2[df2$only_species_level %not% c("Z_other_foraminifera","Z_not_foraminifera"),]
df3$foram <- NA
df3$foram[df3$only_species_level %in% lbf] <- "LBF"
df3$foram[df3$only_species_level %not% lbf] <- "small"

df$foram <- NA
df$foram[df$only_species_level %in% lbf] <- "LBF"
df$foram[df$only_species_level %not% lbf1] <- "small"
df$foram[df$only_species_level == "Z_other_foraminifera"] <- "other foram"
df$foram[df$only_species_level == "Z_not_foraminifera"] <- "not foram"

#number of ESVs assigned to LBF species level
sum(df3$counts[df3$foram == "LBF"])
#percentage of ESVs assigned to LBF species level
sum(df3$counts[df3$foram == "LBF"])/sum(df2$counts) *100
#percentage of reads of ESVs assigned to LBF species level
sum(df$reads[df$foram == "LBF"])/sum(df$reads) *100

#percentage of reads of EVS assigned to none-LBF species level
sum(df$reads[df$foram == "small"])/sum(df$reads) *100



#______Bulk samples only________________________
df <- read.csv("~/Downloads/Dataset/workingdataset-quantitativemetabarcoding-994_new.csv")
df <- df[df$type == "bulk",]

df1 <- df[,c(1,14)]
df1$counts <- 1
df2 <- unique(df1)

#total number of ESVs
sum(df2$counts)

#number of ESVs assigned to species level
sum(df2$counts[df2$only_species_level %not% c("Z_other_foraminifera","Z_not_foraminifera")])
#percentage of ESVs assigned to species level
sum(df2$counts[df2$only_species_level %not% c("Z_other_foraminifera","Z_not_foraminifera")])/sum(df2$counts) *100
#percentage of reads ESVs assigned to species level
sum(df$reads[df$only_species_level %not% c("Z_other_foraminifera","Z_not_foraminifera")])/sum(df$reads) *100

#percentage of reads of ESVs assigned to foraminifera
sum(df$reads[df$only_species_level != "Z_not_foraminifera"])/sum(df$reads) *100
#percentage of ESVs assigned to foraminifera
sum(df2$counts[df2$only_species_level != "Z_not_foraminifera"])/sum(df2$counts) *100
#number of ESVs assigned to foraminifera
sum(df2$counts[df2$only_species_level != "Z_not_foraminifera"])

#percentage of reads of ESVs not assigned to species level
sum(df$reads[df$only_species_level == "Z_other_foraminifera"])/sum(df$reads) *100

#to calculate percentage of ESVs assigned to none-LBF

lbf <- c("Sorites_sp2_Spermonde","Sorites_sp1_Spermonde","Peneroplis_sp2_&_Peneroplis_pertusus_5117_&_Dendritina_ambigua_Spermonde",
         "Peneroplis_sp1_Spermonde","Parasorites_sp","Operculina_complanata","Operculina_ammonoides","Nummulites_venosus",                                                                
         "Neorotalia_gaimardi_&_Baculogypsina_sphaerulata_Spermonde","Neorotalia_calcar_Spermonde","Heterostegina_depressa_sp2_Spermonde",                                              
         "Heterostegina_depressa_sp1_Spermonde","Elphidium_Rik_9327_9337_9350_Spermonde","Calcarina_spengleri_Rik_9352_consensus0_1732_Spermonde",                            
         "Calcarina_sp._Spermonde","Calcarina_sp._5164_Spermonde","Calcarina_hispida_Ambon_&_Calcarina_spengleri_&_Calcarina_sp._5247_Spermonde",      
         "Calcarina_hispida_&_Calcarina_sp_5163_Spermonde","Baculogypsinoides_spinosus_8078_8079_8080",                                         
         "Amphistegina_radiata_Spermonde","Amphistegina_papillosa_Spermonde_3741","Amphistegina_lobifera","Amphistegina_lessonii",                                                             
         "Amphistegina_bicirculata_Rik_9379_9380_9381_9382_Eilat","Amphisorus_&_Amphisorus_SpL_&_Amphisorus_SpS_Spermonde")

lbf1 <- c(lbf, "Z_other_foraminifera")
lbf2 <- c(lbf, "Z_other_foraminifera","Z_not_foraminifera")

df3 <- df2[df2$only_species_level %not% c("Z_other_foraminifera","Z_not_foraminifera"),]
df3$foram <- NA
df3$foram[df3$only_species_level %in% lbf] <- "LBF"
df3$foram[df3$only_species_level %not% lbf] <- "small"

df$foram <- NA
df$foram[df$only_species_level %in% lbf] <- "LBF"
df$foram[df$only_species_level %not% lbf1] <- "small"
df$foram[df$only_species_level == "Z_other_foraminifera"] <- "other foram"
df$foram[df$only_species_level == "Z_not_foraminifera"] <- "not foram"

#number of ESVs assigned to LBF species level
sum(df3$counts[df3$foram == "LBF"])
#percentage of ESVs assigned to LBF species level
sum(df3$counts[df3$foram == "LBF"])/sum(df2$counts) *100
#percentage of reads of ESVs assigned to LBF species level
sum(df$reads[df$foram == "LBF"])/sum(df$reads) *100

#percentage of reads of EVS assigned to none-LBF species level
sum(df$reads[df$foram == "small"])/sum(df$reads) *100

#______MOCK samples only________________________
df <- read.csv("~/Downloads/Dataset/workingdataset-quantitativemetabarcoding-994_new.csv")
df <- df[df$type == "picked",]

df1 <- df[,c(1,14)]
df1$counts <- 1
df2 <- unique(df1)

#total number of ESVs
sum(df2$counts)

#number of ESVs assigned to species level
sum(df2$counts[df2$only_species_level %not% c("Z_other_foraminifera","Z_not_foraminifera")])
#percentage of ESVs assigned to species level
sum(df2$counts[df2$only_species_level %not% c("Z_other_foraminifera","Z_not_foraminifera")])/sum(df2$counts) *100
#percentage of reads ESVs assigned to species level
sum(df$reads[df$only_species_level %not% c("Z_other_foraminifera","Z_not_foraminifera")])/sum(df$reads) *100

#percentage of reads of ESVs assigned to foraminifera
sum(df$reads[df$only_species_level != "Z_not_foraminifera"])/sum(df$reads) *100
#percentage of ESVs assigned to foraminifera
sum(df2$counts[df2$only_species_level != "Z_not_foraminifera"])/sum(df2$counts) *100
#number of ESVs assigned to foraminifera
sum(df2$counts[df2$only_species_level != "Z_not_foraminifera"])

#percentage of reads of ESVs not assigned to species level
sum(df$reads[df$only_species_level == "Z_other_foraminifera"])/sum(df$reads) *100

#to calculate percentage of ESVs assigned to none-LBF

lbf <- c("Sorites_sp2_Spermonde","Sorites_sp1_Spermonde","Peneroplis_sp2_&_Peneroplis_pertusus_5117_&_Dendritina_ambigua_Spermonde",
         "Peneroplis_sp1_Spermonde","Parasorites_sp","Operculina_complanata","Operculina_ammonoides","Nummulites_venosus",                                                                
         "Neorotalia_gaimardi_&_Baculogypsina_sphaerulata_Spermonde","Neorotalia_calcar_Spermonde","Heterostegina_depressa_sp2_Spermonde",                                              
         "Heterostegina_depressa_sp1_Spermonde","Elphidium_Rik_9327_9337_9350_Spermonde","Calcarina_spengleri_Rik_9352_consensus0_1732_Spermonde",                            
         "Calcarina_sp._Spermonde","Calcarina_sp._5164_Spermonde","Calcarina_hispida_Ambon_&_Calcarina_spengleri_&_Calcarina_sp._5247_Spermonde",      
         "Calcarina_hispida_&_Calcarina_sp_5163_Spermonde","Baculogypsinoides_spinosus_8078_8079_8080",                                         
         "Amphistegina_radiata_Spermonde","Amphistegina_papillosa_Spermonde_3741","Amphistegina_lobifera","Amphistegina_lessonii",                                                             
         "Amphistegina_bicirculata_Rik_9379_9380_9381_9382_Eilat","Amphisorus_&_Amphisorus_SpL_&_Amphisorus_SpS_Spermonde")

lbf1 <- c(lbf, "Z_other_foraminifera")
lbf2 <- c(lbf, "Z_other_foraminifera","Z_not_foraminifera")

df3 <- df2[df2$only_species_level %not% c("Z_other_foraminifera","Z_not_foraminifera"),]
df3$foram <- NA
df3$foram[df3$only_species_level %in% lbf] <- "LBF"
df3$foram[df3$only_species_level %not% lbf] <- "small"

df$foram <- NA
df$foram[df$only_species_level %in% lbf] <- "LBF"
df$foram[df$only_species_level %not% lbf1] <- "small"
df$foram[df$only_species_level == "Z_other_foraminifera"] <- "other foram"
df$foram[df$only_species_level == "Z_not_foraminifera"] <- "not foram"

#number of ESVs assigned to LBF species level
sum(df3$counts[df3$foram == "LBF"])
#percentage of ESVs assigned to LBF species level
sum(df3$counts[df3$foram == "LBF"])/sum(df2$counts) *100
#percentage of reads of ESVs assigned to LBF species level
sum(df$reads[df$foram == "LBF"])/sum(df$reads) *100

#percentage of reads of EVS assigned to none-LBF species level
sum(df$reads[df$foram == "small"])/sum(df$reads) *100

#--------NEW - prep data - Keep ASVs that are at least in 2 of the 3 replicates-------------------------------

df <- read.csv("~/Downloads/Dataset/workingdataset-quantitativemetabarcoding-994_new.csv")

df1 <- df
df1$pa[df1$reads > 0] <- 1
df1$pa[df1$reads == 0] <- 0

df2 <- df1 %>% group_by(ESVs, field.nmbr., type, species, island) %>% summarize(counts = sum(pa)) 

df3 <- df2[df2$counts > 1, ]

df4 <- merge(df3[,1:3], df1, by = c("ESVs", "field.nmbr.","type"))
df5 <- df4[,1:(ncol(df4)-1)]

new_df <- df5

new_df$depth <- NA

new_df$depth[new_df$field.nmbr. == "UPG-EG22-021"] <- 17.4
new_df$depth[new_df$field.nmbr. == "UPG-EG22-022"] <- 8.1
new_df$depth[new_df$field.nmbr. == "UPG-EG22-023"] <- 14.6
new_df$depth[new_df$field.nmbr. == "UPG-EG22-024"] <- 27.4
new_df$depth[new_df$field.nmbr. == "UPG-EG22-025"] <- 5.5
new_df$depth[new_df$field.nmbr. == "UPG-EG22-026"] <- 30.5
new_df$depth[new_df$field.nmbr. == "UPG-EG22-027"] <- 1.5
new_df$depth[new_df$field.nmbr. == "UPG-EG22-028"] <- 1
new_df$depth[new_df$field.nmbr. == "UPG-EG22-029"] <- 11.3
new_df$depth[new_df$field.nmbr. == "UPG-EG22-030"] <- 21.7
new_df$depth[new_df$field.nmbr. == "UPG-EG22-031"] <- 2
new_df$depth[new_df$field.nmbr. == "UPG-EG22-032"] <- 2.5
new_df$depth[new_df$field.nmbr. == "UPG-EG22-033"] <- 24.5

new_df$depth[new_df$field.nmbr. == "UPG-EG22-081"] <- 8.7
new_df$depth[new_df$field.nmbr. == "UPG-EG22-082"] <- 1
new_df$depth[new_df$field.nmbr. == "UPG-EG22-083"] <- 6.6
new_df$depth[new_df$field.nmbr. == "UPG-EG22-084"] <- 2
new_df$depth[new_df$field.nmbr. == "UPG-EG22-085"] <- 11.4
new_df$depth[new_df$field.nmbr. == "UPG-EG22-086"] <- 17.6
new_df$depth[new_df$field.nmbr. == "UPG-EG22-087"] <- 1.5
new_df$depth[new_df$field.nmbr. == "UPG-EG22-088"] <- 14.8
new_df$depth[new_df$field.nmbr. == "UPG-EG22-089"] <- 23.2
new_df$depth[new_df$field.nmbr. == "UPG-EG22-090"] <- 20.2
new_df$depth[new_df$field.nmbr. == "UPG-EG22-091"] <- 25.9

new_df$depth[new_df$field.nmbr. == "UPG-EG22-098"] <- 7.1
new_df$depth[new_df$field.nmbr. == "UPG-EG22-099"] <- 16.9
new_df$depth[new_df$field.nmbr. == "UPG-EG22-100"] <- 19.5
new_df$depth[new_df$field.nmbr. == "UPG-EG22-101"] <- 2
new_df$depth[new_df$field.nmbr. == "UPG-EG22-102"] <- 13.6
new_df$depth[new_df$field.nmbr. == "UPG-EG22-103"] <- 1
new_df$depth[new_df$field.nmbr. == "UPG-EG22-104"] <- 25.6
new_df$depth[new_df$field.nmbr. == "UPG-EG22-106"] <- 1.5
new_df$depth[new_df$field.nmbr. == "UPG-EG22-107"] <- 10
new_df$depth[new_df$field.nmbr. == "UPG-EG22-108"] <- 3.5
new_df$depth[new_df$field.nmbr. == "UPG-EG22-109"] <- 22.6

write.csv(new_df, "~/Downloads/Dataset/workingdataset-quantitativemetabarcoding-994_new_replicatefiltered.csv")


#--------NEW - prep data - Create genus level dataset with replicates separated---------------

df <- read.csv("~/Downloads/Dataset/workingdataset-quantitativemetabarcoding-994_new_replicatefiltered.csv")
morpho <- read.csv("~/Downloads/Dataset/melted_count_surface_picked_forams_corrected.csv")

df1 <- df

lbf <- c("Sorites_sp2_Spermonde","Sorites_sp1_Spermonde","Peneroplis_sp2_&_Peneroplis_pertusus_5117_&_Dendritina_ambigua_Spermonde",
         "Peneroplis_sp1_Spermonde","Parasorites_sp","Operculina_complanata","Operculina_ammonoides","Nummulites_venosus",                                                                
         "Neorotalia_gaimardi_&_Baculogypsina_sphaerulata_Spermonde","Neorotalia_calcar_Spermonde","Heterostegina_depressa_sp2_Spermonde",                                              
         "Heterostegina_depressa_sp1_Spermonde","Elphidium_Rik_9327_9337_9350_Spermonde","Calcarina_spengleri_Rik_9352_consensus0_1732_Spermonde",                            
         "Calcarina_sp._Spermonde","Calcarina_sp._5164_Spermonde","Calcarina_hispida_Ambon_&_Calcarina_spengleri_&_Calcarina_sp._5247_Spermonde",      
         "Calcarina_hispida_&_Calcarina_sp_5163_Spermonde","Baculogypsinoides_spinosus_8078_8079_8080",                                         
         "Amphistegina_radiata_Spermonde","Amphistegina_papillosa_Spermonde_3741","Amphistegina_lobifera","Amphistegina_lessonii",                                                             
         "Amphistegina_bicirculata_Rik_9379_9380_9381_9382_Eilat","Amphisorus_&_Amphisorus_SpL_&_Amphisorus_SpS_Spermonde")

df1$foram <- NA
df1$foram[df1$only_species_level %in% lbf] <- "LBF"
df1 <- df1[df1$foram == "LBF",] 
df1 <- as.data.frame(na.omit(df1))

df1$morpho_species <- NA

df1$morpho_species[df1$only_species_level == "Sorites_sp2_Spermonde"] <- "Sorites sp."
df1$morpho_species[df1$only_species_level == "Sorites_sp1_Spermonde"] <- "Sorites sp."
df1$morpho_species[df1$only_species_level == "Peneroplis_sp2_&_Peneroplis_pertusus_5117_&_Dendritina_ambigua_Spermonde"] <- "Peneroplis sp."
df1$morpho_species[df1$only_species_level == "Peneroplis_sp1_Spermonde"] <- "Peneroplis sp."
df1$morpho_species[df1$only_species_level == "Parasorites_sp"] <- "Parasorites sp."
df1$morpho_species[df1$only_species_level == "Operculina_complanata"] <- "Operculina complanata"
df1$morpho_species[df1$only_species_level == "Operculina_ammonoides"] <- "Operculina ammonoides"
df1$morpho_species[df1$only_species_level == "Nummulites_venosus"] <- "Nummulites venosus"
df1$morpho_species[df1$only_species_level == "Neorotalia_gaimardi_&_Baculogypsina_sphaerulata_Spermonde"] <- "Neorotalia gaimardi"
df1$morpho_species[df1$only_species_level == "Neorotalia_calcar_Spermonde"] <- "Neorotalia calcar"
df1$morpho_species[df1$only_species_level == "Heterostegina_depressa_sp2_Spermonde"] <- "Heterostegina depressa"
df1$morpho_species[df1$only_species_level == "Heterostegina_depressa_sp1_Spermonde"] <- "Heterostegina depressa"
df1$morpho_species[df1$only_species_level == "Elphidium_Rik_9327_9337_9350_Spermonde"] <- "Elphidium sp."
df1$morpho_species[df1$only_species_level == "Calcarina_spengleri_Rik_9352_consensus0_1732_Spermonde"] <- "Calcarina spengleri"
df1$morpho_species[df1$only_species_level == "Calcarina_sp._Spermonde"] <- "Calcarina sp."
df1$morpho_species[df1$only_species_level == "Calcarina_sp._5164_Spermonde"] <- "Calcarina sp."
df1$morpho_species[df1$only_species_level == "Calcarina_hispida_Ambon_&_Calcarina_spengleri_&_Calcarina_sp._5247_Spermonde"] <- "Calcarina sp."
df1$morpho_species[df1$only_species_level == "Calcarina_hispida_&_Calcarina_sp_5163_Spermonde"] <- "Calcarina hispida"
df1$morpho_species[df1$only_species_level == "Borelis_schlumbergeri_3737_Spermonde"] <- "Borelis schlumbergeri"
df1$morpho_species[df1$only_species_level == "Baculogypsinoides_spinosus_8078_8079_8080"] <- "Baculogypsinoides spinosus"
df1$morpho_species[df1$only_species_level == "Amphistegina_radiata_Spermonde"] <- "Amphistegina radiata"
df1$morpho_species[df1$only_species_level == "Amphistegina_papillosa_Spermonde_3741"] <- "Amphistegina papillosa"
df1$morpho_species[df1$only_species_level == "Amphistegina_lobifera"] <- "Amphistegina lobifera"
df1$morpho_species[df1$only_species_level == "Amphistegina_lessonii"] <- "Amphistegina lessonii"
df1$morpho_species[df1$only_species_level == "Amphistegina_bicirculata_Rik_9379_9380_9381_9382_Eilat"] <- "Amphistegina bicirculata"
df1$morpho_species[df1$only_species_level == "Amphisorus_&_Amphisorus_SpL_&_Amphisorus_SpS_Spermonde"] <- "Amphisorus sp."

splitspecies <- str_split_fixed(df1$morpho_species, " ", 2)
colnames(splitspecies) <- c("Genus", "Species_alone")
df2 <- cbind(df1, splitspecies)

splitname <- as.data.frame(str_split_fixed(df2$new.name, "_", 3))
colnames(splitname) <- c("Island_letter", "field", "type_replicate")
df2 <- cbind(df2, splitname)

df3 <- df2 %>% group_by(type, island, field.nmbr., Genus, morpho_species, species, type_replicate) %>% summarise(addreads = sum(reads))
df4 <- df3 %>% group_by(type, island, field.nmbr., Genus, type_replicate) %>% summarise(sumreads = sum(addreads))
df5 <- dcast(df4, island + field.nmbr. + Genus ~ 
               type_replicate, fun.aggregate = sum, value.var = "sumreads")
df5$merge <- paste(df5$field.nmbr.,"_",df5$Genus)

morpho$type <- "picked"
colnames(morpho) <- c("field.nmbr","registr.nmbr","Genus","morpho_species", "surface_area_um2","counts", "type")
morpho1 <- morpho %>% group_by(type, field.nmbr, Genus) %>% summarise(totalarea = sum(surface_area_um2),totalcounts = sum(counts))
morpho1 <- na.omit(morpho1)
morpho1$merge <- paste(morpho1$field.nmbr,"_",morpho1$Genus)

new_df <- merge(df5, morpho1[,4:6], by = "merge")

new_df$depth <- NA

new_df$depth[new_df$field.nmbr. == "UPG-EG22-021"] <- 17.4
new_df$depth[new_df$field.nmbr. == "UPG-EG22-022"] <- 8.1
new_df$depth[new_df$field.nmbr. == "UPG-EG22-023"] <- 14.6
new_df$depth[new_df$field.nmbr. == "UPG-EG22-024"] <- 27.4
new_df$depth[new_df$field.nmbr. == "UPG-EG22-025"] <- 5.5
new_df$depth[new_df$field.nmbr. == "UPG-EG22-026"] <- 30.5
new_df$depth[new_df$field.nmbr. == "UPG-EG22-027"] <- 1.5
new_df$depth[new_df$field.nmbr. == "UPG-EG22-028"] <- 1
new_df$depth[new_df$field.nmbr. == "UPG-EG22-029"] <- 11.3
new_df$depth[new_df$field.nmbr. == "UPG-EG22-030"] <- 21.7
new_df$depth[new_df$field.nmbr. == "UPG-EG22-031"] <- 2
new_df$depth[new_df$field.nmbr. == "UPG-EG22-032"] <- 2.5
new_df$depth[new_df$field.nmbr. == "UPG-EG22-033"] <- 24.5

new_df$depth[new_df$field.nmbr. == "UPG-EG22-081"] <- 8.7
new_df$depth[new_df$field.nmbr. == "UPG-EG22-082"] <- 1
new_df$depth[new_df$field.nmbr. == "UPG-EG22-083"] <- 6.6
new_df$depth[new_df$field.nmbr. == "UPG-EG22-084"] <- 2
new_df$depth[new_df$field.nmbr. == "UPG-EG22-085"] <- 11.4
new_df$depth[new_df$field.nmbr. == "UPG-EG22-086"] <- 17.6
new_df$depth[new_df$field.nmbr. == "UPG-EG22-087"] <- 1.5
new_df$depth[new_df$field.nmbr. == "UPG-EG22-088"] <- 14.8
new_df$depth[new_df$field.nmbr. == "UPG-EG22-089"] <- 23.2
new_df$depth[new_df$field.nmbr. == "UPG-EG22-090"] <- 20.2
new_df$depth[new_df$field.nmbr. == "UPG-EG22-091"] <- 25.9

new_df$depth[new_df$field.nmbr. == "UPG-EG22-098"] <- 7.1
new_df$depth[new_df$field.nmbr. == "UPG-EG22-099"] <- 16.9
new_df$depth[new_df$field.nmbr. == "UPG-EG22-100"] <- 19.5
new_df$depth[new_df$field.nmbr. == "UPG-EG22-101"] <- 2
new_df$depth[new_df$field.nmbr. == "UPG-EG22-102"] <- 13.6
new_df$depth[new_df$field.nmbr. == "UPG-EG22-103"] <- 1
new_df$depth[new_df$field.nmbr. == "UPG-EG22-104"] <- 25.6
new_df$depth[new_df$field.nmbr. == "UPG-EG22-106"] <- 1.5
new_df$depth[new_df$field.nmbr. == "UPG-EG22-107"] <- 10
new_df$depth[new_df$field.nmbr. == "UPG-EG22-108"] <- 3.5
new_df$depth[new_df$field.nmbr. == "UPG-EG22-109"] <- 22.6

new_df2 <- new_df %>% group_by(field.nmbr.) %>% 
  mutate(proparea = totalarea/sum(totalarea),propcount = totalcounts/sum(totalcounts), 
         propbulk_A = bulk_A/sum(bulk_A), 
         propbulk_B = bulk_B/sum(bulk_B),
         propbulk_C = bulk_C/sum(bulk_C),
         proppick_A = picked_A/sum(picked_A),
         proppick_B = picked_B/sum(picked_B),
         proppick_C = picked_C/sum(picked_C))

new_df2$island[new_df2$island == "Lumulumu"] <- "Lumu-lumu"
new_df2$island[new_df2$island == "Padjenekang"] <- "Pajenekang"

write.csv(new_df2, "~/Downloads/Dataset/Genuslevel_counts_reads_area_replicate_new.csv", row.names=FALSE)


#--------NEW - prep data - Create genus level dataset with replicates merged---------------

df <- read.csv("~/Downloads/Dataset/workingdataset-quantitativemetabarcoding-994_new_replicatefiltered.csv")
morpho <- read.csv("~/Downloads/Dataset/melted_count_surface_picked_forams_corrected.csv")

df1 <- df

lbf <- c("Sorites_sp2_Spermonde","Sorites_sp1_Spermonde","Peneroplis_sp2_&_Peneroplis_pertusus_5117_&_Dendritina_ambigua_Spermonde",
         "Peneroplis_sp1_Spermonde","Parasorites_sp","Operculina_complanata","Operculina_ammonoides","Nummulites_venosus",                                                                
         "Neorotalia_gaimardi_&_Baculogypsina_sphaerulata_Spermonde","Neorotalia_calcar_Spermonde","Heterostegina_depressa_sp2_Spermonde",                                              
         "Heterostegina_depressa_sp1_Spermonde","Elphidium_Rik_9327_9337_9350_Spermonde","Calcarina_spengleri_Rik_9352_consensus0_1732_Spermonde",                            
         "Calcarina_sp._Spermonde","Calcarina_sp._5164_Spermonde","Calcarina_hispida_Ambon_&_Calcarina_spengleri_&_Calcarina_sp._5247_Spermonde",      
         "Calcarina_hispida_&_Calcarina_sp_5163_Spermonde","Baculogypsinoides_spinosus_8078_8079_8080",                                         
         "Amphistegina_radiata_Spermonde","Amphistegina_papillosa_Spermonde_3741","Amphistegina_lobifera","Amphistegina_lessonii",                                                             
         "Amphistegina_bicirculata_Rik_9379_9380_9381_9382_Eilat","Amphisorus_&_Amphisorus_SpL_&_Amphisorus_SpS_Spermonde")

df1$foram <- NA
df1$foram[df1$only_species_level %in% lbf] <- "LBF"
df1 <- df1[df1$foram == "LBF",] 
df1 <- as.data.frame(na.omit(df1))

df1$morpho_species <- NA

df1$morpho_species[df1$only_species_level == "Sorites_sp2_Spermonde"] <- "Sorites sp."
df1$morpho_species[df1$only_species_level == "Sorites_sp1_Spermonde"] <- "Sorites sp."
df1$morpho_species[df1$only_species_level == "Peneroplis_sp2_&_Peneroplis_pertusus_5117_&_Dendritina_ambigua_Spermonde"] <- "Peneroplis sp."
df1$morpho_species[df1$only_species_level == "Peneroplis_sp1_Spermonde"] <- "Peneroplis sp."
df1$morpho_species[df1$only_species_level == "Parasorites_sp"] <- "Parasorites sp."
df1$morpho_species[df1$only_species_level == "Operculina_complanata"] <- "Operculina complanata"
df1$morpho_species[df1$only_species_level == "Operculina_ammonoides"] <- "Operculina ammonoides"
df1$morpho_species[df1$only_species_level == "Nummulites_venosus"] <- "Nummulites venosus"
df1$morpho_species[df1$only_species_level == "Neorotalia_gaimardi_&_Baculogypsina_sphaerulata_Spermonde"] <- "Neorotalia gaimardi"
df1$morpho_species[df1$only_species_level == "Neorotalia_calcar_Spermonde"] <- "Neorotalia calcar"
df1$morpho_species[df1$only_species_level == "Heterostegina_depressa_sp2_Spermonde"] <- "Heterostegina depressa"
df1$morpho_species[df1$only_species_level == "Heterostegina_depressa_sp1_Spermonde"] <- "Heterostegina depressa"
df1$morpho_species[df1$only_species_level == "Elphidium_Rik_9327_9337_9350_Spermonde"] <- "Elphidium sp."
df1$morpho_species[df1$only_species_level == "Calcarina_spengleri_Rik_9352_consensus0_1732_Spermonde"] <- "Calcarina spengleri"
df1$morpho_species[df1$only_species_level == "Calcarina_sp._Spermonde"] <- "Calcarina sp."
df1$morpho_species[df1$only_species_level == "Calcarina_sp._5164_Spermonde"] <- "Calcarina sp."
df1$morpho_species[df1$only_species_level == "Calcarina_hispida_Ambon_&_Calcarina_spengleri_&_Calcarina_sp._5247_Spermonde"] <- "Calcarina sp."
df1$morpho_species[df1$only_species_level == "Calcarina_hispida_&_Calcarina_sp_5163_Spermonde"] <- "Calcarina hispida"
df1$morpho_species[df1$only_species_level == "Borelis_schlumbergeri_3737_Spermonde"] <- "Borelis schlumbergeri"
df1$morpho_species[df1$only_species_level == "Baculogypsinoides_spinosus_8078_8079_8080"] <- "Baculogypsinoides spinosus"
df1$morpho_species[df1$only_species_level == "Amphistegina_radiata_Spermonde"] <- "Amphistegina radiata"
df1$morpho_species[df1$only_species_level == "Amphistegina_papillosa_Spermonde_3741"] <- "Amphistegina papillosa"
df1$morpho_species[df1$only_species_level == "Amphistegina_lobifera"] <- "Amphistegina lobifera"
df1$morpho_species[df1$only_species_level == "Amphistegina_lessonii"] <- "Amphistegina lessonii"
df1$morpho_species[df1$only_species_level == "Amphistegina_bicirculata_Rik_9379_9380_9381_9382_Eilat"] <- "Amphistegina bicirculata"
df1$morpho_species[df1$only_species_level == "Amphisorus_&_Amphisorus_SpL_&_Amphisorus_SpS_Spermonde"] <- "Amphisorus sp."

splitspecies <- str_split_fixed(df1$morpho_species, " ", 2)
colnames(splitspecies) <- c("Genus", "Species_alone")
df2 <- cbind(df1, splitspecies)

splitname <- as.data.frame(str_split_fixed(df2$new.name, "_", 3))
colnames(splitname) <- c("Island_letter", "field", "type_replicate")
df2 <- cbind(df2, splitname)

df3 <- df2 %>% group_by(type, island, field.nmbr., Genus, morpho_species, species, type_replicate) %>% summarise(addreads = sum(reads))
df4 <- df3 %>% group_by(type, island, field.nmbr., Genus, type_replicate) %>% summarise(sumreads = sum(addreads))
df4b <- df4 %>% group_by(type, island, field.nmbr., Genus) %>% summarise(meanreads = mean(sumreads))
df5 <- dcast(df4b, island + field.nmbr. + Genus ~ 
               type, fun.aggregate = sum, value.var = "meanreads")
df5$merge <- paste(df5$field.nmbr.,"_",df5$Genus)

morpho$type <- "picked"
colnames(morpho) <- c("field.nmbr","registr.nmbr","Genus","morpho_species", "surface_area_um2","counts", "type")
morpho1 <- morpho %>% group_by(type, field.nmbr, Genus) %>% summarise(totalarea = sum(surface_area_um2),totalcounts = sum(counts))
morpho1 <- na.omit(morpho1)
morpho1$merge <- paste(morpho1$field.nmbr,"_",morpho1$Genus)

new_df <- merge(df5, morpho1[,4:6], by = "merge")

#add sampling depth to the data
new_df$depth <- NA

new_df$depth[new_df$field.nmbr. == "UPG-EG22-021"] <- 17.4
new_df$depth[new_df$field.nmbr. == "UPG-EG22-022"] <- 8.1
new_df$depth[new_df$field.nmbr. == "UPG-EG22-023"] <- 14.6
new_df$depth[new_df$field.nmbr. == "UPG-EG22-024"] <- 27.4
new_df$depth[new_df$field.nmbr. == "UPG-EG22-025"] <- 5.5
new_df$depth[new_df$field.nmbr. == "UPG-EG22-026"] <- 30.5
new_df$depth[new_df$field.nmbr. == "UPG-EG22-027"] <- 1.5
new_df$depth[new_df$field.nmbr. == "UPG-EG22-028"] <- 1
new_df$depth[new_df$field.nmbr. == "UPG-EG22-029"] <- 11.3
new_df$depth[new_df$field.nmbr. == "UPG-EG22-030"] <- 21.7
new_df$depth[new_df$field.nmbr. == "UPG-EG22-031"] <- 2
new_df$depth[new_df$field.nmbr. == "UPG-EG22-032"] <- 2.5
new_df$depth[new_df$field.nmbr. == "UPG-EG22-033"] <- 24.5

new_df$depth[new_df$field.nmbr. == "UPG-EG22-081"] <- 8.7
new_df$depth[new_df$field.nmbr. == "UPG-EG22-082"] <- 1
new_df$depth[new_df$field.nmbr. == "UPG-EG22-083"] <- 6.6
new_df$depth[new_df$field.nmbr. == "UPG-EG22-084"] <- 2
new_df$depth[new_df$field.nmbr. == "UPG-EG22-085"] <- 11.4
new_df$depth[new_df$field.nmbr. == "UPG-EG22-086"] <- 17.6
new_df$depth[new_df$field.nmbr. == "UPG-EG22-087"] <- 1.5
new_df$depth[new_df$field.nmbr. == "UPG-EG22-088"] <- 14.8
new_df$depth[new_df$field.nmbr. == "UPG-EG22-089"] <- 23.2
new_df$depth[new_df$field.nmbr. == "UPG-EG22-090"] <- 20.2
new_df$depth[new_df$field.nmbr. == "UPG-EG22-091"] <- 25.9

new_df$depth[new_df$field.nmbr. == "UPG-EG22-098"] <- 7.1
new_df$depth[new_df$field.nmbr. == "UPG-EG22-099"] <- 16.9
new_df$depth[new_df$field.nmbr. == "UPG-EG22-100"] <- 19.5
new_df$depth[new_df$field.nmbr. == "UPG-EG22-101"] <- 2
new_df$depth[new_df$field.nmbr. == "UPG-EG22-102"] <- 13.6
new_df$depth[new_df$field.nmbr. == "UPG-EG22-103"] <- 1
new_df$depth[new_df$field.nmbr. == "UPG-EG22-104"] <- 25.6
new_df$depth[new_df$field.nmbr. == "UPG-EG22-106"] <- 1.5
new_df$depth[new_df$field.nmbr. == "UPG-EG22-107"] <- 10
new_df$depth[new_df$field.nmbr. == "UPG-EG22-108"] <- 3.5
new_df$depth[new_df$field.nmbr. == "UPG-EG22-109"] <- 22.6

#create proportional/relative values for all data types
new_df2 <- new_df %>% group_by(field.nmbr.) %>% 
  mutate(proparea = totalarea/sum(totalarea), propbulk = bulk/sum(bulk), proppick = picked/sum(picked), propcount = totalcounts/sum(totalcounts))

#correct island names
new_df2$island[new_df2$island == "Lumulumu"] <- "Lumu-lumu"
new_df2$island[new_df2$island == "Padjenekang"] <- "Pajenekang"

write.csv(new_df2, "~/Downloads/Dataset/Genuslevel_counts_reads_area_merged_new.csv", row.names=FALSE)

#--------fig5 - Proportion correction between real surface area vs reads (linear, power law) on merged data--------------------
df <- read.csv("~/Downloads/Dataset/Genuslevel_counts_reads_area_merged_new.csv")

df <- df[, c(1:9)]
df$picked[df$totalarea == 0] <- 0

df$area_mm2 <- df$totalarea/1000000
df$picked_million <- df$picked/1000000

#from linear regression of ddPCR gene copy number by surface area (lm) - mean from 999 permutations
permut <- read.csv("~/Downloads/Dataset/ddpcr_powerlaw_linear_coefficients_permutations999.csv")

df$density_factor_ddPCR[df$Genus == "Amphisorus"] <- permut$mean_lma[permut$species == "Amphisorus SpL"]
df$density_factor_ddPCR[df$Genus == "Amphistegina"] <- permut$mean_lma[permut$species == "Amphistegina lessonii"]
df$density_factor_ddPCR[df$Genus == "Baculogypsinoides"] <- permut$mean_lma[permut$species == "Baculogypsinoides spinosus"]
df$density_factor_ddPCR[df$Genus == "Calcarina"] <- permut$mean_lma[permut$species == "Calcarina spengleri"]
df$density_factor_ddPCR[df$Genus == "Heterostegina"] <- permut$mean_lma[permut$species == "Heterostegina depressa"]
df$density_factor_ddPCR[df$Genus == "Neorotalia"] <- permut$mean_lma[permut$species == "Neorotalia gaimardi"]
df$density_factor_ddPCR[df$Genus == "Operculina"] <- permut$mean_lma[permut$species == "Operculina ammonoides"]
df$density_factor_ddPCR[df$Genus == "Elphidium"] <- permut$mean_lma[permut$species == "Calcarina spengleri"]
df$density_factor_ddPCR[df$Genus == "Nummulites"] <- permut$mean_lma[permut$species == "Heterostegina depressa"]
df$density_factor_ddPCR[df$Genus == "Peneroplis"] <- permut$mean_lma[permut$species == "Amphisorus SpL"]
df$density_factor_ddPCR[df$Genus == "Sorites"] <- permut$mean_lma[permut$species == "Amphisorus SpL"]

#from powerlaw regression of ddPCR gene copy number by surface area (log) - mean from 999 permutations

df$power_a_ddPCR[df$Genus == "Amphisorus"] <- permut$mean_a[permut$species == "Amphisorus SpL"]
df$power_a_ddPCR[df$Genus == "Amphistegina"] <- permut$mean_a[permut$species == "Amphistegina lessonii"]
df$power_a_ddPCR[df$Genus == "Baculogypsinoides"] <- permut$mean_a[permut$species == "Baculogypsinoides spinosus"]
df$power_a_ddPCR[df$Genus == "Calcarina"] <- permut$mean_a[permut$species == "Calcarina spengleri"]
df$power_a_ddPCR[df$Genus == "Heterostegina"] <- permut$mean_a[permut$species == "Heterostegina depressa"]
df$power_a_ddPCR[df$Genus == "Neorotalia"] <- permut$mean_a[permut$species == "Neorotalia gaimardi"]
df$power_a_ddPCR[df$Genus == "Operculina"] <- permut$mean_a[permut$species == "Operculina ammonoides"]
df$power_a_ddPCR[df$Genus == "Elphidium"] <- permut$mean_a[permut$species == "Calcarina spengleri"]
df$power_a_ddPCR[df$Genus == "Nummulites"] <- permut$mean_a[permut$species == "Heterostegina depressa"]
df$power_a_ddPCR[df$Genus == "Peneroplis"] <- permut$mean_a[permut$species == "Amphisorus SpL"]
df$power_a_ddPCR[df$Genus == "Sorites"] <- permut$mean_a[permut$species == "Amphisorus SpL"]

df$power_b_ddPCR[df$Genus == "Amphisorus"] <- permut$mean_b[permut$species == "Amphisorus SpL"]
df$power_b_ddPCR[df$Genus == "Amphistegina"] <- permut$mean_b[permut$species == "Amphistegina lessonii"]
df$power_b_ddPCR[df$Genus == "Baculogypsinoides"] <- permut$mean_b[permut$species == "Baculogypsinoides spinosus"]
df$power_b_ddPCR[df$Genus == "Calcarina"] <- permut$mean_b[permut$species == "Calcarina spengleri"]
df$power_b_ddPCR[df$Genus == "Heterostegina"] <- permut$mean_b[permut$species == "Heterostegina depressa"]
df$power_b_ddPCR[df$Genus == "Neorotalia"] <- permut$mean_b[permut$species == "Neorotalia gaimardi"]
df$power_b_ddPCR[df$Genus == "Operculina"] <- permut$mean_b[permut$species == "Operculina ammonoides"]
df$power_b_ddPCR[df$Genus == "Elphidium"] <- permut$mean_b[permut$species == "Calcarina spengleri"]
df$power_b_ddPCR[df$Genus == "Nummulites"] <- permut$mean_b[permut$species == "Heterostegina depressa"]
df$power_b_ddPCR[df$Genus == "Peneroplis"] <- permut$mean_b[permut$species == "Amphisorus SpL"]
df$power_b_ddPCR[df$Genus == "Sorites"] <- permut$mean_b[permut$species == "Amphisorus SpL"]

#optimized correction factor by iteration

df$optimized_factor[df$Genus == "Amphisorus"] <- 0.0968
df$optimized_factor[df$Genus == "Amphistegina"] <- 0.8061
df$optimized_factor[df$Genus == "Baculogypsinoides"] <- 0.0193
df$optimized_factor[df$Genus == "Calcarina"] <- 2.2747
df$optimized_factor[df$Genus == "Elphidium"] <- 1.0388
df$optimized_factor[df$Genus == "Heterostegina"] <- 2.8996
df$optimized_factor[df$Genus == "Neorotalia"] <- 1.9349
df$optimized_factor[df$Genus == "Nummulites"] <- 7.2478
df$optimized_factor[df$Genus == "Operculina"] <- 0.2284
df$optimized_factor[df$Genus == "Peneroplis"] <- 2.5285
df$optimized_factor[df$Genus == "Sorites"] <- 0.1521

#estimated surface area based on linear and log regression from ddPCR results

#linear : reads = density * area
df$estimated_area_lm <- df$picked_million / df$density_factor_ddPCR

#powerlaw: read = a*area^b  ---> log(read) = a*log(area)+ b ---> area = exp((log(read) - b)/a)
df$estimated_area_log <- exp((log(df$picked_million) - df$power_b_ddPCR)/df$power_a_ddPCR)

#optimized by iteration
df$estimated_area_op <- df$picked_million / df$optimized_factor

new_df2 <- df %>% group_by(field.nmbr.) %>% 
  mutate(proparea = area_mm2/sum(area_mm2),
         propcounts = totalcounts/sum(totalcounts),
         propreads_raw = picked_million/sum(picked_million),
         propreads_correctlm = estimated_area_lm/sum(estimated_area_lm),
         propreads_correctlog = estimated_area_log/sum(estimated_area_log),
         propreads_correctop = estimated_area_op/sum(estimated_area_op))

new_df2$raw_area <- new_df2$propreads_raw - new_df2$proparea
new_df2$lm_area <- new_df2$propreads_correctlm - new_df2$proparea
new_df2$log_area <- new_df2$propreads_correctlog - new_df2$proparea
new_df2$op_area <- new_df2$propreads_correctop - new_df2$proparea

new_df2$raw_count <- new_df2$propreads_raw - new_df2$propcounts
new_df2$lm_count <- new_df2$propreads_correctlm - new_df2$propcounts
new_df2$log_count <- new_df2$propreads_correctlog - new_df2$propcounts
new_df2$op_count <- new_df2$propreads_correctop - new_df2$propcounts

#with proportions (for bar charts)
df1 <- new_df2[,c(1:4, 9, 19:23)]
df2 <- melt(df1, id.vars = c("Genus","merge","field.nmbr.","island","depth"), 
            value.name = "value", variable.name = "category")

#with differences (for boxplot)
df3area <- new_df2[,c(1:4, 9, 25:27)]
df4area <- melt(df3area, id.vars = c("Genus","merge","field.nmbr.","island","depth"), 
            value.name = "value", variable.name = "category")
df3count <- new_df2[,c(1:4, 9, 29:31)]
df4count <- melt(df3count, id.vars = c("Genus","merge","field.nmbr.","island","depth"), 
                value.name = "value", variable.name = "category")

#to calculate mean and sd of the proportion difference between (corrected) reads and space
df5 <- df4area 
df5$value[df5$value == 0] <- NA
df6 <- df5 %>% group_by(Genus, category) %>% summarize(meanvalue = round(mean(value*100,na.rm=T), 2), sdvalue = round(sd(value*100,na.rm=T),2))

df7 <- df4count 
df7$value[df7$value == 0] <- NA
df8 <- df7 %>% group_by(Genus, category) %>% summarize(meanvalue = round(mean(value*100,na.rm=T), 2), sdvalue = round(sd(value*100,na.rm=T),2))

df_stat <- cbind(df6, df8)


#different in proportions between reads and area, including first correction from density

p1 <- ggplot(df4area[df4area$category != "log_area",], aes(x=category, y=value*100, col = Genus)) +
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = -5, ymax = 5,
           fill = "black", alpha = 0.1) +
  geom_abline(slope = 0, col = "gray") +
  geom_half_boxplot() +
  geom_half_point(aes(fill = Genus), alpha = 0.3) +
  geom_signif(col = "black", comparisons = list(c("raw_area", "lm_area")), 
              map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05), test = "t.test", 
              y_position = 50, tip_length = 0.01, vjust = .1) +
  facet_grid(.~Genus, scales = "free") + 
  scale_color_manual(values =c("Amphisorus" = "#00B6EB","Amphistegina" = "#FB61D7","Baculogypsinoides" = "#F8766D",
                               "Calcarina" = "#E6194B","Elphidium" = "#00C094","Heterostegina" = "#000075",
                               "Neorotalia" = "#53B400","Nummulites" = "#006400","Operculina" = "#A58AFF",
                               "Peneroplis" = "#911EB4", "Sorites" = "#FFE119")) +
  scale_fill_manual(values =c("Amphisorus" = "#00B6EB","Amphistegina" = "#FB61D7","Baculogypsinoides" = "#F8766D",
                               "Calcarina" = "#E6194B","Elphidium" = "#00C094","Heterostegina" = "#000075",
                               "Neorotalia" = "#53B400","Nummulites" = "#006400","Operculina" = "#A58AFF",
                               "Peneroplis" = "#911EB4", "Sorites" = "#FFE119")) +
  ylab("Reads - Area (%)") + 
  xlab("Data types") +
  theme_classic2() + 
  theme(axis.text.x=element_blank(),axis.title.x=element_blank(),legend.position = "none",
        strip.background = element_blank())

p2 <- ggplot(df4count[df4count$category != "log_count",], aes(x=category, y=value*100, col = Genus)) +
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = -5, ymax = 5,
           fill = "black", alpha = 0.1) +
  geom_abline(slope = 0, col = "gray") +
  geom_half_boxplot() +
  geom_half_point(aes(fill = Genus), alpha = 0.3) +
  geom_signif(col = "black", comparisons = list(c("raw_count", "lm_count")), 
              map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05), test = "t.test", 
              y_position = 50, tip_length = 0.01, vjust = .1) +
  facet_grid(.~Genus, scales = "free") + 
  scale_color_manual(values =c("Amphisorus" = "#00B6EB","Amphistegina" = "#FB61D7","Baculogypsinoides" = "#F8766D",
                               "Calcarina" = "#E6194B","Elphidium" = "#00C094","Heterostegina" = "#000075",
                               "Neorotalia" = "#53B400","Nummulites" = "#006400","Operculina" = "#A58AFF",
                               "Peneroplis" = "#911EB4", "Sorites" = "#FFE119")) +
  scale_fill_manual(values =c("Amphisorus" = "#00B6EB","Amphistegina" = "#FB61D7","Baculogypsinoides" = "#F8766D",
                              "Calcarina" = "#E6194B","Elphidium" = "#00C094","Heterostegina" = "#000075",
                              "Neorotalia" = "#53B400","Nummulites" = "#006400","Operculina" = "#A58AFF",
                              "Peneroplis" = "#911EB4", "Sorites" = "#FFE119")) +
  ylab("Reads - Counts (%)") + 
  xlab("Data types") +
  theme_classic2() + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5),legend.position = "none",
        strip.background = element_blank(), strip.text = element_blank())

layout <- "
A
B
"
p1 + p2 + plot_layout(design = layout)



 
#--------fig6 - Depth gradient picked, area, counts - LM - with correction on proportion reads + t-test-----------------
df <- read.csv("~/Downloads/Dataset/Genuslevel_counts_reads_area_merged_new.csv")

df <- df[, c(1:9)]
df$picked[df$totalarea == 0] <- 0

df$area_mm2 <- df$totalarea/1000000
df$picked_million <- df$picked/1000000
#from linear regression of ddPCR gene copy number by surface area (lm) - mean from 999 permutations
permut <- read.csv("~/Downloads/Dataset/ddpcr_powerlaw_linear_coefficients_permutations999.csv")

df$density_factor_ddPCR[df$Genus == "Amphisorus"] <- permut$mean_lma[permut$species == "Amphisorus SpL"]
df$density_factor_ddPCR[df$Genus == "Amphistegina"] <- permut$mean_lma[permut$species == "Amphistegina lessonii"]
df$density_factor_ddPCR[df$Genus == "Baculogypsinoides"] <- permut$mean_lma[permut$species == "Baculogypsinoides spinosus"]
df$density_factor_ddPCR[df$Genus == "Calcarina"] <- permut$mean_lma[permut$species == "Calcarina spengleri"]
df$density_factor_ddPCR[df$Genus == "Heterostegina"] <- permut$mean_lma[permut$species == "Heterostegina depressa"]
df$density_factor_ddPCR[df$Genus == "Neorotalia"] <- permut$mean_lma[permut$species == "Neorotalia gaimardi"]
df$density_factor_ddPCR[df$Genus == "Operculina"] <- permut$mean_lma[permut$species == "Operculina ammonoides"]
df$density_factor_ddPCR[df$Genus == "Elphidium"] <- permut$mean_lma[permut$species == "Calcarina spengleri"]
df$density_factor_ddPCR[df$Genus == "Nummulites"] <- permut$mean_lma[permut$species == "Heterostegina depressa"]
df$density_factor_ddPCR[df$Genus == "Peneroplis"] <- permut$mean_lma[permut$species == "Amphisorus SpL"]
df$density_factor_ddPCR[df$Genus == "Sorites"] <- permut$mean_lma[permut$species == "Amphisorus SpL"]

#from powerlaw regression of ddPCR gene copy number by surface area (log) - mean from 999 permutations

df$power_a_ddPCR[df$Genus == "Amphisorus"] <- permut$mean_a[permut$species == "Amphisorus SpL"]
df$power_a_ddPCR[df$Genus == "Amphistegina"] <- permut$mean_a[permut$species == "Amphistegina lessonii"]
df$power_a_ddPCR[df$Genus == "Baculogypsinoides"] <- permut$mean_a[permut$species == "Baculogypsinoides spinosus"]
df$power_a_ddPCR[df$Genus == "Calcarina"] <- permut$mean_a[permut$species == "Calcarina spengleri"]
df$power_a_ddPCR[df$Genus == "Heterostegina"] <- permut$mean_a[permut$species == "Heterostegina depressa"]
df$power_a_ddPCR[df$Genus == "Neorotalia"] <- permut$mean_a[permut$species == "Neorotalia gaimardi"]
df$power_a_ddPCR[df$Genus == "Operculina"] <- permut$mean_a[permut$species == "Operculina ammonoides"]
df$power_a_ddPCR[df$Genus == "Elphidium"] <- permut$mean_a[permut$species == "Calcarina spengleri"]
df$power_a_ddPCR[df$Genus == "Nummulites"] <- permut$mean_a[permut$species == "Heterostegina depressa"]
df$power_a_ddPCR[df$Genus == "Peneroplis"] <- permut$mean_a[permut$species == "Amphisorus SpL"]
df$power_a_ddPCR[df$Genus == "Sorites"] <- permut$mean_a[permut$species == "Amphisorus SpL"]

df$power_b_ddPCR[df$Genus == "Amphisorus"] <- permut$mean_b[permut$species == "Amphisorus SpL"]
df$power_b_ddPCR[df$Genus == "Amphistegina"] <- permut$mean_b[permut$species == "Amphistegina lessonii"]
df$power_b_ddPCR[df$Genus == "Baculogypsinoides"] <- permut$mean_b[permut$species == "Baculogypsinoides spinosus"]
df$power_b_ddPCR[df$Genus == "Calcarina"] <- permut$mean_b[permut$species == "Calcarina spengleri"]
df$power_b_ddPCR[df$Genus == "Heterostegina"] <- permut$mean_b[permut$species == "Heterostegina depressa"]
df$power_b_ddPCR[df$Genus == "Neorotalia"] <- permut$mean_b[permut$species == "Neorotalia gaimardi"]
df$power_b_ddPCR[df$Genus == "Operculina"] <- permut$mean_b[permut$species == "Operculina ammonoides"]
df$power_b_ddPCR[df$Genus == "Elphidium"] <- permut$mean_b[permut$species == "Calcarina spengleri"]
df$power_b_ddPCR[df$Genus == "Nummulites"] <- permut$mean_b[permut$species == "Heterostegina depressa"]
df$power_b_ddPCR[df$Genus == "Peneroplis"] <- permut$mean_b[permut$species == "Amphisorus SpL"]
df$power_b_ddPCR[df$Genus == "Sorites"] <- permut$mean_b[permut$species == "Amphisorus SpL"]

#optimized correction factor by iteration (Fridi)

# df$optimized_factor[df$Genus == "Amphisorus"] <- 0.0968
# df$optimized_factor[df$Genus == "Amphistegina"] <- 0.8061
# df$optimized_factor[df$Genus == "Baculogypsinoides"] <- 0.0193
# df$optimized_factor[df$Genus == "Calcarina"] <- 2.2747
# df$optimized_factor[df$Genus == "Elphidium"] <- 1.0388
# df$optimized_factor[df$Genus == "Heterostegina"] <- 2.8996
# df$optimized_factor[df$Genus == "Neorotalia"] <- 1.9349
# df$optimized_factor[df$Genus == "Nummulites"] <- 7.2478
# df$optimized_factor[df$Genus == "Operculina"] <- 0.2284
# df$optimized_factor[df$Genus == "Peneroplis"] <- 2.5285
# df$optimized_factor[df$Genus == "Sorites"] <- 0.1521

#estimated surface area based on linear and log regression from ddPCR results

#linear : reads = density * area
df$estimated_area_lm <- df$picked_million / df$density_factor_ddPCR

#powerlaw: read = a*area^b  ---> log(read) = a*log(area)+ b ---> area = exp((log(read) - b)/a)
df$estimated_area_log <- exp((log(df$picked_million) - df$power_b_ddPCR)/df$power_a_ddPCR)

#optimized by iteration
#df$estimated_area_op <- df$picked_million / df$optimized_factor

new_df2 <- df %>% group_by(field.nmbr.) %>% 
  mutate(prop_area = area_mm2/sum(area_mm2)*100,
         prop_count = totalcounts/sum(totalcounts)*100,
         prop_raw = picked_million/sum(picked_million)*100,
         prop_lm = estimated_area_lm/sum(estimated_area_lm)*100,
         prop_log = estimated_area_log/sum(estimated_area_log)*100)
         #prop_op = estimated_area_op/sum(estimated_area_op))

df2 <- as.data.frame(new_df2[,c(1:4,9,17:21)])


#_____________first plot______________________
area_raw <- lm(prop_area ~ prop_raw + 0, data = df2)
area_lm <- lm(prop_area ~ prop_lm + 0, data = df2)

p1 <- ggplot(df2, aes(y = prop_area)) +
  stat_smooth(aes(x = prop_raw), method = "lm", formula = y~x+0, se = TRUE, col = "#E6194B", fill = "#E6194B", alpha = 0.1) + 
  stat_cor(aes(x = prop_raw), col = "#E6194B", label.y = 95) +
  geom_point(aes(x = prop_raw), col = "#E6194B", shape = 16, size = 2, alpha = 0.5) +
  annotate(geom = "text",x = 0, y = 100, label = "pre-corrected number of reads", col = "#E6194B") +
  annotate(geom = "text", x = 60, y = 95, label = paste("slope = ", round(area_raw$coefficients[[1]], digits = 3), sep = ""), col = "#E6194B") +
  stat_smooth(aes(x = prop_lm), method = "lm", formula = y~x+0, se = TRUE, col = "#000075", fill = "#000075", alpha = 0.1) + 
  stat_cor(aes(x = prop_lm), col = "#000075", label.y = 80) +
  geom_point(aes(x = prop_lm), col = "#000075", shape = 17, size = 2, alpha = 0.5) +
  annotate(geom = "text",x = 0, y = 85, label = "post-corrected number of reads", col = "#000075") +
  annotate(geom = "text",x = 60, y = 80, label = paste("slope = ", round(area_lm$coefficients[[1]], digits = 3), sep = ""), col = "#000075") +
  geom_abline(linetype = "dashed") +
  theme_classic()+
  ylab("Proportional surface area (%)") +
  xlab("Relative number of reads (pre- and post-corrected) (%)")

count_raw <- lm(prop_count ~ prop_raw + 0, data = df2)
count_lm <- lm(prop_count ~ prop_lm + 0, data = df2)

p2 <- ggplot(df2, aes(y = prop_count)) +
  stat_smooth(aes(x = prop_raw), method = "lm", formula = y~x+0, se = TRUE, col = "#E6194B", fill = "#E6194B", alpha = 0.1) + 
  stat_cor(aes(x = prop_raw), col = "#E6194B", label.y = 95) +
  geom_point(aes(x = prop_raw), col = "#E6194B", shape = 16, size = 2, alpha = 0.5) +
  annotate(geom = "text",x = 0, y = 100, label = "pre-corrected number of reads", col = "#E6194B") +
  annotate(geom = "text", x = 60, y = 95, label = paste("slope = ", round(count_raw$coefficients[[1]], digits = 3), sep = ""), col = "#E6194B") +
  stat_smooth(aes(x = prop_lm), method = "lm", formula = y~x+0, se = TRUE, col = "#000075", fill = "#000075", alpha = 0.1) + 
  stat_cor(aes(x = prop_lm), col = "#000075", label.y = 80) +
  geom_point(aes(x = prop_lm), col = "#000075", shape = 17, size = 2, alpha = 0.5) +
  annotate(geom = "text",x = 0, y = 85, label = "post-corrected number of reads", col = "#000075") +
  annotate(geom = "text",x = 60, y = 80, label = paste("slope = ", round(count_lm$coefficients[[1]], digits = 3), sep = ""), col = "#000075") +
  geom_abline(linetype = "dashed") +
  theme_classic()+
  ylab("Relative abundance (%)") +
  xlab("Relative number of reads (pre- and post-corrected) (%)")

layout <- "
AB
"
p1 + p2 + plot_layout(design = layout)



df2 <- melt(df2, id.vars = c("merge","island","field.nmbr.","depth","Genus"), 
            value.name = "value", variable.name = "type")

orderline <- c("prop_raw", "prop_lm", "prop_log","prop_area", "prop_count")
df2 <- df2 %>% mutate(type = factor(type, levels = orderline)) %>%
  arrange(type)

df2a <- dcast(df2, island + depth + type + field.nmbr. ~ Genus, 
              fun.aggregate = sum, value.var = "value")
df2a[is.na(df2a)] <- 0
df2a <- melt(df2a, id.vars = c("field.nmbr.","island","depth","type"), 
             value.name = "value", variable.name = "Genus")

df2a$type_genus <- paste(df2a$type,df2a$Genus, sep = "_")


#Depth gradient mean reads merged triplicates, with correction, without bulk
df20 <- df2a[df2a$type %not% "prop_bulk",]

df20$depth <- as.character(df20$depth)

depthord <- c("1","1.5","2","2.5","3.5","5.5","6.6","7.1","8.1","8.7","10", "11.3","11.4","13.6","14.6","14.8","16.9", 
              "17.4","17.6","19.5","20.2","21.7","22.6","23.2","24.5","25.6","25.9","27.4","30.5")
depthord <- rev(depthord)

ggplot(df20, aes(y = value, x = factor(depth, depthord), fill = Genus)) +
  geom_bar(stat = 'identity', position = 'fill', width=1, alpha = 0.7)+
  scale_fill_manual(values = c("Amphisorus" = "#00B6EB","Amphistegina" = "#FB61D7","Baculogypsinoides" = "#F8766D",
                               "Calcarina" = "#E6194B","Elphidium" = "#00C094","Heterostegina" = "#000075",
                               "Neorotalia" = "#53B400","Nummulites" = "#006400","Operculina" = "#A58AFF",
                               "Peneroplis" = "#911EB4", "Sorites" = "#FFE119")) +
  theme_classic()+
  facet_grid(island~factor(type, levels=c("prop_count","prop_area","prop_raw", "prop_lm")), scales = "free") +
  theme(legend.position = "bottom", legend.title = element_blank(),strip.background = element_blank()) +
  ylab("Proportion") +
  xlab("Water depth (m)")+
  coord_flip() 


#___________paired pairwise t-test between proportion of methods_______________________

#transform categorical variables into factors
df20 <- as.data.frame(df20)

#pairwise_t_test

stat_test_type <- pairwise_t_test(df20, value ~ type, 
                                     p.adjust.method = "bonferroni")

write.csv(stat_test_type, "~/Downloads/Dataset/pairwise_t_test_results_pickedvspickedcorrvsareavscount_NEWdata.csv", row.names=FALSE)

#___________paired pairwise t-test between proportion of methods with Euclidean distances_________________________

counts <- df2a[df2a$type %in% "prop_count",]
counts1 <- dcast(counts, field.nmbr. ~ Genus,  fun.aggregate = sum, value.var = "value")
rownames(counts1) <- counts1$field.nmbr.
counts1 <- counts1[,-1]

measured <- df2a[df2a$type %in% "prop_area",]
measured1 <- dcast(measured, field.nmbr. ~ Genus,  fun.aggregate = sum, value.var = "value")
rownames(measured1) <- measured1$field.nmbr.
measured1 <- measured1[,-1]

raw <- df2a[df2a$type %in% "prop_raw",]
raw1 <- dcast(raw, field.nmbr. ~ Genus,  fun.aggregate = sum, value.var = "value")
rownames(raw1) <- raw1$field.nmbr.
raw1 <- raw1[,-1]

corrected_lm <- df2a[df2a$type %in% "prop_lm",]
corrected_lm1 <- dcast(corrected_lm, field.nmbr. ~ Genus,  fun.aggregate = sum, value.var = "value")
rownames(corrected_lm1) <- corrected_lm1$field.nmbr.
corrected_lm1 <- corrected_lm1[,-1]

corrected_log <- df2a[df2a$type %in% "prop_log",]
corrected_log1 <- dcast(corrected_log, field.nmbr. ~ Genus,  fun.aggregate = sum, value.var = "value")
rownames(corrected_log1) <- corrected_log1$field.nmbr.
corrected_log1 <- corrected_log1[,-1]

# compute two distance variable using Euclidean distance
raw_meas.distance <- sqrt(apply((measured1 - raw1)^2,1,sum))
lm_meas.distance <- sqrt(apply((measured1 - corrected_lm1)^2,1,sum))
log_meas.distance <- sqrt(apply((measured1 - corrected_log1)^2,1,sum))

raw_count.distance <- sqrt(apply((counts1 - raw1)^2,1,sum))
lm_count.distance <- sqrt(apply((counts1 - corrected_lm1)^2,1,sum))
log_count.distance <- sqrt(apply((counts1 - corrected_log1)^2,1,sum))

result_ttest <- data.frame()
result_ttest[1:16,] <- NA

# run ordinary t-test (means)
result_ttest$name[1] <- t.test(lm_meas.distance, raw_meas.distance, alternative = "less")$data.name
result_ttest$alternative[1] <- t.test(lm_meas.distance, raw_meas.distance, alternative = "less")$alternative
result_ttest$p_value[1] <- t.test(lm_meas.distance, raw_meas.distance, alternative = "less")$p.value
result_ttest$mean_x[1] <- t.test(lm_meas.distance, raw_meas.distance, alternative = "less")$estimate[[1]]
result_ttest$mean_y[1] <- t.test(lm_meas.distance, raw_meas.distance, alternative = "less")$estimate[[2]]

result_ttest$name[2] <- t.test(log_meas.distance, raw_meas.distance, alternative = "less")$data.name
result_ttest$alternative[2] <- t.test(log_meas.distance, raw_meas.distance, alternative = "less")$alternative
result_ttest$p_value[2] <- t.test(log_meas.distance, raw_meas.distance, alternative = "less")$p.value
result_ttest$mean_x[2] <- t.test(log_meas.distance, raw_meas.distance, alternative = "less")$estimate[[1]]
result_ttest$mean_y[2] <- t.test(log_meas.distance, raw_meas.distance, alternative = "less")$estimate[[2]]

result_ttest$name[3] <- t.test(log_meas.distance, lm_meas.distance, alternative = "less")$data.name
result_ttest$alternative[3] <- t.test(log_meas.distance, lm_meas.distance, alternative = "less")$alternative
result_ttest$p_value[3] <- t.test(log_meas.distance, lm_meas.distance, alternative = "less")$p.value
result_ttest$mean_x[3] <- t.test(log_meas.distance, lm_meas.distance, alternative = "less")$estimate[[1]]
result_ttest$mean_y[3] <- t.test(log_meas.distance, lm_meas.distance, alternative = "less")$estimate[[2]]

result_ttest$name[4] <- t.test(lm_meas.distance, raw_meas.distance, alternative = "greater")$data.name
result_ttest$alternative[4] <- t.test(lm_meas.distance, raw_meas.distance, alternative = "greater")$alternative
result_ttest$p_value[4] <- t.test(lm_meas.distance, raw_meas.distance, alternative = "greater")$p.value
result_ttest$mean_x[4] <- t.test(lm_meas.distance, raw_meas.distance, alternative = "greater")$estimate[[1]]
result_ttest$mean_y[4] <- t.test(lm_meas.distance, raw_meas.distance, alternative = "greater")$estimate[[2]]

result_ttest$name[5] <- t.test(log_meas.distance, raw_meas.distance, alternative = "greater")$data.name
result_ttest$alternative[5] <- t.test(log_meas.distance, raw_meas.distance, alternative = "greater")$alternative
result_ttest$p_value[5] <- t.test(log_meas.distance, raw_meas.distance, alternative = "greater")$p.value
result_ttest$mean_x[5] <- t.test(log_meas.distance, raw_meas.distance, alternative = "greater")$estimate[[1]]
result_ttest$mean_y[5] <- t.test(log_meas.distance, raw_meas.distance, alternative = "greater")$estimate[[2]]

result_ttest$name[6] <- t.test(log_meas.distance, lm_meas.distance, alternative = "greater")$data.name
result_ttest$alternative[6] <- t.test(log_meas.distance, lm_meas.distance, alternative = "greater")$alternative
result_ttest$p_value[6] <- t.test(log_meas.distance, lm_meas.distance, alternative = "greater")$p.value
result_ttest$mean_x[6] <- t.test(log_meas.distance, lm_meas.distance, alternative = "greater")$estimate[[1]]
result_ttest$mean_y[6] <- t.test(log_meas.distance, lm_meas.distance, alternative = "greater")$estimate[[2]]

result_ttest$name[7] <- t.test(lm_count.distance, raw_count.distance, alternative = "less")$data.name
result_ttest$alternative[7] <- t.test(lm_count.distance, raw_count.distance, alternative = "less")$alternative
result_ttest$p_value[7] <- t.test(lm_count.distance, raw_count.distance, alternative = "less")$p.value
result_ttest$mean_x[7] <- t.test(lm_count.distance, raw_count.distance, alternative = "less")$estimate[[1]]
result_ttest$mean_y[7] <- t.test(lm_count.distance, raw_count.distance, alternative = "less")$estimate[[2]]

result_ttest$name[8] <- t.test(log_count.distance, raw_count.distance, alternative = "less")$data.name
result_ttest$alternative[8] <- t.test(log_count.distance, raw_count.distance, alternative = "less")$alternative
result_ttest$p_value[8] <- t.test(log_count.distance, raw_count.distance, alternative = "less")$p.value
result_ttest$mean_x[8] <- t.test(log_count.distance, raw_count.distance, alternative = "less")$estimate[[1]]
result_ttest$mean_y[8] <- t.test(log_count.distance, raw_count.distance, alternative = "less")$estimate[[2]]

result_ttest$name[9] <- t.test(log_count.distance, lm_count.distance, alternative = "less")$data.name
result_ttest$alternative[9] <- t.test(log_count.distance, lm_count.distance, alternative = "less")$alternative
result_ttest$p_value[9] <- t.test(log_count.distance, lm_count.distance, alternative = "less")$p.value
result_ttest$mean_x[9] <- t.test(log_count.distance, lm_count.distance, alternative = "less")$estimate[[1]]
result_ttest$mean_y[9] <- t.test(log_count.distance, lm_count.distance, alternative = "less")$estimate[[2]]

result_ttest$name[10] <- t.test(lm_count.distance, raw_count.distance, alternative = "greater")$data.name
result_ttest$alternative[10] <- t.test(lm_count.distance, raw_count.distance, alternative = "greater")$alternative
result_ttest$p_value[10] <- t.test(lm_count.distance, raw_count.distance, alternative = "greater")$p.value
result_ttest$mean_x[10] <- t.test(lm_count.distance, raw_count.distance, alternative = "greater")$estimate[[1]]
result_ttest$mean_y[10] <- t.test(lm_count.distance, raw_count.distance, alternative = "greater")$estimate[[2]]

result_ttest$name[11] <- t.test(log_count.distance, raw_count.distance, alternative = "greater")$data.name
result_ttest$alternative[11] <- t.test(log_count.distance, raw_count.distance, alternative = "greater")$alternative
result_ttest$p_value[11] <- t.test(log_count.distance, raw_count.distance, alternative = "greater")$p.value
result_ttest$mean_x[11] <- t.test(log_count.distance, raw_count.distance, alternative = "greater")$estimate[[1]]
result_ttest$mean_y[11] <- t.test(log_count.distance, raw_count.distance, alternative = "greater")$estimate[[2]]

result_ttest$name[12] <- t.test(log_count.distance, lm_count.distance, alternative = "greater")$data.name
result_ttest$alternative[12] <- t.test(log_count.distance, lm_count.distance, alternative = "greater")$alternative
result_ttest$p_value[12] <- t.test(log_count.distance, lm_count.distance, alternative = "greater")$p.value
result_ttest$mean_x[12] <- t.test(log_count.distance, lm_count.distance, alternative = "greater")$estimate[[1]]
result_ttest$mean_y[12] <- t.test(log_count.distance, lm_count.distance, alternative = "greater")$estimate[[2]]

result_ttest$name[13] <- t.test(log_count.distance, log_meas.distance, alternative = "less")$data.name
result_ttest$alternative[13] <- t.test(log_count.distance, log_meas.distance, alternative = "less")$alternative
result_ttest$p_value[13] <- t.test(log_count.distance, log_meas.distance, alternative = "less")$p.value
result_ttest$mean_x[13] <- t.test(log_count.distance, log_meas.distance, alternative = "less")$estimate[[1]]
result_ttest$mean_y[13] <- t.test(log_count.distance, log_meas.distance, alternative = "less")$estimate[[2]]

result_ttest$name[14] <- t.test(log_count.distance, log_meas.distance, alternative = "greater")$data.name
result_ttest$alternative[14] <- t.test(log_count.distance, log_meas.distance, alternative = "greater")$alternative
result_ttest$p_value[14] <- t.test(log_count.distance, log_meas.distance, alternative = "greater")$p.value
result_ttest$mean_x[14] <- t.test(log_count.distance, log_meas.distance, alternative = "greater")$estimate[[1]]
result_ttest$mean_y[14] <- t.test(log_count.distance, log_meas.distance, alternative = "greater")$estimate[[2]]

result_ttest$name[15] <- t.test(lm_count.distance, lm_meas.distance, alternative = "less")$data.name
result_ttest$alternative[15] <- t.test(lm_count.distance, lm_meas.distance, alternative = "less")$alternative
result_ttest$p_value[15] <- t.test(lm_count.distance, lm_meas.distance, alternative = "less")$p.value
result_ttest$mean_x[15] <- t.test(lm_count.distance, lm_meas.distance, alternative = "less")$estimate[[1]]
result_ttest$mean_y[15] <- t.test(lm_count.distance, lm_meas.distance, alternative = "less")$estimate[[2]]

result_ttest$name[16] <- t.test(lm_count.distance, lm_meas.distance, alternative = "greater")$data.name
result_ttest$alternative[16] <- t.test(lm_count.distance, lm_meas.distance, alternative = "greater")$alternative
result_ttest$p_value[16] <- t.test(lm_count.distance, lm_meas.distance, alternative = "greater")$p.value
result_ttest$mean_x[16] <- t.test(lm_count.distance, lm_meas.distance, alternative = "greater")$estimate[[1]]
result_ttest$mean_y[16] <- t.test(lm_count.distance, lm_meas.distance, alternative = "greater")$estimate[[2]]

write.csv(result_ttest, "~/Downloads/Dataset/welchtwosample_t_test_results_new.csv", row.names=FALSE)

result_vtest <- data.frame()
result_vtest[1:16,] <- NA

# run ordinary F-test (variances)
result_vtest$name[1] <- var.test(lm_meas.distance, raw_meas.distance, alternative = "less")$data.name
result_vtest$alternative[1] <- var.test(lm_meas.distance, raw_meas.distance, alternative = "less")$alternative
result_vtest$p_value[1] <- var.test(lm_meas.distance, raw_meas.distance, alternative = "less")$p.value
result_vtest$ratio_variance[1] <- var.test(lm_meas.distance, raw_meas.distance, alternative = "less")$estimate[[1]]

result_vtest$name[2] <- var.test(log_meas.distance, raw_meas.distance, alternative = "less")$data.name
result_vtest$alternative[2] <- var.test(log_meas.distance, raw_meas.distance, alternative = "less")$alternative
result_vtest$p_value[2] <- var.test(log_meas.distance, raw_meas.distance, alternative = "less")$p.value
result_vtest$ratio_variance[2] <- var.test(log_meas.distance, raw_meas.distance, alternative = "less")$estimate[[1]]

result_vtest$name[3] <- var.test(log_meas.distance, lm_meas.distance, alternative = "less")$data.name
result_vtest$alternative[3] <- var.test(log_meas.distance, lm_meas.distance, alternative = "less")$alternative
result_vtest$p_value[3] <- var.test(log_meas.distance, lm_meas.distance, alternative = "less")$p.value
result_vtest$ratio_variance[3] <- var.test(log_meas.distance, lm_meas.distance, alternative = "less")$estimate[[1]]

result_vtest$name[4] <- var.test(lm_meas.distance, raw_meas.distance, alternative = "greater")$data.name
result_vtest$alternative[4] <- var.test(lm_meas.distance, raw_meas.distance, alternative = "greater")$alternative
result_vtest$p_value[4] <- var.test(lm_meas.distance, raw_meas.distance, alternative = "greater")$p.value
result_vtest$ratio_variance[4] <- var.test(lm_meas.distance, raw_meas.distance, alternative = "greater")$estimate[[1]]

result_vtest$name[5] <- var.test(log_meas.distance, raw_meas.distance, alternative = "greater")$data.name
result_vtest$alternative[5] <- var.test(log_meas.distance, raw_meas.distance, alternative = "greater")$alternative
result_vtest$p_value[5] <- var.test(log_meas.distance, raw_meas.distance, alternative = "greater")$p.value
result_vtest$ratio_variance[5] <- var.test(log_meas.distance, raw_meas.distance, alternative = "greater")$estimate[[1]]

result_vtest$name[6] <- var.test(log_meas.distance, lm_meas.distance, alternative = "greater")$data.name
result_vtest$alternative[6] <- var.test(log_meas.distance, lm_meas.distance, alternative = "greater")$alternative
result_vtest$p_value[6] <- var.test(log_meas.distance, lm_meas.distance, alternative = "greater")$p.value
result_vtest$ratio_variance[6] <- var.test(log_meas.distance, lm_meas.distance, alternative = "greater")$estimate[[1]]

result_vtest$name[7] <- var.test(lm_count.distance, raw_count.distance, alternative = "less")$data.name
result_vtest$alternative[7] <- var.test(lm_count.distance, raw_count.distance, alternative = "less")$alternative
result_vtest$p_value[7] <- var.test(lm_count.distance, raw_count.distance, alternative = "less")$p.value
result_vtest$ratio_variance[7] <- var.test(lm_count.distance, raw_count.distance, alternative = "less")$estimate[[1]]

result_vtest$name[8] <- var.test(log_count.distance, raw_count.distance, alternative = "less")$data.name
result_vtest$alternative[8] <- var.test(log_count.distance, raw_count.distance, alternative = "less")$alternative
result_vtest$p_value[8] <- var.test(log_count.distance, raw_count.distance, alternative = "less")$p.value
result_vtest$ratio_variance[8] <- var.test(log_count.distance, raw_count.distance, alternative = "less")$estimate[[1]]

result_vtest$name[9] <- var.test(log_count.distance, lm_count.distance, alternative = "less")$data.name
result_vtest$alternative[9] <- var.test(log_count.distance, lm_count.distance, alternative = "less")$alternative
result_vtest$p_value[9] <- var.test(log_count.distance, lm_count.distance, alternative = "less")$p.value
result_vtest$ratio_variance[9] <- var.test(log_count.distance, lm_count.distance, alternative = "less")$estimate[[1]]

result_vtest$name[10] <- var.test(lm_count.distance, raw_count.distance, alternative = "greater")$data.name
result_vtest$alternative[10] <- var.test(lm_count.distance, raw_count.distance, alternative = "greater")$alternative
result_vtest$p_value[10] <- var.test(lm_count.distance, raw_count.distance, alternative = "greater")$p.value
result_vtest$ratio_variance[10] <- var.test(lm_count.distance, raw_count.distance, alternative = "greater")$estimate[[1]]

result_vtest$name[11] <- var.test(log_count.distance, raw_count.distance, alternative = "greater")$data.name
result_vtest$alternative[11] <- var.test(log_count.distance, raw_count.distance, alternative = "greater")$alternative
result_vtest$p_value[11] <- var.test(log_count.distance, raw_count.distance, alternative = "greater")$p.value
result_vtest$ratio_variance[11] <- var.test(log_count.distance, raw_count.distance, alternative = "greater")$estimate[[1]]

result_vtest$name[12] <- var.test(log_count.distance, lm_count.distance, alternative = "greater")$data.name
result_vtest$alternative[12] <- var.test(log_count.distance, lm_count.distance, alternative = "greater")$alternative
result_vtest$p_value[12] <- var.test(log_count.distance, lm_count.distance, alternative = "greater")$p.value
result_vtest$ratio_variance[12] <- var.test(log_count.distance, lm_count.distance, alternative = "greater")$estimate[[1]]

result_vtest$name[13] <- var.test(log_count.distance, log_meas.distance, alternative = "less")$data.name
result_vtest$alternative[13] <- var.test(log_count.distance, log_meas.distance, alternative = "less")$alternative
result_vtest$p_value[13] <- var.test(log_count.distance, log_meas.distance, alternative = "less")$p.value
result_vtest$ratio_variance[13] <- var.test(log_count.distance, log_meas.distance, alternative = "less")$estimate[[1]]

result_vtest$name[14] <- var.test(log_count.distance, log_meas.distance, alternative = "greater")$data.name
result_vtest$alternative[14] <- var.test(log_count.distance, log_meas.distance, alternative = "greater")$alternative
result_vtest$p_value[14] <- var.test(log_count.distance, log_meas.distance, alternative = "greater")$p.value
result_vtest$ratio_variance[14] <- var.test(log_count.distance, log_meas.distance, alternative = "greater")$estimate[[1]]

result_vtest$name[15] <- var.test(lm_count.distance, lm_meas.distance, alternative = "less")$data.name
result_vtest$alternative[15] <- var.test(lm_count.distance, lm_meas.distance, alternative = "less")$alternative
result_vtest$p_value[15] <- var.test(lm_count.distance, lm_meas.distance, alternative = "less")$p.value
result_vtest$ratio_variance[15] <- var.test(lm_count.distance, lm_meas.distance, alternative = "less")$estimate[[1]]

result_vtest$name[16] <- var.test(lm_count.distance, lm_meas.distance, alternative = "greater")$data.name
result_vtest$alternative[16] <- var.test(lm_count.distance, lm_meas.distance, alternative = "greater")$alternative
result_vtest$p_value[16] <- var.test(lm_count.distance, lm_meas.distance, alternative = "greater")$p.value
result_vtest$ratio_variance[16] <- var.test(lm_count.distance, lm_meas.distance, alternative = "greater")$estimate[[1]]

write.csv(result_vtest, "~/Downloads/Dataset/f_test_results.csv", row.names=FALSE)




#___________ANOSIM_________________________________________

df20$value_ID <- paste(df20$field.nmbr., "_", df20$type, sep = "")

df20 <- df20[df20$type != "prop_log",]

df5 <- dcast(df20, island + depth + type + field.nmbr. + value_ID ~ Genus, 
             fun.aggregate = sum, value.var = "value")
df5$island <- as.factor(as.character(df5$island))
df5$field.nmbr. <- as.factor(as.character(df5$field.nmbr.))
df5$type <- as.factor(as.character(df5$type))

par(mfrow = c(2, 1))

#I need to isolate values (N), Species code (x) and Sample codes (y) with col and row names
df_NMDS_vegan <- dcast(df20, value_ID ~ Genus, fun.aggregate = sum, value.var = "value")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$value_ID
df_NMDS_vegan <- df_NMDS_vegan[,-1]


#transform --> inversion columns and rows only if species (x) and samples (y) - Samples must be the rows and species the columns
#dt_t <- t(df)
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns
ord <- metaMDS(df_NMDS_vegan, distance = "bray")
env <- df5[,6:ncol(df5)]

#to calculate the species vectors going over the NMDS
en <- envfit(ord, env, permutations = 999, na.rm = TRUE)

#groups for plots
island <- droplevels(df5$island)
type <- droplevels(df5$type)
field <- droplevels(df5$field.nmbr.)

pchtype <- c(0,21,15,19, 20) 
colfield <- topo.colors(35)

#Significance is the p value
ano <- anosim(df_NMDS_vegan, type)

#for plot
plot(ord, disp="sites", type = "n", main = paste ("Similarity between data types - ANOSIM: p =", ano$signif, " and R =", round(ano$statistic, 3)))
#ordiellipse(ord, location, col= colils, cex=2, conf = 0.99)
points(ord, disp="sites", pch = pchtype[type], cex=1.5, col = colfield[field])
with(ord, legend(x = "topright", legend = levels(field), col = colfield, pch = 15))
with(ord, legend(x = "topleft", legend = levels(type), col = "black", pch = pchtype))

ano <- anosim(df_NMDS_vegan, field)

#for plot
plot(ord, disp="sites", type = "n", main = paste ("Similarity between sampling sites - ANOSIM: p =", ano$signif, " and R =", round(ano$statistic, 3)))
#ordiellipse(ord, location, col= colils, cex=2, conf = 0.99)
points(ord, disp="sites", pch = pchtype[type], cex=1.5, col = colfield[field])
with(ord, legend(x = "topright", legend = levels(field), col = colfield, pch = 15))
with(ord, legend(x = "topleft", legend = levels(type), col = "black", pch = pchtype))


#--------fig7a,b,c - NMDS bulk vs sorted - t-test and ANOSIM genus level--------

df <- read.csv("~/Downloads/Dataset/workingdataset-quantitativemetabarcoding-994_new_replicatefiltered.csv")

df1 <- df

lbf <- c("Sorites_sp2_Spermonde","Sorites_sp1_Spermonde","Peneroplis_sp2_&_Peneroplis_pertusus_5117_&_Dendritina_ambigua_Spermonde",
         "Peneroplis_sp1_Spermonde","Parasorites_sp","Operculina_complanata","Operculina_ammonoides","Nummulites_venosus",                                                                
         "Neorotalia_gaimardi_&_Baculogypsina_sphaerulata_Spermonde","Neorotalia_calcar_Spermonde","Heterostegina_depressa_sp2_Spermonde",                                              
         "Heterostegina_depressa_sp1_Spermonde","Elphidium_Rik_9327_9337_9350_Spermonde","Calcarina_spengleri_Rik_9352_consensus0_1732_Spermonde",                            
         "Calcarina_sp._Spermonde","Calcarina_sp._5164_Spermonde","Calcarina_hispida_Ambon_&_Calcarina_spengleri_&_Calcarina_sp._5247_Spermonde",      
         "Calcarina_hispida_&_Calcarina_sp_5163_Spermonde","Baculogypsinoides_spinosus_8078_8079_8080",                                         
         "Amphistegina_radiata_Spermonde","Amphistegina_papillosa_Spermonde_3741","Amphistegina_lobifera","Amphistegina_lessonii",                                                             
         "Amphistegina_bicirculata_Rik_9379_9380_9381_9382_Eilat","Amphisorus_&_Amphisorus_SpL_&_Amphisorus_SpS_Spermonde")

df1 <- df1[df1$only_species_level %in% lbf,]

#LBF
df1$species[df1$only_species_level == "Sorites_sp2_Spermonde"] <- "Sorites sp1"
df1$species[df1$only_species_level == "Sorites_sp1_Spermonde"] <- "Sorites sp2"
df1$species[df1$only_species_level == "Peneroplis_sp2_&_Peneroplis_pertusus_5117_&_Dendritina_ambigua_Spermonde"] <- "Peneroplis sp2/Dendritina ambigua"
df1$species[df1$only_species_level == "Peneroplis_sp1_Spermonde"] <- "Peneroplis sp1"
df1$species[df1$only_species_level == "Parasorites_sp"] <- "Parasorites sp."
df1$species[df1$only_species_level == "Operculina_complanata"] <- "Operculina complanata"
df1$species[df1$only_species_level == "Operculina_ammonoides"] <- "Operculina ammonoides"
df1$species[df1$only_species_level == "Nummulites_venosus"] <- "Nummulites venosus"
df1$species[df1$only_species_level == "Neorotalia_gaimardi_&_Baculogypsina_sphaerulata_Spermonde"] <- "Neorotalia gaimardi/Baculogypsina sphearulata"
df1$species[df1$only_species_level == "Neorotalia_calcar_Spermonde"] <- "Neorotalia calcar"
df1$species[df1$only_species_level == "Heterostegina_depressa_sp2_Spermonde"] <- "Heterostegina depressa sp2"
df1$species[df1$only_species_level == "Heterostegina_depressa_sp1_Spermonde"] <- "Heterostegina depressa sp1"
df1$species[df1$only_species_level == "Elphidium_Rik_9327_9337_9350_Spermonde"] <- "Elphidium sp."
df1$species[df1$only_species_level == "Calcarina_spengleri_Rik_9352_consensus0_1732_Spermonde"] <- "Calcarina spengleri"
df1$species[df1$only_species_level == "Calcarina_sp._Spermonde"] <- "Calcarina sp1"
df1$species[df1$only_species_level == "Calcarina_sp._5164_Spermonde"] <- "Calcarina sp2"
df1$species[df1$only_species_level == "Calcarina_hispida_Ambon_&_Calcarina_spengleri_&_Calcarina_sp._5247_Spermonde"] <- "Calcarina sp3"
df1$species[df1$only_species_level == "Calcarina_hispida_&_Calcarina_sp_5163_Spermonde"] <- "Calcarina hispida"
df1$species[df1$only_species_level == "Baculogypsinoides_spinosus_8078_8079_8080"] <- "Baculogypsinoides spinosus"
df1$species[df1$only_species_level == "Amphistegina_radiata_Spermonde"] <- "Amphistegina radiata"
df1$species[df1$only_species_level == "Amphistegina_papillosa_Spermonde_3741"] <- "Amphistegina papillosa"
df1$species[df1$only_species_level == "Amphistegina_lobifera"] <- "Amphistegina lobifera"
df1$species[df1$only_species_level == "Amphistegina_lessonii"] <- "Amphistegina lessonii"
df1$species[df1$only_species_level == "Amphistegina_bicirculata_Rik_9379_9380_9381_9382_Eilat"] <- "Amphistegina bicirculata"
df1$species[df1$only_species_level == "Amphisorus_&_Amphisorus_SpL_&_Amphisorus_SpS_Spermonde"] <- "Amphisorus sp."


splitspecies <- str_split_fixed(df1$species, " ", 2)
colnames(splitspecies) <- c("Genus", "Species_alone")
df2 <- cbind(df1, splitspecies)

splitrep <- str_split_fixed(df2$new.name, "_", 3)
colnames(splitrep) <- c("letter_island", "fieldsample", "type_rep")
df2 <- cbind(df2, splitrep)

unique(df2$Genus)

new_df2 <- df2 %>% group_by(field.nmbr., depth, island, type, type_rep, Genus) %>% 
  summarize(sumreads = sum(reads))


#___________paired pairwise t-test between replicates

#compare the proportion of genus, because picked C had small error in pooling, so a factor different
new_df2 <- new_df2 %>% group_by(type_rep, field.nmbr.) %>% 
  mutate(propreads = sumreads/sum(sumreads))

#transform categorical variables into factors
new_df2 <- as.data.frame(new_df2)

new_df2$type <- as.factor(as.character(new_df2$type))
new_df2$type_rep <- as.factor(as.character(new_df2$type_rep))
new_df2$sumreads <- as.numeric(new_df2$sumreads)

#pairwise_t_test

stat_test_typerep <- pairwise_t_test(new_df2, propreads ~ type_rep, 
                                  p.adjust.method = "bonferroni")

stat_test_type <- pairwise_t_test(new_df2, propreads ~ type, 
                                 p.adjust.method = "bonferroni")

t_test_results <- rbind(stat_test_typerep, stat_test_type)

write.csv(t_test_results, "~/Downloads/Dataset/pairwise_t_test_results_bulkvspicked_includingrep_genus_NEW.csv", row.names=FALSE)

#___________ANOSIM between sample genus level


#To build the NMDS, I have to "unmelt" the dataframe with reshape
new_df2$ID <- paste(new_df2$field.nmbr.,"_",new_df2$type_rep, sep = "")

#___________NMDS BADI_______________________________
df_NMDS <- dcast(new_df2[new_df2$island == "Badi",], field.nmbr. + island + depth + type + type_rep + ID ~ Genus, fun.aggregate = sum, value.var = "sumreads")

df_NMDS$type <- as.factor(as.character(df_NMDS$type))

#I need to isolate values (N), Species code (x) and Sample codes (y) with col and row names
df_NMDS_vegan <- dcast(new_df2[new_df2$island == "Badi",], ID ~ Genus, fun.aggregate = sum, value.var = "sumreads")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$ID
df_NMDS_vegan <- df_NMDS_vegan[,-1]


#transform --> inversion columns and rows only if species (x) and samples (y) - Samples must be the rows and species the columns
#dt_t <- t(df)
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns
ord <- metaMDS(df_NMDS_vegan, distance = "bray")
env <- df_NMDS[,7:ncol(df_NMDS)]

#to calculate the species vectors going over the NMDS
en <- envfit(ord, env, permutations = 999, na.rm = TRUE)

#groups for plots
site <- droplevels(as.factor(df_NMDS$field.nmbr.))
island <- droplevels(as.factor(df_NMDS$island))
type <- droplevels(df_NMDS$type)
type_rep <- droplevels(df_NMDS$type_rep)

colils <- c('#E6194B','#3399FF','#660066')
field <- palette(rainbow(35)) 
pch <- c(1,19)
pchils <- c(0,8,19)
coltype <- c('#E6194B','#3399FF')

#Significance is the p value
ano <- anosim(df_NMDS_vegan, type)
ano2 <- anosim(df_NMDS_vegan, site)

#for plot
plot(ord, disp="sites", type = "n", main = paste ("Difference between data type - ANOSIM: p =", ano$signif, " and R =", round(ano$statistic, 3), "Difference between sampling sites - ANOSIM: p =", ano2$signif, " and R =", round(ano2$statistic, 3)))
points(ord, disp="sites", pch = pchils[island], cex=1.5, 
       col = coltype[type])
with(ord, legend(x = "topright", legend = levels(island), pch = pchils))
with(ord, legend(x = "bottomright", legend = levels(type), col = coltype, pch = 15))


#second option NMDS

plot(ord, disp="sites", type = "n", main = paste ("Badi - ANOSIM (sample types): p =", ano$signif, " and R =", round(ano$statistic, 3), "ANOSIM (sampling site): p =", ano2$signif, " and R =", round(ano2$statistic, 3)))
points(ord, disp="sites", pch = pch[type], cex=1.5, col = colils[island])
ordihull(ord, groups=site,  draw = "polygon", col = "black", alpha = 0.1)
with(ord, legend(x = "topright", legend = levels(island), col = colils, pch = 15))
with(ord, legend(x = "bottomright", legend = levels(type), pch = pch, col = "black"))

#___________NMDS PADJENEKANG_______________________________
df_NMDS <- dcast(new_df2[new_df2$island == "Padjenekang",], field.nmbr. + island + depth + type + type_rep + ID ~ Genus, fun.aggregate = sum, value.var = "sumreads")

df_NMDS$type <- as.factor(as.character(df_NMDS$type))

#I need to isolate values (N), Species code (x) and Sample codes (y) with col and row names
df_NMDS_vegan <- dcast(new_df2[new_df2$island == "Padjenekang",], ID ~ Genus, fun.aggregate = sum, value.var = "sumreads")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$ID
df_NMDS_vegan <- df_NMDS_vegan[,-1]


#transform --> inversion columns and rows only if species (x) and samples (y) - Samples must be the rows and species the columns
#dt_t <- t(df)
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns
ord <- metaMDS(df_NMDS_vegan, distance = "bray")
env <- df_NMDS[,7:ncol(df_NMDS)]

#to calculate the species vectors going over the NMDS
en <- envfit(ord, env, permutations = 999, na.rm = TRUE)

#groups for plots
site <- droplevels(as.factor(df_NMDS$field.nmbr.))
island <- droplevels(as.factor(df_NMDS$island))
type <- droplevels(df_NMDS$type)
type_rep <- droplevels(df_NMDS$type_rep)

colils <- c('#660066')
pch <- c(1,19)

#Significance is the p value
ano <- anosim(df_NMDS_vegan, type)
ano2 <- anosim(df_NMDS_vegan, site)

#second option NMDS

plot(ord, disp="sites", type = "n", main = paste ("Pajenekang - ANOSIM (sample types): p =", ano$signif, " and R =", round(ano$statistic, 3), "ANOSIM (sampling site): p =", ano2$signif, " and R =", round(ano2$statistic, 3)))
points(ord, disp="sites", pch = pch[type], cex=1.5, col = colils[island])
ordihull(ord, groups=site,  draw = "polygon", col = "black", alpha = 0.1)
with(ord, legend(x = "topright", legend = levels(island), col = colils, pch = 15))
with(ord, legend(x = "bottomright", legend = levels(type), pch = pch, col = "black"))

#___________NMDS LUMULUMU_______________________________
df_NMDS <- dcast(new_df2[new_df2$island == "Lumulumu",], field.nmbr. + island + depth + type + type_rep + ID ~ Genus, fun.aggregate = sum, value.var = "sumreads")

df_NMDS$type <- as.factor(as.character(df_NMDS$type))

#I need to isolate values (N), Species code (x) and Sample codes (y) with col and row names
df_NMDS_vegan <- dcast(new_df2[new_df2$island == "Lumulumu",], ID ~ Genus, fun.aggregate = sum, value.var = "sumreads")
rownames(df_NMDS_vegan) <- df_NMDS_vegan$ID
df_NMDS_vegan <- df_NMDS_vegan[,-1]


#transform --> inversion columns and rows only if species (x) and samples (y) - Samples must be the rows and species the columns
#dt_t <- t(df)
df_NMDS_vegan <- decostand(df_NMDS_vegan, "total") #Make sure Samples are in rows and Species in columns
ord <- metaMDS(df_NMDS_vegan, distance = "bray")
env <- df_NMDS[,7:ncol(df_NMDS)]

#to calculate the species vectors going over the NMDS
en <- envfit(ord, env, permutations = 999, na.rm = TRUE)

#groups for plots
site <- droplevels(as.factor(df_NMDS$field.nmbr.))
island <- droplevels(as.factor(df_NMDS$island))
type <- droplevels(df_NMDS$type)
type_rep <- droplevels(df_NMDS$type_rep)

colils <- c('#3399FF','#660066')
pch <- c(1,19)

#Significance is the p value
ano <- anosim(df_NMDS_vegan, type)
ano2 <- anosim(df_NMDS_vegan, site)

#second option NMDS

plot(ord, disp="sites", type = "n", main = paste ("Lumu-lumu - ANOSIM (sample types): p =", ano$signif, " and R =", round(ano$statistic, 3), "ANOSIM (sampling site): p =", ano2$signif, " and R =", round(ano2$statistic, 3)))
points(ord, disp="sites", pch = pch[type], cex=1.5, col = colils[island])
ordihull(ord, groups=site,  draw = "polygon", col = "black", alpha = 0.1)
with(ord, legend(x = "topright", legend = levels(island), col = colils, pch = 15))
with(ord, legend(x = "bottomright", legend = levels(type), pch = pch, col = "black"))

#--------fig7d - Diversity index bulk vs sorted----------------------------
df <- read.csv("~/Downloads/Dataset/workingdataset-quantitativemetabarcoding-994_new_replicatefiltered.csv")

df1 <- df

lbf <- c("Sorites_sp2_Spermonde","Sorites_sp1_Spermonde","Peneroplis_sp2_&_Peneroplis_pertusus_5117_&_Dendritina_ambigua_Spermonde",
         "Peneroplis_sp1_Spermonde","Parasorites_sp","Operculina_complanata","Operculina_ammonoides","Nummulites_venosus",                                                                
         "Neorotalia_gaimardi_&_Baculogypsina_sphaerulata_Spermonde","Neorotalia_calcar_Spermonde","Heterostegina_depressa_sp2_Spermonde",                                              
         "Heterostegina_depressa_sp1_Spermonde","Elphidium_Rik_9327_9337_9350_Spermonde","Calcarina_spengleri_Rik_9352_consensus0_1732_Spermonde",                            
         "Calcarina_sp._Spermonde","Calcarina_sp._5164_Spermonde","Calcarina_hispida_Ambon_&_Calcarina_spengleri_&_Calcarina_sp._5247_Spermonde",      
         "Calcarina_hispida_&_Calcarina_sp_5163_Spermonde","Baculogypsinoides_spinosus_8078_8079_8080",                                         
         "Amphistegina_radiata_Spermonde","Amphistegina_papillosa_Spermonde_3741","Amphistegina_lobifera","Amphistegina_lessonii",                                                             
         "Amphistegina_bicirculata_Rik_9379_9380_9381_9382_Eilat","Amphisorus_&_Amphisorus_SpL_&_Amphisorus_SpS_Spermonde")

df1 <- df1[df1$only_species_level %in% lbf,]

#LBF
df1$species[df1$only_species_level == "Sorites_sp2_Spermonde"] <- "Sorites sp1"
df1$species[df1$only_species_level == "Sorites_sp1_Spermonde"] <- "Sorites sp2"
df1$species[df1$only_species_level == "Peneroplis_sp2_&_Peneroplis_pertusus_5117_&_Dendritina_ambigua_Spermonde"] <- "Peneroplis sp2/Dendritina ambigua"
df1$species[df1$only_species_level == "Peneroplis_sp1_Spermonde"] <- "Peneroplis sp1"
df1$species[df1$only_species_level == "Parasorites_sp"] <- "Parasorites sp."
df1$species[df1$only_species_level == "Operculina_complanata"] <- "Operculina complanata"
df1$species[df1$only_species_level == "Operculina_ammonoides"] <- "Operculina ammonoides"
df1$species[df1$only_species_level == "Nummulites_venosus"] <- "Nummulites venosus"
df1$species[df1$only_species_level == "Neorotalia_gaimardi_&_Baculogypsina_sphaerulata_Spermonde"] <- "Neorotalia gaimardi/Baculogypsina sphearulata"
df1$species[df1$only_species_level == "Neorotalia_calcar_Spermonde"] <- "Neorotalia calcar"
df1$species[df1$only_species_level == "Heterostegina_depressa_sp2_Spermonde"] <- "Heterostegina depressa sp2"
df1$species[df1$only_species_level == "Heterostegina_depressa_sp1_Spermonde"] <- "Heterostegina depressa sp1"
df1$species[df1$only_species_level == "Elphidium_Rik_9327_9337_9350_Spermonde"] <- "Elphidium sp."
df1$species[df1$only_species_level == "Calcarina_spengleri_Rik_9352_consensus0_1732_Spermonde"] <- "Calcarina spengleri"
df1$species[df1$only_species_level == "Calcarina_sp._Spermonde"] <- "Calcarina sp1"
df1$species[df1$only_species_level == "Calcarina_sp._5164_Spermonde"] <- "Calcarina sp2"
df1$species[df1$only_species_level == "Calcarina_hispida_Ambon_&_Calcarina_spengleri_&_Calcarina_sp._5247_Spermonde"] <- "Calcarina sp3"
df1$species[df1$only_species_level == "Calcarina_hispida_&_Calcarina_sp_5163_Spermonde"] <- "Calcarina hispida"
df1$species[df1$only_species_level == "Baculogypsinoides_spinosus_8078_8079_8080"] <- "Baculogypsinoides spinosus"
df1$species[df1$only_species_level == "Amphistegina_radiata_Spermonde"] <- "Amphistegina radiata"
df1$species[df1$only_species_level == "Amphistegina_papillosa_Spermonde_3741"] <- "Amphistegina papillosa"
df1$species[df1$only_species_level == "Amphistegina_lobifera"] <- "Amphistegina lobifera"
df1$species[df1$only_species_level == "Amphistegina_lessonii"] <- "Amphistegina lessonii"
df1$species[df1$only_species_level == "Amphistegina_bicirculata_Rik_9379_9380_9381_9382_Eilat"] <- "Amphistegina bicirculata"
df1$species[df1$only_species_level == "Amphisorus_&_Amphisorus_SpL_&_Amphisorus_SpS_Spermonde"] <- "Amphisorus sp."


splitspecies <- str_split_fixed(df1$species, " ", 2)
colnames(splitspecies) <- c("Genus", "Species_alone")
df2 <- cbind(df1, splitspecies)

splitrep <- str_split_fixed(df2$new.name, "_", 3)
colnames(splitrep) <- c("letter_island", "fieldsample", "type_rep")
df2 <- cbind(df2, splitrep)

unique(df2$Genus)

new_df2 <- df2 %>% group_by(field.nmbr., depth, island, type, type_rep, Genus) %>% 
  summarize(sumreads = sum(reads))

new_df2$counts <- 1

df3 <- new_df2 %>% group_by(field.nmbr., depth, island, type, type_rep) %>% 
  summarize(sumcounts = sum(counts))

ggplot(df3, aes(type, sumcounts, shape = type)) +
  geom_half_boxplot() +
  geom_half_point_panel(aes(col = island), size = 2) +
  geom_signif(col = "black", comparisons = list(c("bulk", "picked")), 
              map_signif_level=FALSE, test = "t.test", 
              y_position = 10, tip_length = 0.01, vjust = .1) +
  scale_shape_manual(values = c(1,19)) +
  scale_color_manual(values = c('#E6194B','#3399FF','#660066')) +
  theme_classic() + 
  theme(strip.text.y = element_text(angle = 0), axis.title.x = element_blank(),
        legend.position = "bottom") +
  ylab("Genus richness")
