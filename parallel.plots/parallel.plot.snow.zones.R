#Script for making parallel plots (significant) using info taken from the R 
#Graph Gallery
#https://www.data-to-viz.com/graph/parallel.html

install.packages("hrbrthemes")
install.packages("patchwork")
install.packages("GGally")
install.packages("viridis")
install.packages("dpylr")
install.packages("ggsignif")
install.packages("ggstance")
install.packages("gridExtra")

library(dplyr)
library(tidyverse)
library(hrbrthemes)
library(patchwork)
library(GGally)
library(viridis)
library(ggplot2)
library(ggsignif)
library(ggstance)
library(readxl)
library(gridExtra)
library(scales)
setwd("/Volumes/Avery.Snow/snow.bloom.2018/stats.analysis/parallel.plots")
abc.plot.original <- read_excel("JUMP.format.paired.wilcoxon.xlsx")

#-----------------------------------------------------------------------------------------------------------------------------

#Paired Rank Sum Values were analyzed in JUMP and are stored in a separate file ("Paired.Stats.Wilcoxon.Results"). The t-values used here are simply to 
  #generate stars and are not representative of the actual analysis

#Richness is the number of species per sample
#Evenness is the measure of OTU relative abundance per sample
#Simpson's index Simpson's Index of Diversity 1 - D 
  #between 0 and 1 the greater the value, the greater the sample diversity
  #the index represents the probability that two individuals randomly selected
  #from a sample will belong to different species.

#creating a parallel plot for snow sample sites across bloom sections A, B, and C
#we will look at algae and bacterial richness
#and RedOx and Conductivity
#and all other environmental factors found significant by our Paired Wilcoxon-Rank Sum test

#%>% is a infix operator, not part of base R. It passess the argument on the left, in this case
#abc.plot to the arguments of the right. In essence it is equivalent in the following way
#head(abc.plot) and abc.plot %>% head()
#possibly an easier way to visualize elements in longer functions

#-----------------------------------------------------------------------------------------------------------------------------

#introduce jitter to all our continuous data in order to visualize overlapping lines

abc.plot.jittered <- data.frame(lapply(abc.plot.original[,c(3:77)], jitter))
abc.plot.cat <- abc.plot.original[,c(1:2)]
abc.plot <- cbind(abc.plot.cat,abc.plot.jittered)

#-----------------------------------------------------------------------------------------------------------------------------

#Nitrate Concentration

nitrate.abc <- data.frame(abc.plot[,c(1,2,23,48,73)])
nitrate.names <- c("Sample_ID","Location","M","P","A")
colnames(nitrate.abc) <- nitrate.names

nitrate <- nitrate.abc %>%
  ggparcoord(
    columns = 3:5, groupColumn = 2,
    showPoints = FALSE, 
    title = "Nitrate",
    alphaLines = 0,
    scale = "globalminmax"
  ) + 
  theme_classic()+
  geom_line(aes( linetype=Location, color=Location), size = 0.25)+
  geom_point(size = 0.65)+
  scale_linetype_manual(values=c("solid", "dashed"))+
  scale_color_manual(values=c('#000000','#000000'))+
  ylab("mV")+
  scale_x_discrete(expand = c(0.04, 0.04))+
  
  theme(
    text = element_text(size=15),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size=18.5, hjust = 0.5),
    legend.position = "none",
    axis.line.x.bottom = element_line(size = 1),
    axis.title.x = element_blank(),
    axis.ticks.x.bottom = element_line(size = 0.90),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 15),
    legend.title = element_blank(),
    line = element_blank(),
    legend.key.size = unit(38,'points'),
    plot.margin = margin(10,2.5,5,17.5)
  )
nitrate
  
#-----------------------------------------------------------------------------------------------------------------------------

#Potassium Concentration

potassium.abc <- data.frame(abc.plot[,c(1,2,25,50,75)])
potassium.names <- c("Sample_ID","Location","M","P","A")
colnames(potassium.abc) <- potassium.names

potassium <- potassium.abc %>%
  ggparcoord(
    columns = 3:5, groupColumn = 2,
    showPoints = FALSE, 
    title = "Potassium",
    alphaLines = 0, 
    scale = 'globalminmax' # 'centerObs'; must be one of 'std', 'robust', 'uniminmax', 'globalminmax', 'center', or 'centerObs'
  ) + 
  theme_classic()+
  geom_line(aes(linetype=Location, color=Location), size= 0.25)+
  geom_point(size = 0.65)+
  geom_signif(comparisons = list(c("M", "A")), test = 't.test', map_signif_level = c(" "=0.25), textsize = 8, tip_length = 0.02, size = 0.075, y_position = -152)+
  geom_signif(comparisons = list(c("M", "P")), test = 't.test', map_signif_level = c(" "=0.25), textsize = 8, tip_length = 0.02, size = 0.075)+
  scale_linetype_manual(values=c("solid", "dashed"))+
  scale_color_manual(values=c('#000000','#000000'))+
  ylab("mV")+
  scale_x_discrete(expand = c(0.04, 0.04))+
  
  theme(
    text = element_text(size=15),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size=18.5, hjust = 0.5),
    legend.position = "none",
    axis.line.x.bottom = element_line(size = 1),
    axis.title.x = element_blank(),
    axis.ticks.x.bottom = element_line(size = 0.90),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 15),
    legend.title = element_blank(),
    line = element_blank(),
    legend.key.size = unit(38,'points'),
    plot.margin = margin(10,2.5,5,17.5)
  )

#-----------------------------------------------------------------------------------------------------------------------------

#Resistivity

resist.abc <- data.frame(abc.plot[,c(1,2,20,45,70)])
resist.names <- c("Sample_ID","Location","M","P","A")
colnames(resist.abc) <- resist.names

res <- resist.abc %>%
  ggparcoord(
    columns = 3:5, groupColumn = 2,
    showPoints = FALSE, 
    title = "Resistivity",
    alphaLines = 0, 
    scale = 'globalminmax' # 'centerObs'; must be one of 'std', 'robust', 'uniminmax', 'globalminmax', 'center', or 'centerObs'
  ) + 
  theme_classic()+
  geom_line(aes(linetype=Location, color=Location), size= 0.25)+
  geom_point(size = 0.65)+
  geom_signif(comparisons = list(c("M", "P")), test = 't.test', map_signif_level = c(" "=0.50), textsize = 8, tip_length = 0.02, size = 0.075)+
  geom_signif(comparisons = list(c("M", "A")), test = 't.test', map_signif_level = c(" "=0.50), y_position = 1.120, textsize = 8, tip_length = 0.02, size = 0.075)+
  scale_linetype_manual(values=c("solid", "dashed"))+
  scale_color_manual(values=c('#000000','#000000'))+
  ylab(expression("M"*Omega~"cm"))+
  scale_x_discrete(expand = c(0.04, 0.04))+
  
  theme(
    text = element_text(size=15),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size=18.5, hjust = 0.5),
    legend.position = "none",
    axis.line.x.bottom = element_line(size = 1),
    axis.title.x = element_blank(),
    axis.ticks.x.bottom = element_line(size = 0.90),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 15),
    legend.title = element_blank(),
    line = element_blank(),
    legend.key.size = unit(38,'points'),
    plot.margin = margin(10,2.5,5,17.5)
  )

#-----------------------------------------------------------------------------------------------------------------------------

#ORP

ORP.abc <- data.frame(abc.plot[,c(1,2,16,41,66)])
ORP.names <- c("Sample_ID","Location","M","P","A")
colnames(ORP.abc) <- ORP.names

orp <- ORP.abc %>%
  ggparcoord(
    columns = 3:5, groupColumn = 2,
    showPoints = FALSE, 
    title = "ORP",
    alphaLines = 0,
    scale = "globalminmax"
  ) + 
  theme_classic()+
  geom_line(aes( linetype=Location, color=Location), size = 0.25)+
  geom_point(size = 0.65)+
  geom_signif(comparisons = list(c("M", "P")), map_signif_level = c(" "=0.50), textsize = 8, tip_length = 0.02, size = 0.075)+
  scale_linetype_manual(values=c("solid", "dashed"))+
  scale_color_manual(values=c('#000000','#000000'))+
  ylab("mV")+
  scale_x_discrete(expand = c(0.04, 0.04))+
  
  theme(
    text = element_text(size=15),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size=18.5, hjust = 0.5),
    legend.position = "none",
    axis.line.x.bottom = element_line(size = 1),
    axis.title.x = element_blank(),
    axis.ticks.x.bottom = element_line(size = 0.90),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 15),
    legend.title = element_blank(),
    line = element_blank(),
    legend.key.size = unit(38,'points'),
    plot.margin = margin(10,2.5,5,17.5)
  )

#-----------------------------------------------------------------------------------------------------------------------------

#Conductivity

conductivity.abc <- data.frame(abc.plot[,c(1,2,19,44,69)])
conductivity.names <- c("Sample_ID","Location","M","P","A")
colnames(conductivity.abc) <- conductivity.names
conductivity.abc.rm <- conductivity.abc #removal of outlier sample from rockies
conductivity.abc.rm[2,] <- NA

con <- conductivity.abc.rm %>%
  ggparcoord(
    columns = 3:5, groupColumn = 2,
    showPoints = FALSE, 
    title = "Conductivity",
    alphaLines = 0, 
    missing = 'exclude',
    scale = 'globalminmax' #must be one of 'std', 'robust', 'uniminmax', 'globalminmax', 'center', or 'centerObs',
  ) + 
  theme_classic()+
  geom_line(aes(linetype=Location, color=Location), size= 0.25)+
  geom_point(size = 0.65)+
  geom_signif(comparisons = list(c("M", "A")), test ='t.test', map_signif_level = c(" "=0.50), y_position = 12.25 ,textsize = 8, tip_length = 0.02, size = 0.075)+
  scale_linetype_manual(values=c("solid", "dashed"))+
  scale_color_manual(values=c('#000000','#000000'))+
  ylab("µS/cm")+
  scale_x_discrete(expand = c(0.04, 0.04))+
  
  theme(
    text = element_text(size=15),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size=18.5, hjust = 0.5),
    legend.position = "none",
    axis.line.x.bottom = element_line(size = 1),
    axis.title.x = element_blank(),
    axis.ticks.x.bottom = element_line(size = 0.90),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 15),
    legend.title = element_blank(),
    line = element_blank(),
    legend.key.size = unit(38,'points'),
    plot.margin = margin(10,2.5,5,17.5)
  )

#-----------------------------------------------------------------------------------------------------------------------------

#Algal Richness (Sobs.Algae)

algae.rich.abc <- data.frame(abc.plot[,c(1,2,9,34,59)])
algae.rich.names <- c("Sample_ID","Location","M","P","A")
colnames(algae.rich.abc) <- algae.rich.names

algae.rich <- algae.rich.abc %>%
  ggparcoord(
    columns = 3:5, groupColumn = 2,
    showPoints = FALSE, 
    title = "Richness (Algae)",
    alphaLines = 0, 
    scale = 'globalminmax' # 'centerObs'; must be one of 'std', 'robust', 'uniminmax', 'globalminmax', 'center', or 'centerObs'
  ) + 
  theme_classic()+
  geom_line(aes(linetype=Location, color=Location), size= 0.25)+
  geom_point(size = 0.65)+
  geom_signif(comparisons = list(c("P", "A")), test = 't.test', map_signif_level = c(" "=0.50), textsize = 8, tip_length = 0.02, size = 0.075)+
  scale_linetype_manual(values=c("solid", "dashed"))+
  scale_color_manual(values=c('#000000','#000000'))+
  ylab("# of OTU")+
  scale_x_discrete(expand = c(0.04, 0.04))+
  
  theme(
    text = element_text(size=15),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 15, hjust = 0.5),
    legend.position = "none",
    axis.line.x.bottom = element_line(size = 1),
    axis.title.x = element_blank(),
    axis.ticks.x.bottom = element_line(size = 0.90),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 15),
    legend.title = element_blank(),
    line = element_blank(),
    legend.key.size = unit(38,'points'),
    plot.margin = margin(10,2.5,5,17.5)
  )

#-----------------------------------------------------------------------------------------------------------------------------

#Bacterial Richness (Sobs.Bac)

bac.rich.abc <- data.frame(abc.plot[,c(1,2,12,37,62)])
bac.rich.names <- c("Sample_ID","Location","M","P","A")
colnames(bac.rich.abc) <- bac.rich.names

bac.rich <- bac.rich.abc %>%
  ggparcoord(
    columns = 3:5, groupColumn = 2,
    showPoints = FALSE, 
    title = "Richness (Bacteria)",
    alphaLines = 0, 
    scale = 'globalminmax' # 'centerObs'; must be one of 'std', 'robust', 'uniminmax', 'globalminmax', 'center', or 'centerObs'
  ) + 
  theme_classic()+
  geom_line(aes(linetype=Location, color=Location), size= 0.25)+
  geom_point(size = 0.65)+
  geom_signif(comparisons = list(c("P", "A")), test = 't.test', map_signif_level = c(" "=0.50), textsize = 8, tip_length = 0.02, size = 0.075)+
  scale_linetype_manual(values=c("solid", "dashed"))+
  scale_color_manual(values=c('#000000','#000000'))+
  ylab("# of OTU")+
  scale_x_discrete(expand = c(0.04, 0.04))+
  
  theme(
    text = element_text(size=15),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 15, hjust = 0.5),
    legend.position = "none",
    axis.line.x.bottom = element_line(size = 1),
    axis.title.x = element_blank(),
    axis.ticks.x.bottom = element_line(size = 0.90),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 15),
    legend.title = element_blank(),
    line = element_blank(),
    legend.key.size = unit(38,'points'),
    plot.margin = margin(10,2.5,5,17.5)
  )

#-----------------------------------------------------------------------------------------------------------------------------

# 1-Simpson's Diversity (Algae)

simp.div.alg.abc <- data.frame(abc.plot[,c(1,2,11,36,61)])
simp.div.alg.names <- c("Sample_ID","Location","M","P","A")
colnames(simp.div.alg.abc) <- simp.div.alg.names

simp.div.alg <- simp.div.alg.abc %>%
  ggparcoord(
    columns = 3:5, groupColumn = 2,
    showPoints = FALSE, 
    title = "1-Simpson's Diversity (Algae)",
    alphaLines = 0, 
    scale = 'globalminmax' # 'centerObs'; must be one of 'std', 'robust', 'uniminmax', 'globalminmax', 'center', or 'centerObs'
  ) + 
  theme_classic()+
  geom_line(aes(linetype=Location, color=Location), size= 0.25)+
  geom_point(size = 0.65)+
  geom_signif(comparisons = list(c("P", "A")), test = 't.test', map_signif_level = c(" "=0.50), textsize = 8, tip_length = 0.02, size = 0.075)+
  scale_linetype_manual(values=c("solid", "dashed"))+
  scale_color_manual(values=c('#000000','#000000'))+
  ylab("1 - D")+
  scale_x_discrete(expand = c(0.04, 0.04))+
  
  theme(
    text = element_text(size=15),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size= 15, hjust = 0.5),
    legend.position = "none",
    axis.line.x.bottom = element_line(size = 1),
    axis.title.x = element_blank(),
    axis.ticks.x.bottom = element_line(size = 0.90),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 15),
    legend.title = element_blank(),
    line = element_blank(),
    legend.key.size = unit(38,'points'),
    plot.margin = margin(10,2.5,5,17.5)
  )

#-----------------------------------------------------------------------------------------------------------------------------

#Simpson's Evenness Community

simp.even.com.abc <- data.frame(abc.plot[,c(1,2,4,29,55)])
simp.even.com.names <- c("Sample_ID","Location","M","P","A")
colnames(simp.even.com.abc) <- simp.even.com.names

simp.even.com <- simp.even.com.abc %>%
  ggparcoord(
    columns = 3:5, groupColumn = 2,
    showPoints = FALSE, 
    title = "Simpson Evenness (Community)",
    alphaLines = 0, 
    scale = 'globalminmax' # 'centerObs'; must be one of 'std', 'robust', 'uniminmax', 'globalminmax', 'center', or 'centerObs'
  ) + 
  theme_classic()+
  geom_line(aes(linetype=Location, color=Location), size= 0.25)+
  geom_point(size = 0.65)+
  geom_signif(comparisons = list(c("P", "A")), test = 't.test', map_signif_level = c(" "=0.50), textsize = 8, tip_length = 0.02, size = 0.075)+
  scale_linetype_manual(values=c("solid", "dashed"))+
  scale_color_manual(values=c('#000000','#000000'))+
  ylab("Relative Abundance")+
  scale_x_discrete(expand = c(0.04, 0.04))+
  
  theme(
    text = element_text(size=15),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size=18.5, hjust = 0.5),
    legend.position = "none",
    axis.line.x.bottom = element_line(size = 1),
    axis.title.x = element_blank(),
    axis.ticks.x.bottom = element_line(size = 0.90),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 15),
    legend.title = element_blank(),
    line = element_blank(),
    legend.key.size = unit(38,'points'),
    plot.margin = margin(10,2.5,5,17.5)
  )

#-----------------------------------------------------------------------------------------------------------------------------

#Pollen Concentration

pollen.abc <- data.frame(abc.plot[,c(1,2,27,52,77)])
pollen.names <- c("Sample_ID","Location","M","P","A")
colnames(pollen.abc) <- pollen.names

pol <- pollen.abc %>%
  ggparcoord(
    columns = 3:5, groupColumn = 2,
    showPoints = FALSE, 
    title = "[Pollen]",
    alphaLines = 0, 
    scale = 'globalminmax' # 'centerObs'; must be one of 'std', 'robust', 'uniminmax', 'globalminmax', 'center', or 'centerObs'
  ) + 
  theme_classic()+
  geom_line(aes(linetype=Location, color=Location), size= 0.25)+
  geom_point(size = 0.65)+
  geom_signif(comparisons = list(c("M", "P")), test = 't.test', map_signif_level = c(" "=0.50), textsize = 8, tip_length = 0.02, size = 0.075)+
  geom_signif(comparisons = list(c("M", "A")), test = 't.test', map_signif_level = c(" "=0.50), textsize = 8, tip_length = 0.02, size = 0.075, y_position = 4000)+
  scale_linetype_manual(values=c("solid", "dashed"))+
  scale_color_manual(values=c('#000000','#000000'))+
  ylab("Pollen Cells/mL")+
  scale_x_discrete(expand = c(0.04, 0.04))+
  
  theme(
    text = element_text(size=15),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size=18.5, hjust = 0.5),
    legend.position = "none",
    axis.line.x.bottom = element_line(size = 1),
    axis.title.x = element_blank(),
    axis.ticks.x.bottom = element_line(size = 0.90),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 15),
    legend.title = element_blank(),
    line = element_blank(),
    legend.key.size = unit(38,'points'),
    plot.margin = margin(10,2.5,5,17.5)
  )

#-----------------------------------------------------------------------------------------------------------------------------

#Algal Cell Concentration

cell.abc <- data.frame(abc.plot[,c(1,2,26,51,76)])
cell.names <- c("Sample_ID","Location","M","P","A")
colnames(cell.abc) <- cell.names
cell.abc.rm <- cell.abc #removal of outlier sample 
cell.abc.rm[7,] <- NA

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE, digits = 0)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}

cell <- cell.abc.rm %>%
  ggparcoord(
    columns = 3:5, groupColumn = 2,
    showPoints = FALSE, 
    title = "[Snow Algae]",
    alphaLines = 0, 
    scale = 'globalminmax' # 'centerObs'; must be one of 'std', 'robust', 'uniminmax', 'globalminmax', 'center', or 'centerObs'
  ) + 
  theme_classic()+
  geom_line(aes(linetype=Location, color=Location), size= 0.25)+
  geom_point(size = 0.65)+
  geom_signif(comparisons = list(c("M", "A")), test = 't.test', map_signif_level = c(" "=0.50), textsize = 8, tip_length = 0.02, size = 0.075, y_position = 203500)+
  geom_signif(comparisons = list(c("P", "A")), test = 't.test', map_signif_level = c(" "=0.50), textsize = 8, tip_length = 0.02, size = 0.075)+
  geom_signif(comparisons = list(c("M", "P")), test = 't.test', map_signif_level = c(" "=0.50), textsize = 8, tip_length = 0.02, size = 0.075)+
  scale_linetype_manual(values=c("solid", "dashed"))+
  scale_color_manual(values=c('#000000','#000000'))+
  ylab("Snow Algae Cells/mL")+
  scale_x_discrete(expand = c(0.04, 0.04))+
  scale_y_continuous(labels = fancy_scientific)+
  
  theme(
    text = element_text(size=15),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size=18.5, hjust = 0.5),
    legend.position = "none",
    axis.line.x.bottom = element_line(size = 1),
    axis.title.x = element_blank(),
    axis.ticks.x.bottom = element_line(size = 0.90),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 13),
    legend.title = element_blank(),
    line = element_blank(),
    legend.key.size = unit(38,'points'),
    plot.margin = margin(10,2.5,5,17.5)
  )

#-----------------------------------------------------------------------------------------------------------------------------

#Dissolved Oxygen (DO%)
DO.abc <- data.frame(abc.plot[,c(1,2,18,42,67)])
DO.names <- c("Sample_ID","Location","M","P","A")
colnames(DO.abc) <- DO.names

DO <- DO.abc %>%
  ggparcoord(
    columns = 3:5, groupColumn = 2,
    showPoints = FALSE, 
    title = "DO",
    alphaLines = 0,
    scale = "globalminmax"
  ) + 
  theme_classic()+
  geom_line(aes( linetype=Location, color=Location), size = 0.25)+
  geom_point(size = 0.65)+
  
  scale_linetype_manual(values=c("solid", "dashed"))+
  scale_color_manual(values=c('#000000','#000000'))+
  ylab("%")+
  scale_x_discrete(expand = c(0.04, 0.04))+
  
  theme(
    text = element_text(size=15),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size=18.5, hjust = 0.5),
    legend.position = "none",
    axis.line.x.bottom = element_line(size = 1),
    axis.title.x = element_blank(),
    axis.ticks.x.bottom = element_line(size = 0.90),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 15),
    legend.title = element_blank(),
    line = element_blank(),
    legend.key.size = unit(38,'points'),
    plot.margin = margin(10,2.5,5,17.5)
  )
DO

#Dissolved Oxygen (DO ppm)
DO.abc <- data.frame(abc.plot[,c(1,2,17,43,68)])
DO.names <- c("Sample_ID","Location","M","P","A")
colnames(DO.abc) <- DO.names

DO <- DO.abc %>%
  ggparcoord(
    columns = 3:5, groupColumn = 2,
    showPoints = FALSE, 
    title = "DO",
    alphaLines = 0,
    scale = "globalminmax"
  ) + 
  theme_classic()+
  geom_line(aes( linetype=Location, color=Location), size = 0.25)+
  geom_point(size = 0.65)+
  
  scale_linetype_manual(values=c("solid", "dashed"))+
  scale_color_manual(values=c('#000000','#000000'))+
  ylab("ppm")+
  scale_x_discrete(expand = c(0.04, 0.04))+
  
  theme(
    text = element_text(size=15),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size=18.5, hjust = 0.5),
    legend.position = "none",
    axis.line.x.bottom = element_line(size = 1),
    axis.title.x = element_blank(),
    axis.ticks.x.bottom = element_line(size = 0.90),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 15),
    legend.title = element_blank(),
    line = element_blank(),
    legend.key.size = unit(38,'points'),
    plot.margin = margin(10,2.5,5,17.5)
  )
DO

#-----------------------------------------------------------------------------------------------------------------------------
#while nitrate and simp.div.alg were found to be statistically significant in our Wilcoxon-Rank Sum Test
  #I chose to exclude them from the final graphic as a true trend was not clear and may have been influenced
  #by outliers. Resistivity was also excluded as it can be easily derived/inferred from conductivity measures.
  #DO was included as it shows a clear trend, and barely makes the cutoff for significance

#For cell concentration and conductivity a single sample was removed from each as it was considered an outlier.
  #This was done to better express the trends of the remaining data points.


grid <- grid.arrange(orp, DO, con, potassium, pol, cell, algae.rich, bac.rich, simp.div.alg, ncol = 3, nrow = 3)
ggsave(filename = "full.parallel2.png", plot = grid, height = 12.5, width = 12, dpi = 1500)
