#Taxon Abundance of Algae, Bacteria, and Fungi for Snow Community Analysis Using Group Abundance
#The only issue with using group.abundance is there is not a way to filter out <1% relative OTU
#abundance samples and add these back as bars, it may be better to plot in ggplot using phyloseq
#asterisks and circles denoting significance for rockies or cascades were added manually in PowerPoint after images were generated

setwd("/Users/BrownLab/Desktop/Students/Avery/Snow_Community_Project/Stats/Abundance")
install.packages("RAM")
library(RAM)
library(readxl)
library(tidyverse)

Algae <- read_excel("algae.FINAL.shared.METACODER.xlsx", 
                          sheet = "Algae Total")
Bacteria <- read_excel("bac.FINAL.shared.METACODER.xlsx", 
                          sheet = "Sheet2")
Fungi <- read_excel("fungi.FINAL.shared.METACODER.xlsx", 
                          sheet = "Sheet3")

########################################################################################################################
########################################################################################################################

#Algae Relative Abundance (Percent-Species)
algae.graph.percent <- group.abundance(data=list(Algae=Algae), 
                               rank="s", 
                               drop.unclassified = FALSE,
                               count = FALSE,
                               main = "",
                               ggplot2 = TRUE,
                               top = 13) +
  
  scale_y_continuous(expand = c(0,0), labels = function(x) paste0(x*100, "%")) +
  ylab("OTU Abundance (%)") +
  xlab("Sample") +
  geom_bar(width = 1, stat = "identity", show.legend = TRUE) +
  geom_vline(xintercept = c(3.5, 6.5, 9.5, 15.5, 18.5, 21.5, 24.5, 27.5), linetype = "dashed", size = 0.1, color = "black") +
  geom_vline(xintercept = 12.5, size = 0.3, color = "black") +
  
  theme(
  panel.background = element_rect(fill = "transparent") # bg of the panel
  , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
  , panel.grid.major = element_blank() # get rid of major grid
  , panel.grid.minor = element_blank() # get rid of minor grid
  , legend.background = element_rect(fill = "transparent") # get rid of legend bg
  , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  , strip.background = element_blank()
  , strip.text = element_blank()
)

algae.graph.percent
ggsave(algae.graph.percent, filename = "algae.abund.percent.species.1pecent.png", width=12, dpi = 1200)

########################################################################################################################
########################################################################################################################


########################################################################################################################
########################################################################################################################

#Bacteria Relative Abundance (Percent-Family)
bac.graph.percent <- group.abundance(data=list(Bacteria=Bacteria), 
                             rank="f", 
                             drop.unclassified = FALSE,
                             count = FALSE,
                             top = 11,
                             main = "",
                             ggplot2 = TRUE) +
  scale_y_continuous(expand = c(0,0), labels = function(x) paste0(x*100, "%")) +
  ylab("OTU Abundance (%)") +
  xlab("Sample") +
  geom_bar(width = 1, stat = "identity") +
  geom_vline(xintercept = c(3.5, 6.5, 9.5, 15.5, 18.5, 21.5, 24.5, 27.5), linetype = "dashed", size = 0.1, color = "black") +
  geom_vline(xintercept = 12.5, size = 0.3, color = "black") +
  
  theme(
    panel.background = element_rect(fill = "transparent") # bg of the panel
    , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
    , panel.grid.major = element_blank() # get rid of major grid
    , panel.grid.minor = element_blank() # get rid of minor grid
    , legend.background = element_rect(fill = "transparent") # get rid of legend bg
    , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    , strip.background = element_blank()
    , strip.text = element_blank()
  )

bac.graph.percent    
ggsave(bac.graph.percent, filename = "bac.abund.percent.family.1pecent.png", width = 12, dpi = 1200)

########################################################################################################################
########################################################################################################################

#Fungi Relative Abundance (Percent)
fungi.graph.percent <- group.abundance(data=list(Fungi=Fungi), 
                               rank="f", 
                               drop.unclassified = FALSE,
                               count = FALSE,
                               main = "",
                               ggplot2 = TRUE) +
  scale_y_continuous(expand = c(0,0), labels = function(x) paste0(x*100, "%")) +
  ylab("OTU Abundance (%)") +
  xlab("Sample") +
  geom_bar(width = 1, stat = "identity") +
  geom_vline(xintercept = c(3.5, 6.5, 9.5, 15.5, 18.5, 21.5, 24.5, 27.5), linetype = "dashed", size = 0.1, color = "black") +
  geom_vline(xintercept = 12.5, size = 0.3, color = "black") +
  
  theme(
    panel.background = element_rect(fill = "transparent") # bg of the panel
    , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
    , panel.grid.major = element_blank() # get rid of major grid
    , panel.grid.minor = element_blank() # get rid of minor grid
    , legend.background = element_rect(fill = "transparent") # get rid of legend bg
    , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    , strip.background = element_blank()
    , strip.text = element_blank()
  )

fungi.graph.percent     
ggsave(fungi.graph.percent, filename = "fungi.abund.percent.family.1percent.jpg", width = 12, dpi = 1500)
