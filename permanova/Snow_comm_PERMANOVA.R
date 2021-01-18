install.packages("vegan")
install.packages("RVAideMemoire")
library(vegan, RVAideMemoire)

#Universal Map File
setwd("/Users/BrownLab/Desktop/Students/Avery/Snow_Community_Project/Community")
map<-read.csv(file="com_Snow_Community_PERMANOVA_MAP.csv", header=TRUE, row.names=c(1))

#Algae Bray-Curtis
setwd("/Users/BrownLab/Desktop/Students/Avery/Snow_Community_Project/Algae")
dist<-read.csv(file="algae_bray_square.csv", header=TRUE, row.names = c(1))
dist<-as.dist(dist)

#Fungi Bray-Curtis
setwd("/Users/BrownLab/Desktop/Students/Avery/Snow_Community_Project/Fungi")
dist<-read.csv(file="fungi_bray_square.csv", header=TRUE, row.names = c(1))
dist<-as.dist(dist)

#Bacteria Bray-Curtis
setwd("/Users/BrownLab/Desktop/Students/Avery/Snow_Community_Project/Bac")
dist<-read.csv(file="bac_bray_square.csv", header=TRUE, row.names = c(1))
dist<-as.dist(dist)

#Community Bray-Curtis
setwd("/Users/BrownLab/Desktop/Students/Avery/Snow_Community_Project/Community")
dist<-read.csv(file="com_bray_square.csv", header=TRUE, row.names = c(1))
dist<-as.dist(dist)

#Calculate PERMANOVA Values for Bray-Curtis (take the top value of each set, order matters here)
adonis(dist ~ Location + pH + ORP + DO + Conductivity + Resistivity + TDS + Salinity + Nitrate + Ammonium + Potassium + Cell_Count + Pollen, data=map)
adonis(dist ~ pH + Location + ORP + DO + Conductivity + Resistivity + TDS + Salinity + Nitrate + Ammonium + Potassium + Cell_Count + Pollen, data=map)
adonis(dist ~ ORP + pH + Location + DO + Conductivity + Resistivity + TDS + Salinity + Nitrate + Ammonium + Potassium + Cell_Count + Pollen, data=map)
adonis(dist ~ DO + ORP + pH + Location + Conductivity + Resistivity + TDS + Salinity + Nitrate + Ammonium + Potassium + Cell_Count + Pollen, data=map)
adonis(dist ~ Conductivity + DO + ORP + pH + Location + Resistivity + TDS + Salinity + Nitrate + Ammonium + Potassium + Cell_Count + Pollen, data=map)
adonis(dist ~ Resistivity + Conductivity + DO + ORP + pH + Location + TDS + Salinity + Nitrate + Ammonium + Potassium + Cell_Count + Pollen, data=map)
adonis(dist ~ TDS + Resistivity + Conductivity + DO + ORP + pH + Location + Salinity + Nitrate + Ammonium + Potassium + Cell_Count + Pollen, data=map)
adonis(dist ~ Salinity + TDS + Resistivity + Conductivity + DO + ORP + pH + Location + Nitrate + Ammonium + Potassium + Cell_Count + Pollen, data=map)
adonis(dist ~ Nitrate + Salinity + TDS + Resistivity + Conductivity + DO + ORP + pH + Location + Ammonium + Potassium + Cell_Count + Pollen, data=map)
adonis(dist ~ Ammonium + Nitrate + Salinity + TDS + Resistivity + Conductivity + DO + ORP + pH + Location + Potassium + Cell_Count + Pollen, data=map)
adonis(dist ~ Potassium + Ammonium + Nitrate + Salinity + TDS + Resistivity + Conductivity + DO + ORP + pH + Location + Cell_Count + Pollen, data=map)
adonis(dist ~ Cell_Count + Potassium + Ammonium + Nitrate + Salinity + TDS + Resistivity + Conductivity + DO + ORP + pH + Location + Pollen, data=map)
adonis(dist ~ Pollen + Cell_Count + Potassium + Ammonium + Nitrate + Salinity + TDS + Resistivity + Conductivity + DO + ORP + pH + Location, data=map)


#install John Aitchison's log-ratio approach to dissimilarity measures
install.packages("coda.base")
install.packages("compositions")
library(compositions)
library(coda.base)

#Algae Aitchison Values (better suited for compositional data; see Gloor et al. 2017)
setwd("/Users/BrownLab/Desktop/Students/Avery/Snow_Community_Project/Algae")
algae <- read.delim("algae.FINAL.shared")
algae <- as.matrix(algae[1:31,4:94])
algae.clr <- clr(algae)
dist.aitch <- coda.base::dist(x = algae.clr, method = "aitchison", diag = TRUE)
dist.aitch

#Bacteria Aitchison Values
#Fungi Aitchison Values
#Community Aitchison Values



adonis(dist.aitch ~ Location + pH + ORP + DO + Conductivity + Resistivity + TDS + Salinity + Nitrate + Ammonium + Potassium + Cell_Count + Pollen, data=map)
adonis(dist.aitch ~ pH + Location + ORP + DO + Conductivity + Resistivity + TDS + Salinity + Nitrate + Ammonium + Potassium + Cell_Count + Pollen, data=map)
adonis(dist.aitch ~ ORP + pH + Location + DO + Conductivity + Resistivity + TDS + Salinity + Nitrate + Ammonium + Potassium + Cell_Count + Pollen, data=map)
adonis(dist.aitch ~ DO + ORP + pH + Location + Conductivity + Resistivity + TDS + Salinity + Nitrate + Ammonium + Potassium + Cell_Count + Pollen, data=map)
adonis(dist.aitch ~ Conductivity + DO + ORP + pH + Location + Resistivity + TDS + Salinity + Nitrate + Ammonium + Potassium + Cell_Count + Pollen, data=map)
adonis(dist.aitch ~ Resistivity + Conductivity + DO + ORP + pH + Location + TDS + Salinity + Nitrate + Ammonium + Potassium + Cell_Count + Pollen, data=map)
adonis(dist.aitch ~ TDS + Resistivity + Conductivity + DO + ORP + pH + Location + Salinity + Nitrate + Ammonium + Potassium + Cell_Count + Pollen, data=map)
adonis(dist.aitch ~ Salinity + TDS + Resistivity + Conductivity + DO + ORP + pH + Location + Nitrate + Ammonium + Potassium + Cell_Count + Pollen, data=map)
adonis(dist.aitch ~ Nitrate + Salinity + TDS + Resistivity + Conductivity + DO + ORP + pH + Location + Ammonium + Potassium + Cell_Count + Pollen, data=map)
adonis(dist.aitch ~ Ammonium + Nitrate + Salinity + TDS + Resistivity + Conductivity + DO + ORP + pH + Location + Potassium + Cell_Count + Pollen, data=map)
adonis(dist.aitch ~ Potassium + Ammonium + Nitrate + Salinity + TDS + Resistivity + Conductivity + DO + ORP + pH + Location + Cell_Count + Pollen, data=map)
adonis(dist.aitch ~ Cell_Count + Potassium + Ammonium + Nitrate + Salinity + TDS + Resistivity + Conductivity + DO + ORP + pH + Location + Pollen, data=map)
adonis(dist.aitch ~ Pollen + Cell_Count + Potassium + Ammonium + Nitrate + Salinity + TDS + Resistivity + Conductivity + DO + ORP + pH + Location, data=map)
