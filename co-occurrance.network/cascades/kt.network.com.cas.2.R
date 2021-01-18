################################################################################################################################

#Kendall-Tau Microbial Correlation Network Analysis
#Shawn Brown Lab University of Memphis
#Modified and Created by Avery Tucker

################################################################################################################################
#Set Directory

setwd("/Users/BrownLab/Desktop/Students/Avery/projects/snow.bloom.2018/stats.analysis/kendall.tau.network/cascades")

################################################################################################################################
#Install Packages

#Install Bioconductor and any associated packages in developer mode
BiocManager::install(pkgs = c("phyloseq"))

install.packages("devtools")
install.packages("GGally")
install.packages("ggraph")
install.packages("igraph")
install.packages("intergraph")
install.packages("Matrix")
install.packages("network")
install_github("zdk123/SpiecEasi")
install.packages("sna")
install.packages("stringr")
install.packages("WGCNA")

library(devtools)
library(GGally)
library(ggraph)
library(igraph)
library(intergraph)
library(Matrix)
library(network)
library(phyloseq)
library(SpiecEasi)
library(sna)
library(stringr)
library(WGCNA)

################################################################################################################################
#Read in taxa names and make a phyloseq object
taxa <- as.matrix(read.csv("all.lineage.cascade.csv", row.names = 1))
phylo.taxa <- tax_table(taxa)
phyloseq <- phyloseq(phylo.taxa)

#Read in snow community data matrix
abund <- as.matrix(read.csv("all.abund.cascade.csv", row.names = 1))
data <- t(abund)

#Check that all OTU names are represented in each column (should equal 0)
a <- row.names(taxa)
b <- row.names(abund)
a[! (a %in% b)]

################################################################################################################################
#igraph Network - Full Community

#Generate kt correlation matrix
kt.cor <- as.matrix(cor(data, use = "p", method = c("kendall")))

#Select and parse correlation value threshold
correlation.threshold <- 0.4
kt.thresh <- matrix(data = 0, nrow = length(kt.cor[,1]), 
                    ncol = length(kt.cor[1,]), 
                    dimnames = list(row.names(kt.cor), colnames(kt.cor)))

for(i in 1:nrow(kt.cor)) {
  for(j in 1:ncol(kt.cor)) {
    if(abs(kt.cor[i,j]) >= correlation.threshold) {
      kt.thresh[i,j] <- (kt.cor[i,j])
    }
  }
}

#Generate kt p-value matrix
nSamples = nrow(data)
kt.p = as.matrix(corPvalueStudent(as.matrix(kt.cor), nSamples))

#Select and parse p-value threshold
p.value.threshold <- 0.1 # set the p-value threshold

for(i in 1:nrow(kt.p)) {
  for(j in 1:ncol(kt.p)) {
    if(abs(kt.p[i,j]) >= p.value.threshold) {
      kt.thresh[i,j] <- 0
    }
  }
}

#Generate sparse matrix thresholded for correlation score
diag(kt.thresh) <- 0
kt.thresh <- Matrix(kt.thresh, sparse = TRUE)

### Create igraph objects
ig.thresh <- adj2igraph(kt.thresh,
                        vertex.attr = list(name = taxa_names(phyloseq), 
                                           rmEmptyNodes = FALSE, 
                                           diag = FALSE))

################################################################################################################################
#igraph Full Network Visualizations (basic)

#Visualize using igraph plotting
plot(ig.thresh, layout = layout.fruchterman.reingold, vertex.label=NA, main="Community Correlation Network")

#Histogram visualization of edge weights
edge.list <- summary(kt.thresh*kt.cor)
hist(edge.list[,3], main='Edge Weight Distribution', xlab = 'Edge Weights', xlim = range(-1:1))

#Degree(s) statistics
dd.spar <- degree.distribution(ig.thresh)
plot(0:(length(dd.spar)-1), dd.spar, ylim=c(0,0.1), type='b', ylab="Frequency", xlab="Degree(s)", main="Degree Distributions")

################################################################################################################################
#igraph to GEPHI workflow
#Adapted from https://gist.github.com/Vessy/6047440

#Start with our community igraph network
gD <- ig.thresh

#Create the following categories for GEPHI: nodes, edges, nodesAtt, edgesAtt, nodesVizAtt, edgesVizAtt, edgesWeight
#Node and edge attributes may be useful for other network manipulations in Gephi
#Node/edge visual attributes are used for network visualization

#Nodes
nodes_df <- data.frame(ID = c(1:vcount(gD)), NAME = V(gD)$name)

#Edges
edges_df <- as.data.frame(get.edges(gD, c(1:ecount(gD))))

#NodesAtt
V(gD)$degree <- degree(gD)
nodes_att <- data.frame(DEG = V(gD)$degree)#, BET = V(gD)$betweenness) 

#EdgesAtt
edges_att <- data.frame(WGH = E(gD)$weight) 

#NodesVizAtt
#Node Coordinates
#2D coordinates result in a better (2D) plot than 3D coordinates
nodes_coord <- as.data.frame(layout.fruchterman.reingold(gD, weights = abs(E(gD)$weight), dim = 2, niter = 10000))
nodes_coord <- cbind(nodes_coord, rep(0, times = nrow(nodes_coord)))

#Calculate node size
#And we will assign a node size for each node based on its community relative abundance
nodes_size <- rowSums(abund)/sum(abund)

#Define node color
#We'll interpolate node colors based on the node degree using the "colorRampPalette" function from the "grDevices" library
library("grDevices")

#Determine node color - here I use red for bacteria, blue for fungi, and green for algae (hexadecimal based)
v <- as_ids(V(gD))
nodes_col <- vector(length = length(v))
for(i in 1:length(nodes_col)) {
  if(str_detect(v[i], "B")) {
    nodes_col[i] <- c("#ff0000")
  }
  else if(str_detect(v[i], "A")) {
    nodes_col[i] <- c("#2fff00")
  }
  else {
    nodes_col[i] <- c("#0033ff")
  }
}
#Transform it into a data frame (we have to transpose it first)
nodes_col_df <- as.data.frame(t(col2rgb(nodes_col, alpha = FALSE)))
#And add alpha (between 0 and 1). The alpha from "col2rgb" function takes values from 0-255, so we cannot use it
nodes_col_df <- cbind(nodes_col_df, alpha = rep(1, times = nrow(nodes_col_df)))
#Assign visual attributes to nodes (colors have to be 4dimensional - RGBA)
nodes_att_viz <- list(color = nodes_col_df, position = nodes_coord, size = nodes_size)


#EdgesVizAtt
#Assign visual attributes to edges
edges_col <- vector(length = length(E(gD)$weight))
for(i in 1:length(E(gD)$weight)) {
  if(E(gD)$weight[i] > 0) {
    edges_col[i] <- c("#000000")
  }
  else {
    edges_col[i] <- c("#de1f12")
  }
}
edges_col_df <- as.data.frame(t(col2rgb(edges_col, alpha = FALSE)))
edges_col_df <- cbind(edges_col_df, alpha = rep(1, times = nrow(edges_col_df)))
edges_att_viz <-list(color = edges_col_df)


# Write the network into a gexf (Gephi) file
write.gexf(nodes = nodes_df, edges = edges_df, edgesWeight = abs(E(gD)$weight), 
           nodesVizAtt = nodes_att_viz, edgesVizAtt = edges_att_viz, edgesAtt = edges_att, nodesAtt = nodes_att, 
           defaultedgetype = "undirected", 
           output = "snow.com.rockies.gexf")

################################################################################################################################
#Generation of Network Subplots

#Plot Sanguina nivaloides network (AOtu0001)
otu <- c("AOtu0001")
title.name <- "Sanguina nivaloides (AOtu0001) (Cascades)"
otu.sub <- induced_subgraph(ig.thresh, vids = unlist(neighborhood(graph=ig.thresh, order=1, nodes=otu)))
transitivity(otu.sub, type = "localaverage", isolates = "zero")
average.path.length(otu.sub)
E(otu.sub)
V(otu.sub)

pdf(file = paste(title.name, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 10) # The height of the plot in inches

plot_network_2(otu.sub, physeq=phyloseq, color='Class', shape=NULL, title= title.name, abundance=abund,
               type="taxa", label="value", layout.method=layout.star,
               line_alpha=0.4, line_color=color, line_weight=TRUE, hjust = 1, 
               point_alpha=0.95, point_label_size=5, otu.center = otu)
dev.off()

#Plot S. nivaloides network (AOtu0002)
otu <- c("AOtu0002")
title.name <- "Sanguina nivaloides (AOtu0002) (Cascades)"
otu.sub <- induced_subgraph(ig.thresh, vids = unlist(neighborhood(graph=ig.thresh, order=1, nodes=otu)))
transitivity(otu.sub, type = "localaverage", isolates = "zero")
average.path.length(otu.sub)
E(otu.sub)
V(otu.sub)

pdf(file = paste(title.name, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 10) # The height of the plot in inches

plot_network_2(otu.sub, physeq=phyloseq, color='Class', shape=NULL, title= title.name, abundance=abund,
               type="taxa", label="value", layout.method=layout.star,
               line_alpha=0.4, line_color=color, line_weight=TRUE, hjust = 1, 
               point_alpha=0.95, point_label_size=5, otu.center = otu)
dev.off()

#Plot Soletalia (BOtu0001)
otu <- c("BOtu0001")
title.name <- "Soletalia sp. (BOtu0001) (Cascades)"
otu.sub <- induced_subgraph(ig.thresh, vids = unlist(neighborhood(graph=ig.thresh, order=1, nodes=otu)))
transitivity(otu.sub, type = "localaverage", isolates = "zero")
average.path.length(otu.sub)
E(otu.sub)
V(otu.sub)

pdf(file = paste(title.name, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 14, # The width of the plot in inches
    height = 11) # The height of the plot in inches

plot_network_2(otu.sub, physeq=phyloseq, color='Class', shape=NULL, title= title.name, abundance=abund,
               type="taxa", label="value", layout.method=layout.graphopt,
               line_alpha=0.4, line_color=color, line_weight=TRUE, hjust = 1, 
               point_alpha=0.95, point_label_size=2.5)
dev.off()

#Plot Hymnobacter Network (BOtu0002)
otu <- c("BOtu0002")
title.name <- "Hymnobacter sp. (BOtu0002) OTU (Cascades)"
otu.sub <- induced_subgraph(ig.thresh, vids = unlist(neighborhood(graph=ig.thresh, order=1, nodes=otu)))
transitivity(otu.sub, type = "localaverage", isolates = "zero")
average.path.length(otu.sub)
E(otu.sub)
V(otu.sub)

pdf(file = paste(title.name, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 10) # The height of the plot in inches

plot_network_2(otu.sub, physeq=phyloseq, color='Class', shape=NULL, title= title.name, abundance=abund,
               type="taxa", label="value", layout.method=layout.star,
               line_alpha=0.4, line_color=color, line_weight=TRUE, hjust = 1, 
               point_alpha=0.95, point_label_size=5, otu.center = otu)
dev.off()

#Plot Ferruginibacter Network (BOtu0003)
otu <- c("BOtu0003")
title.name <- "Ferruginibacter sp. (BOtu0003) (Cascades)"
otu.sub <- induced_subgraph(ig.thresh, vids = unlist(neighborhood(graph=ig.thresh, order=1, nodes=otu)))
transitivity(otu.sub, type = "localaverage", isolates = "zero")
average.path.length(otu.sub)
E(otu.sub)
V(otu.sub)

pdf(file = paste(title.name, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 10) # The height of the plot in inches

plot_network_2(otu.sub, physeq=phyloseq, color='Class', shape=NULL, title= title.name, abundance=abund,
               type="taxa", label="value", layout.method=layout.star,
               line_alpha=0.4, line_color=color, line_weight=TRUE, hjust = 1, 
               point_alpha=0.95, point_label_size=5, otu.center = otu)
dev.off()

#Plot Polaromonas Network (BOtu0004)
otu <- c("BOtu0004")
title.name <- "Polaromonas sp. (BOtu0004) (Cascades)"
otu.sub <- induced_subgraph(ig.thresh, vids = unlist(neighborhood(graph=ig.thresh, order=1, nodes=otu)))
transitivity(otu.sub, type = "localaverage", isolates = "zero")
average.path.length(otu.sub)
E(otu.sub)
V(otu.sub)

pdf(file = paste(title.name, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 10) # The height of the plot in inches

plot_network_2(otu.sub, physeq=phyloseq, color='Class', shape=NULL, title= title.name, abundance=abund,
               type="taxa", label="value", layout.method=layout.star,
               line_alpha=0.4, line_color=color, line_weight=TRUE, hjust = 1, 
               point_alpha=0.95, point_label_size=5, otu.center = otu)
dev.off()

#Plot unclassified Chytridiomycota (FOtu0001)
otu <- c("FOtu0001")
title.name <- "unclassified Chytridiomycota (FOtu0001) (Cascades)"
otu.sub <- induced_subgraph(ig.thresh, vids = unlist(neighborhood(graph=ig.thresh, order=1, nodes=otu)))
transitivity(otu.sub, type = "localaverage", isolates = "zero")
average.path.length(otu.sub)
E(otu.sub)
V(otu.sub)

pdf(file = paste(title.name, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 14.5, # The width of the plot in inches
    height = 10) # The height of the plot in inches

plot_network_2(otu.sub, physeq=phyloseq, color='Class', shape=NULL, title= title.name, abundance=abund,
               type="taxa", label="value", layout.method=layout.davidson.harel,
               line_alpha=0.4, line_color=color, line_weight=TRUE, hjust = 1, 
               point_alpha=0.95, point_label_size=4, otu.center = otu)
dev.off()

#Plot unclassified Microbotryomycete (FOtu0002)
otu <- c("FOtu0002")
title.name <- "unclassified Microbotryomycete \n(FOtu0002) (Cascades)"
otu.sub <- induced_subgraph(ig.thresh, vids = unlist(neighborhood(graph=ig.thresh, order=1, nodes=otu)))
transitivity(otu.sub, type = "localaverage", isolates = "zero")
average.path.length(otu.sub)
E(otu.sub)
V(otu.sub)

pdf(file = paste(title.name, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 10) # The height of the plot in inches

plot_network_2(otu.sub, physeq=phyloseq, color='Class', shape=NULL, title= title.name, abundance=abund,
               type="taxa", label="value", layout.method=layout.star,
               line_alpha=0.4, line_color=color, line_weight=TRUE, hjust = 1, 
               point_alpha=0.95, point_label_size=5, otu.center = otu)
dev.off()

#Plot Leucosporidiales unclassified (FOtu0003) Network
otu <- c("FOtu0003")
title.name <- "Leucosporidiaceae unclassified (FOtu0003) (Cascades)"
g <- induced_subgraph(ig.thresh, vids = unlist(neighborhood(graph=ig.thresh, order=1, nodes=otu)))
transitivity(otu.sub, type = "localaverage", isolates = "zero")
average.path.length(otu.sub)
E(otu.sub)
V(otu.sub)

pdf(file = paste(title.name, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 22, # The width of the plot in inches
    height = 12) # The height of the plot in inches

plot_network_2(g, physeq=phyloseq, color='Class', shape=NULL, title= title.name, abundance=abund,
               type="taxa", label="value", layout.method=layout.graphopt,
               line_alpha=0.4, line_color=color, line_weight=TRUE, hjust = 1, 
               point_alpha=0.95, point_label_size=5)
dev.off()

#Plot Leucosporidiales unclassified (FOtu0005) Network
otu <- c("FOtu0005")
title.name <- "Leucosporidiales unclassified \n(FOtu0005) (Cascades)"
otu.sub <- induced_subgraph(ig.thresh, vids = unlist(neighborhood(graph=ig.thresh, order=1, nodes=otu)))
transitivity(otu.sub, type = "localaverage", isolates = "zero")
average.path.length(otu.sub)
E(otu.sub)
V(otu.sub)

pdf(file = paste(title.name, ".pdf", sep = ""),   # The directory you want to save the file in
    width = 14, # The width of the plot in inches
    height = 10) # The height of the plot in inches

plot_network_2(otu.sub, physeq=phyloseq, color='Class', shape=NULL, title= title.name, abundance=abund,
               type="taxa", label="value", layout.method=layout.star,
               line_alpha=0.4, line_color=color, line_weight=TRUE, hjust = 1, 
               point_alpha=0.95, point_label_size=5, otu.center = otu)
dev.off()
