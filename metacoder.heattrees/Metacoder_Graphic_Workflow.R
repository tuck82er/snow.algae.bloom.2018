# Heat Trees and Stats for Snow Community Samples 03-20-2019 - Avery Ezra Tucker
# Created by Zachary Foster during time at Grunwald Lab

# reference website (Grunwald Lab)
# https://grunwaldlab.github.io/metacoder_documentation/example.html


# websites with good discussion
# adding colors https://github.com/grunwaldlab/metacoder/issues/210
# great tool for finding good four color gradients http://tristen.ca/hcl-picker/#/hlc/6/1/AE6942/55F4FA
  # all gradients should have four colors with the final color ending in grey
# heat tree arguments https://rdrr.io/cran/metacoder/man/heat_tree.html
# hexadecimal color picker https://www.google.com/search?client=safari&rls=en&q=hexadecimal+color+wheel&ie=UTF-8&oe=UTF-8
# github https://github.com/cran/metacoder

# packages and directory
install.packages("taxa")
install.packages("vegan")
install.packages("tidyverse") # to create tibbled data frames
install.packages("metacoder")
library(metacoder)
library(tidyverse)
library(taxa)
library(vegan)
library(readxl)
library(tibble)
setwd("/Users/BrownLab/Desktop/Students/Avery/Snow_Community_Project/Stats/Metacoder/Algae")

######################################################################################################################################################################################
######################################################################################################################################################################################
# Individual algae taxonomic heat tree

# load a map file (modified PERMANOVA) and the OTU/lineage_string/sample file (Modified .shared from Mothur)
algae.total <- read_excel("algae.FINAL.shared.METACODER.xlsx", sheet = "Algae All Sites (No H2O)")
algae.total.map <- read_excel("algae_Snow_Community_PERMANOVA_MAP.xlsx", sheet = "Algae_Map (Numeric Names)")

obj.algae.total <- parse_tax_data(algae.total,
                      class_cols = "Lineage", # the column that contains taxonomic information
                      class_sep = ";", # The character used to separate taxa in the classification
                      class_regex = "^(.+)__(.+)$", # Regex identifying where the data for each taxon is
                      class_key = c(tax_rank = "info", # A key describing each regex capture group
                                    tax_name = "taxon_name"))

# if error occurs with obj not showing OTU column or not reading "col", restart R
# can omit class_regex and class_key if using plain taxon names
             
# classify otus in terms of proporotion making up the taxon
obj.algae.total$data$tax_data <- calc_obs_props(obj.algae.total, "tax_data")

# find the abundance of OTUs for each taxon
# for cols do not include the header, therefore use 2:32
obj.algae.total$data$tax_abund <- calc_taxon_abund(obj.algae.total, "tax_data", cols = 2:31)

# calculate the number of samples that have reads for each taxon
obj.algae.total$data$tax_occ <- calc_n_samples(obj.algae.total, "tax_abund", groups = algae.total.map$Location)

set.seed(1)

# algae generation, color and labeling of individual taxonomic trees
#c an add multiple node_color groups together for complete and subgrouped trees

# Algal Total OTU Diversity
heat_tree(obj.algae.total,
          title =  sprintf("Algal Total OTU Diversity"),
          node_label = taxon_names,
          node_size = n_obs,
          node_color = Rockies+Cascades,
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads",
          layout = "davidson-harel",
          initial_layout = "reingold-tilford")

# Algal Diversity Rockies (individual plot vs total plot)
heat_tree(obj.algae.total,
          title =  sprintf("Algal OTU Diversity: Rockies"),
          node_label = taxon_names,
          node_size = n_obs,
          node_color = Rockies,
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads",
          layout = "davidson-harel",
          initial_layout = "reingold-tilford")

# Algal Diversity Cascades (individual plot vs total plot)
heat_tree(obj.algae.total,
          title =  sprintf("Algal OTU Diversity: Cascades"),
          node_label = taxon_names,
          node_size = n_obs,
          node_color = Cascades,
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads",
          layout = "davidson-harel",
          initial_layout = "reingold-tilford")

######################################################################################################################################################################################
# Algae Heat Tree Matrix

# create a differential heat tree comparing different groupings
# in this case we will look at the differences in a, b, and c sites 
# between the Cascades and Rockies

# add a "taxon_id" column to the first column of our existing data frame (algae.total.map) in order
# for the compare_groups() to properly recognize our data map
# we can use our tibble package to accomplish this

algae.double <- sapply(algae.total[,c(3:32)], as.integer)
algae.names <- algae.total[,c(1,2)]
algae.total.numeric <- cbind(algae.names,algae.double)
algae.total.tibble <- as.tibble(algae.total.numeric)

# SAMPLE IDS MUST BE REPRESENTED BY NUMERIC CHARACTERS UNITS OR ELSE IT WILL NOT BE RECOGNIZED!!!!!!
algae.total.map$Sample_ID <- sapply(algae.total.map$Sample_ID, as.character)

obj.algae.total.heat.matrix <- parse_tax_data(algae.total.tibble,
                                  class_cols = "Lineage", # the column that contains taxonomic information
                                  class_sep = ";", # The character used to separate taxa in the classification
                                  class_regex = "^(.+)__(.+)$", # Regex identifying where the data for each taxon is
                                  class_key = c(tax_rank = "info", # A key describing each regex capture group
                                                tax_name = "taxon_name"))

# Convert counts to proportions
obj.algae.total.heat.matrix$data$otu_table <- calc_obs_props(obj.algae.total.heat.matrix, data = "tax_data", cols = algae.total.map$Sample_ID)

# Get per-taxon counts
obj.algae.total.heat.matrix$data$tax_table <- calc_taxon_abund(obj.algae.total.heat.matrix, data = "otu_table", cols = algae.total.map$Sample_ID)

# Calculate difference between groups
obj.algae.total.heat.matrix$data$diff_table <- compare_groups(obj.algae.total.heat.matrix, data = 'tax_table',
                                      cols = algae.total.map$Sample_ID, # What columns of sample data to use
                                      groups = algae.total.map$Snow) # What category each sample is assigned to
set.seed(1)

print(obj.algae.total.heat.matrix$data$diff_table)

# log2 median ratio
heat_tree_matrix(obj.algae.total.heat.matrix,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_range = diverging_palette(),
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 node_size_axis_label = "Number of OTUs",
                 node_color_axis_label = "Log2 ratio median proportions",
                 layout = "davidson-harel", # The primary layout algorithm
                 initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                 output_file = "algae_log2_differential_heat_matrix.jpg") # Saves the plot as a pdf file

heat_tree_matrix(obj.algae.total.heat.matrix,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_range = diverging_palette(),
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 node_size_axis_label = "Number of OTUs",
                 node_color_axis_label = "Log2 ratio median proportions",
                 layout = "davidson-harel", # The primary layout algorithm
                 initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                 output_file = "algae_log2_differential_heat_matrix.jpg") # Saves the plot as a pdf file


######################################################################################################################################################################################
# Individual Bacteria taxonomic heat tree

# load a map file (modified PERMANOVA) and the OTU/lineage_string/sample file (Modified .shared from Mothur)
bac.total <- read_excel("~/Desktop/Students/Avery/SnowCommSeqProcessing/Graphics/Metacoder/bac.FINAL.shared.METACODER.xlsx", sheet = "bac All Sites (No H2O)")
bac.total.map <- read.csv("~/Desktop/Students/Avery/SnowCommSeqProcessing/Graphics/Metacoder/bac_Snow_Community_PERMANOVA_MAP.csv")

obj.bac.total <- parse_tax_data(bac.total,
                                  class_cols = "Lineage", # the column that contains taxonomic information
                                  class_sep = ";") # The character used to separate taxa in the classification
                                  #class_regex = "^(.+)__(.+)$", # Regex identifying where the data for each taxon is
                                  #class_key = c(tax_rank = "info", # A key describing each regex capture group
                                                tax_name = "taxon_name"))

# if error occurs with obj not showing OTU column or not reading "col", restart R
# can omit class_regex and class_key if using plain taxon names

# classify otus in terms of proporotion making up the taxon
obj.bac.total$data$tax_data <- calc_obs_props(obj.bac.total, "tax_data")

# find the abundance of OTUs for each taxon
# for cols do not include the header, therefore use 2:32
obj.bac.total$data$tax_abund <- calc_taxon_abund(obj.bac.total, "tax_data", cols = 2:31)

# calculate the number of samples that have reads for each taxon
obj.bac.total$data$tax_occ <- calc_n_samples(obj.bac.total, "tax_abund", groups = bac.total.map$Location)

set.seed(1)

# bac generation, color and labeling of individual taxonomic trees
# can add multiple node_color groups together for complete and subgrouped trees

# alternative coloring set:yellow,magenta,blue,grey  c("#E5E5E5","#004cff","#dc5ef2","#f4be41")
# coloring set: brown, green, teal, grey c("#E5E5E5","#40BEA3","#499FB2","#895875")

# Bacterial Total OTU Diversity
heat_tree(obj.bac.total,
          title =  sprintf("Bacterial Total OTU Diversity"),
          node_color_range = c("#E5E5E5","#004cff","#dc5ef2","#f4be41"), 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = Rockies+Cascades,
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads",
          layout = "davidson-harel",
          initial_layout = "reingold-tilford")

# Bacterial Diversity Rockies (individual plot vs total plot)
heat_tree(obj.bac.total,
          title =  sprintf("Bacterial OTU Diversity: Rockies"),
          node_color_range = c("#E5E5E5","#004cff","#dc5ef2","#f4be41"),
          node_label = taxon_names,
          node_size = n_obs,
          node_color = Rockies,
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads",
          layout = "davidson-harel",
          initial_layout = "reingold-tilford")

# Bacterial Diversity Cascades (individual plot vs total plot)
heat_tree(obj.bac.total,
          title =  sprintf("Bacterial OTU Diversity: Cascades"),
          node_color_range = c("#E5E5E5","#004cff","#dc5ef2","#f4be41"),
          node_label = taxon_names,
          node_size = n_obs,
          node_color = Cascades,
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads",
          layout = "davidson-harel",
          initial_layout = "reingold-tilford")

######################################################################################################################################################################################

# Bacteria Heat Tree Matrix

# create a differential heat tree comparing different groupings
# in this case we will look at the differences in a, b, and c sites 
# between the Cascades and Rockies

# add a "taxon_id" column to the first column of our existing data frame (bac.total.map) in order
# for the compare_groups() to properly recognize our data map
# we can use our tibble package to accomplish this

bac.double <- sapply(bac.total[,c(3:32)], as.integer)
bac.names <- bac.total[,c(1,2)]
bac.total.numeric <- cbind(bac.names,bac.double)
bac.total.tibble <- as.tibble(bac.total.numeric)

# SAMPLE IDS MUST BE REPRESENTED BY NUMERIC CHARACTERS UNITS OR ELSE IT WILL NOT BE RECOGNIZED!!!!!!

bac.total.map$Sample_ID <- sapply(bac.total.map$Sample_ID, as.character)

obj.bac.total.heat.matrix <- parse_tax_data(bac.total.tibble,
                                            class_cols = "Lineage", # the column that contains taxonomic information
                                            class_sep = ";") # The character used to separate taxa in the classification
                                            class_regex = "^(.+)__(.+)$", # Regex identifying where the data for each taxon is
                                            class_key = c(tax_rank = "info", # A key describing each regex capture group
                                                          tax_name = "taxon_name"))

# Convert counts to proportions
obj.bac.total.heat.matrix$data$otu_table <- calc_obs_props(obj.bac.total.heat.matrix, data = "tax_data", cols = bac.total.map$Sample_ID)

# Get per-taxon counts
obj.bac.total.heat.matrix$data$tax_table <- calc_taxon_abund(obj.bac.total.heat.matrix, data = "otu_table", cols = bac.total.map$Sample_ID)

# Calculate difference between groups
obj.bac.total.heat.matrix$data$diff_table <- compare_groups(obj.bac.total.heat.matrix, data = 'tax_table',
                                                            cols = bac.total.map$Sample_ID, # What columns of sample data to use
                                                            groups = bac.total.map$Snow) # What category each sample is assigned to

print(obj.bac.total.heat.matrix$data$diff_table)

heat_tree_matrix(obj.bac.total.heat.matrix,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_range = diverging_palette(),
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 node_size_axis_label = "Number of OTUs",
                 node_color_axis_label = "Log2 ratio median proportions",
                 layout = "davidson-harel", # The primary layout algorithm
                 initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                 output_file = "bac_log2_differential_heat_matrix.jpg") # Saves the plot as a pdf file

######################################################################################################################################################################################

# Individual Fungi taxonomic heat tree

# load a map file (modified PERMANOVA) and the OTU/lineage_string/sample file (Modified .shared from Mothur)
fungi.total <- read_excel("/Users/BrownLab/Desktop/Students/Avery/Snow_Community_Project/Stats/Metacoder/Fungal/fungi.FINAL.shared.METACODERs.xlsx", sheet = "Fungi Total (family level)")
fungi.total.map <- read_excel("/Users/BrownLab/Desktop/Students/Avery/Snow_Community_Project/Stats/Metacoder/Fungal/Fungi_Snow_Community_PERMANOVA_MAP.xlsx", sheet = "Numeric")

obj.fungi.total <- parse_tax_data(fungi.total,
                                class_cols = "Lineage", # the column that contains taxonomic information
                                class_sep = ";") # The character used to separate taxa in the classification
                                #class_regex = "^(.+)__(.+)$", # Regex identifying where the data for each taxon is
                                #class_key = c(tax_rank = "info", # A key describing each regex capture group
                                #tax_name = "taxon_name"))
                                
# if error occurs with obj not showing OTU column or not reading "col", restart R
# can omit class_regex and class_key if using plain taxon names

# classify otus in terms of proporotion making up the taxon
obj.fungi.total$data$tax_data <- calc_obs_props(obj.fungi.total, "tax_data")

# find the abundance of OTUs for each taxon
# for cols do not include the header, therefore use 2:32
obj.fungi.total$data$tax_abund <- calc_taxon_abund(obj.fungi.total, "tax_data", cols = 2:31)

# calculate the number of samples that have reads for each taxon
obj.fungi.total$data$tax_occ <- calc_n_samples(obj.fungi.total, "tax_abund", groups = fungi.total.map$Location)

set.seed(1)

# fungi generation, color and labeling of individual taxonomic trees
# can add multiple node_color groups together for complete and subgrouped trees

# Fungal Total OTU Diversity
heat_tree(obj.fungi.total,
          title =  sprintf("Fungal Total OTU Diversity (Family)"),
          node_color_range = c("#E5E5E5","#58E0B4","#9EA84A","#AE6942"), #f4be41"
          node_label = taxon_names,
          node_size = n_obs,
          node_color = Rockies+Cascades,
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads",
          layout = "davidson-harel",
          initial_layout = "reingold-tilford")

# Fungal Diversity Rockies (individual plot vs total plot)
heat_tree(obj.fungi.total,
          title =  sprintf("Fungal OTU Diversity: Rockies"),
          node_color_range = c("#E5E5E5","#58E0B4","#9EA84A","#AE6942"),
          node_label = taxon_names,
          node_size = n_obs,
          node_color = Rockies,
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads",
          layout = "davidson-harel",
          initial_layout = "reingold-tilford")

# Fungal Diversity Cascades (individual plot vs total plot)
heat_tree(obj.fungi.total,
          title =  sprintf("Fungal OTU Diversity: Cascades"),
          node_color_range = c("#E5E5E5","#58E0B4","#9EA84A","#AE6942"),
          node_label = taxon_names,
          node_size = n_obs,
          node_color = Cascades,
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads",
          layout = "davidson-harel",
          initial_layout = "reingold-tilford")

######################################################################################################################################################################################
######################################################################################################################################################################################

# Fungi Heat Tree Matrix

# create a differential heat tree comparing different groupings
# in this case we will look at the differences in a, b, and c sites 
# between the Cascades and Rockies

# add a "taxon_id" column to the first column of our existing data frame (fungi.total.map) in order
# for the compare_groups() to properly recognize our data map
# we can use our tibble package to accomplish this

fungi.double <- sapply(fungi.total[,c(3:32)], as.integer)
fungi.names <- fungi.total[,c(1,2)]
fungi.total.numeric <- cbind(fungi.names,fungi.double)
fungi.total.tibble <- as.tibble(fungi.total.numeric)

# SAMPLE IDS MUST BE REPRESENTED BY NUMERIC CHARACTERS UNITS OR ELSE IT WILL NOT BE RECOGNIZED!!!!!!

fungi.total.map$Sample_ID <- sapply(fungi.total.map$Sample_ID, as.character)

obj.fungi.total.heat.matrix <- parse_tax_data(fungi.total.tibble,
                                              class_cols = "Lineage", # the column that contains taxonomic information
                                              class_sep = ";") # The character used to separate taxa in the classification
                                              #class_regex = "^(.+)__(.+)$") # Regex identifying where the data for each taxon is
                                              #class_key = c(tax_rank = "info", # A key describing each regex capture group
                                                            #tax_name = "taxon_name"))

# Convert counts to proportions
obj.fungi.total.heat.matrix$data$otu_table <- calc_obs_props(obj.fungi.total.heat.matrix, data = "tax_data", cols = fungi.total.map$Sample_ID)

# Get per-taxon counts
obj.fungi.total.heat.matrix$data$tax_table <- calc_taxon_abund(obj.fungi.total.heat.matrix, data = "otu_table", cols = fungi.total.map$Sample_ID)

# Calculate difference between groups
obj.fungi.total.heat.matrix$data$diff_table <- compare_groups(obj.fungi.total.heat.matrix, data = 'tax_table',
                                                              cols = fungi.total.map$Sample_ID, # What columns of sample data to use
                                                              groups = fungi.total.map$Snow) # What category each sample is assigned to

set.seed(100000000)

print(obj.fungi.total.heat.matrix$data$diff_table)

heat_tree_matrix(obj.fungi.total.heat.matrix,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_range = diverging_palette(),
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 node_size_axis_label = "Number of OTUs",
                 node_color_axis_label = "Log2 ratio median proportions",
                 layout = "davidson-harel", # The primary layout algorithm
                 initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                 output_file = "fungi_log2_differential_heat_matrix.jpg") # Saves the plot as a pdf file


######################################################################################################################################################################################
######################################################################################################################################################################################

# Individual Community taxonomic heat tree

setwd("/Users/BrownLab/Desktop/Students/Avery/Snow_Community_Project/Stats/Metacoder/Community")
# load a map file (modified PERMANOVA) and the OTU/lineage_string/sample file (Modified .shared from Mothur)
community.total <- read_excel("/Users/BrownLab/Desktop/Students/Avery/Snow_Community_Project/Stats/Metacoder/Fungal/community.FINAL.shared.METACODERs.xlsx")
community.total.map <- read_excel("/Users/BrownLab/Desktop/Students/Avery/Snow_Community_Project/Stats/Metacoder/Fungal/fungi_Snow_Community_PERMANOVA_MAP.xlsx", sheet = "Numeric")

obj.community.total <- parse_tax_data(community.total,
                                  class_cols = "Lineage", # the column that contains taxonomic information
                                  class_sep = ";") # The character used to separate taxa in the classification
# class_regex = "^(.+)__(.+)$", # Regex identifying where the data for each taxon is
# class_key = c(tax_rank = "info", # A key describing each regex capture group
# tax_name = "taxon_name"))

# if error occurs with obj not showing OTU column or not reading "col", restart R
# can omit class_regex and class_key if using plain taxon names

# classify otus in terms of proporotion making up the taxon
obj.community.total$data$tax_data <- calc_obs_props(obj.community.total, "tax_data")

# find the abundance of OTUs for each taxon
# for cols do not include the header, therefore use 2:32
obj.community.total$data$tax_abund <- calc_taxon_abund(obj.community.total, "tax_data", cols = 2:31)

# calculate the number of samples that have reads for each taxon
obj.community.total$data$tax_occ <- calc_n_samples(obj.community.total, "tax_abund", groups = community.total.map$Location)

set.seed(1)

# community generation, color and labeling of individual taxonomic trees
# can add multiple node_color groups together for complete and subgrouped trees

# Community Total OTU Diversity
heat_tree(obj.community.total,
          title =  sprintf("Fungal Total OTU Diversity (Family)"),
          node_color_range = c("#E5E5E5","#58E0B4","#9EA84A","#AE6942"), #f4be41"
          node_label = taxon_names,
          node_size = n_obs,
          node_color = Rockies+Cascades,
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads",
          layout = "davidson-harel",
          initial_layout = "reingold-tilford")

# Community Diversity Rockies (individual plot vs total plot)
heat_tree(obj.community.total,
          title =  sprintf("Fungal OTU Diversity: Rockies"),
          node_color_range = c("#E5E5E5","#58E0B4","#9EA84A","#AE6942"),
          node_label = taxon_names,
          node_size = n_obs,
          node_color = Rockies,
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads",
          layout = "davidson-harel",
          initial_layout = "reingold-tilford")

# Communityl Diversity Cascades (individual plot vs total plot)
heat_tree(obj.community.total,
          title =  sprintf("Fungal OTU Diversity: Cascades"),
          node_color_range = c("#E5E5E5","#58E0B4","#9EA84A","#AE6942"),
          node_label = taxon_names,
          node_size = n_obs,
          node_color = Cascades,
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads",
          layout = "davidson-harel",
          initial_layout = "reingold-tilford")

######################################################################################################################################################################################

# Community Heat Tree Matrix

# create a differential heat tree comparing different groupings
# in this case we will look at the differences in a, b, and c sites 
# between the Cascades and Rockies

# add a "taxon_id" column to the first column of our existing data frame (community.total.map) in order
# for the compare_groups() to properly recognize our data map
# we can use our tibble package to accomplish this

community.double <- sapply(community.total[,c(4:33)], as.integer)
community.names <- community.total[,c(1,2)]
community.total.numeric <- cbind(community.names,community.double)
community.total.tibble <- as.tibble(community.total.numeric)

# SAMPLE IDS MUST BE REPRESENTED BY NUMERIC CHARACTERS UNITS OR ELSE IT WILL NOT BE RECOGNIZED!!!!!!

community.total.map$Sample_ID <- sapply(community.total.map$Sample_ID, as.character)

obj.community.total.heat.matrix <- parse_tax_data(community.total.tibble,
                                              class_cols = "Lineage", # the column that contains taxonomic information
                                              class_sep = ";") # The character used to separate taxa in the classification
# class_regex = "^(.+)__(.+)$") # Regex identifying where the data for each taxon is
# class_key = c(tax_rank = "info", # A key describing each regex capture group
# tax_name = "taxon_name"))

# Convert counts to proportions
obj.community.total.heat.matrix$data$otu_table <- calc_obs_props(obj.community.total.heat.matrix, data = "tax_data", cols = community.total.map$Sample_ID)

# Get per-taxon counts
obj.community.total.heat.matrix$data$tax_table <- calc_taxon_abund(obj.community.total.heat.matrix, data = "otu_table", cols = community.total.map$Sample_ID)

# Calculate difference between groups
obj.community.total.heat.matrix$data$diff_table <- compare_groups(obj.community.total.heat.matrix, data = 'tax_table',
                                                              cols = community.total.map$Sample_ID, # What columns of sample data to use
                                                              groups = community.total.map$Snow) # What category each sample is assigned to

# save file of data matrix including taxanomic data (names) and calculated statistical values
obj.com.difftable <- as.data.frame(obj.community.total.heat.matrix$data$diff_table)
write.csv(obj.com.difftable, file = "com.abc.pairwise.diff_table.csv")

obj.com.taxatable <- obj.community.total.heat.matrix$data$tax_table
write.csv(obj.com.taxatable, file = "com.abc.pairwise.taxatable.csv")
View(obj.community.total.heat.matrix$)

taxonnames <- taxon_names(obj.community.total.heat.matrix)
write.csv(taxonnames, "com.abc.pairwise.names.csv")

# generate tree matrix
set.seed(100000000)
heat_tree_matrix(obj.community.total.heat.matrix,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_range = diverging_palette(),
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 node_size_axis_label = "Number of OTUs",
                 node_color_axis_label = "Log2 ratio median proportions",
                 layout = "davidson-harel", # The primary layout algorithm
                 initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                 output_file = "community_log2_differential_heat_matrix.jpg",
                 node_size_axis_label = "NA") # Saves the plot as a pdf file
make_plot_legend(100, 100, tick_size = 0, title = NA)

######################################################################################################################################################################################

# Cascade Community Heat Tree Matrix

# create a differential heat tree comparing different groupings
# in this case we will look at the differences in a, b, and c sites 
# between the Cascades and Rockies

# add a "taxon_id" column to the first column of our existing data frame (community.total.map) in order
# for the compare_groups() to properly recognize our data map
# we can use our tibble package to accomplish this

community.double.cas <- sapply(community.total[,c(13:30)], as.integer)
community.names.cas <- community.total[,c(1,2)]
community.total.numeric.cas <- cbind(community.names.cas,community.double.cas)
community.total.tibble.cas <- as.tibble(community.total.numeric.cas)
community.total.map.cas <- community.total.map[10:27,]

# SAMPLE IDS MUST BE REPRESENTED BY NUMERIC CHARACTERS UNITS OR ELSE IT WILL NOT BE RECOGNIZED!!!!!!

community.total.map.cas$Sample_ID <- sapply(community.total.map.cas$Sample_ID, as.character)



obj.community.total.heat.matrix.cas <- parse_tax_data(community.total.tibble.cas,
                                                  class_cols = "Lineage", # the column that contains taxonomic information
                                                  class_sep = ";") # The character used to separate taxa in the classification
# class_regex = "^(.+)__(.+)$") # Regex identifying where the data for each taxon is
# class_key = c(tax_rank = "info", # A key describing each regex capture group
# tax_name = "taxon_name"))

# Convert counts to proportions
obj.community.total.heat.matrix.cas$data$otu_table <- calc_obs_props(obj.community.total.heat.matrix.cas, data = "tax_data", cols = community.total.map.cas$Sample_ID)

# Get per-taxon counts
obj.community.total.heat.matrix.cas$data$tax_table <- calc_taxon_abund(obj.community.total.heat.matrix.cas, data = "otu_table", cols = community.total.map.cas$Sample_ID)

# Calculate difference between groups
obj.community.total.heat.matrix.cas$data$diff_table <- compare_groups(obj.community.total.heat.matrix.cas, data = 'tax_table',
                                                                  cols = community.total.map.cas$Sample_ID, # What columns of sample data to use
                                                                  groups = community.total.map.cas$Snow) # What category each sample is assigned to

# save file of data matrix including taxanomic data (names) and calculated statistical values
obj.com.difftable.cas <- as.data.frame(obj.community.total.heat.matrix.cas$data$diff_table)
write.csv(obj.com.difftable.cas, file = "com.cas.abc.pairwise.diff_table.csv")

obj.com.taxatable.cas <- obj.community.total.heat.matrix.cas$data$tax_table
write.csv(obj.com.taxatable.cas, file = "com.cas.abc.pairwise.taxatable.csv")

taxonnames <- taxon_names(obj.community.total.heat.matrix.cas)
write.csv(taxonnames, "com.cas.abc.pairwise.names.csv")

# generate tree matrix (Cascades)
set.seed(1)
heat_tree_matrix(obj.community.total.heat.matrix.cas,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_range = diverging_palette(),
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 node_color_axis_label = "Log2 ratio median proportions",
                 layout = "davidson-harel", # The primary layout algorithm
                 initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                 output_file = "community_cas_log2_differential_heat_matrix.pdf",
                 node_size_axis_label = "NA") # Saves the plot as a pdf file
make_plot_legend(100, 100, tick_size = 0, title = NA)

######################################################################################################################################################################################

# Rockies Community Heat Tree Matrix

# create a differential heat tree comparing different groupings
# in this case we will look at the differences in a, b, and c sites 
# between the Rockiess and Rockies

# add a "taxon_id" column to the first column of our existing data frame (community.total.map) in order
# for the compare_groups() to properly recognize our data map
# we can use our tibble package to accomplish this

community.double.roc <- sapply(community.total[,c(4:12, 31:33)], as.integer)
community.names.roc <- community.total[,c(1,2)]
community.total.numeric.roc <- cbind(community.names.roc,community.double.roc)
community.total.tibble.roc <- as.tibble(community.total.numeric.roc)
community.total.map.roc <- community.total.map[c(1:9, 28:30),]

# SAMPLE IDS MUST BE REPRESENTED BY NUMERIC CHARACTERS UNITS OR ELSE IT WILL NOT BE RECOGNIZED!!!!!!

community.total.map.roc$Sample_ID <- sapply(community.total.map.roc$Sample_ID, as.character)



obj.community.total.heat.matrix.roc <- parse_tax_data(community.total.tibble.roc,
                                                      class_cols = "Lineage", # the column that contains taxonomic information
                                                      class_sep = ";") # The character used to separate taxa in the classification
# class_regex = "^(.+)__(.+)$") # Regex identifying where the data for each taxon is
# class_key = c(tax_rank = "info", # A key describing each regex capture group
# tax_name = "taxon_name"))

# Convert counts to proportions
obj.community.total.heat.matrix.roc$data$otu_table <- calc_obs_props(obj.community.total.heat.matrix.roc, data = "tax_data", cols = community.total.map.roc$Sample_ID)

# Get per-taxon counts
obj.community.total.heat.matrix.roc$data$tax_table <- calc_taxon_abund(obj.community.total.heat.matrix.roc, data = "otu_table", cols = community.total.map.roc$Sample_ID)

# Calculate difference between groups
obj.community.total.heat.matrix.roc$data$diff_table <- compare_groups(obj.community.total.heat.matrix.roc, data = 'tax_table',
                                                                      cols = community.total.map.roc$Sample_ID, # What columns of sample data to use
                                                                      groups = community.total.map.roc$Snow) # What category each sample is assigned to

# save file of data matrix including taxanomic data (names) and calculated statistical values
obj.com.difftable.roc <- as.data.frame(obj.community.total.heat.matrix.roc$data$diff_table)
write.csv(obj.com.difftable.roc, file = "com.roc.abc.pairwise.diff_table.csv")

obj.com.taxatable.roc <- obj.community.total.heat.matrix.roc$data$tax_table
write.csv(obj.com.taxatable.roc, file = "com.roc.abc.pairwise.taxatable.csv")

taxonnames <- taxon_names(obj.community.total.heat.matrix.roc)
write.csv(taxonnames, "com.roc.abc.pairwise.names.csv")

# generate tree matrix (Rockiess)
set.seed(1)
heat_tree_matrix(obj.community.total.heat.matrix.roc,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_range = diverging_palette(),
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 node_color_axis_label = "Log2 ratio median proportions",
                 layout = "davidson-harel", # The primary layout algorithm
                 initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                 output_file = "community_roc_log2_differential_heat_matrix.pdf",
                 node_size_axis_label = "NA") # Saves the plot as a pdf file
make_plot_legend(100, 100, tick_size = 0, title = NA)