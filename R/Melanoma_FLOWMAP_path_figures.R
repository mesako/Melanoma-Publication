# rm(list = ls())

#######################################
######## LOAD AND SET-UP FILES ########
#######################################

library(igraph)
library(scales)

# markers with highest variance
clustering.var <- c("Ki67", "Mart1", "HIF1a", "LDH", "AMPKA",
                    "p_ERK1", "PFK", "p_ACAC", "Slug", "p_LKB")

# file.name <- "data_with_diff_clusters_added_values_from_allmarkers.rds"
# df.with.vals <- readRDS(file = file.name)

################################
######## WRITE CSV FILE ########
################################

# which cell belong to which cluster (in original excel table)
write.csv(x = df.with.vals, file = "single_cell_with_added_values.csv")

#######################################################
######## MODIFY GRAPHML FILE TO DIM MAJOR PATH ########
#######################################################

graph.file <- "melanoma_latest_graph.graphml"
latest.graph <- read.graph(file = graph.file, format = "graphml")

# Rclusterpp cluster order
# major path: 1 => 7 => 6 => 8 => 11 => 10 => 12
# minor path: 3 => 2 => 9 => 13 => 14
# in the middle: 5 => 4

# set SNAI and Ic values of one or the other to -1000

# major path figure
# set 3, 2, 9, 13, 14, 5, 4 to -1000
remove.set <- c(3, 2, 9, 13, 14, 5, 4)
node.clusters <- get.vertex.attribute(latest.graph, "Rclusterpp_Cluster")
node.remove <- which(node.clusters %in% remove.set)

major.only.graph <- latest.graph
major.only.graph <- set.vertex.attribute(major.only.graph, "Rclusterpp_SNAI",
                                         index = node.remove,
                                         value = rep(-1000, times = length(node.remove)))
major.only.graph <- set.vertex.attribute(major.only.graph, "Rclusterpp_Ic",
                                         index = node.remove,
                                         value = rep(-1000, times = length(node.remove)))
major.graph.file <- "melanoma_latest_graph_major_only.graphml"
write.graph(graph = major.only.graph, file = major.graph.file, format = "graphml")

#######################################################
######## MODIFY GRAPHML FILE TO DIM MINOR PATH ########
#######################################################

# minor path figure
# set 7, 6, 8, 11, 10, 12, 5, 4 to -1000
remove.set <- c(1, 7, 6, 8, 11, 10, 12, 5, 4)
node.clusters <- get.vertex.attribute(latest.graph, "Rclusterpp_Cluster")
node.remove <- which(node.clusters %in% remove.set)

minor.only.graph <- latest.graph
minor.only.graph <- set.vertex.attribute(minor.only.graph, "Rclusterpp_SNAI",
                                         index = node.remove,
                                         value = rep(-1000, times = length(node.remove)))
minor.only.graph <- set.vertex.attribute(minor.only.graph, "Rclusterpp_Ic",
                                         index = node.remove,
                                         value = rep(-1000, times = length(node.remove)))
minor.graph.file <- "melanoma_latest_graph_amended_minor_only.graphml"
write.graph(graph = minor.only.graph, file = minor.graph.file, format = "graphml")


#####################################################
######## EXPORT COLOR PALETTE FOR BOTH PATHS ########
#####################################################

# major path palette
colors <- c("#1500BF", "#C3C7C3", "#D10100")
color.palette <- colorRampPalette(colors)
use.palette <- color.palette(100)

major.SNAI.vals <- sort(unique(get.vertex.attribute(major.only.graph, "Rclusterpp_SNAI")))
major.SNAI.vals <- major.SNAI.vals[-1]

major.Ic.vals <- sort(unique(get.vertex.attribute(major.only.graph, "Rclusterpp_Ic")))
major.Ic.vals <- major.Ic.vals[-1]

x <- rescale(major.SNAI.vals, to = c(0, 100))
v <- seq(0, 100, by = 1)
distrib <- findInterval(x, v, all.inside = TRUE)
major.SNAI.colors <- use.palette[distrib]

x <- rescale(major.Ic.vals, to = c(0, 100))
v <- seq(0, 100, by = 1)
distrib <- findInterval(x, v, all.inside = TRUE)
major.Ic.colors <- use.palette[distrib]

major.export <- cbind.data.frame(major.SNAI.vals, major.SNAI.colors,
                                 major.Ic.vals, major.Ic.colors)
write.csv(major.export, file = "major_path_vals_colors.csv", quote = FALSE)

# minor path palette
minor.SNAI.vals <- sort(unique(get.vertex.attribute(minor.only.graph, "Rclusterpp_SNAI")))
minor.SNAI.vals <- minor.SNAI.vals[-1]

minor.Ic.vals <- sort(unique(get.vertex.attribute(minor.only.graph, "Rclusterpp_Ic")))
minor.Ic.vals <- minor.Ic.vals[-1]

x <- rescale(minor.SNAI.vals, to = c(0, 100))
v <- seq(0, 100, by = 1)
distrib <- findInterval(x, v, all.inside = TRUE)
minor.SNAI.colors <- use.palette[distrib]

x <- rescale(minor.Ic.vals, to = c(0, 100))
v <- seq(0, 100, by = 1)
distrib <- findInterval(x, v, all.inside = TRUE)
minor.Ic.colors <- use.palette[distrib]

minor.export <- cbind.data.frame(minor.SNAI.vals, minor.SNAI.colors,
                                 minor.Ic.vals, minor.Ic.colors)
write.csv(minor.export, file = "amended_minor_path_vals_colors.csv", quote = FALSE)
