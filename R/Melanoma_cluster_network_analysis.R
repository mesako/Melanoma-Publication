rm(list = ls())

#######################################
######## LOAD AND SET-UP FILES ########
#######################################

library(igraph)
library(scales)
library(Hmisc)
library(reshape2)
library(ggplot2)

# markers with highest variance
clustering.var <- c("Ki67", "Mart1", "HIF1a", "LDH", "AMPKA",
                    "p_ERK1", "PFK", "p_ACAC", "Slug", "p_LKB")

file.name <- "data_with_diff_clusters_added_values_from_allmarkers.rds"
df.with.vals <- readRDS(file = file.name)

non.markers <- c("day", "Cell_Order_Number", "lambda0", "lambda1",
                 "lambda2", "lambda3", colnames(df.with.vals)[27:ncol(df.with.vals)])

AssignMarkerClasses <- function(graph, marker.class.key) {
  new.graph <- graph
  for (x in 1:length(marker.class.key)) {
    V(new.graph)$class[get.vertex.attribute(new.graph, "name") %in% marker.class.key[[x]]] <- names(marker.class.key)[x]
  }
  return(new.graph)
}

ConverttoColor <- function(class.label, color.code) {
  new.color <- color.code[[class.label]]
  return(new.color)
}

ClusterCorrelationNetworks <- function(temp.cluster, file.prefix,
                                       marker.classes, color.code) {
  temp.pearson.corr <- rcorr(as.matrix(temp.cluster), type = "pearson")
  file.name <- paste(file.prefix, "_pearson_corr.csv", sep = "")
  temp.pearson.corr$n <- NULL
  lapply(temp.pearson.corr, function(x) write.table(data.frame(x), file = file.name, append = TRUE, sep = ","))
  
  graph <- graph.adjacency(temp.pearson.corr$r, weighted = TRUE, mode = "undirected", diag = FALSE)
  temp.graph <- graph.adjacency(temp.pearson.corr$P, weighted = TRUE, mode = "undirected", diag = FALSE)
  graph <- set.edge.attribute(graph, "p_value", value = get.edge.attribute(temp.graph, "weight"))
  
  graph.sig.0.05 <- delete.edges(graph, E(graph)[p_value > 0.05])
  E(graph.sig.0.05)$color[E(graph.sig.0.05)$weight > 0] <- "red"
  E(graph.sig.0.05)$color[E(graph.sig.0.05)$weight < 0] <- "blue"
  graph.sig.0.05 <- AssignMarkerClasses(graph.sig.0.05, marker.classes)
  V(graph.sig.0.05)$color <- unlist(lapply(V(graph.sig.0.05)$class, ConverttoColor, color.code = color.code))
  
  plot.name <- paste(file.prefix, "_pearson_correlation_network_p_value_0.05.png", sep = "")
  png(plot.name, width = 5, height = 5, units = "in", res = 600)
  plot(graph.sig.0.05, layout = layout_in_circle, vertex.label.family = "Arial",
       edge.width = abs(E(graph.sig.0.05)$weight * 3), vertex.label.color = "black")
  dev.off()

  graph.sig.0.01 <- delete.edges(graph.sig.0.05, E(graph.sig.0.05)[p_value > 0.01])
  
  plot.name <- paste(file.prefix, "_pearson_correlation_network_p_value_0.01.png", sep = "")
  png(plot.name, width = 5, height = 5, units = "in", res = 600)
  plot(graph.sig.0.01, layout = layout_in_circle, vertex.label.family = "Arial",
       edge.width = abs(E(graph.sig.0.01)$weight * 3), vertex.label.color = "black")
  dev.off()
  
  graph.sig.0.001 <- delete.edges(graph.sig.0.01, E(graph.sig.0.01)[p_value > 0.001])
  
  plot.name <- paste(file.prefix, "_pearson_correlation_network_p_value_0.001.png", sep = "")
  png(plot.name, width = 5, height = 5, units = "in", res = 600)
  plot(graph.sig.0.001, layout = layout_in_circle, vertex.label.family = "Arial",
       edge.width = abs(E(graph.sig.0.001)$weight * 3), vertex.label.color = "black")
  dev.off()
  
  graph.sig.0.001.r.0.5 <- delete.edges(graph.sig.0.001, E(graph.sig.0.001)[abs(weight) < 0.5])
  
  plot.name <- paste(file.prefix, "_pearson_correlation_network_p_value_0.001_r_0.5.png", sep = "")
  png(plot.name, width = 5, height = 5, units = "in", res = 600)
  plot(graph.sig.0.001.r.0.5, layout = layout_in_circle, vertex.label.family = "Arial",
       edge.width = abs(E(graph.sig.0.001.r.0.5)$weight * 3), vertex.label.color = "black")
  dev.off()
  
  file.name <- paste(file.prefix, "_pearson_correlation_network_p_value_0.05", sep = "")
  file.name <- paste(file.name, ".graphml", sep = "")
  V(graph.sig.0.05)$label <- V(graph.sig.0.05)$name
  E(graph.sig.0.05)$weight <- abs(E(graph.sig.0.05)$weight)
  write.graph(graph.sig.0.05, file.name, format = "graphml")
  
  file.name <- paste(file.prefix, "_pearson_correlation_network_p_value_0.01", sep = "")
  file.name <- paste(file.name, ".graphml", sep = "")
  V(graph.sig.0.01)$label <- V(graph.sig.0.01)$name
  E(graph.sig.0.01)$weight <- abs(E(graph.sig.0.01)$weight)
  write.graph(graph.sig.0.01, file.name, format = "graphml")
  
  file.name <- paste(file.prefix, "_pearson_correlation_network_p_value_0.001", sep = "")
  file.name <- paste(file.name, ".graphml", sep = "")
  V(graph.sig.0.001)$label <- V(graph.sig.0.001)$name
  E(graph.sig.0.001)$weight <- abs(E(graph.sig.0.001)$weight)
  write.graph(graph.sig.0.001, file.name, format = "graphml")
  
  file.name <- paste(file.prefix, "_pearson_correlation_network_p_value_0.001_r_0.5", sep = "")
  file.name <- paste(file.name, ".graphml", sep = "")
  V(graph.sig.0.001.r.0.5)$label <- V(graph.sig.0.001.r.0.5)$name
  E(graph.sig.0.001.r.0.5)$weight <- abs(E(graph.sig.0.001.r.0.5)$weight)
  write.graph(graph.sig.0.001.r.0.5, file.name, format = "graphml")
  
  temp.spearman.corr <- rcorr(as.matrix(temp.cluster), type = "spearman")
  file.name <- paste(file.prefix, "_spearman_corr.csv", sep = "")
  temp.spearman.corr$n <- NULL
  lapply(temp.pearson.corr, function(x) write.table(data.frame(x), file = file.name, append = TRUE, sep = ","))
  
  temp2 <- temp.spearman.corr$r
  temp2[which(temp2 == 0)] <- 0.0000001
  graph.spearman <- graph.adjacency(temp2, weighted = TRUE, mode = "undirected", diag = FALSE)
  temp.graph.spearman <- graph.adjacency(temp.spearman.corr$P, weighted = TRUE, mode = "undirected", diag = FALSE)
  graph.spearman <- set.edge.attribute(graph.spearman, "p_value", value = get.edge.attribute(temp.graph.spearman, "weight"))
  
  graph.spearman.sig.0.05 <- delete.edges(graph.spearman, E(graph.spearman)[p_value > 0.05])
  E(graph.spearman.sig.0.05)$color[E(graph.spearman.sig.0.05)$weight > 0] <- "red"
  E(graph.spearman.sig.0.05)$color[E(graph.spearman.sig.0.05)$weight < 0] <- "blue"
  V(graph.spearman.sig.0.05)$color <- V(graph.sig.0.05)$color
  
  plot.name <- paste(file.prefix, "_spearman_correlation_network_p_value_0.05.png", sep = "")
  png(plot.name, width = 5, height = 5, units = "in", res = 600)
  plot(graph.spearman.sig.0.05, layout = layout_in_circle, vertex.label.family = "Arial",
       edge.width = abs(E(graph.spearman.sig.0.05)$weight * 3), vertex.label.color = "black")
  dev.off()
  
  graph.spearman.sig.0.01 <- delete.edges(graph.spearman.sig.0.05, E(graph.spearman.sig.0.05)[p_value > 0.01])
  
  plot.name <- paste(file.prefix, "_spearman_correlation_network_p_value_0.01.png", sep = "")
  png(plot.name, width = 5, height = 5, units = "in", res = 600)
  plot(graph.spearman.sig.0.01, layout = layout_in_circle, vertex.label.family = "Arial",
       edge.width = abs(E(graph.spearman.sig.0.01)$weight * 3), vertex.label.color = "black")
  dev.off()
  
  graph.spearman.sig.0.001 <- delete.edges(graph.spearman.sig.0.01, E(graph.spearman.sig.0.01)[p_value > 0.001])
  
  plot.name <- paste(file.prefix, "_spearman_correlation_network_p_value_0.001.png", sep = "")
  png(plot.name, width = 5, height = 5, units = "in", res = 600)
  plot(graph.spearman.sig.0.001, layout = layout_in_circle, vertex.label.family = "Arial",
       edge.width = abs(E(graph.spearman.sig.0.001)$weight * 3), vertex.label.color = "black")
  dev.off()
  
  file.name <- paste(file.prefix, "_spearman_correlation_network_p_value_0.05", sep = "")
  file.name <- paste(file.name, ".graphml", sep = "")
  V(graph.spearman.sig.0.05)$label <- V(graph.spearman.sig.0.05)$name
  E(graph.spearman.sig.0.05)$weight <- abs(E(graph.spearman.sig.0.05)$weight)
  write.graph(graph.spearman.sig.0.05, file.name, format = "graphml")
  
  file.name <- paste(file.prefix, "_spearman_correlation_network_p_value_0.01", sep = "")
  file.name <- paste(file.name, ".graphml", sep = "")
  V(graph.spearman.sig.0.01)$label <- V(graph.spearman.sig.0.01)$name
  E(graph.spearman.sig.0.01)$weight <- abs(E(graph.spearman.sig.0.01)$weight)
  write.graph(graph.spearman.sig.0.01, file.name, format = "graphml")
  
  file.name <- paste(file.prefix, "_spearman_correlation_network_p_value_0.001", sep = "")
  file.name <- paste(file.name, ".graphml", sep = "")
  V(graph.spearman.sig.0.001)$label <- V(graph.spearman.sig.0.001)$name
  E(graph.spearman.sig.0.001)$weight <- abs(E(graph.spearman.sig.0.001)$weight)
  write.graph(graph.spearman.sig.0.001, file.name, format = "graphml")
  
  all.graphs <- list(graph.sig.0.05 = graph.sig.0.05,
                     graph.sig.0.01 = graph.sig.0.01,
                     graph.sig.0.001 = graph.sig.0.001,
                     graph.sig.0.001.r.0.5 = graph.sig.0.001.r.0.5,
                     graph.spearman.sig.0.05 = graph.spearman.sig.0.05,
                     graph.spearman.sig.0.01 = graph.spearman.sig.0.01,
                     graph.spearman.sig.0.001 = graph.spearman.sig.0.001)
  names(all.graphs) <- c("Pearson_p0.05", "Pearson_p0.01", "Pearson_p0.001",
                         "Pearson_p0.001_r0.5", "Spearman_p0.05", "Spearman_p0.01",
                         "Spearman_p0.001")
  return(all.graphs)
}

metabolite.class <- c("LDH", "PFK", "Glucose", "PKM2")
proliferation.class <- c("Ki67")
resistant.class <- c("AXL", "NCadherin", "NGFR", "TNFR", "Slug")
metabolic.reg.class <- c("pLKB", "pACAC", "PDK1", "HIF1a", "AMPKA")
sensitive.class <- c("Mart1", "MITF")
signaling.class <- c("pERK1", "pNFkB", "pSrc")

marker.classes <- list("metabolite" = metabolite.class, "proliferation" = proliferation.class,
                       "resistant" = resistant.class, "metabolic_regulator" = metabolic.reg.class,
                       "sensitive" = sensitive.class, "signaling" = signaling.class)

color.code <- list("metabolic_regulator" = "red", "resistant" = "orange",
                   "proliferation" = "cyan", "metabolite" = "green",
                   "sensitive" = "yellow", "signaling" = "pink")

non.markers <- c("day", "Cell_Order_Number", "lambda0", "lambda1",
                 "lambda2", "lambda3", colnames(df.with.vals)[27:ncol(df.with.vals)])

# Export the correlation matrix for Rclusterpp clusters
# Spearman correlation network of major path in Rclusterpp
# (cluster 7) with p-value cut-off 0.001
these.graphs <- list()
i <- 1
for (x in unique(df.with.vals$Rclusterpp_Cluster)) {
  file.prefix <- paste("Rclusterpp_cluster_", x, sep = "")
  temp.cluster <- df.with.vals[which(df.with.vals$Rclusterpp_Cluster == x), ]
  temp.cluster <- temp.cluster[, setdiff(colnames(temp.cluster), non.markers)]
  colnames(temp.cluster) <- c("Ki67", "Glucose", "MITF", "Mart1", "PFK",
                              "pACAC", "pLKB", "PDK1", "PKM2", "LDH",
                              "NGFR", "HIF1a", "TNFR", "NCadherin", "AXL",
                              "pERK1", "pNFkB", "pSrc", "Slug", "AMPKA")
  these.graphs[[i]] <- ClusterCorrelationNetworks(temp.cluster, file.prefix, marker.classes, color.code)
  i <- i + 1
}

IdentifyNodeHubs <- function(all.graphs, file.prefix) {
  # Node degree
  all.degrees <- c()
  for (i in 1:length(all.graphs)) {
    this.degree <- degree(all.graphs[[i]], mode = "all")
    all.degrees <- cbind(all.degrees, this.degree)
  }
  colnames(all.degrees) <- names(all.graphs)
  
  all.degrees.m <- melt(all.degrees)
  plot.name <- paste(file.prefix, "_correlation_network_node_degree.png", sep = "")
  node.degree.barplot <- ggplot(data = all.degrees.m, aes(x = Var1, y = value, fill = Var2)) +
    geom_bar(stat = "identity", position = "dodge")
  print(node.degree.barplot)
  ggsave(filename = plot.name, width = 12, height = 4, units = c("in"), dpi = 400)
  
  all.degrees.rescaled <- apply(all.degrees, 2, rescale)
  
  file.name <- paste(file.prefix, "_all_marker_degrees_rescaled.csv", sep = "")
  write.csv(all.degrees.rescaled, file = file.name)
  
  all.degrees.rescaled.m <- melt(all.degrees.rescaled)
  plot.name <- paste(file.prefix, "_correlation_network_node_degree_rescaled.png", sep = "")
  node.degree.rescaled.barplot <- ggplot(data = all.degrees.rescaled.m, aes(x = Var1, y = value, fill = Var2)) +
    geom_bar(stat = "identity", position = "dodge")
  print(node.degree.rescaled.barplot)
  ggsave(filename = plot.name, width = 12, height = 4, units = c("in"), dpi = 400)
  
  # Kleinberg's hub centrality scores
  all.hub.scores <- c()
  for (i in 1:length(all.graphs)) {
    this.hub.score <- hub_score(all.graphs[[i]])$vector
    all.hub.scores <- cbind(all.hub.scores, this.hub.score)
  }
  colnames(all.hub.scores) <- names(all.graphs)
  
  file.name <- paste(file.prefix, "_all_marker_hub_scores.csv", sep = "")
  write.csv(all.hub.scores, file = file.name)
  
  all.hs.m <- melt(all.hub.scores)
  plot.name <- paste(file.prefix, "_correlation_network_hub_score.png", sep = "")
  hub.scores.barplot <- ggplot(data = all.hs.m, aes(x = Var1, y = value, fill = Var2)) +
    geom_bar(stat = "identity", position = "dodge")
  print(hub.scores.barplot)
  ggsave(filename = plot.name, width = 12, height = 4, units = c("in"), dpi = 400)
}

IdentifyNetworkCommunities <- function(all.graphs, file.prefix) {
  for (i in 1:length(all.graphs)) {
    temp <- all.graphs[[i]]
    if (length(E(temp)) != 0) {
      E(temp)$weight <- abs(E(all.graphs[[i]])$weight)
      # Community detection based on edge betweenness (Girvan-Newman)
      temp.ceb <- cluster_edge_betweenness(temp) 
      this.graph.name <- names(all.graphs)[i]
      print(this.graph.name)
      print(temp)
      plot.name <- paste(file.prefix, "_correlation_network_", this.graph.name,
                         "_GN_communities.png", sep = "")
      png(plot.name, width = 4, height = 4, units = "in", res = 600)
      dendPlot(temp.ceb, mode = "hclust")
      dev.off()
      
      # Community detection based on greedy optimization of modularity
      temp.cfg <- cluster_fast_greedy(temp)
      plot.name <- paste(file.prefix, "_correlation_network_", this.graph.name,
                         "_FG_communities.png", sep = "")
      png(plot.name, width = 4, height = 4, units = "in", res = 600)
      plot(temp.cfg, temp)
      dev.off() 
    }
  }
}

# Do community analysis on this network
for (x in 1:length(these.graphs)) {
  print(x)
  file.prefix <- paste("Rclusterpp_cluster_", x, sep = "")
  this.graph <- these.graphs[[x]]
  IdentifyNodeHubs(this.graph, file.prefix)
  IdentifyNetworkCommunities(this.graph, file.prefix)
}

# Export the values and error bar in the line plots from
# path_comparison_marker_change_by_day
cutoff <- 0.3
pos.lambda2 <- df.with.vals[which(df.with.vals$lambda2 >= cutoff), ]
neg.lambda2 <- df.with.vals[which(df.with.vals$lambda2 < -cutoff), ]
temp.pos.lambda2 <- cbind(pos.lambda2, path = rep("minor", times = nrow(pos.lambda2)))
temp.neg.lambda2 <- cbind(neg.lambda2, path = rep("major", times = nrow(neg.lambda2)))
combined <- rbind(temp.pos.lambda2, temp.neg.lambda2)

skip.channels <- c("Cell_Order_Number", "lambda0",
                   "lambda3", "day", "Seurat_Cluster",
                   "Rclusterpp_Cluster", "flowMeans_Cluster",
                   "SC3_Cluster", "Seurat_SNAI", "Rclusterpp_Ic",
                   "Rclusterpp_SNAI", "Seurat_Ic", "flowMeans_Ic",
                   "flowMeans_SNAI", "SC3_Ic", "SC3_SNAI")

sem <- function(x) {
  y <- sd(x) / sqrt(length(x))
  return(y)
}

all.minor.means <- c()
all.minor.sems <- c()
all.major.means <- c()
all.major.sems <- c()
new.col.names <- c()

minor.data <- combined[which(combined$path == "minor"), ]
major.data <- combined[which(combined$path == "major"), ]
all.days <- unique(combined$day)

for (marker in colnames(pos.lambda2)) {
  if (marker %in% skip.channels) {
    next
  }
  for (x in all.days) {
    all.minor.means <- c(all.minor.means, mean(minor.data[which(minor.data$day == x), marker]))
    all.minor.sems <- c(all.minor.sems, sem(minor.data[which(minor.data$day == x), marker]))
    all.major.means <- c(all.major.means, mean(major.data[which(major.data$day == x), marker]))
    all.major.sems <- c(all.major.sems, sem(major.data[which(major.data$day == x), marker]))
  }
  new.col.names <- c(new.col.names, marker)
}

minor.data <- rbind(all.minor.means[seq(from = 1, by = 4, to = length(all.minor.means))],
                    all.minor.sems[seq(from = 1, by = 4, to = length(all.minor.sems))],
                    all.minor.means[seq(from = 2, by = 4, to = length(all.minor.means))],
                    all.minor.sems[seq(from = 2, by = 4, to = length(all.minor.sems))],
                    all.minor.means[seq(from = 3, by = 4, to = length(all.minor.means))],
                    all.minor.sems[seq(from = 3, by = 4, to = length(all.minor.sems))],
                    all.minor.means[seq(from = 4, by = 4, to = length(all.minor.means))],
                    all.minor.sems[seq(from = 4, by = 4, to = length(all.minor.sems))])
rownames(minor.data) <- as.vector(rbind(paste("Day", all.days, "Mean"),
                                       paste("Day", all.days, "SEM")))
colnames(minor.data) <- new.col.names
write.csv(minor.data, file = "minor_path_mean_SEM_data.csv")

major.data <- rbind(all.major.means[seq(from = 1, by = 4, to = length(all.major.means))],
                    all.major.sems[seq(from = 1, by = 4, to = length(all.major.sems))],
                    all.major.means[seq(from = 2, by = 4, to = length(all.major.means))],
                    all.major.sems[seq(from = 2, by = 4, to = length(all.major.sems))],
                    all.major.means[seq(from = 3, by = 4, to = length(all.major.means))],
                    all.major.sems[seq(from = 3, by = 4, to = length(all.major.sems))],
                    all.major.means[seq(from = 4, by = 4, to = length(all.major.means))],
                    all.major.sems[seq(from = 4, by = 4, to = length(all.major.sems))])
rownames(major.data) <- as.vector(rbind(paste("Day", all.days, "Mean"),
                                        paste("Day", all.days, "SEM")))
colnames(major.data) <- new.col.names
write.csv(major.data, file = "major_path_mean_SEM_data.csv")
