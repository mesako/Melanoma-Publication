# rm(list = ls())

##################################################
######## PACKAGE INSTALLATION AND LOADING ########
##################################################

library(FLOWMAPR)
library(readxl)
library(Seurat)

################################################
######## LOAD SINGLE-CELL MELANOMA DATA ########
################################################

file <- "/Users/mesako/Downloads/melanoma-dataset.xlsx"
df <- read_excel(file)
df <- as.data.frame(df)
df.keep <- subset.data.frame(df, select = c("day", "Cell Order Number"))
df.transform <- subset.data.frame(df, select = setdiff(colnames(df), c("day", "Cell Order Number")))
df.transform <- apply(df.transform, 2, log)
final.df <- cbind(df.keep, df.transform)
colnames(final.df) <- gsub(colnames(final.df), pattern = " ", replacement = "_")
colnames(final.df) <- gsub(colnames(final.df), pattern = "-", replacement = "_")

df.by.time <- list()
for (i in 1:length(unique(final.df$day))) {
  this.time <- unique(final.df$day)[i]
  temp.data <- final.df[final.df$day == this.time, ]
  temp.data <- temp.data[, setdiff(colnames(final.df), c("day", "Cell_Order_Number"))]
  df.by.time[[i]] <- temp.data
}
orig.times <- unique(final.df$day)

#################################################################
######## SEURAT CLUSTERING WITHIN AND BETWEEN TIMEPOINTS ########
#################################################################

melanoma.data <- as.data.frame(t(final.df))
melanoma.data <- CreateSeuratObject(raw.data = melanoma.data, min.cells = 3, min.genes = 20, 
                                    project = "Melanoma")
clustering.var <- c("Ki67", "Mart1", "HIF1a", "LDH", "AMPKA",
                    "p_ERK1", "PFK", "p_ACAC", "Slug", "p_LKB")
melanoma.data@scale.data <- melanoma.data@data
melanoma.data@var.genes <- clustering.var
melanoma.data <- RunPCA(object = melanoma.data, do.print = FALSE)

clustering.res1 <- FindClusters(object = melanoma.data, reduction.type = "pca", dims.use = 1:6, 
                                resolution = 1, print.output = 0, save.SNN = TRUE)
seurat.cluster.assign.1 <- clustering.res1@ident

melanoma.data.by.time <- list()
clustering.by.time <- list()
for (i in 1:length(df.by.time)) {
  melanoma.data.by.time[[i]] <- as.data.frame(t(df.by.time[[i]]))
  melanoma.data.by.time[[i]] <- CreateSeuratObject(raw.data = melanoma.data.by.time[[i]], min.cells = 3, min.genes = 19, 
                                                   project = paste("Melanoma_timepoint", orig.times[i], sep = ""))
  melanoma.data.by.time[[i]]@scale.data <- melanoma.data.by.time[[i]]@data
  melanoma.data.by.time[[i]]@var.genes <- clustering.var
  melanoma.data.by.time[[i]] <- RunPCA(object = melanoma.data.by.time[[i]], do.print = FALSE)
  melanoma.data.by.time[[i]] <- FindClusters(object = melanoma.data.by.time[[i]], reduction.type = "pca", dims.use = 1:6, 
                                             resolution = 1, print.output = 0, save.SNN = TRUE)
  clustering.by.time[[i]] <- melanoma.data.by.time[[i]]@ident
}

####################################################
######## MAP CLUSTER VALUES TO SINGLE CELLS ########
####################################################

# USING SEURAT CLUSTERS FOR NOW
# Seurat clusters with all timepoints pooled, resolution = 1
seurat.cluster.assign.1
df.with.clusters <- cbind(final.df, Seurat_Pooled_Clusters = as.numeric(unname(seurat.cluster.assign.1)))

# Seurat clusters with each timepoint separately, resolution = 1
j <- 0
cluster.vector <- c()
for (i in 1:length(clustering.by.time)) {
  temp.vector <- as.numeric(unname(clustering.by.time[[i]]))
  temp.vector <- temp.vector + j
  print(unique(temp.vector))
  cluster.vector <- c(cluster.vector, temp.vector)
  j <- max(temp.vector)
  print(j)
}

for (i in 1:j) {
  print(sum(cluster.vector == i))
}

df.with.clusters <- cbind(df.with.clusters, Seurat_Separate_Cluster = cluster.vector)

######################################################################
######## CALCULATE AND SELECT TOP WITHIN-CLUSTER CORRELATIONS ########
######################################################################

library(Hmisc)

non.markers <- c("day", "Cell_Order_Number", "Seurat_Pooled_Clusters", "Seurat_Separate_Cluster")

cluster.correlations <- list()
for (i in unique(df.with.clusters$Seurat_Pooled_Clusters)) {
  this.cluster <- df.with.clusters[df.with.clusters$Seurat_Pooled_Clusters == i, ]
  this.cluster.subset <- this.cluster[, setdiff(colnames(this.cluster), non.markers)]
  cor(this.cluster.subset, method = "pearson")
  temp.pearson.corr <- rcorr(as.matrix(this.cluster.subset), type = "pearson")
  temp.pearson.corr$P[lower.tri(temp.pearson.corr$P, diag = FALSE)] <- NA
  signif.corr <- which(temp.pearson.corr$P < 0.01, arr.ind = TRUE)
  param.names <- c()
  correlation.values <- c()
  for (j in 1:nrow(signif.corr)) {
    this.name <- colnames(temp.pearson.corr$P)[unname(signif.corr[j, ])]
    this.name <- paste(this.name, collapse = "vs")
    param.names <- c(param.names, this.name)
    correlation.values <- c(correlation.values,
                            temp.pearson.corr$r[unname(signif.corr[j, ])[1],
                                                unname(signif.corr[j, ])[2]])
  }
  this.clus.corr <- data.frame(correlation.values)
  rownames(this.clus.corr) <- param.names
  cluster.correlations[[i]] <- this.clus.corr
}

all.signif.corr <- c()
for (i in 1:length(cluster.correlations)) {
  all.signif.corr <- c(all.signif.corr, rownames(cluster.correlations[[i]]))
}
pooled.signif.corr.freq <- table(all.signif.corr)
pooled.signif.corr.freq <- sort(pooled.signif.corr.freq)

cluster.correlations.by.time <- list()
for (i in unique(df.with.clusters$Seurat_Separate_Cluster)) {
  print(sum(df.with.clusters$Seurat_Separate_Cluster == i))
  this.cluster <- df.with.clusters[df.with.clusters$Seurat_Separate_Cluster == i, ]
  this.cluster.subset <- this.cluster[, setdiff(colnames(this.cluster), non.markers)]
  cor(this.cluster.subset, method = "pearson")
  temp.pearson.corr <- rcorr(as.matrix(this.cluster.subset), type = "pearson")
  temp.pearson.corr$P[lower.tri(temp.pearson.corr$P, diag = FALSE)] <- NA
  signif.corr <- which(temp.pearson.corr$P < 0.01, arr.ind = TRUE)
  param.names <- c()
  correlation.values <- c()
  for (j in 1:nrow(signif.corr)) {
    this.name <- colnames(temp.pearson.corr$P)[unname(signif.corr[j, ])]
    this.name <- paste(this.name, collapse = "vs")
    param.names <- c(param.names, this.name)
    correlation.values <- c(correlation.values,
                            temp.pearson.corr$r[unname(signif.corr[j, ])[1],
                                                unname(signif.corr[j, ])[2]])
  }
  this.clus.corr <- data.frame(correlation.values)
  rownames(this.clus.corr) <- param.names
  cluster.correlations[[i]] <- this.clus.corr
}

all.signif.corr <- c()
for (i in 1:length(cluster.correlations)) {
  all.signif.corr <- c(all.signif.corr, rownames(cluster.correlations[[i]]))
}
sep.signif.corr.freq <- table(all.signif.corr)
sep.signif.corr.freq <- sort(sep.signif.corr.freq)

top.n.per <- 0.20
top.corr.separate <- tail(names(sep.signif.corr.freq), n = (length(sep.signif.corr.freq) * top.n.per))
top.corr.pooled <- tail(names(pooled.signif.corr.freq), n = (length(pooled.signif.corr.freq) * top.n.per))

setdiff(top.corr.separate, top.corr.pooled)
setdiff(top.corr.pooled, top.corr.separate)
chosen.corr <- intersect(top.corr.separate, top.corr.pooled)

#####################################################################
######## ADD CHOSEN CLUSTER CORRELATIONS TO INDIVIDUAL CELLS ########
#####################################################################

# calculating correlations within time-separate Seurat clusters
# add as marker ("MITFvsMart1") in df.with.clusters to feed into FLOW-MAP
df.with.corr <- df.with.clusters
chosen.corr.list <- strsplit(chosen.corr, split = "vs")

for (i in 1:length(chosen.corr.list)) {
  this.corr <- chosen.corr.list[[i]]
  corr.vals <- c()
  for (j in unique(df.with.clusters$Seurat_Separate_Cluster)) {
    this.cluster <- df.with.clusters[df.with.clusters$Seurat_Separate_Cluster == j, ]
    corr.val <- cor(this.cluster[, this.corr[1]], this.cluster[, this.corr[2]], method = "pearson")
    corr.vals <- c(corr.vals, rep(corr.val, times = nrow(this.cluster)))
  }
  df.with.corr <- cbind(df.with.corr, corr.vals)
  colnames(df.with.corr)[ncol(df.with.corr)] <- chosen.corr[i]
}

#######################################################
######## GENERATION OF FLOW-MAPS WITH CLUSTERS ########
#######################################################

mode <- "single"
save.folder <- "/Users/mesako/Desktop/20180113_Melanoma_Project"
per <- 1
seed.X <- 3
# markers with highest variance
clustering.var <- c("Ki67", "Mart1", "HIF1a", "LDH", "AMPKA",
                    "p_ERK1", "PFK", "p_ACAC", "Slug", "p_LKB")
project.name <- "Melanoma_cluster_by_topvar_logtransform_clusters_topcorrelations"

clustering <- FALSE
name.sort <- FALSE
downsample <- FALSE
savePDFs <- TRUE
which.palette <- "bluered"
time.col.label <- "day"
condition.col.label <- NULL
minimum <- 2
maximum <- 20
distance.metric <- "euclidean"

df <- FLOWMAPR::RestructureDF(df.with.corr, time.col.label = time.col.label, 
                              condition.col.label = condition.col.label)$new.df
graph <- FLOWMAPR::FLOWMAPfromDF(mode = mode, df = df, project.name = project.name,
                                 time.col.label = time.col.label, condition.col.label = condition.col.label,
                                 clustering.var = clustering.var, distance.metric = distance.metric,
                                 minimum = minimum, maximum = maximum, per = per, save.folder = save.folder,
                                 name.sort = name.sort, clustering = clustering, seed.X = seed.X,
                                 savePDFs = savePDFs, which.palette = which.palette)
