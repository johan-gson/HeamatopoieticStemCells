library(tidyverse)
library(monocle3)

#ok, so install monocle3 instead
#BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                       'limma', 'S4Vectors', 'SingleCellExperiment',
#                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))

#install.packages("devtools")
#devtools::install_github('cole-trapnell-lab/leidenbase')
#devtools::install_github('cole-trapnell-lab/monocle3')

#BiocManager::install("slingshot")

baseDir = "C:/work/R/MonocleStemCells/HeamatopoieticStemCells/" #path to the repo
fig____path = paste0(baseDir, "figures/") #make sure to create this folder

#All cells
#####################

cmTibb = read_csv(paste0(baseDir, "data/All cellsSCC_mod.csv"))
#cmTibb = read_csv(paste0(baseDir, "data/CD34+cd38- most immature cells_one+ or two+Seurat_mod.csv"))

startMeta = which(colnames(cmTibb) == "SampleID")

cmDf = as.data.frame(cmTibb[,2:(startMeta-1)])
rownames(cmDf) = cmTibb[[1]]

#metaDf = as.data.frame(cmTibb[,c(startMeta, dim(cmTibb)[2])])
metaDf = as.data.frame(cmTibb[,-1])
#log transform the expression data in meta
metaDf[,1:(startMeta-2)] = log2(metaDf[,1:(startMeta-2)] + 1)

rownames(metaDf) = cmTibb[[1]]

#not sure what this is good for
metaGenes = as.data.frame(colnames(cmDf))
colnames(metaGenes) = "gene_short_name"
rownames(metaGenes) = colnames(cmDf)
dd = t(as.matrix(cmDf))

cds = new_cell_data_set(as(dd, "sparseMatrix"), cell_metadata = metaDf, gene_metadata = metaGenes)

## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 15)

## Step 2: Remove batch effects with cell alignment
cds <- align_cds(cds, alignment_group = "SampleID")

## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)

## Step 4: Cluster the cells
cds <- cluster_cells(cds)

#plot_cells(cds)

## Step 5: Learn a graph
cds <- learn_graph(cds)

## Step 6: Order cells
cds <- order_cells(cds)

plot_cells(cds)
colData(cds)
plot_cells(cds, color="CD38 (Ab)")
plot_cells(cds, color="CD34")
plot_cells(cds, color="CD38")

plot_cells(cds, color="HBB")
plot_cells(cds, color="CXCR4")
plot_cells(cds, color="CD19")

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

plot_cells(cds, color="SampleID")

plot_cells(cds, color="Seurat_Clusters_VOEB")


#So, for some reason this call fails here, but not for the immature cells below.
#I get "Error: 'rBind' is defunct.". This is described here: 
# https://github.com/cole-trapnell-lab/monocle3/issues/509
#maybe it helps to install the develop branch of monocle3, may be worth a try
cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
#extract significant genes (using multiple testing, FDR) and sort
signGenes <- subset(cds_pr_test_res, q_value < 0.05)
srt = sort(signGenes$q_value, index.return = TRUE)
signGenesSort = signGenes[srt$ix,]

signGenesSort

write_tsv(signGenesSort, file=paste0(baseDir, "results/DEGenesAll.txt"))

#the most immature cells
##############################

set.seed(1)

cmTibb2 = read_csv(paste0(baseDir, "data/CD34+cd38- most immature cells_one+ or two+Seurat_mod.csv"))

startMeta2 = which(colnames(cmTibb2) == "SampleID")

cmDf2 = as.data.frame(cmTibb2[,2:(startMeta2-1)])
rownames(cmDf2) = cmTibb2[[1]]

metaDf2 = as.data.frame(cmTibb2[,-1])
#log transform the expression data in meta
metaDf2[,1:(startMeta2-2)] = log2(metaDf2[,1:(startMeta2-2)] + 1)

rownames(metaDf2) = cmTibb2[[1]]

#not sure what this is good for
metaGenes2 = as.data.frame(colnames(cmDf2))
colnames(metaGenes2) = "gene_short_name"
rownames(metaGenes2) = colnames(cmDf2)
dd2 = t(as.matrix(cmDf2))

cds2 = new_cell_data_set(as(dd2, "sparseMatrix"), cell_metadata = metaDf2, gene_metadata = metaGenes2)

## Step 1: Normalize and pre-process the data
cds2 <- preprocess_cds(cds2, num_dim = 15)

## Step 2: Remove batch effects with cell alignment
cds2 <- align_cds(cds2, alignment_group = "SampleID")

## Step 3: Reduce the dimensions using UMAP
cds2 <- reduce_dimension(cds2)

## Step 4: Cluster the cells
cds2 <- cluster_cells(cds2)

plot_cells(cds2)


## Step 5: Learn a graph
cds2 <- learn_graph(cds2)

## Step 6: Order cells
cds2 <- order_cells(cds2)

plot_cells(cds2)


plot_cells(cds2,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

plot_cells(cds2, color="CD38 (Ab)")

plot_cells(cds2, color="SampleID")

plot_cells(cds2, color="Seurat_Clusters_QE4W")


cds_pr_test_res2 <- graph_test(cds2, neighbor_graph="principal_graph", cores=4)
#extract significant genes (using multiple testing, FDR) and sort
signGenes2 <- subset(cds_pr_test_res2, q_value < 0.05)
srt2 = sort(signGenes2$q_value, index.return = TRUE)
signGenesSort2 = signGenes2[srt2$ix,]

signGenesSort2

write_tsv(signGenesSort2, file=paste0(baseDir, "results/DEGenesCD38-.txt"))


#An experiment with a straight line and development
cmTibb2 = read_csv(paste0(baseDir, "data/CD34+cd38- most immature cells_one+ or two+Seurat_mod.csv"))

startMeta2 = which(colnames(cmTibb2) == "SampleID")

cmDf2 = as.data.frame(cmTibb2[,2:(startMeta2-1)])
rownames(cmDf2) = cmTibb2[[1]]

metaDf2 = as.data.frame(cmTibb2[,-1])
#log transform the expression data in meta
metaDf2[,1:(startMeta2-2)] = log2(metaDf2[,1:(startMeta2-2)] + 1)

rownames(metaDf2) = cmTibb2[[1]]


#make a small test with a linear development trajectory
#library(Seurat)

#d = CreateSeuratObject(counts = t(cmDf2), project = "test", min.cells = 0, min.features = 0)
#d = NormalizeData(d, normalization.method = "LogNormalize", scale.factor = 10000)
#d = FindVariableFeatures(d, selection.method = "vst", nfeatures = 2000)
#d = ScaleData(d)
#d = RunPCA(d, features = VariableFeatures(object = d))
#d = FindNeighbors(d, dims = 1:10)
#d = FindClusters(d, resolution = 0.5)
#d <- RunUMAP(d, dims = 1:10)
#DimPlot(d, reduction = "umap")
#FeaturePlot(d, "CD38 (Ab)")
#Hmm, this is not batch corrected....




#try using the existing UMAP from monocle3 
##############################################

umapData = cds2@int_colData@listData$reducedDims@listData$UMAP
plot(umapData[,1],umapData[,2] ) #looks good

startPoint = c(-5,0)
endPoint = c(5,0)
dirVector = endPoint - startPoint
unitVector = dirVector /sqrt(dirVector[1]^2 + dirVector[2]^2)

#project all points onto the vector using scalar product
scProd = rep(NA, nrow(umapData))
for (i in 1:nrow(umapData)) {
  scProd[i] = umapData[i,1]* unitVector[1] + umapData[i,2]* unitVector[2]
}

#discard the cells that are too far to the right (a small cluster that will mess things up)
cellFilt = scProd > 5
sum(cellFilt)  #41 cells, discard them

cmDf2Filt = cmDf2[!cellFilt,]
scProdFilt = scProd[!cellFilt]
umapDataFilt = umapData[!cellFilt,]
plot(umapDataFilt[,1],umapDataFilt[,2] ) #looks good

srt = sort(scProdFilt, index.return = TRUE)

#make a sliding window with 300 cells that sums up the average gene expression as function of scProdFilt
#check that the order is the same:
#all(rownames(umapDataFilt) == rownames(cmDf2Filt))#TRUE, ok

windowSize = 300
nGenes = ncol(cmDf2Filt)
numCells = nrow(cmDf2Filt)
#this is a bit confusing. The genes are columns in the cmDf2Filt matrix - in the 
#expr matrix I'm changing that to rows
numWindows = numCells-windowSize + 1
expr = matrix(NA, nrow = nGenes, ncol = numWindows)
cdProdVals = rep(NA, numWindows)
sortedCmDf2 = cmDf2Filt[srt$ix,]

#This takes a minute or so to run
for(windowStart in 1:ncol(expr)) {
  expr[,windowStart] = as.numeric(colSums(sortedCmDf2[windowStart:(windowStart+windowSize-1),]))
  cdProdVals[windowStart] = mean(srt$x[windowStart:(windowStart+windowSize-1)])
}

exprCPM = expr
#normalize the expression to CPM
for(windowStart in 1:ncol(expr)) {
  exprCPM[,windowStart] = expr[,windowStart] * 10^6 / sum(expr[,windowStart])
}

#test: 
#colSums(exprCPM) # ok

colnames(cmDf2Filt)
#"CD38 (Ab)" 1-51 are Ab, 30 is CD38
#which(colnames(cmDf2Filt) == "CD38")
#plot(cdProdVals, exprCPM[30,])
#plot(cdProdVals, exprCPM[colnames(cmDf2Filt) == "CD38",])
#plot(cdProdVals, exprCPM[1,]) #nice
#plot(cdProdVals, exprCPM[2,]) #perhaps
#plot(cdProdVals, exprCPM[3,]) #perhaps
#plot(cdProdVals, exprCPM[4,]) #perhaps
#plot(cdProdVals, exprCPM[12,]) #perhaps


plotFacetWrap = function(indices) {
  numGenes = length(indices)
  dd = rep(NA, numWindows*length(indices))
  for (i in 1:numGenes) {
    offs = numWindows*(i-1)
    dd[(offs+1):(offs + numWindows)] = exprCPM[indices[i],]
  }
  
  fac = factor(rep(indices, each=numWindows), indices, colnames(cmDf2Filt)[indices])
  x = rep(cdProdVals, numGenes)
  
  
  df = tibble(x = x, y = dd, gene=fac)
  
  p <- ggplot(df, aes(x=x, y=y, )) + geom_line() + facet_wrap(vars(gene), scales="free")
  
  return(p)  
  
}






#plot all Ab using a facet grid in ggplot ( gene 1-51)
numGenes = 51
dd = rep(NA, numWindows*numGenes)
for (i in 1:numGenes) {
  offs = numWindows*(i-1)
  dd[(offs+1):(offs + numWindows)] = exprCPM[i,]
}

fac = factor(rep(1:numGenes, each=numWindows), 1:numGenes, colnames(cmDf2Filt)[1:numGenes])
x = rep(cdProdVals, numGenes)


df = tibble(x = x, y = dd, gene=fac)

p <- ggplot(df, aes(x=x, y=y, )) + geom_line() + facet_wrap(vars(gene), scales="free")
p

p1_51 = plotFacetWrap(1:51)
ggsave(
  paste0(fig____path, "FW_1_51.png"),
  plot = p1_51,
  width = 14, height = 12, dpi = 300)

p52_100 = plotFacetWrap(52:100)
ggsave(
  paste0(fig____path, "FW_52_100.png"),
  plot = p52_100,
  width = 14, height = 12, dpi = 300)

p101_200 = plotFacetWrap(101:200)
ggsave(
  paste0(fig____path, "FW_101_200.png"),
  plot = p101_200,
  width = 14, height = 20, dpi = 300)

p201_300 = plotFacetWrap(201:300)
ggsave(
  paste0(fig____path, "FW_201_300.png"),
  plot = p201_300,
  width = 14, height = 20, dpi = 300)

p301_400 = plotFacetWrap(301:400)
ggsave(
  paste0(fig____path, "FW_301_400.png"),
  plot = p301_400,
  width = 14, height = 20, dpi = 300)

p401_500 = plotFacetWrap(401:500)
ggsave(
  paste0(fig____path, "FW_401_500.png"),
  plot = p401_500,
  width = 14, height = 20, dpi = 300)

p501_600 = plotFacetWrap(501:600)
ggsave(
  paste0(fig____path, "FW_501_600.png"),
  plot = p501_600,
  width = 14, height = 20, dpi = 300)

p601_648 = plotFacetWrap(601:648)
ggsave(
  paste0(fig____path, "FW_601_648.png"),
  plot = p601_648,
  width = 14, height = 20, dpi = 300)





print(colnames(metaDf2)[500:length(metaDf2)],max.print=2000)
#UMAP_X_TZ8V is the one I should use!

#Seurat_Clusters_CDOJ
#Seurat_Clusters_C3N0
#Seurat_Clusters_3X55
#Seurat_Clusters_ZNRE
#Seurat_Clusters_VOEB
#Seurat_Clusters_AX5O
#Seurat_Clusters_TZ8V
#Seurat_Clusters_6QLE
#Seurat_Clusters_QE4W

#Test slingshot


