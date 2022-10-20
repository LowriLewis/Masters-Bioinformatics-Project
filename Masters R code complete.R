# lIBRARIES
library(pryr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(clustree)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
library(patchwork)
library(remotes)
library(limma)
library(stringr)
library(tidyr)
library(Matrix)
library(ggpubr)
library(viridis)
library(ggpubr)
library(BiocManager)
library(slingshot)
library(tidyverse)
library(tidymodels)
library(scales)
library(SeuratObject)
library(slingshot)
library(Matrix)



#CODE
setwd("........")

###Reading in the Ramos Data Set
ramosData <- read.delim("...", sep =";", header = T)

###Renaming Column X to "Gene_Name"
ramosData <- rename(ramosData, Gene_Name=X)

###Setting Gene Name as rownames and removing duplicate column
rownames(ramosDataFile) <- ramosDataFile$Gene_Name
ramosDataFile <- ramosDataFile[,-1]

###Assessing counts
counts_per_cell <- Matrix::colSums(ramosData) 
counts_per_gene <- Matrix::rowSums(ramosData) 
genes_per_cell <- Matrix::rowSums(ramosData>0)

###Creating a Seurat Object 
ramosSeurat <- CreateSeuratObject(counts=ramosData)

###Reading in the 2 sets of Metadata
Meta <- read.delim("GSE192935_Metadata_Human_SC_RNAseq_Fig1_and_2.csv.gz", sep=";", header=TRUE, row.names = "X")
Meta2 <- read.csv("Metadata_PatientNo.Realloc.csv", header = T) ####Patient numbers were changed from original input of "1-6" to "1-4, 6-7" as patients in this study were referred to as 1,2,3,4,6,7

###Separating Metadata column orig.ident to obtain patient number alone (to bind metadata together) and removing redundant information generated
Meta <- cbind(Meta, strsplit2(Meta$orig.ident, split="_"))
drop <- c("1", "2")
Meta <- Meta[,!(names(Meta) %in% drop)]
#Renaming the patient column
names(Meta)[7] <- "Patient"
#Keeping only the last character of the Patient column (i.e. the number)
Meta$Patient <- as.integer(str_sub(Meta$Patient, -1))
#Renaming the Patient column of metadata 2 to match the other metadata
names(Meta2)[1] <- "Patient"

###Merging Metadata
allMeta <- left_join(Meta, Meta2, by="Patient")
#Setting rownames again (as they were lost in merging)
rownames(allMeta) <- rownames(Meta)

###Adding Metadata to Seurat Object
ramosSeurat <- AddMetaData(ramosSeurat, allMeta)

###Identifying percent Mitochondrial and Ribosomal genes - they have already worked this out (for mito) but better for us to determine it ourselves too
ramosSeurat[["percent.mt.me"]] <- PercentageFeatureSet(ramosSeurat, pattern = "^MT-")
ramosSeurat[["percent.ribo"]] <- PercentageFeatureSet(ramosSeurat, pattern = "^RPL|^RPS")

###Checking mt.me and their percent mito match
sum(ramosSeurat[["percent.mt.me"]]) == sum(ramosSeurat[["percent.mito"]]) #Error as not both numerically alike
head(ramosSeurat[["percent.mt.me"]])
head(ramosSeurat[["percent.mito"]]) #can see now that their decimal places are commas, and also 100x less than my values.

#Turning their comma decimal separators for %mito to full stops in their percent.mito metadata (to turn it into a numeric vector)
ramosSeurat$percent.mito <- as.numeric(gsub(",", ".", ramosSeurat$percent.mito))     
#Multipllying the entire column by 100 to get them to the same values 
ramosSeurat$percent.mito <- ramosSeurat$percent.mito *100

sum(ramosSeurat[["percent.mt.me"]]) == sum(ramosSeurat[["percent.mito"]]) #Returns true -> can see after reading paper that they already removed mito%>10-20% and also <200genes per cell -> won't filter further for m=now but may end up doing it depending on QC

###Visualising QC metrics as a violin plot
VlnPlot(ramosSeurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
###Checking for batch effects - there are some so need to normalise and scale separately then integrate 
VlnPlot(ramosSeurat, features = c("nFeature_RNA"), group.by = "orig.ident") + ggtitle("nFeature_RNA, grouped by orig.ident")
VlnPlot(ramosSeurat, features = c("nCount_RNA"), group.by = "orig.ident") + ggtitle("nCount_RNA, grouped by orig.ident")
VlnPlot(ramosSeurat, features = c("percent.mt.me"), group.by = "orig.ident") + ggtitle("PercentMT, grouped by orig.ident")
VlnPlot(ramosSeurat, features = c("nFeature_RNA"), group.by = "Patient") + ggtitle("nFeature_RNA, grouped by Patient")
VlnPlot(ramosSeurat, features = c("nCount_RNA"), group.by = "Patient") + ggtitle("nCount_RNA, grouped by Patient")
VlnPlot(ramosSeurat, features = c("percent.mt.me"), group.by = "Patient") + ggtitle("PercentMT, grouped by Patient")
VlnPlot(ramosSeurat, features = c("nFeature_RNA"), group.by = "Tissue") + ggtitle("nFeature_RNA, grouped by Tissue")
VlnPlot(ramosSeurat, features = c("nCount_RNA"), group.by = "Tissue") + ggtitle("nCount_RNA, grouped by Tissue")
VlnPlot(ramosSeurat, features = c("percent.mt.me"), group.by = "Tissue") + ggtitle("PercentMT, grouped by Tissue")

plot1 <- FeatureScatter(ramosSeurat, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(ramosSeurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2 #plots show that we could maybe cut down on mito% further, but in doing so would risk cutting out large amount of blood sample...

###Identifying the most highly variable features
ramosSeurat <- FindVariableFeatures(ramosSeurat, nfeatures = 2000)
#Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(ramosSeurat), 10)
# plot variable features with top10 most variable genes labelled
plot1 <- VariableFeaturePlot(ramosSeurat, cols = c("lightgrey", "navyblue"))
LabelPoints(plot = plot1, points = top10, repel = TRUE) +ggtitle("2000 most variable genes pre-normalisation")+theme(axis.line = element_line(size = 1), plot.title = element_text(hjust=0.5))

###Normalising and Scaling RNA data slot for plotting dimensionality reduction differential expression analysis
ramosSeurat <- NormalizeData(ramosSeurat, scale.factor=10000)
ramosSeurat <- ScaleData(ramosSeurat)
###Elbow plot to identify the most descriptive PCs
ramosSeurat <- RunPCA(object = ramosSeurat, npcs = 100)
PCAPlot(ramosSeurat, group.by="Tissue")
ElbowPlot(ramosSeurat, ndims=100) +ggtitle("Elbow plot for dimensionality reduction") +theme(plot.title = element_text(hjust=0.5)) + geom_segment(aes(x=50,xend=50,yend=1.5,y=5), arrow=arrow(length = unit(0.5, "cm")))
#looks like most information is captured by first 50 dimensions- will use dimensions 1:50 from now on

###Subsetting the Seurat into samples for normalisation
BLD3 <- subset(ramosSeurat, subset = (Tissue=="Blood" & Patient=="3"))
LN1 <- subset(ramosSeurat, subset = (Tissue=="Lymph_Node" & Patient=="1"))
LN3 <- subset(ramosSeurat, subset = (Tissue=="Lymph_Node" & Patient=="3"))
LN5 <- subset(ramosSeurat, subset = (Tissue=="Lymph_Node" & Patient=="5"))
LN6 <- subset(ramosSeurat, subset = (Tissue=="Lymph_Node" & Patient=="6"))
LN7 <- subset(ramosSeurat, subset = (Tissue=="Lymph_Node" & Patient=="7"))
TM2 <- subset(ramosSeurat, subset = (Tissue=="Tumor" & Patient=="2"))
TM3 <- subset(ramosSeurat, subset = (Tissue=="Tumor" & Patient=="3"))
TM5 <- subset(ramosSeurat, subset = (Tissue=="Tumor" & Patient=="5"))
TM6 <- subset(ramosSeurat, subset = (Tissue=="Tumor" & Patient=="6"))

###Creating a list of the separated Seurat samples
SC.list = list(BLD3, LN1, LN3, LN5, LN6, LN7, TM2, TM3, TM5, TM6)
names(SC.list) = c("10X_BLD_Pt3","10X_LN_Pt1", "10X_LN_Pt3", "10X_LN_Pt5", "10X_LN_Pt6", "10X_LN_Pt7",
                   "10X_TM_Pt2", "10X_TM_Pt3", "10X_TM_Pt5", "10X_TM_Pt6")

###Normalisation and scaling with SCTransform separately for cluster analysis (on a loop so all is done in one function)
for (i in 1:length(x = SC.list)) {
  SC.list[[i]] <- SCTransform(object = SC.list[[i]])
}

###Finding integration anchors
SC.Features <- SelectIntegrationFeatures(object.list = SC.list)
SC.list <- PrepSCTIntegration(SC.list, anchor.features = SC.Features)
SC.list <- lapply(SC.list, FUN=RunPCA, features=SC.Features)
SC.anchors <- FindIntegrationAnchors(object.list = SC.list, normalization.method="SCT", reduction="rpca", anchor.features=SC.Features)
###Integrating the samples
SC.integrated <- IntegrateData(anchorset = SC.anchors, normalization.method = "SCT")

###UMAP plotting
DefaultAssay(SC.integrated)
MaxDim = 50

SC.integrated <- RunPCA(object = SC.integrated)
SC.integrated <- RunUMAP(object = SC.integrated, reduction = "pca", dims = 1:MaxDim)
DimPlot(object = SC.integrated, reduction = "umap", group.by = "orig.ident", split.by="Tissue") + ggtitle("UMAP of all samples across tissues") +theme(plot.title = element_text(hjust=0.5, size=15), axis.title = element_text(size=12), axis.line = element_line(size=1.2))

DefaultAssay(SC.integrated) <- "integrated"
clusterRes <- c(seq(0,0.15, 0.05), seq(0.2,1,0.1))

###Finding out clusters through running Find neighbours and find clusters
SC.integrated <- FindNeighbors(object = SC.integrated, dims = 1:MaxDim)
SC.integrated <- FindClusters(object = SC.integrated, resolution= clusterRes)
clustree(SC.integrated, prefix = "integrated_snn_res.")  #Shows how the clusters split over different resolutions

###Plotting all of the different clusters at resolutions 0-1.0
plotList<-list()
for (j in 1:length(clusterRes)) {
  Idents(SC.integrated) <- paste0("RNA_snn_res.", clusterRes[j])
  plotList[[j]]<-DimPlot(object=SC.integrated, reduction = "umap", label=TRUE, group.by = paste0("integrated_snn_res.", clusterRes[j])) +
    ggtitle(paste0("Resolution ", clusterRes[j])) + theme(legend.position = "null", plot.title=element_text(vjust=0.5))
}

grid.arrange(grobs=plotList) ###Resolution of 0.1 or 0.15 look good? population highlighted stemming off main bulk? ---- This population likely to be non-classical monos looking at FCGR3A expression in dotplots below



#Assessing % Mito and Ribo gene clustering
FeaturePlot(SC.integrated, features = "percent.mito", split.by = "Tissue", pt.size=0.2) +theme(legend.position = "right") + labs(col="Expression (%)")
FeaturePlot(SC.integrated, features = "percent.ribo", split.by = "Tissue", pt.size=0.2) + theme(legend.position = "right") + labs(col="Expression (%)")
###May be some high mito% RNA clustering (see x=2,y=-3 cluster in blood and tumour... will keep in mind downstream!)

###Dotplots for cluster identification - determining the optimal resolution to use for clustering
DefaultAssay(SC.integrated) <- "RNA"
stress <- c("HSPA1A", "HSPB1")
hypoxia <- c("SIGLEC1","TNFRSF1B","ADORA2B","HCAR3","IL6ST","LILRA2","PLAUR","TNFRSF10D","CLEC5A", "TREM1")
markers <- c('FCER1A', 'CST3', "FCGR3A", "CD14", "MS4A7", "LYZ", "APOE", "TREM2","FOLR2", "APOC1", "CD163")
cycling <- c("MKI67", "TOP2A", "CDC20")
FCGRs <- c("FCGR1A", "FCGR1B", "FCGR2A", "FCGR2B", "FCGR3A", "FCGR3B")

plotList <- list()
plotList[[1]] <- DotPlot(SC.integrated, features = FCGRs) + ggtitle("FCGR expression") + labs(y="Cluster", x="Gene")+ theme(plot.title = element_text(hjust=0.5, size=25), axis.text.x = element_text(size=20, angle=60, vjust = 1, hjust=1), axis.title = element_text(size=20), axis.text.y = element_text(size=15), axis.line = element_line(size=1.5))
plotList[[2]] <- DotPlot(SC.integrated,features=stress) + ggtitle("Markers of stress") + labs(y="Cluster", x="Gene")+ theme(plot.title = element_text(hjust=0.5, size=25), axis.text.x = element_text(size=20),  axis.title = element_text(size=20), axis.text.y = element_text(size=15), axis.line = element_line(size=1.5))
plotList[[3]] <- DotPlot(SC.integrated,features=markers) + ggtitle("Markers of interest") + labs(y="Cluster", x="Gene")+ theme(plot.title = element_text(hjust=0.5, size=25), axis.text.x = element_text(size=20, angle=60, vjust = 1, hjust=1),  axis.title = element_text(size=20), axis.text.y = element_text(size=15), axis.line = element_line(size=1.5))
plotList[[4]] <- DotPlot(SC.integrated,features=hypoxia) + ggtitle("Markers of hypoxia") + labs(y="Cluster", x="Gene")+ theme(plot.title = element_text(hjust=0.5, size=25), axis.text.x = element_text(size=20, angle=60, vjust = 1, hjust=1),  axis.title = element_text(size=20), axis.text.y = element_text(size=15), axis.line = element_line(size=1.5))
plotList[[5]] <- DotPlot(SC.integrated,features=cycling) + ggtitle("Markers of cell cycling") + labs(y="Cluster", x="Gene")+ theme(plot.title = element_text(hjust=0.5, size=25), axis.text.x = element_text(size=20, angle=60, vjust = 1, hjust=1),  axis.title = element_text(size=20), axis.text.y = element_text(size=15), axis.line = element_line(size=1.5))

grid.arrange(grobs=plotList)


###Feature Plots for identifying clusters --NOTE - all legends have different scales!
###FCGRs
plotList <- list()
for(i in 1:length(FCGRs)) {
  plotList[[i]] <- FeaturePlot(SC.integrated, features=FCGRs[[i]], pt.size=0.2)  + ggtitle(FCGRs[i])
}
grid.arrange(grobs=plotList)

###HYPOXIA
plotList <- list()
for(i in 1:length(hypoxia)) {
  plotList[[i]] <- FeaturePlot(SC.integrated, features=hypoxia[[i]], pt.size=0.2)  + ggtitle(hypoxia[i])
}
grid.arrange(grobs=plotList)

###STRESS
plotList <- list()
for(i in 1:length(stress)) {
  plotList[[i]] <- FeaturePlot(SC.integrated, features=stress[[i]], pt.size=0.3)  + ggtitle(stress[i])
}
grid.arrange(grobs=plotList)

###CYCLING
plotList <- list()
for(i in 1:length(cycling)) {
  plotList[[i]] <- FeaturePlot(SC.integrated, features=cycling[[i]], pt.size=0.3)  + ggtitle(cycling[i])
}
grid.arrange(grobs=plotList)

###MARKERS
plotList <- list()
for(i in 1:length(markers)) {
  plotList[[i]] <- FeaturePlot(SC.integrated, features=markers[[i]], pt.size=0.3)  + ggtitle(markers[i])
}
grid.arrange(grobs=plotList)
#######Looking at results (cycling cells especially) - resolution 0.1 may be better - defines the population of cycling cells nicely and also nice showing of macrophage cluster

###Violinplots for identifying clusters
VlnPlot(SC.integrated, assay="RNA", features=stress, slot = "data")
VlnPlot(SC.integrated, assay="RNA", features=hypoxia, slot = "data")
VlnPlot(SC.integrated, assay="RNA", features=markers, slot= "data") 
VlnPlot(SC.integrated, assay="RNA", features=cycling, slot="data")
VlnPlot(SC.integrated, assay="RNA", features=FCGRs, slot="data")

#Feature plot splitting by tissue 
DefaultAssay(SC.integrated) <- "RNA"
FeaturePlot(object = SC.integrated, features=FCGRs, reduction = "umap", split.by = "Tissue", slot = "data", pt.size = 0.3)

##Assessing differences in hypoxia based on patient
plotList <- list()
for(i in 1:length(hypoxia)) {
  plotList[[i]] <- VlnPlot(SC.integrated, features = hypoxia[[i]], assay="RNA", split.by = "Patient" , cols = viridis(6)) + ggtitle(hypoxia[[i]]) + labs(x="Cluster", y="Expression") + theme(plot.title=element_text(size=25),axis.text = element_text(size=15), axis.title = element_text(size=20), axis.line = element_line(size = 1.5))
}

grid.arrange(grobs=plotList[])

###Determining cluster identities using Res of 0.1
DefaultAssay(object = SC.integrated) <- "integrated"
Idents(SC.integrated) <- SC.integrated@meta.data[["integrated_snn_res.0.1"]]

#Finding Markers between clusters and selecting those with adjusted p value of <5%
SC.integrated <- PrepSCTFindMarkers(SC.integrated)
allMarkers <- FindAllMarkers(SC.integrated, test.use="wilcox", only.pos = T, min.pct = 0.1, logfc.threshold = 0.1)
SC.markers <- allMarkers[which(allMarkers$p_val_adj<0.05),]

#Extracting the defining cell markers from those that are differentially expressed between them
top10_markers <- as.data.frame(SC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC))
top10_markers

###Asssigning and extracting markers 
marker_genes_pan_T=c('CD3D','CD3G')
marker_genes_naive_CD4_T = c('IL7R','CCR7')
marker_genes_CD8_T = c('CD8A')
marker_genes_Classical_Mono=c('CD14', 'LYZ')
marker_genes_Non_Classical_Mono =c('MS4A7') #Including FCGR3A separately in dotplot, not as
marker_genes_NK = c('GNLY', 'NKG7')
marker_genes_DC=c('FCER1A', 'CST3')
marker_genes_B=c('MS4A1')
marker_genes_Platelets=c('PPBP')
marker_genes_stressed=c("HSPA1A", "HSPB1")
marker_genes_cycling=c("MKI67", "TOP2A", "CDC20")

df.markers<-data.frame("celltype"=c(rep('pan_T', length(marker_genes_pan_T)),
                                    rep('naive_CD4_T', length(marker_genes_naive_CD4_T)),
                                    rep('CD8_T', length(marker_genes_CD8_T)), 
                                    rep('Classical_Mono',length(marker_genes_Classical_Mono)), 
                                    rep('Non_Classical_Mono', length(marker_genes_Non_Classical_Mono)), 
                                    rep('NK',length(marker_genes_NK)), 
                                    rep('B-cell',length(marker_genes_B)), 
                                    rep('Platelets', length(marker_genes_Platelets)),
                                    'NK/Non_Classical_Mono', rep('Stressed cells', length(marker_genes_stressed)), rep('Cycling cells', length(marker_genes_cycling))),
                       "gene"=c(marker_genes_pan_T,
                                marker_genes_naive_CD4_T, 
                                marker_genes_CD8_T,
                                marker_genes_Classical_Mono,
                                marker_genes_Non_Classical_Mono,
                                marker_genes_NK,
                                marker_genes_B,
                                marker_genes_Platelets, 
                                'FCGR3A', marker_genes_stressed, marker_genes_cycling))

markerNames <- allMarkers[allMarkers$gene %in% df.markers$gene,]
markerNames

####CLUSTERS - 0 = Classical monos, 1 = macrophages (and TAMS), 2 = Macrophages, 3 = non-classical monos, 4 = cycling cells

###Want to subcluster the blood population to ensure we start trajectory analysis from the most circulating-like monocyte population

####subsetting the blood population
SCBlood <- subset(SC.integrated, subset=(Tissue=="Blood"))

###Reclustering the blood population to find where CD14 expression is highest...
plotList <- list()
for(i in 1:length(idents)){
  Idents(SCBlood) <- paste0("integrated_snn_res.", idents[[i]])
  plotList[[i]] <- DimPlot(SCBlood, reduction="umap", label=T) + ggtitle(paste0("Blood cluster Res", idents[i]))
}

grid.arrange(grobs=plotList) ###Resolution 0.7 looks best... clearer classical monocyte population and cycling cell population is separated


### SLINGSHOT TRAJECTORY ANALYSIS - on resolution 0.7 on the integrated data set again
Idents(macroMono) <- macroMono$integrated_snn_res.0.7

#Extracting data into format required
dimRed <- SC.integrated@reductions$umap@cell.embeddings
clustering <- SC.integrated$integrated_snn_res.0.7

#Conducting trajectory analysis
sds <- slingshot(dimRed, clusterLabels = clustering, start.clus=5)
lineages <- getLineages(dimRed, clustering, start.clus=5)
lineages <- as.SlingshotDataSet(lineages)

#Plotting the lineage
pal <- viridis(15)
plot(dimRed, col=pal[clustering], asp=0.5, cex=0.2, pch=16)
lines(lineages, lwd=3, col='black')
#Applying curves to trajectory analysis
curves <- getCurves(lineages)
curves <- as.SlingshotDataSet(curves)

#Plotting with curves
plot(dimRed, col=pal[clustering], asp=0.5, cex=0.2, pch=16)
lines(curves, lwd=3, col='black')


####Determining what the new different clusters are
allMarkersRes0.7 <- FindAllMarkers(SC.integrated, logfc.threshold = 0.1, min.pct = 0.1, only.pos = T)
allMarkersRes0.7 <- allMarkersRes0.7[which(allMarkersRes0.7$p_val_adj<0.05),]

top5_markers <- as.data.frame(allMarkersRes0.7 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))
top5_markers


###Assigning new cluster names
Idents(SC.integrated) <- "integrated_snn_res.0.7"
levels(SC.integrated)

SC.integrated <- RenameIdents(SC.integrated, "0" = "Neutrophils (MXD1+)" = "0", 
                              "1" = "Myeloid DCs" = "1", 
                              "2" = "Classical monos" = "2", 
                              "3" =  "Classical monos (THBS1+)" = "3",
                              "4" = "Intermediate monos" = "4", 
                              "5" = "Myeloid DCs" = "5", 
                              "6" = "TAMS (APOE+)" = "6", 
                              "7" = "Myeloid DCs" = "7", 
                              "8" = "Neutrophils (CXCL8)" = "8", 
                              "9" = "TAMS (APOC1+)"= "9", 
                              "10" = "Macrophages (SPP1+)"= "10",  
                              "11" = "Non-classical monos" = "11", 
                              "12" = "Macrophages/TAMS (IDO1+)" = "12", 
                              "13" = "Differentiating monos"= "13", 
                              "14" = "Cycling cells"= "14")

###UMAP including all cell cluster names 
DimPlot(SC.integrated, reduction="umap", label=T, label.size = 3, repel = T) + ggtitle("All samples cell clusters (res0.7)") +theme(plot.title = element_text(hjust = 0.5, size=15), axis.line = element_line(size=1.2), axis.title = element_text(size=12), axis.text = element_text(size = 10))

###Subsetting the macromono population - excluding the cycling cell population too though
levels(SC.integrated)
MacroMono2 <- subset(SC.integrated, idents=c("Classical monos",
                                             "Classical monos (THBS1+)",
                                             "Intermediate monos", 
                                             "TAMS (APOE+)", 
                                             "TAMS (APOC1+)", 
                                             "Macrophages (SPP1+)",  
                                             "Non-classical monos", 
                                             "Macrophages (S100A10+)"))

###UMAP plot of the macromono population and splitting by tissue
DimPlot(MacroMono2, reduction="umap", label=T, label.size = 3, repel = T) + ggtitle("Monocytes/macrophage clustering (res0.7)") +theme(plot.title = element_text(hjust = 0.5, size=15), axis.line = element_line(size=1.2), axis.title = element_text(size=12), axis.text = element_text(size = 10))
DimPlot(MacroMono2, reduction="umap", label=T, label.size = 3, repel = T, split.by="Tissue") + ggtitle("Monocytes/macrophage clustering (res0.7)") +theme(plot.title = element_text(hjust = 0.5, size=15), axis.line = element_line(size=1.2), axis.title = element_text(size=12), axis.text = element_text(size = 10))

### SLINGSHOT TRAJECTORY ANALYSIS - on resolution 0.7 on the integrated, macromono-extracted data set - NOT WHOLE DATASET THIS TIME
levels(MacroMono2)
pal <-brewer.pal(n=9, "Set2")

#Assigning base colour to current UMAP plot to match future trajectory analysis plot
DimPlot(MacroMono2, reduction="umap", label=T, label.size = 3, repel = T, cols=pal) + ggtitle("Monocyte/Macrophage UMAP clustering (res0.7)") +theme(plot.title = element_text(hjust = 0.5, size=15), axis.line = element_line(size=1.2), axis.title = element_text(size=12), axis.text = element_text(size = 10), legend.position = "null")

DefaultAssay(MacroMono2) <- "RNA"
FeaturePlot(MacroMono2, features=c("APOC1", "APOE", "TREM2", "FOLR2")) ##TAMs and TR/infilt markers - for interest
#Extracting data into format required
dimRed <- MacroMono2@reductions$umap@cell.embeddings
clustering <- MacroMono2@active.ident
levels(clustering)
#Conducting trajectory analysis
lineages <- getLineages(dimRed, clustering, start.clus=3)
lineages <- as.SlingshotDataSet(lineages)
lineages

#Plotting the lineage
par(mfrow = c(1, 2))
plot(dimRed, col=pal[clustering], asp=1, cex=0.1, pch=16)
lines(lineages, lwd=3, col='black')

#Applying curves to trajectory analysis
curves <- getCurves(lineages)
curves <- as.SlingshotDataSet(curves)

#Plotting with curves
plot(dimRed, col=pal[clustering], asp=1, cex=0.1, pch=16)
lines(curves, lwd=3, col='black')

###Looking at monocyte markers found at ref belwo to check mono clustering - looks good! start of trajectory analysis looks good for classical monos.
###markers from https://ashpublications.org/blood/article/118/5/e16/29016/Gene-expression-profiling-reveals-the-defining
DefaultAssay(SCBlood) <- "RNA"
SCBlood <-NormalizeData(SCBlood)
SCBlood <-ScaleData(SCBlood)

IMMarkers <- c("HLA-DRA",  "CLEC10A", "CD86") #"GFRA2",
classMarkers <- c("S100A12", "CD14", "CCR2", "CD99")
nonClassMarkers <- c("SH2D1B", "FCGR3A", "VMO1", "CX3CR1")

FeaturePlot(SCBlood, features=IMMarkers, slot="data", pt.size = 0.3)

DefaultAssay(MacroMono2) <- "RNA"
DotPlot(MacroMono2, features= FCGRs) + ggtitle("FCGR expression") + labs(y="Cluster", x="Gene")+ theme(plot.title = element_text(hjust=0.5, size=25), axis.text.x = element_text(size=20, angle=60, vjust = 1, hjust=1), axis.title = element_text(size=20), axis.text.y = element_text(size=15), axis.line = element_line(size=1.5))
DotPlot(MacroMono2,features=stress) + ggtitle("Markers of stress") + labs(y="Cluster", x="Gene")+ theme(plot.title = element_text(hjust=0.5, size=25), axis.text.x = element_text(size=20),  axis.title = element_text(size=20), axis.text.y = element_text(size=15), axis.line = element_line(size=1.5))
DotPlot(MacroMono2,features=markers) + ggtitle("Cell markers") + labs(y="Cluster", x="Gene")+ theme(plot.title = element_text(hjust=0.5, size=25), axis.text.x = element_text(size=20, angle=60, vjust = 1, hjust=1),  axis.title = element_text(size=20), axis.text.y = element_text(size=15), axis.line = element_line(size=1.5))
DotPlot(MacroMono2,features=hypoxia) + ggtitle("Markers of hypoxia") + labs(y="Cluster", x="Gene")+ theme(plot.title = element_text(hjust=0.5, size=25), axis.text.x = element_text(size=20, angle=60, vjust = 1, hjust=1),  axis.title = element_text(size=20), axis.text.y = element_text(size=15), axis.line = element_line(size=1.5))
DotPlot(MacroMono2,features=cycling) + ggtitle("Markers of cell cycling") + labs(y="Cluster", x="Gene")+ theme(plot.title = element_text(hjust=0.5, size=25), axis.text.x = element_text(size=20, angle=60, vjust = 1, hjust=1),  axis.title = element_text(size=20), axis.text.y = element_text(size=15), axis.line = element_line(size=1.5))


#### Rosanna's code for professional dotplotting
marker_genes_Classical_Mono=c('CD14', 'LYZ')
marker_genes_Non_Classical_Mono =c('MS4A7', 'FCGR3A')
marker_genes_NK = c('GNLY', 'NKG7')
marker_genes_DC=c('FCER1A', 'CST3')
marker_genes_MacroTAMS=c("APOE", "TREM2","FOLR2", "APOC1", "CD163")


df.markers<-data.frame("celltype"=c(rep('Classical_Mono',length(marker_genes_Classical_Mono)), 
                                    rep('Non_Classical_Mono', length(marker_genes_Non_Classical_Mono)), 
                                    rep('NK',length(marker_genes_NK)), 
                                    rep('DC',length(marker_genes_DC)),
                                    rep('Macro_TAMS',length(marker_genes_MacroTAMS))), 
                       "gene"=c(marker_genes_Classical_Mono,
                                marker_genes_Non_Classical_Mono,
                                marker_genes_NK,
                                marker_genes_DC,
                                marker_genes_MacroTAMS))



DotPlot(SC.integrated, assay="SCT", features = df.markers$gene, scale = FALSE) + ggtitle("Markers of Interest") + 
  xlab("") +
  ylab("Cluster") +
  annotate("text", x = 1:dim(df.markers)[1], y = length(unique(SC.integrated@active.ident)) + 1, label = df.markers$celltype, angle=90, hjust=0, vjust=0.5) +
  coord_cartesian(clip="off", ylim=c(1,length(unique(SC.integrated@active.ident)) + 5 )) +
  theme(axis.text.x =element_text(angle=90, vjust = 0.5, hjust=1),
        axis.text.y =element_text(face="bold"),
        plot.title = element_text(size=20, hjust=0.5))



###CODE I DIDNT GET ROUND TO USING
####Finding markers between diifferent cell types - will use against trajectories plotted by slingshot.
classVAPOE <- FindMarkers(MacroMono2, ident.1="Classical monos", ident.2 ="TAMS (APOC1+)")
classVAPOE <- classVAPOE[which(classVAPOE$p_val_adj<0.05),]

APOEVAPOC <- FindMarkers(MacroMono2, ident.1="TAMS (APOC1+)", ident.2 ="TAMS (APOE+)")
APOEVAPOC <- APOEVAPOC[which(APOEVAPOC$p_val_adj<0.05),]

classVAPOE <- FindMarkers(MacroMono2, ident.1="Classical monos", ident.2 ="TAMS (APOE+)")
classVAPOE <- classVAPOE[which(classVAPOE$p_val_adj<0.05),]
