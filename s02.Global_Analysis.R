library(dplyr)
  library(Seurat)
  
  # Load the Human Gonod dataset
  pbmc.data <- Read10X(data.dir = "../data/filtered_gene_bc_matrices_mex/hg19/")
  # Initialize the Seurat object with the raw (non-normalized data).
  pbmc <- CreateSeuratObject(counts = pbmc.data, project = "humanGonod", min.cells = 3, min.features = 200)
  
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  
  # Visualize QC metrics as a violin plot
  VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  
  plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  CombinePlots(plots = list(plot1, plot2))
  
  pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  
  #pbmc[["RNA"]]@data
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(pbmc), 10)
  
  # plot variable features with and without labels
  plot1 <- VariableFeaturePlot(pbmc)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  CombinePlots(plots = list(plot1, plot2))
  
  pbmc <- ScaleData(pbmc)
  
  pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
  
  # Examine and visualize PCA results a few different ways
  print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
  
  VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
  DimPlot(pbmc, reduction = "pca",dims = c(1,3))
  
  PCA_Result_coord <- as.data.frame(pbmc@reductions$pca@cell.embeddings)
  PCA_Result_coord$Type <- gsub(".*-","",rownames(PCA_Result_coord))
  Type2Weeks <-c("6W","7W","7W","7W","8W","9W","10W","13W","15W","15W","16W","18W","19W","22W","22W","23W","23W")
  names(Type2Weeks) <-as.character(c(1:17))
  PCA_Result_coord$Week <- Type2Weeks[PCA_Result_coord$Type]
  Type2Sex <- c("Male","Female","Male","Female","Male","Male","Male","Female","Male","Male","Female","Female","Male","Female","Female","Male","Male")
  names(Type2Sex) <-as.character(c(1:17))
  PCA_Result_coord$Sex <- Type2Sex[PCA_Result_coord$Type]
  write.table(file="PCA_Result_coord.txt",PCA_Result_coord,quote = F,sep="\t")
  
  tiff("PCA_Result_week_PC1_2.tiff",antialias ="cleartype",height = 400,width = 500 )
  ggplot(data=PCA_Result_coord,aes(x=PC_1,y=PC_2,color=Week))+geom_point()+scale_color_d3(palette = "category20")+theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
  dev.off()
  
  tiff("PCA_Result_week_PC1_3.tiff",antialias ="cleartype",height = 400,width = 500 )
  ggplot(data=PCA_Result_coord,aes(x=PC_1,y=PC_3,color=Week))+geom_point()+scale_color_d3(palette = "category20")+theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
  dev.off()
  
  tiff("PCA_Result_week_PC2_3.tiff",antialias ="cleartype",height = 400,width = 500 )
  ggplot(data=PCA_Result_coord,aes(x=PC_2,y=PC_3,color=Week))+geom_point()+scale_color_d3(palette = "category20")+theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
  dev.off()
  
  tiff("PCA_Result_Sex_PC1_2.tiff",antialias ="cleartype",height = 400,width = 500 )
  ggplot(data=PCA_Result_coord,aes(x=PC_1,y=PC_2,color=Sex))+geom_point()+scale_color_d3(palette = "category20")+theme_bw()+theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
  dev.off()
  
  tiff("PCA_Result_Sex_PC1_3.tiff",antialias ="cleartype",height = 400,width = 500 )
  ggplot(data=PCA_Result_coord,aes(x=PC_1,y=PC_3,color=Sex))+geom_point()+scale_color_d3(palette = "category20")+theme_bw()+theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
  dev.off()
  
  tiff("PCA_Result_Sex_PC2_3.tiff",antialias ="cleartype",height = 400,width = 500 )
  ggplot(data=PCA_Result_coord,aes(x=PC_2,y=PC_3,color=Sex))+geom_point()+scale_color_d3(palette = "category20")+theme_bw()+theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
  dev.off()
  
  PCA_Result_loading <- as.data.frame(pbmc@reductions$pca@feature.loadings)
  pdf("Dim1_PCA_Result_heatmap.pdf",height = 6,width = 8)
  DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
  dev.off()
  # computation time
  pbmc <- JackStraw(pbmc, num.replicate = 20)
  pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
  
  JackStrawPlot(pbmc, dims = 1:15)
  ElbowPlot(pbmc)
  
  pbmc <- FindNeighbors(pbmc, dims = 1:10)
  pbmc <- FindClusters(pbmc, resolution = 0.5)
  head(Idents(pbmc), 5)
  # If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
  # 'umap-learn')
  pbmc <- RunUMAP(pbmc, dims = 1:10)
  save(file="pbmc.Robj",pbmc)
   # note that you can set `label = TRUE` or use the LabelClusters function to help label
  # individual clusters
  
  DimPlot(pbmc, reduction = "umap")
  
  umap_result_coord <- as.data.frame(pbmc@reductions$umap@cell.embeddings)
  umap_result_loading <- as.data.frame(pbmc@reductions$umap@feature.loadings)
  umap_result_coord<- cbind(umap_result_coord,PCA_Result_coord[rownames(umap_result_coord),c("Type","Week","Sex")])
  umap_result_coord$Cluster <- pbmc@meta.data[rownames(umap_result_coord),"RNA_snn_res.0.5"]
  write.table(file="umap_result_coord.txt",umap_result_coord,quote = F,sep = "\t")
  
  tiff("UMAP_Result_Sex_dim_1_10.tiff",antialias ="cleartype",height = 400,width = 500 )
  ggplot(data=umap_result_coord,aes(x=UMAP_1,y=UMAP_2,color=Sex))+geom_point(size=0.5)+scale_color_manual(values = c("deeppink","skyblue"))+theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
  dev.off()
  
  
  tiff("UMAP_Result_Weeks_dim_1_10.tiff",antialias ="cleartype",height = 400,width = 500 )
  ggplot(data=umap_result_coord,aes(x=UMAP_1,y=UMAP_2,color=Week))+geom_point(size=0.5)+scale_color_d3(palette = "category20")+theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
  dev.off()
  
 umap_result_coord$CellType <- "Soma"
  umap_result_coord[umap_result_coord$Cluster%in%c(5,16,18),"CellType"] <- "PGC"
  umap_result_coord[umap_result_coord$Cluster%in%c(14,11),"CellType"] <- "Immune"
  umap_result_coord[umap_result_coord$Cluster%in%c(15,21),"CellType"] <- "Blood"
  umap_result_coord[umap_result_coord$Cluster%in%c(17),"CellType"] <- "Endo"
  write.table(file="umap_result_coord_withCellType.txt",umap_result_coord,quote = F,sep="\t")
  
  tiff("UMAP_Result_CellType_dim_1_10.tiff",antialias ="cleartype",height = 400,width = 500 )
  ggplot(data=umap_result_coord,aes(x=UMAP_1,y=UMAP_2,color=CellType))+geom_point(size=0.5)+scale_color_manual(values = brewer.pal(8,"Set1")[c(1,2,5,4,3)])+theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
  dev.off()


  pdf("Embryo_CellTypeRatio_barplot.pdf",height = 6,width = 10)
  ggplot(data=Embryo_CellType_Ratio,aes(x=Embryo,y=Ratio,fill=CellType))+geom_bar(stat = "identity")+theme_bw()+scale_fill_locuszoom()+theme(axis.text.x = element_text(angle = 90,hjust = 1))
  dev.off()
  
  pdf("Male_Embryo_PGC_Ratio_dotplot.pdf",height = 6,width = 8)
  ggplot(data=Embryo_CellType_Ratio[Embryo_CellType_Ratio$Sex=="Male"&Embryo_CellType_Ratio$CellType=="PGC",],aes(x=Embryo,y=Ratio,group=1))+geom_point(color="skyblue",size=4)+geom_line(lty="longdash",color="black")+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust = 1))
  dev.off()
  
  pdf("Female_Embryo_PGC_Ratio_dotplot.pdf",height = 6,width = 4)
  ggplot(data=Embryo_CellType_Ratio[Embryo_CellType_Ratio$Sex=="Female"&Embryo_CellType_Ratio$CellType=="PGC",],aes(x=Embryo,y=Ratio,group=1))+geom_point(color="deeppink",size=4)+geom_line(lty="longdash",color="black")+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust = 1))
  dev.off() 
