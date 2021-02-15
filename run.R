library(Seurat)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggrepel)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DESeq2)
library(SCENIC)
library(SCopeLoomR)
library(AUCell)
library(SingleCellExperiment)
library(dplyr)
library(patchwork)




##### Functions required to create stacked violin plot ######

modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}



### setting the working directory
setwd("/home/diagnostics/Desktop/PhD_Projects/2020-02-24-Thomas_Pietchman")
### Creating required directories
dir.create("figures")
dir.create("tables")



#######################################################################################################
###################################### INTEGRATED DATA ################################################
#######################################################################################################

data <- Read10X(data.dir = "./data/D10_2_D11_2_D12_2/filtered_feature_bc_matrix")
data.hto <- as.sparse(data$`Antibody Capture`)
data.rna <- as.sparse(data$`Gene Expression`)[1:33694,]
data.vir <- as.sparse(data$`Gene Expression`)[33695:33706,]
rm(data)

### Creating a new seurat object with the new tables
data <- CreateSeuratObject(counts = data.rna, project = "data", min.cells = 3)
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000) ### increase the number of variable features
data <- ScaleData(data, features = rownames(data))


data[["HTO"]] <- CreateAssayObject(counts = data.hto)
data[["VIRUS"]] <- CreateAssayObject(counts = data.vir)

rm(data.hto, data.rna, data.vir)

### adding the sample name to meta data table
data@meta.data$Infection <- gsub("3", "ctrl" ,gsub("2", "bystdr",gsub("1", "inf", gsub(".*-", "", colnames(data@assays$RNA)))))
data@meta.data$Infection2 <- gsub("3", "mock" ,gsub("2", "bystander",gsub("1", "infected", gsub(".*-", "", colnames(data@assays$RNA)))))




data <- NormalizeData(data, assay = "HTO", normalization.method = "CLR")
data <- HTODemux(data, assay = "HTO", positive.quantile = 0.998)
table(data$HTO_classification.global)


Idents(data) <- "HTO_maxID"

a <- ggplot(data.frame(t(data@assays$HTO@data), data@meta.data), aes(x= Hashtag1.TotalA, y=HTO_maxID)) + 
  ggridges::geom_density_ridges2(aes(x = Hashtag1.TotalA, fill = HTO_maxID), alpha = 0.5) +
  theme_classic() +
  scale_fill_brewer(palette = "Set2") +
  xlab("Log Normalized Expression") +
  ylab("") +
  theme(legend.position = "none") +
  ggtitle("Hashtag1")

b <- ggplot(data.frame(t(data@assays$HTO@data), data@meta.data), aes(x= Hashtag2.TotalA, y=HTO_maxID)) + 
  ggridges::geom_density_ridges2(aes(x = Hashtag2.TotalA, fill = HTO_maxID), alpha = 0.5) +
  theme_classic() +
  scale_fill_brewer(palette = "Set2") +
  xlab("Log Normalized Expression") +
  ylab("") +
  theme(legend.position = "none") +
  ggtitle("Hashtag2")

c <- ggplot(data.frame(t(data@assays$HTO@data), data@meta.data), aes(x= Hashtag3.TotalA, y=HTO_maxID)) + 
  ggridges::geom_density_ridges2(aes(x = Hashtag3.TotalA, fill = HTO_maxID), alpha = 0.5) +
  theme_classic() +
  scale_fill_brewer(palette = "Set2") +
  xlab("Log Normalized Expression") +
  ylab("") +
  theme(legend.position = "none") +
  ggtitle("Hashtag3")

d <- ggplot(data.frame(t(data@assays$HTO@data), data@meta.data), aes(x= Hashtag4.TotalA, y=HTO_maxID)) + 
  ggridges::geom_density_ridges2(aes(x = Hashtag4.TotalA, fill = HTO_maxID), alpha = 0.5) +
  theme_classic() +
  scale_fill_brewer(palette = "Set2") +
  xlab("Log Normalized Expression") +
  ylab("") +
  theme(legend.position = "none") +
  ggtitle("Hashtag4")  

ggarrange(a,b,c,d, ncol = 1, nrow = 4)
ggsave("./figures/ridges_plot.svg", width = 5, height = 8)
rm(a,b,c,d)




Idents(data) <- "HTO_classification.global"
ggplot(data@meta.data, aes(x = HTO_classification.global, y = nCount_RNA, color = HTO_classification.global)) +
  geom_jitter(size = 0.1) +
  geom_violin(alpha = 0.3) +
  scale_color_brewer(palette = "Set1") +
  theme_classic() +
  xlab("") +
  ylab("nUMI (log transformed)") +
  scale_y_continuous(trans='log10') +
  theme(legend.position = "none")

ggsave("./figures/HTO_nUMI.svg")

ggplot(data@meta.data, aes(x = HTO_classification.global, y = nFeature_RNA, color = HTO_classification.global)) +
  geom_jitter(size = 0.1) +
  geom_violin(alpha = 0.3) +
  scale_color_brewer(palette = "Set1") +
  theme_classic() +
  xlab("") +
  ylab("nGene (log transformed)") +
  scale_y_continuous(trans='log10') +
  theme(legend.position = "none")

ggsave("./figures/HTO_nGene.svg")



meta2 <- data@meta.data
# we remove negative cells from the object
data <- subset(data, idents = "Negative", invert = TRUE)
data <- subset(data, idents = "Doublet", invert = TRUE)
###########################################
########### Normal QC Procedure ###########
###########################################
### Finding mt genes
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
a <- ggplot(data@meta.data, aes(x = factor(data@meta.data$Infection, levels = c("ctrl", "bystdr", "inf")), y = percent.mt, color = Infection)) +
  geom_violin() +
  geom_jitter(size = 0.1) +
  scale_color_manual(values = c('#999999','#E69F00', '#56B4E9')) +
  scale_x_discrete(labels=c("ctrl" = "mock", "bystdr" = "bystander", "inf" = "infected")) +
  theme_classic() +
  xlab("") +
  geom_hline(yintercept = 30, linetype="dashed") +
  ylab("Mitochondrial Percentage") +
  theme(legend.position = "none")


b <- ggplot(data@meta.data, aes(x = factor(data@meta.data$Infection, levels = c("ctrl", "bystdr", "inf")), y = nCount_RNA, color = Infection)) +
  geom_violin() +
  geom_jitter(size = 0.1) +
  scale_color_manual(values = c('#999999','#E69F00', '#56B4E9')) +
  scale_x_discrete(labels=c("ctrl" = "mock", "bystdr" = "bystander", "inf" = "infected")) +
  theme_classic() +
  xlab("") +
  geom_hline(yintercept = 5000, linetype="dashed") +
  geom_hline(yintercept = 100000, linetype="dashed") +
  ylab("nUMI (log scale)") +
  scale_y_continuous(trans='log10') +
  theme(legend.position = "none")


c <- ggplot(data@meta.data, aes(x = factor(data@meta.data$Infection, levels = c("ctrl", "bystdr", "inf")), y = nFeature_RNA, color = Infection)) +
  geom_violin() +
  geom_jitter(size = 0.1) +
  scale_color_manual(values = c('#999999','#E69F00', '#56B4E9')) +
  scale_x_discrete(labels=c("ctrl" = "mock", "bystdr" = "bystander", "inf" = "infected")) +
  theme_classic() +
  xlab("") +
  geom_hline(yintercept = 1500, linetype="dashed") +
  geom_hline(yintercept = 10000, linetype="dashed") +
  ylab("nGene (log scale)") +
  scale_y_continuous(trans='log10') +
  theme(legend.position = "none")

ggarrange(a,b,c, ncol = 1, nrow = 3)
ggsave("./figures/QC_violinPlot.svg", height = 13, width = 5)
rm(a,b,c)



a <- ggplot(data@meta.data[data@meta.data$Infection == "bystdr", ], aes(x = Infection, y = (nCount_VIRUS + 1), color = Infection)) +
  geom_violin() +
  geom_jitter(size = 0.1) +
  scale_color_manual(values = "#999999") +
  scale_x_discrete(labels=c("bystdr" = "bystander")) +
  theme_classic() +
  xlab("") +
  geom_hline(yintercept = 50, linetype="dashed") +
  ylab("") +
  scale_y_continuous(trans='log10') +
  theme(legend.position = "none")


b <- ggplot(data@meta.data[data@meta.data$Infection == "inf", ], aes(x = Infection, y = (nCount_VIRUS + 1), color = Infection)) +
  geom_violin() +
  geom_jitter(size = 0.1) +
  scale_color_manual(values = '#56B4E9') +
  scale_x_discrete(labels=c("inf" = "infected")) +
  theme_classic() +
  xlab("") +
  geom_hline(yintercept = 100, linetype="dashed") +
  ylab("") +
  scale_y_continuous(trans='log10') +
  theme(legend.position = "none")

c <- ggplot(data@meta.data[data@meta.data$Infection == "ctrl", ], aes(x = Infection, y = (nCount_VIRUS + 1), color = Infection)) +
  geom_violin() +
  geom_jitter(size = 0.1) +
  scale_color_manual(values = '#E69F00') +
  scale_x_discrete(labels=c("ctrl" = "mock")) +
  theme_classic() +
  xlab("") +
  #geom_hline(yintercept = 100, linetype="dashed") +
  ylab("Number of viral UMIs") +
  scale_y_continuous(trans='log10') +
  theme(legend.position = "none")

ggarrange(c,a,b, ncol = 3, nrow = 1)
ggsave("./figures/viral_umis_per_condition.svg", height = 5, width = 8)
rm(a,b,c)

ggplot(data@meta.data, aes(x = factor(data@meta.data$Infection, levels = c("ctrl", "bystdr", "inf")), y = (nCount_VIRUS + 1), color = Infection)) +
  geom_violin(scale = "width") +
  geom_jitter(size = 0.1) +
  scale_color_manual(values = c('#999999','#E69F00', '#56B4E9')) +
  theme_classic() +
  xlab("") +
  ylab("number of viral UMIs") +
  scale_y_continuous(trans='log10') +
  theme(legend.position = "none")

ggsave("./figures/viral_umis_per_condition_2.svg", height = 4, width = 7)

ggplot(data@meta.data, aes(x = factor(data@meta.data$Infection, levels = c("ctrl", "bystdr", "inf")), y = (nCount_VIRUS + 1), color = Infection)) +
  geom_violin(scale = "width") +
  geom_jitter(size = 0.1) +
  scale_color_manual(values = c('#999999','#E69F00', '#56B4E9')) +
  theme_classic() +
  xlab("") +
  geom_hline(yintercept = 50, linetype="dashed") +
  geom_hline(yintercept = 100, linetype="dashed") +
  ylab("number of viral UMIs") +
  scale_y_continuous(trans='log10') +
  theme(legend.position = "none")

ggsave("./figures/viral_umis_per_condition_3.svg", height = 4, width = 7)

### Selecting the QC thresholds
nGene_cells <- rownames(data@meta.data[data@meta.data$nFeature_RNA > 1500 & data@meta.data$nFeature_RNA < 10000, ])
nUMI_cells <- rownames(data@meta.data[data@meta.data$nCount_RNA > 5000 & data@meta.data$nCount_RNA < 100000, ])
mito_cells <- rownames(data@meta.data[data@meta.data$percent.mt < 30, ])
bystdr_cells <- rownames(data@meta.data[data@meta.data$Infection == "bystdr" & data@meta.data$nCount_VIRUS < 50, ])
inf_cells <- rownames(data@meta.data[data@meta.data$Infection == "inf" & data@meta.data$nCount_VIRUS > 100, ])
ctrl_cells <- rownames(data@meta.data[data@meta.data$Infection == "ctrl", ])
infection_cells <- c(inf_cells,ctrl_cells,bystdr_cells)
cells <- intersect(intersect(intersect(nGene_cells,nUMI_cells), mito_cells),infection_cells)

### Subsetting cells
data <- subset(data, cells = cells)

### dimension reduction
data <- RunPCA(data, features = VariableFeatures(object = data), npcs = 50)
DimPlot(data, reduction = "pca", dims = c(1, 2))

data <- FindNeighbors(data, dims = 1:35)

data <- FindClusters(data, resolution = 0.7)
data <- RunUMAP(data, dims = 1:35, umap.method = "umap-learn")

DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()



#########################################



### Showing the expression of Hashtag antibodies
df <- data.frame(t(data@assays$HTO@data), UMAP_1 = data@reductions$umap@cell.embeddings[,1], UMAP_2 = data@reductions$umap@cell.embeddings[,2])

Hashtag1 <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=Hashtag1.TotalA)) + 
  geom_jitter(size = 1) +
  scale_color_gradientn(colours = rev(brewer.pal(11, "Spectral")), limits=c(quantile(df$Hashtag1, probs = c(0.05, 0.95), names = FALSE)), oob = scales::squish) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  ggtitle("Hashtag1")

Hashtag2 <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=Hashtag2.TotalA)) + 
  geom_jitter(size = 1) +
  scale_color_gradientn(colours = rev(brewer.pal(11, "Spectral")), limits=c(quantile(df$Hashtag2, probs = c(0.05, 0.95), names = FALSE)), oob = scales::squish) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  ggtitle("Hashtag2")

Hashtag3<- ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=Hashtag3.TotalA)) + 
  geom_jitter(size = 1) +
  scale_color_gradientn(colours = rev(brewer.pal(11, "Spectral")), limits=c(quantile(df$Hashtag3, probs = c(0.05, 0.95), names = FALSE)), oob = scales::squish) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  ggtitle("Hashtag3")

Hashtag4 <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=Hashtag4.TotalA)) + 
  geom_jitter(size = 1) +
  scale_color_gradientn(colours = rev(brewer.pal(11, "Spectral")), limits=c(quantile(df$Hashtag4, probs = c(0.05, 0.95), names = FALSE)), oob = scales::squish) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  ggtitle("Hashtag4")

ggarrange(Hashtag1, Hashtag2, Hashtag3, Hashtag4,  ncol = 2, nrow = 2)
ggsave("./figures/12_Hashtag expression.png", height = 10, width = 14)
rm(Hashtag1,Hashtag2,Hashtag3, Hashtag4)


df <- cbind(df,HTO =  data$HTO_classification)

ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=HTO)) +
  geom_jitter() +
  scale_color_brewer(palette = "Set2") +
  theme_classic() +
  theme(legend.title = element_blank()) 

ggsave("./figures/HTO_classification.svg", width = 7, height = 6)

##########################################
##########################################

cluster_renamed <- gsub("0", "cilliated cells (bystander)",
                        gsub("9","cilliated cells (mock)",
                             gsub("8","intermediate",
                                  gsub("7","unassigned (A)",
                                       gsub("6","cilliated cells (mock)",
                                            gsub("5", "cilliated cells (infected)",
                                                 gsub("4", "basal cells",
                                                      gsub("3", "mixed goblet and club cells (B)",
                                                           gsub("2", "mixed goblet and club cells (B)",
                                                                gsub("1","mixed goblet and club cells (A)",
                                                                     gsub("16", "ionocyte",
                                                                          gsub("15","cilliated cells (bystander)",
                                                                               gsub("14", "cilliated cells two",
                                                                                    gsub("13", "cilliated cells (infected B)",
                                                                                         gsub("12","moucus cells",
                                                                                              gsub("11", "unassigned (B)",
                                                                                                   gsub("10", "proliferating basal cells",as.character(data@meta.data$seurat_clusters))))))))))))))))))


data$seurat_clusters <- factor(cluster_renamed)
Idents(data) <- cluster_renamed

df <- data.frame(data@reductions$umap@cell.embeddings, cluster = data$seurat_clusters)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
  geom_point(size = 2) +
  scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(14)) +
  theme_classic() 



DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 3) + 
  NoLegend() +
  scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(14))



ggsave("./figures/UMAP.svg", width = 6, height = 5)


ggplot(data@meta.data, aes(x = factor(data@meta.data$Infection2, levels = c("mock", "bystander", "infected")), y = (nCount_VIRUS + 1), color = Infection)) +
  geom_violin(scale = "width") +
  geom_jitter(size = 0.1) +
  scale_color_manual(values=c('#999999','#E69F00', '#56B4E9'))+
  theme_classic() +
  xlab("") +
  ylab("number of viral UMIs") +
  scale_y_continuous(trans='log10') +
  theme(legend.position = "none") 
  
  

ggsave("./figures/viral_umis_per_condition_4.svg", height = 4, width = 4)


#####################################
#####################################
#####################################
# DE expressionon of all clusters 
Idents(data) <- data$seurat_clusters
markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers2 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

Marker.Genes <- markers2$gene
Marker.Genes <- unique(Marker.Genes)
DotPlot(data, features = Marker.Genes) + RotatedAxis() +
  xlab("") +
  ylab("") +
  theme() +
  scale_color_viridis_c(direction = 1, option = "plasma") +
  theme(axis.text.y = element_text(face = "italic")) +
  coord_flip() 

ggsave("./figures/DotPlot2.svg", height = 30, width = 8)


features<- Marker.Genes
StackedVlnPlot(obj = data, features = features) 
ggsave("./figures/Vln2.svg", height = 49, width = 8)

## Comparision of mock, bystander and infected
data_cilliated <- subset(data, idents = c("cilliated cells (mock)", "cilliated cells (bystander)", "cilliated cells (infected)"))
Idents(data_cilliated) <- factor(Idents(data_cilliated), levels= c("cilliated cells (mock)", "cilliated cells (bystander)", "cilliated cells (infected)"))
data_cilliated$seurat_clusters <- factor(x = data_cilliated$seurat_clusters, levels = c("cilliated cells (mock)", "cilliated cells (bystander)", "cilliated cells (infected)"))
DimPlot(data_cilliated) 

markers <- FindAllMarkers(data_cilliated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file = "cilliated_marker.csv")
markers2 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
StackedVlnPlot(obj = data_cilliated, features = markers2$gene[1:20], cols = c('#E69F00','#999999', '#56B4E9')) 
ggsave("./figures/bystdr_marker.svg", height = 12, width = 3.5)

StackedVlnPlot(obj = data_cilliated, features = markers2$gene[21:40], cols = c('#E69F00','#999999', '#56B4E9')) 
ggsave("./figures/infected_marker.svg", height = 12, width = 3.5)

StackedVlnPlot(obj = data_cilliated, features = markers2$gene[41:60], cols = c('#E69F00','#999999', '#56B4E9')) 
ggsave("./figures/mock_marker.svg", height = 12, width = 3.5)


write.csv(markers, file = "cilliated_marker.csv")


############################################################################
## Here the goal is to find genes that are specific to infected condition ##
############################################################################

markers_infected <- markers[markers$cluster == "cilliated cells (infected)", ]
markers_infected <- markers_infected$gene
### we calculate the mean of each DE gene in each condition and then select the genes with 2 times higher in infected

cil_inf <- subset(data_cilliated, idents = c("cilliated cells (infected)"))
cil_bystdr <- subset(data_cilliated, idents = c( "cilliated cells (bystander)"))
cil_mock <- subset(data_cilliated, idents = c("cilliated cells (mock)"))

med_inf <- rowMeans(as.matrix(cil_inf@assays$RNA@data[markers_infected,]))
med_bystdr <- rowMeans(as.matrix(cil_bystdr@assays$RNA@data[markers_infected,]))
med_mock <- rowMeans(as.matrix(cil_mock@assays$RNA@data[markers_infected,]))

df_med <- data.frame(med_inf, med_bystdr, med_mock)
rownames(df_med) <- markers_infected

int <- df_med[df_med$med_inf > df_med$med_bystdr * 2 & df_med$med_inf > df_med$med_mock * 2,]
int2 <- data.frame(med_inf = as.character(round(int[,1],3)), med_bystdr = as.character(round(int[,2],3)), med_mock = as.character(round(int[,3],3)))
rownames(int2) <- rownames(int)
write.table(int2, file = "med_infected_specific.csv", sep = "\t")


StackedVlnPlot(obj = data_cilliated, features = rownames(int), cols = c('#E69F00','#999999', '#56B4E9')) 
ggsave("./figures/infected_specific.svg", height = 30, width = 3.5)



## GO term analysis aon the infected specific genes with clusterprofiler
genes <- rownames(int)
gene.df <- bitr(genes, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene         = gene.df$ENSEMBL,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05)



res <- ego@result
pos <- data.frame(res[1:40,], wave = rep("positive",40))

ggplot(pos, aes(x = factor(Description, levels = rev(pos$Description)), y = p.adjust, fill = wave)) +
  geom_bar(stat="identity") +
  scale_y_continuous(trans="log10") +
  scale_x_discrete(position = "top") +
  geom_hline(yintercept = 0.05, linetype="dashed") + 
  coord_flip() +
  theme_classic() +
  xlab("") +
  ylab("adjusted p-value") +
  scale_fill_brewer(palette="Dark2") +
  theme(legend.position = "none")

ggsave("./figures/infected_specific_GSE.svg", width = 8, height = 5)






## Here ISGs are defined as genes which are up-regulated in bystander condition
data_cilliated <- subset(data, idents = c("cilliated cells (mock)", "cilliated cells (bystander)"))
Idents(data_cilliated) <- factor(Idents(data_cilliated), levels= c("cilliated cells (mock)", "cilliated cells (bystander)"))
data_cilliated$seurat_clusters <- factor(x = data_cilliated$seurat_clusters, levels = c("cilliated cells (mock)", "cilliated cells (bystander)", "cilliated cells (infected)"))
DimPlot(data_cilliated) 

markers <- FindAllMarkers(data_cilliated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ISGs <- markers[markers$cluster == "cilliated cells (bystander)",]

write.csv(ISGs, file = "ISGs.csv")


#######################################
#######################################
###  nGene and nUMI projected on umap
df <- data.frame(data@reductions$umap@cell.embeddings, nUMI = data$nCount_RNA, nGene = data$nFeature_RNA)

ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = nUMI)) +
  geom_point(size = 1) +
  theme_classic() +
  scale_color_gradientn(colours = rev(brewer.pal(11, "Spectral"))) +
  theme(legend.title = element_blank(), legend.position = "none") 
ggsave("./figures/nUMI_umap.svg", width = 6, height = 5)

ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = nGene)) +
  geom_point(size = 1) +
  theme_classic() +
  scale_color_gradientn(colours = rev(brewer.pal(11, "Spectral"))) +
  theme(legend.title = element_blank(), legend.position = "none") 
ggsave("./figures/nGene_umap.svg", width = 6, height = 5)


### Dotplot
Marker.Genes <- c("FOXJ1", "TP73", "CCDC78", "MUC5B", "SPDEF", "MUC5AC", "SCGB3A1", "SCGB1A1", "FOXI1", "CFTR", "KRT5", "TP63", "NOTCH3")
DotPlot(data, features = Marker.Genes) + RotatedAxis() +
  xlab("") +
  ylab("") +
  theme(legend.position = "none") +
  scale_color_viridis_c(direction = -1, option = "plasma")

ggsave("./figures/DotPlot.svg", height = 6, width = 6)



## UMAP with infected, bystander and mock seperated by color 
df <- data.frame(data@reductions$umap@cell.embeddings, cluster = data$seurat_clusters)
df <- data.frame(df, infection = data@meta.data$Infection)
a <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, color=infection)) +
    geom_jitter(size = 0.8) +
    #scale_color_brewer(palette = "Dark2") +
  scale_color_manual(values=c('#999999','#E69F00', '#56B4E9'))+
  theme_classic() +
    theme(legend.position = "none") 
ggsave("./figures/Infection.svg", a, height = 5, width = 6)
  

### cell cycle assignment
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
  
############################################################################
############################################################################
############################################################################
  
  
  ### doing some measures on viral content of infected cells
data_infected <- subset(data, cells = rownames(data@meta.data[data@meta.data$Infection == "inf",]))
df <- as.data.frame(data_infected@assays$VIRUS@counts)
df <- df[,names(sort(colSums(df), decreasing = FALSE))]
vir.gene <- rownames(df)
df <- data.frame(lapply(df, function(x) x/sum(x)))
rownames(df) <- vir.gene
df <- reshape2::melt(as.matrix(df))
  
  
df$Var1 <- gsub("virus", "", df$Var1)
ggplot(df, aes(x = factor(Var1,levels = c("NS1", "NS2", "N", "P", "mGFP", "M", "SH", "G", "F", "M2-1", "M2-2", "L")), y = value)) +
    geom_jitter(size = 0.2, color = '#999999') +
    geom_boxplot( outlier.size=0.5, notch=TRUE, notchwidth = 0.2, width = 0.6, outlier.shape = NA) +
    #scale_y_continuous(trans='log10') +
    #scale_color_brewer(palette = "Paired") +
    theme_classic() +
    xlab("Genes ordered from 5'to 3' end") +
    ylab("Ratio") +
    theme(legend.position = "none")
ggsave("./figures/viral_genes_boxplot.svg", height = 4, width = 8)
  
### Showing the expression of canonical markers
dir.create("markers")
df <- data.frame(data@reductions$umap@cell.embeddings, cluster = data$seurat_clusters)
GENE <- c("CFTR", "FOXI1","FOXJ1", "TP73", "MUC5B","KRT5","TP63", "NOTCH3", "KRT14","MUC5B","SPDEF", "SCGB3A1", "SCGB1A1","DAPL1","MUC5AC")
  
  for (i in GENE) {
    if (i %in% rownames(data@assays$RNA@data)) {
      df$gene <- data@assays$RNA@data[i, ]
      ggplot(df, aes(x=UMAP_1, y=UMAP_2, color = gene)) +
        geom_jitter(size = 1) + 
        scale_color_viridis(discrete=FALSE, option = "plasma") +
        theme_classic() + 
        ggtitle(i) +
        theme(legend.title = element_blank())
      
      
      ggsave(paste("./markers/", i,".svg", sep = ""), width = 6, height = 5)
      
    } else {
      print("Gene not found")
    }
  }
  
  
  ### Showing infection load on the umap
  data <- NormalizeData(data, assay = "VIRUS", normalization.method = "CLR")
  colSums(data@assays$VIRUS@data)
  Infection_load <- exp(data@assays$VIRUS@data)
  Infection_load <- log(colSums(Infection_load))
  
  data$Infection_load <- Infection_load#data.frame(load = colSums(as.data.frame(data@assays$VIRUS@data)))
  
  ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = data$Infection_load)) +
    geom_point(size = 1) +
    theme_classic() +
    scale_color_gradientn(colours = rev(brewer.pal(11, "Spectral"))) +
      theme(legend.title = element_blank(), legend.position = "none") 
  ggsave("./figures/viral_load_integrated.svg", width = 6, height = 5)
  
  
  
## The focus of the rest of the analysis is on infected data
#######################################################################################################
###################################### INFECTED DATA ##################################################
#######################################################################################################
meta.data <- read.csv("/home/diagnostics/Desktop/PhD_Projects/2020-02-24-Thomas_Pietchman/D10_2_D11_2_D12_2_analysis/meta.data.csv", row.names = 1)
pos <- grep("-1", rownames(meta.data))
meta.data <- meta.data[pos,]
rownames(meta.data) <- gsub("-1","",rownames(meta.data))
ISG <- read.csv("ISG.csv")
ISG <- as.character(ISG$x)
  
### loading the data 
data <- Read10X(data.dir = "./data/D10_2/filtered_feature_bc_matrix")
data.hto <- as.sparse(data$`Antibody Capture`)
data.hto <- data.hto[,rownames(meta.data)]
data.rna <- as.sparse(data$`Gene Expression`)[1:33694,]
data.rna <- data.rna[,rownames(meta.data)]
data.vir <- as.sparse(data$`Gene Expression`)[33695:33706,]
data.vir <- data.vir[,rownames(meta.data)]
rm(data)
  
### Creating a new seurat object with the new tables
data <- CreateSeuratObject(counts = data.rna, project = "data", min.cells = 3)
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 1000) ### increase the number of variable features
data <- ScaleData(data, features = rownames(data))
  
  
  
data[["HTO"]] <- CreateAssayObject(counts = data.hto)
data[["VIRUS"]] <- CreateAssayObject(counts = data.vir)
rm(data.hto, data.rna, data.vir, pos)
  
  
data <- NormalizeData(data, assay = "HTO", normalization.method = "CLR")
data <- NormalizeData(data, assay = "VIRUS", normalization.method = "CLR")
  
### dimension reduction
data <- RunPCA(data, features = VariableFeatures(object = data), npcs = 50)
DimPlot(data, reduction = "pca", dims = c(1, 2))
  
data <- FindNeighbors(data, dims = 1:25)
data <- FindClusters(data, resolution = 0.2)
  

cluster_renamed <- gsub("4","cilliated cell (C)",
                        gsub("3", "cilliated cells (B)",
                        gsub("2", "mixed club and goblet cells",
                             gsub("1","unassigned",
                                  gsub("0", "cilliated cells (A)",as.character(data@meta.data$seurat_clusters))))))
  
### Creating umap
data <- RunUMAP(data, dims = 1:25, umap.method = "umap-learn")
coord <- data@reductions$umap@cell.embeddings
df <- data.frame(coord, cluster = cluster_renamed)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
    geom_point(size = 2) +
    theme_classic() +
    scale_color_brewer(palette = "Dark2") +
    theme(legend.title = element_blank()) +
    theme(legend.position = "none")
ggsave("./figures/umap_infected.svg", width = 7, height = 5.5)
  
#save(data, file = "infected.RData")
#### creating the violin plot that shows the expression of marker genes
  
Marker.Genes <- c("FOXJ1", "TP73", "CCDC78", "MUC5B", "SPDEF", "MUC5AC", "SCGB3A1", "SCGB1A1", "ZNF467")
Marker.Genes <- t(data@assays$RNA@data[Marker.Genes,])
Marker.Genes <- data.frame(cell = rownames(data@meta.data),
                     cluster = cluster_renamed,
                     Marker.Genes)
 
FOXJ1 <- ggplot(Marker.Genes, aes(x=cluster, y=FOXJ1, fill = cluster)) +
    geom_violin(alpha = 1, lwd = 0.7, scale = "width") + 
    scale_fill_brewer(palette = "Dark2") +
    theme_classic() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) 
TP73 <- ggplot(Marker.Genes, aes(x=cluster, y=TP73, fill = cluster)) +
    geom_violin(alpha = 1, lwd = 0.7, scale = "width") + 
    scale_fill_brewer(palette = "Dark2") +
    theme_classic() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
CCDC78 <- ggplot(Marker.Genes, aes(x=cluster, y=CCDC78, fill = cluster)) +
    geom_violin(alpha = 1, lwd = 0.7, scale = "width") + 
    scale_fill_brewer(palette = "Dark2") +
    theme_classic() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
MUC5B <- ggplot(Marker.Genes, aes(x=cluster, y=MUC5B, fill = cluster)) +
    geom_violin(alpha = 1, lwd = 0.7, scale = "width") + 
    scale_fill_brewer(palette = "Dark2") +
    theme_classic() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
SPDEF <- ggplot(Marker.Genes, aes(x=cluster, y=SPDEF, fill = cluster)) +
    geom_violin(alpha = 1, lwd = 0.7, scale = "width") + 
    scale_fill_brewer(palette = "Dark2") + 
    theme_classic() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
MUC5AC <- ggplot(Marker.Genes, aes(x=cluster, y=MUC5AC, fill = cluster)) +
    geom_violin(alpha = 1, lwd = 0.7, scale = "width") + 
    scale_fill_brewer(palette = "Dark2") +
    theme_classic() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
SCGB3A1 <- ggplot(Marker.Genes, aes(x=cluster, y=SCGB3A1, fill = cluster)) +
    geom_violin(alpha = 1, lwd = 0.7, scale = "width") + 
    scale_fill_brewer(palette = "Dark2") +
    theme_classic() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
SCGB1A1 <- ggplot(Marker.Genes, aes(x=cluster, y=SCGB1A1, fill = cluster)) +
    geom_violin(alpha = 1, lwd = 0.7, scale = "width") + 
    scale_fill_brewer(palette = "Dark2") +
    theme_classic() +
    xlab("")+
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size=0))

ggarrange(FOXJ1, TP73, CCDC78, MUC5B, SPDEF, SCGB3A1, SCGB1A1,
            ncol = 1, nrow = 7) 
  
ggsave("./figures/marker_infected.svg", width = 3, height = 8)
  
 
  
### ZNF467  
ZNF467 <- ggplot(Marker.Genes, aes(x=cluster, y=ZNF467, fill = cluster)) +
    geom_violin(alpha = 1, lwd = 0.7, scale = "width") + 
    scale_fill_brewer(palette = "Dark2") +
    theme_classic() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) 
  
ggsave("./figures/ZNF_vln.svg", width = 5, height = 3)
  
  
  
  
###################################################################
##################################################################
##################################################################
### doing some measures on viral content of infected cells
df <- as.data.frame(data@assays$VIRUS@counts)
df <- df[,names(sort(colSums(df), decreasing = FALSE))]
vir.gene <- rownames(df)
df <- data.frame(lapply(df, function(x) x/sum(x)))
rownames(df) <- vir.gene
df <- reshape2::melt(as.matrix(df))
  
  
df$Var1 <- gsub("virus", "", df$Var1)
ggplot(df, aes(x = factor(Var1,levels = c("NS1", "NS2", "N", "P", "mGFP", "M", "SH", "G", "F", "M2-1", "M2-2", "L")), y = value)) +
    geom_jitter(size = 0.2, color = '#999999') +
    geom_boxplot( outlier.size=0.5, notch=TRUE, notchwidth = 0.2, width = 0.6, outlier.shape = NA) +
    theme_classic() +
    xlab("Genes ordered from 5' to 3' end") +
    ylab("ratio") +
    theme(legend.position = "none")
ggsave("./figures/viral_genes_boxplot2.svg", height = 4, width = 8)



 
#################################################################
############### Assigning patients to each cell #################
#################################################################
### the tables required are available in the clone of the working directory I gave you
  
pers <- read.table("./D10_2_output/clusters.tsv", header = TRUE, row.names = 1)
rownames(pers) <- gsub("-1","",rownames(pers))
pers <- pers[colnames(data),]
table(pers$status)
  
df <- data.frame(coord)
  
  
  
df <- data.frame(df, patient = pers$assignment, status = pers$status)
pos <- grep("doublet", pers$status)
cluster <- as.character(df$patient)
for (i in pos) {
    cluster[i] <- "doublet"
}
df$patient <- cluster
 
 
patient <- cluster
  
ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = patient)) +
    geom_point(size = 2) +
    theme_classic() +
    scale_color_brewer(palette = "Accent") +
    theme(legend.title = element_blank()) +
    ggtitle("patient")
ggsave("./figures/umap_infected_patient.svg", width = 7, height = 5.5)
  
  
#################################################################
#################################################################
#################################################################
### adding meta data information from previous analysis
### Here meta data of integrated data is stored in a csv file and is being used to assign meta data information of the infected data
data$HTO_classification <- meta.data$HTO_classification
data$previous_cluster <- meta.data$RNA_snn_res.0.5
### adding infection load
data$Infection_load <- data.frame(load = log(colSums(exp(data@assays$VIRUS@data))))
  
  
  

cell_order <- names(sort(data$Infection_load, decreasing = FALSE))
  
  
#### distribution of infection in each subset
df <- data.frame(cluster = cluster_renamed, Infection_load = data$Infection_load, day = data$HTO_classification, patient = patient)
ggplot(df, aes(x = cluster, y= Infection_load, color = cluster)) +
    geom_jitter(size = 2) +
    geom_boxplot(color = "black", notch=TRUE, notchwidth = 0.2, width = 0.3, outlier.shape = NA, alpha = 0.5) +
    #geom_violin() +
    scale_y_continuous(trans='log10') +
    theme_classic() +
    ylab("log10(infection load)") +
    xlab("")+
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12))
    #theme(axis.text.x = element_text(angle = 90))
ggsave("./figures/infection_load_distribution2.svg", width = 6, height = 6) 
  
  
ggplot(df, aes(x = cluster, y= Infection_load, color = day)) +
    geom_jitter(size = 2) +
    geom_boxplot(color = "black", notch=TRUE, notchwidth = 0.2, width = 0.3, outlier.shape = NA, alpha = 0.5) +
    #geom_violin() +
    scale_y_continuous(trans='log10') +
    theme_classic() +
    ylab("log10(infection load)") +
    xlab("")+
    scale_color_brewer(palette = "Set2") +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12))
ggsave("./figures/infection_load_distribution_time2.svg", width = 6, height = 6) 
  
  
ggplot(df, aes(x = cluster, y= Infection_load, color = patient)) +
    geom_jitter(size = 2) +
    geom_boxplot(color = "black", notch=TRUE, notchwidth = 0.2, width = 0.3, outlier.shape = NA, alpha = 0.5) +
    #geom_violin() +
    scale_y_continuous(trans='log10') +
    theme_classic() +
    ylab("log10(infection load)") +
    xlab("")+
    scale_color_brewer(palette = "Accent") +
    theme() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12))
ggsave("./figures/infection_load_distribution_patient2.svg", width = 6, height = 6)
  
### Showing infection load on the umap
df <- data.frame(coord,  infection_load = data$Infection_load)
ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = infection_load)) +
    geom_point(size = 2) +
    theme_classic() +
    scale_color_gradientn(colours = rev(brewer.pal(11, "Spectral"))) +
    theme(legend.title = element_blank()) +
    theme(legend.title = element_blank(), legend.position = c(0.9,0.5))
ggsave("./figures/viral_load_infected.svg", width = 7, height = 5.5)
  

  
ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = data@assays$RNA@data["ZNF467",])) +
    geom_point(size = 2) +
    theme_classic() +
    scale_color_gradientn(colours = rev(brewer.pal(11, "Spectral"))) +
    theme(legend.title = element_blank()) +
    theme(legend.title = element_blank(), legend.position = c(0.9,0.5))
ggsave("./figures/ZNF467_scatter.svg", width = 7, height = 5.5)
  
    
### cell cycle assignment
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
ggplot(df, aes(x=UMAP_1, y=UMAP_2, color = data$Phase)) +
    geom_jitter(size = 2) + 
    scale_color_brewer(palette = "Set1") +
    theme_classic() + 
    theme(legend.title = element_blank(), legend.position = c(0.9,0.5))
  
ggsave("./figures/cell_cycle_infected.svg", width = 7, height = 5.5)
  
### Focusing on the analysis of cluster zero (0)
cells <- rownames(data@meta.data[data$seurat_clusters == 0,])
data_subset <- subset(data, cells = cells)
  
  
#####################################################
#####################################################
##################  Interferons #######################
#####################################################
  
  Receptors <- c("IFNB1", "IFNL1", "IFNL2", "IFNL3")
  data_receptor <- data.frame(t(data@assays$RNA@data[Receptors, ]))
  data_receptor <- data.frame(data_receptor, infection = data$Infection_load)
  
  a10 <- ggplot(data_receptor, aes(x = infection, y = IFNB1)) +
    geom_jitter() +
    theme_classic() +
    xlab("Infection load") +
    ggtitle(paste("R = ", as.character(cor(data_receptor[,"IFNB1"], data_receptor[,"infection"]))))
  
  a20 <- ggplot(data_receptor, aes(x = infection, y = IFNL1)) +
    geom_jitter() +
    theme_classic() +
    xlab("Infection load") +
    ggtitle(paste("R = ", as.character(cor(data_receptor[,"IFNL1"], data_receptor[,"infection"]))))
  
  a30 <- ggplot(data_receptor, aes(x = infection, y = IFNL2)) +
    geom_jitter() +
    theme_classic() +
    xlab("Infection load") +
    ggtitle(paste("R = ", as.character(cor(data_receptor[,"IFNL2"], data_receptor[,"infection"]))))
  
  a40 <- ggplot(data_receptor, aes(x = infection, y = IFNL3)) +
    geom_jitter() +
    theme_classic() +
    xlab("Infection load") +
    ggtitle(paste("R = ", as.character(cor(data_receptor[,"IFNL3"], data_receptor[,"infection"]))))
  
  ggarrange(a10,a20,a30,a40, ncol = 4, nrow = 1)
  ggsave("./figures/Receptor.svg", width = 12, height = 3.5)

  
  ### Finding the genes that are correlated with viral load
  
  correlation <- data.frame(cor = apply(data_subset@assays$RNA@data, 1, cor , y = data_subset$Infection_load), ident = rep("ident",21624), label = rownames(data_subset@assays$RNA@data))
  correlation <- correlation[complete.cases(correlation),]
  
  color <- c(1:21349)
  for (i in 1:21349) {
    if (correlation$cor[i] > 0.4 | correlation$cor[i] < -0.4) {
      color[i] <- "strong"
    } else {
      color[i] <- "weak"
    }
  }
  
  correlation <- data.frame(correlation, condition = color)
  
  ggplot(correlation, aes(x=ident, y=cor)) +
    geom_jitter(aes(x=ident, y=cor, color = condition),size = 1, position = position_jitter(seed = 1)) +
    geom_boxplot(outlier.shape = NA) +
    theme_classic() +
    geom_hline(yintercept = 0.4,  linetype="dashed") +
    geom_hline(yintercept = -0.4,  linetype="dashed") +
    xlab("") +
      scale_color_manual(values=c('#F8766D', '#999999')) +
    ylab("coefiicient of correlation") +
    coord_flip() +
    theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  ggsave("./figures/correlation_boxplot3.svg", width = 11, height = 2)
  
  
  
  df2 <- data.frame(data_subset@reductions$umap@cell.embeddings, gene = data_subset@assays$RNA@data["ZNF467",])
  ggplot(df2, aes(x=UMAP_1, y = UMAP_2, color = gene)) +
    geom_jitter(size=2) +
    theme_classic() +
    scale_color_gradientn(colours = rev(brewer.pal(11, "Spectral")))
  
  ggsave("./figures/ZNF476.svg", width = 7, height = 5.5)
  #################################################
  
  
  ### creating a heatmap for highly correlated genes
  cell_order <- names(sort(data_subset$Infection_load, decreasing = FALSE))
  Infection_load_order <- as.numeric(sort(data_subset$Infection_load, decreasing = FALSE))
  time <- data.frame(data$HTO_classification)[cell_order, ]
  
  hm_gene <- correlation[correlation$cor > 0.4 | correlation$cor < -0.4, ]
  norm <- data_subset@assays$RNA@data[rownames(hm_gene), cell_order]
  col <- colnames(norm)
  scaled <- t(apply(norm, 1, scale))
  colnames(scaled) <- col
  
  col_fun = colorRamp2(c(min(data_subset$Infection_load), max(data_subset$Infection_load)), c("white", "red"))
  ha = HeatmapAnnotation(Infection_load = Infection_load_order ,simple_anno_size = unit(1, "cm"), show_legend = FALSE, col = list(Infection_load = col_fun, time = c("Hashtag1-TotalA" = "#e41a1c", "Hashtag2-TotalA" = "#377eb8", "Hashtag3-TotalA" = "#4daf4a", "Hashtag4-TotalA" = "#984ea3")))
  svg("./figures/correlated_gene_heatmap.svg", width = 11, height = 20)
  p <- Heatmap(scaled, cluster_columns = FALSE, cluster_rows = TRUE,top_annotation = ha, col = colorRamp2(seq(-2,2,0.1), viridis(41)), show_column_names = FALSE, show_heatmap_legend = FALSE)
  print(p)
  dev.off()
  

  ### GO term analysis of positive and negative correlated genes
  positive <- rownames(correlation[correlation$cor > 0.3,])
  negative <- rownames(correlation[correlation$cor < -0.3,])
  gene.df <- bitr(positive, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
  ego_positive <- enrichGO(gene         = gene.df$ENSEMBL,
                               OrgDb         = org.Hs.eg.db,
                               keyType       = 'ENSEMBL',
                               ont           = "BP",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.01,
                               qvalueCutoff  = 0.05)
  
  write.csv(ego_positive, file = "./tables/GSE_positive.csv")
  
  
  gene.df <- bitr(negative, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
  ego_negative <- enrichGO(gene         = gene.df$ENSEMBL,
                           OrgDb         = org.Hs.eg.db,
                           keyType       = 'ENSEMBL',
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.05)
  
  write.csv(ego_negative, file = "./tables/GSE_negative.csv")
  
  pos <- data.frame(ego_positive[1:25,], wave = rep("positive",25))
  neg <- data.frame(ego_negative[1:25,], wave = rep("negative",25))
  #neg$Count <- neg$Count * -1
  neg$p.adjust <- log10(neg$p.adjust)
  pos$p.adjust <- log10(pos$p.adjust) * -1
  df <- rbind(pos, neg)
  
  df <- df[complete.cases(df),]
  df$row <- make.unique(paste(df$Description, ":",df$ID))
  df$row <- factor(df$row, levels = rev(df$row))
  
  ggplot(df, aes(x = row, y = p.adjust, fill = wave)) +
    geom_bar(stat="identity") +
    coord_flip() +
    theme_classic() +
    geom_vline(xintercept = 0.5) + 
    xlab("") +
    ylab("log10(p.adjust)") +
    scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
    theme(legend.position = "none")
  
  ggsave("./figures/GSE_pos_neg.svg", width = 10, height = 7) 
  
  
  ########################################
  ######### Pseudobulk analysis ##########
  ########################################
  ### Finding highly variable genes in cluster 0
  data_subset <- FindVariableFeatures(data_subset, selection.method = "vst", nfeatures = 5000)
  var_gene <-VariableFeatures(data_subset)



  ### Pooling and normalizing counts
  count <- data_subset@assays$RNA@counts
  count_order <- count[,cell_order]
  
  pooled <- data.frame(matrix(0, nrow=21624, ncol=20))
  
  for (i in 1:14) {
    columns <- c(1:50) + 50*(i-1)
    print(columns)
    pooled[,i] <- rowSums(as.matrix(count_order[,columns]))
  }
  
  for (i in 15:20) {
    columns <- c(1:51) + (51*(i-1))-14
    print(columns)
    pooled[,i] <- rowSums(as.matrix(count_order[,columns]))
  }
  
  rownames(pooled) <- rownames(count_order)
  
  ### Normalizing data
  sf <- estimateSizeFactorsForMatrix(pooled)
  norm <- t(t(pooled) / sf)
  log_norm <- log(norm + 1)
  log_norm <- log_norm[var_gene,] 
  
  
  res.pca <- prcomp(t(log_norm), scale = FALSE)
  ggplot(data.frame(PC1 = res.pca$x[,"PC1"], PC2 = res.pca$x[,"PC2"], label = colnames(log_norm), level = 1:20), aes(x = PC1, y = PC2)) +
    geom_jitter(size = 7, color = '#56B4E9') +
    theme_classic() +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_vline(xintercept=0, linetype="dashed") +
    geom_text(aes(label = level), color = "black")
  ggsave("./figures/pca_bulk.svg", height = 5, width = 6)
  
  #### Running differential expression analysis between different groups with DESeq2
  condition <- data.frame(time = colnames(pooled), treatment = c("end", "end", "first", "first", "first", "first", "second", "second", "second", "second", "third", "third", "fourth", "fourth", "fourth", "fifth", "fifth", "fifth", "end", "end"), stringsAsFactors = FALSE)
  rownames(condition) <- condition$time
  pooled_de <- pooled[var_gene,]
  
  dds <- DESeqDataSetFromMatrix(countData = pooled_de,
                                colData = condition,
                                design= ~treatment)
  dds <- DESeq(dds)
  resultsNames(dds)
  
  res_first <- results(dds, name = "treatment_first_vs_end")
  res_first <- data.frame(res_first)
  res_first <- res_first[complete.cases(res_first),]
  
  res_second <- results(dds, name = "treatment_second_vs_end")
  res_second <- data.frame(res_second)
  res_second <- res_second[complete.cases(res_second),]
  
  res_third <- results(dds, name = "treatment_third_vs_end")
  res_third <- data.frame(res_third)
  res_third <- res_third[complete.cases(res_third),]
  
  res_fourth <- results(dds, name = "treatment_fourth_vs_end")
  res_fourth <- data.frame(res_fourth)
  res_fourth <- res_fourth[complete.cases(res_fourth),]
  
  res_fifth <- results(dds, name = "treatment_fifth_vs_end")
  res_fifth <- data.frame(res_fifth)
  res_fifth <- res_fifth[complete.cases(res_fifth),]
  
  

  ### upregulated and down regulated genes 
  res_first_hm <- res_first[res_first$pvalue < 0.05,]
  res_first_hm_pos <- res_first_hm[res_first_hm$log2FoldChange > 0,]
  res_first_hm_neg <- res_first_hm[res_first_hm$log2FoldChange < 0,]
  
  res_second_hm <- res_second[res_second$pvalue < 0.05,]
  res_second_hm_pos <- res_second_hm[res_second_hm$log2FoldChange > 0,]
  res_second_hm_neg <- res_second_hm[res_second_hm$log2FoldChange < 0,]
  
  res_third_hm <- res_third[res_third$pvalue < 0.05,]
  res_third_hm_pos <- res_third_hm[res_third_hm$log2FoldChange > 0,]
  res_third_hm_neg <- res_third_hm[res_third_hm$log2FoldChange < 0,]
  
  res_fourth_hm <- res_fourth[res_fourth$pvalue < 0.05,]
  res_fourth_hm_pos <- res_fourth_hm[res_fourth_hm$log2FoldChange > 0,]
  res_fourth_hm_neg <- res_fourth_hm[res_fourth_hm$log2FoldChange < 0,]
  
  res_fifth_hm <- res_fifth[res_fifth$pvalue < 0.05,]
  res_fifth_hm_pos <- res_fifth_hm[res_fifth_hm$log2FoldChange > 0,]
  res_fifth_hm_neg <- res_fifth_hm[res_fifth_hm$log2FoldChange < 0,]
  
  hm_gene_pos <- unique(c(rownames(res_first_hm_pos), rownames(res_second_hm_pos), rownames(res_third_hm_pos), rownames(res_fourth_hm_pos), rownames(res_fifth_hm_pos)))
  hm_gene_neg <- unique(c(rownames(res_first_hm_neg), rownames(res_second_hm_neg), rownames(res_third_hm_neg), rownames(res_fourth_hm_neg), rownames(res_fifth_hm_neg)))
  hm_gene <- unique(c(hm_gene_pos, hm_gene_neg))
  
  remove_gene_neg <- intersect(hm_gene_pos, hm_gene_neg)
  hm_gene_pos <- hm_gene_pos[!hm_gene_pos %in% remove_gene_neg]
  
  ##### heatmap of positive genes

  hetmap_df <- log_norm[hm_gene_pos,]
  col <- colnames(hetmap_df)
  scaled <- t(apply(hetmap_df, 1, scale))
  colnames(scaled) <- col
  

  library(pheatmap)
  set.seed(123)
  out <- pheatmap(scaled,cluster_cols = FALSE, cutree_rows = 3, clustering_method = "ward.D2", filename = "./figures/pheatmap.pdf", width = 15, height = 15, show_rownames = FALSE, show_colnames = FALSE, treeheight_row = 0, treeheight_col = 0, border_color = "NA")
  # set.seed(123)
  # out <- pheatmap(scaled, cluster_cols = FALSE, cutree_rows = 3, clustering_method = "ward.D2", filename = "./figures/pheatmap2.pdf", width = 15, height = 15, show_rownames = TRUE, show_colnames = FALSE, treeheight_row = 0, treeheight_col = 0)
  # set.seed(123)
  # out <- pheatmap(scaled, cluster_cols = FALSE, cutree_rows = 3, clustering_method = "ward.D2", show_rownames = FALSE, show_colnames = FALSE, treeheight_row = 0, treeheight_col = 0)
  # 

  
  
  #### finding gene in each cluster
  set.seed(123)
  clu_heatmap <- sort(cutree(out$tree_row, k=3))
  
  late_intermediate <- names(clu_heatmap[clu_heatmap == 3])
  early_intermediate <- names(clu_heatmap[clu_heatmap == 1])
  intermediate <- names(clu_heatmap[clu_heatmap == 2])


  write.csv(early_intermediate, file = "early_intermediate.csv")
  write.csv(intermediate, file = "intermediate.csv")
  write.csv(late, file = "late_intermediate.csv")

  
  
  ### heatmap of negative genes
  hetmap_df <- log_norm[hm_gene_neg,]
  col <- colnames(hetmap_df)
  scaled <- t(apply(hetmap_df, 1, scale))
  colnames(scaled) <- col

  library(pheatmap)
  set.seed(123)
  out <- pheatmap(scaled,cluster_cols = FALSE, cutree_rows = 3, clustering_method = "ward.D2", filename = "./figures/pheatmap_pos.pdf", width = 15, height = 30, show_rownames = FALSE, show_colnames = FALSE, treeheight_row = 0, treeheight_col = 0, border_color = "NA")
  # set.seed(123)
  # out <- pheatmap(scaled, cluster_cols = FALSE, cutree_rows = 3, clustering_method = "ward.D2", filename = "./figures/pheatmap2_pos.pdf", width = 15, height = 30, show_rownames = TRUE, show_colnames = FALSE, treeheight_row = 0, treeheight_col = 0)
  # set.seed(123)
  # out <- pheatmap(scaled, cluster_cols = FALSE, cutree_rows = 3, clustering_method = "ward.D2", show_rownames = FALSE, show_colnames = FALSE, treeheight_row = 0, treeheight_col = 0)
  # 

  
  
  
  set.seed(123)
  clu_heatmap <- sort(cutree(out$tree_row, k=3))
  
  early <- names(clu_heatmap[clu_heatmap == 3])
  late <- names(clu_heatmap[clu_heatmap == 1])
  up_down_up <- names(clu_heatmap[clu_heatmap == 2])

  write.csv(early, file = "early.csv")
  write.csv(up_down_up, file = "up_down_up.csv")
  write.csv(late, file = "late.csv")
  
  ### positive and negative heatmaps are attached afterwards via an image editing program
  
  
  #######################################################################
  ### Gene Ontology Enrichment analysis of each wave of transcription ###
  #######################################################################
  gene.df <- bitr(early, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
  ego_early <- enrichGO(gene         = gene.df$ENSEMBL,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = 'ENSEMBL',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
  
  write.csv(ego_early, file = "./tables/GSE_early.csv")
  
  gene.df <- bitr(early_intermediate, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
  ego_early_intermediate <- enrichGO(gene         = gene.df$ENSEMBL,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05)
  
  write.csv(ego_early_intermediate, file = "./tables/GSE_early_intermediate.csv")
  
  gene.df <- bitr(intermediate, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
  ego_intermediate <- enrichGO(gene         = gene.df$ENSEMBL,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05)
  
  write.csv(ego_intermediate, file = "./tables/GSE_intermediate.csv")
  
  gene.df <- bitr(late_intermediate, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
  ego_late_intermediate <- enrichGO(gene         = gene.df$ENSEMBL,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.1,
                        qvalueCutoff  = 0.1)
  
  write.csv(ego_late_intermediate, file = "./tables/GSE_late_intermediate.csv")
  
  gene.df <- bitr(late, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
  ego_late <- enrichGO(gene         = gene.df$ENSEMBL,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05)
  
  write.csv(ego_late, file = "./tables/GSE_late.csv")
  
  gene.df <- bitr(up_down_up, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
  ego_up_down_up <- enrichGO(gene         = gene.df$ENSEMBL,
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05)
  
  write.csv(ego_up_down_up, file = "./tables/GSE_up_down_up.csv")
  
  df <- rbind(data.frame(ego_early[1:10,], wave = rep("early",10)),
              data.frame(ego_early_intermediate[1:10,], wave = rep("early_intermediate",10)),
              data.frame(ego_intermediate[1:10,], wave = rep("intermediate",10)),
              data.frame(ego_late_intermediate[1:10,], wave = rep("late_intermediate",10)),
              data.frame(ego_late[1:10,], wave = rep("late",10)),
              data.frame(ego_up_down_up[1:10,], wave = rep("up_down_up",10)))
  
  df <- df[complete.cases(df),]
  df$row <- make.unique(paste(df$ID, ":",df$Description))
  df$row <- factor(df$row, levels = rev(df$row))
  df$log10_pvalue <- log10(df$p.adjust)
  
  ggplot(df, aes(x = row, y = p.adjust, fill = wave)) +
    geom_bar(stat="identity") +
    scale_y_continuous(trans="log10") +
    scale_x_discrete(position = "top") +
    #scale_x_reverse() +
    geom_hline(yintercept = 0.05, linetype="dashed") + 
    coord_flip() +
    theme_classic() +
    xlab("") +
    ylab("adjusted p-value") +
    scale_fill_brewer(palette="Dark2") +
    theme(legend.position = "none")
      
  ggsave("./figures/GSE_six_waves.svg", width = 12, height = 7.5)
  
  

  ######################################################################
  ### Creating the line graphs related to each wave of transcription ###
  ######################################################################
  ### We find mean expression of each set of genes 
  genes <- unique(c(early, early_intermediate, intermediate, late_intermediate,late, up_down_up))
  df <- log_norm[genes,]
  col <- colnames(df)
  scaled <- t(apply(df, 1, scale))
  colnames(scaled) <- col
  
  
  ### Early intermediate
  scaled_early <- data.frame(scaled[early,])
  mean_early <- colMeans(scaled_early)
  
  scaled_early_intermediate <- data.frame(scaled[early_intermediate,])
  mean_early_intermediate <- colMeans(scaled_early_intermediate)
  
  scaled_intermediate <- data.frame(scaled[intermediate,])
  mean_intermediate <- colMeans(scaled_intermediate)
  
  scaled_late_intermediate <- data.frame(scaled[late_intermediate,])
  mean_late_intermediate <- colMeans(scaled_late_intermediate)
  
  scaled_late <- data.frame(scaled[late,])
  mean_late <- colMeans(scaled_late)
  
  scaled_up_down_up <- data.frame(scaled[up_down_up,])
  mean_up_down_up <- colMeans(scaled_up_down_up)
  
  df <- data.frame(early = mean_early,
                   early_intermediate = mean_early_intermediate,
                   intermediate = mean_intermediate,
                   late_intermediate = mean_late_intermediate,
                   late = mean_late,
                   up_down_up = mean_up_down_up,
                   lib = col)
  
  df$lib <- factor(df$lib,levels = col)
  df$time <- c(1:20)
  
  
  ggplot(df, aes(x = time, y = early)) + 
    geom_point() +
    geom_smooth(data = df, aes(x = time, y = early),method = lm, se = FALSE, formula = y ~ splines::bs(x, 3), fullrange=TRUE, color = "#E69F00") +
    theme_classic() +
    xlab("Pseudotime") +
    ylab("Mean expression") 
  ggsave("./figures/wave_early.svg", height = 4, width = 5)
  
  ggplot(df, aes(x = time, y = early_intermediate)) + 
    geom_point() +
    geom_smooth(data = df, aes(x = time, y = early_intermediate),method = lm, se = FALSE, formula = y ~ splines::bs(x, 3), fullrange=TRUE, color = "#E69F00") +
    theme_classic() +
    xlab("Pseudotime") +
    ylab("Mean expression") 
  ggsave("./figures/wave_early_intermediate.svg", height = 4, width = 5)
  
  
  ggplot(df, aes(x = time, y = intermediate)) + 
    geom_point() +
    geom_smooth(data = df, aes(x = time, y = intermediate),method = lm, se = FALSE, formula = y ~ splines::bs(x, 3), fullrange=TRUE, color = "#E69F00") +
    theme_classic() +
    xlab("Pseudotime") +
    ylab("Mean expression") 
  ggsave("./figures/wave_intermediate.svg", height = 4, width = 5)
  
  
  ggplot(df, aes(x = time, y = late_intermediate)) + 
    geom_point() +
    geom_smooth(data = df, aes(x = time, y = late_intermediate),method = lm, se = FALSE, formula = y ~ splines::bs(x, 3), fullrange=TRUE, color = "#E69F00") +
    theme_classic() +
    xlab("Pseudotime") +
    ylab("Mean expression") 
  ggsave("./figures/wave_late_intermediate.svg", height = 4, width = 5)
  
  ggplot(df, aes(x = time, y = late)) + 
    geom_point() +
    geom_smooth(data = df, aes(x = time, y = late),method = lm, se = FALSE, formula = y ~ splines::bs(x, 3), fullrange=TRUE, color = "#E69F00") +
    theme_classic() +
    xlab("Pseudotime") +
    ylab("Mean expression") 
  ggsave("./figures/wave_late.svg", height = 4, width = 5)
  
  ggplot(df, aes(x = time, y = up_down_up)) + 
    geom_point() +
    geom_smooth(data = df, aes(x = time, y = up_down_up),method = lm, se = FALSE, formula = y ~ splines::bs(x, 3), fullrange=TRUE, color = "#E69F00") +
    theme_classic() +
    xlab("Pseudotime") +
    ylab("Mean expression") 
  ggsave("./figures/wave_up_down_up.svg", height = 4, width = 5)
  
  
  
  
  #######################################################################
  ###################### SCENIC on INFECTED cells #######################
  #######################################################################
  ### the SCENIC results are included in the working directory clone I sent you 
  scenic.dir <- file.path("./SCENIC/results/D10_2/")
  regulonsAUC <- readRDS(file.path(scenic.dir, "int/3.4_regulonAUC.Rds"))
  AUC <- getAUC(regulonsAUC)
  
  auc <- AUC[,cell_order]  
  scaled <- t(apply(auc, 1, scale))
  colnames(scaled) <- cell_order
  
  
  ha = HeatmapAnnotation(Infection_load = Infection_load_order)
  pdf("./figures/SCENIC_infected_cilliated_single.pdf", width = 10, height = 20)
  p <- Heatmap(scaled, cluster_columns = FALSE, cluster_rows = TRUE,top_annotation = ha, show_column_names = FALSE, row_names_gp = gpar(fontsize = 3), col = colorRamp2(seq(-2,2,0.1), viridis(41)))
  print(p)
  dev.off()
  
  #### Interesting regulons
  regs_int <- c("JUND (391g)", "TP63 (32g)", "DDIT3 (219g)", "MYC_extended (231g)", "KLF7_extended (278g)", "NFKB2 (1295g)", "NFKB1 (754g)", 
                "IRF7 (1382g)", "IRF2 (1171g)", "STAT1 (1493g)", "STAT2 (2253g)", "STAT2 (2253g)", "ZNF467_extended (237g)",
                "CEBPG (305g)", "BCL3 (327g)", "ETV7_extended (626g)", "ATF4 (549g)", "MXI1_extended (664g)", "BCL6_extended (42g)", "MAX_extended (2186g)")

  auc <- AUC[regs_int,cell_order]  
  scaled <- t(apply(auc, 1, scale))
  colnames(scaled) <- cell_order
  
  
  ha = HeatmapAnnotation(Infection_load = Infection_load_order)
  svg("./figures/SCENIC_infected_cilliated_single_selected.svg", width = 10, height = 7)
  p <- Heatmap(scaled, cluster_columns = FALSE, cluster_rows = TRUE,top_annotation = ha, show_column_names = FALSE, col = colorRamp2(seq(-2,2,0.1), viridis(41)))
  print(p)
  dev.off()
  
  #################
  #################
  ## ZNF476 auc expression
  reg <- c("ZNF467_extended (237g)")
  auc <- AUC[reg,colnames(data_subset)]
    
  
  
  df2 <- data.frame(data_subset@reductions$umap@cell.embeddings, gene = data_subset@assays$RNA@data["ZNF467",], reg = auc)
  ggplot(df2, aes(x=UMAP_1, y = UMAP_2, color = reg)) +
    geom_jitter(size=2) +
    theme_classic() +
    scale_color_gradientn(colours = rev(brewer.pal(11, "RdYlGn")))
  
  ggsave("./figures/ZNF476_reg.svg", width = 7, height = 5.5)
  
  
  ###################
  ###################
  reg <- c("ZNF467_extended (237g)", "ZNF467 (63g)")
  scaled_ZNF <- scaled[reg, ]
  
  ha = HeatmapAnnotation(Infection_load = Infection_load_order)
  svg("./figures/ZNF467.svg", width = 10, height = 0.8)
  p <- Heatmap(scaled_ZNF, cluster_columns = FALSE, cluster_rows = TRUE,top_annotation = ha, show_column_names = FALSE, row_names_gp = gpar(fontsize = 3), col = colorRamp2(seq(-2,2,0.1), viridis(41)))
  print(p)
  dev.off()
  
  
  
  ####################################################################################################################

#### GO analysis of ZNF regulon  
  
  
  regulons <- readRDS("./SCENIC/results/D10_2/int/3.1_regulons_forAUCell.Rds")
  
  reg <- "ZNF467_extended (237g)"
  genes <- regulons[[reg]]
  write.csv(genes, file = "reg_ZNF467_extended.csv")
  
  gene.df <- bitr(genes, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
  ego_ZNF467 <- enrichGO(gene         = gene.df$ENSEMBL,
                             OrgDb         = org.Hs.eg.db,
                             keyType       = 'ENSEMBL',
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.05)
  
  write.csv(ego_ZNF467, file = "ZNF467_GO.csv")
  
  
  res <- ego_ZNF467@result
  pos <- data.frame(res[1:25,], wave = rep("positive",25))

  
  ggplot(pos, aes(x = factor(Description, levels = rev(pos$Description)), y = p.adjust, fill = wave)) +
    geom_bar(stat="identity") +
    scale_y_continuous(trans="log10") +
    scale_x_discrete(position = "top") +
    #scale_x_reverse() +
    geom_hline(yintercept = 0.05, linetype="dashed") + 
    coord_flip() +
    theme_classic() +
    xlab("") +
    ylab("adjusted p-value") +
    scale_fill_brewer(palette="Dark2") +
    theme(legend.position = "none")
  
  ggsave("./figures/ZNF_GSE.svg", width = 8, height = 5)
  
  
  
############################################################################################
################# Signature of the two other cilliated cell populations ####################
############################################################################################
  ### creating a few graphs that we need for the next steps:
  coord <- data.frame(data@reductions$umap@cell.embeddings)
  coord$color <- data@meta.data$seurat_clusters
  ggplot(coord, aes(x = UMAP_1, y = UMAP_2, color = color)) +
    geom_point() +
    theme_classic() +
    scale_color_brewer(palette = "Dark2")
 
  Idents(data) <- data$seurat_clusters
  markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(markers, file = "markers.csv")
  
  
  htmarker <- markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_logFC)
  scaled <- data@assays$RNA@scale.data[htmarker$gene,]
  cluster <- sort(data$seurat_clusters)
  scaled <- scaled[,names(cluster)]
  
  png("./figures/heatmap_DE_infected.png", width = 15, height = 20)
  p1 <- Heatmap(scaled, cluster_columns = FALSE, cluster_rows = FALSE, col = colorRamp2(seq(-2,2,0.1), viridis(41)), show_column_names = FALSE, show_row_names = TRUE, show_heatmap_legend = FALSE)
  print(p1)
  dev.off()
  
  
  htmarker <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
  Marker.Genes <- htmarker$gene
  Marker.Genes <- unique(Marker.Genes)
  DotPlot(data, features = Marker.Genes) + RotatedAxis() +
    xlab("") +
    ylab("") +
    theme() +
    scale_color_viridis_c(direction = 1, option = "plasma") +
    theme(axis.text.y = element_text(face = "italic")) #+
    #coord_flip() 
  
  ggsave("./figures/DotPlot5.svg", height = 3.5, width = 23)
  
  
  
  
  ### 3 vs 0
  clu3 <- rownames(data@meta.data[data@meta.data$seurat_clusters == 3,])
  coord3 <- coord[clu3,]
  
  clu3 <- rownames(data@meta.data[data@meta.data$seurat_clusters == 3,])
  clu3 <- data@meta.data[clu3,]
  
  
  clu0 <- rownames(data@meta.data[data@meta.data$seurat_clusters == 0,])
  clu0 <- data@meta.data[clu0,]
  
  
  cells <- c(rownames(clu3), rownames(clu0))
  data_subset_0vs3 <- subset(data, cells = cells)
  
  coord <- data.frame(data_subset_0vs3@reductions$umap@cell.embeddings)
  coord$color <- data_subset_0vs3@meta.data$seurat_clusters
  ggplot(coord, aes(x = UMAP_1, y = UMAP_2, color = color)) +
    geom_point() +
    theme_classic() +
    scale_color_manual(values=c('#1b9e77','#e7298a')) +
    theme(legend.position = "none")
  ggsave("./figures/cluster0vs3.svg", width = 5, height = 4)
  
    df <- data.frame(cluster = data_subset_0vs3$seurat_clusters, Infection_load = data_subset_0vs3$Infection_load)
    ggplot(df, aes(x = cluster, y= Infection_load, color = cluster)) +
      geom_jitter(size = 2) +
      #geom_violin() +
      scale_y_continuous(trans='log10') +
      theme_classic() +
      ylab("log10 (infection load)") +
      scale_color_manual(values=c('#1b9e77','#66a61e')) +
      theme(legend.position = "none")
ggsave("./figures/infection_load_0vs_3.svg", width = 5, height = 4)
  
  

  Idents(data_subset_0vs3) <- data_subset_0vs3$seurat_clusters
  cluster3.markers <- FindMarkers(data_subset_0vs3, ident.1 = 3, ident.2 = 0, min.pct = 0.25)
  write.csv(cluster3.markers, file = "./tables/cluster0vs3.markers.csv")
  
  
  cluster3_list <- cluster3.markers[cluster3.markers$avg_logFC > 0,]
  cluster0_list <- cluster3.markers[cluster3.markers$avg_logFC < 0,]
  
  
    hm_gene <- c(rownames(cluster3_list)[1:40], rownames(cluster0_list)[1:40])
    scaled <- data_subset_0vs3@assays$RNA@scale.data[hm_gene,]
    cell <- c(rownames(data_subset_0vs3@meta.data[data_subset_0vs3@meta.data$seurat_clusters == 0,]), rownames(data_subset_0vs3@meta.data[data_subset_0vs3@meta.data$seurat_clusters == 3,]))
    scaled <- scaled[,cell]
  
  
  
    svg("./figures/heatmap_0vs3.svg", width = 8, height = 12)
    p2 <- Heatmap(scaled, cluster_columns = FALSE, cluster_rows = TRUE, col = colorRamp2(seq(-2,2,0.1), viridis(41)), show_column_names = FALSE, show_row_names = TRUE, show_heatmap_legend = FALSE)
    print(p2)
    dev.off()
  
  
  
  
  clu3_gene <- rownames(cluster3_list)
  gene.df <- bitr(clu3_gene, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
  ego_clu3 <- enrichGO(gene         = gene.df$ENSEMBL,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = 'ENSEMBL',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
  
  
  clu0_gene <- rownames(cluster0_list)
  gene.df <- bitr(clu0_gene, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
  ego_clu0 <- enrichGO(gene         = gene.df$ENSEMBL,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = 'ENSEMBL',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
  
  
  
  pos <- data.frame(ego_clu3[1:20,], wave = rep("positive",20))
  neg <- data.frame(ego_clu0[1:20,], wave = rep("negative",20))
  #neg$Count <- neg$Count * -1
  neg$p.adjust <- log10(neg$p.adjust)
  pos$p.adjust <- log10(pos$p.adjust) * -1
  df <- rbind(pos, neg)
  
  df <- df[complete.cases(df),]
  df$row <- make.unique(paste(df$ID, ":",df$Description))
  df$row <- factor(df$row, levels = rev(df$row))
  
  ggplot(df, aes(x = row, y = p.adjust, fill = wave)) +
    geom_bar(stat="identity") +
    # scale_y_continuous(trans="log10") +
    scale_x_discrete(position = "top") +
    #geom_hline(yintercept = 0.05, linetype="dashed") + 
    coord_flip() +
    theme_classic() +
    xlab("") +
    ylab("log10(p.adjust)") +
    scale_fill_manual(values=c('#e7298a','#1b9e77')) +
    theme(legend.position = "none") #+
    #ggtitle("clu 0 vs clu 3")
  
  ggsave("./figures/GSE_clu0_clu3.svg", width = 10, height =6)  
  
  
  ggplot(df, aes(x = row, y = p.adjust, color = wave)) + 
    geom_point(aes(size = Count), alpha = 1) +
    scale_size(range = c(0.5, 12)) + # Adjust the range of points size
    theme_classic() +
    coord_flip()
  

  
  pos <- data.frame(ego_positive[1:25,], wave = rep("positive",25))
  neg <- data.frame(ego_negative[1:25,], wave = rep("negative",25))
  #neg$Count <- neg$Count * -1
  neg$p.adjust <- log10(neg$p.adjust)
  pos$p.adjust <- log10(pos$p.adjust) * -1
  df <- rbind(pos, neg)
  
  
  
  
  
  #########################################################
  ########################################################
  #########################################################
  ### 3 vs 0
  clu1 <- rownames(data@meta.data[data@meta.data$seurat_clusters == 1,])
  coord1 <- coord[clu1,]
  
  clu1 <- rownames(data@meta.data[data@meta.data$seurat_clusters == 1,])
  clu1 <- data@meta.data[clu1,]
  
  
  clu0 <- rownames(data@meta.data[data@meta.data$seurat_clusters == 0,])
  clu0 <- data@meta.data[clu0,]
  
  
  cells <- c(rownames(clu1), rownames(clu0))
  data_subset_0vs1 <- subset(data, cells = cells)
  
  coord <- data.frame(data_subset_0vs1@reductions$umap@cell.embeddings)
  coord$color <- data_subset_0vs1@meta.data$seurat_clusters
  ggplot(coord, aes(x = UMAP_1, y = UMAP_2, color = color)) +
    geom_point() +
    theme_classic() +
    scale_color_manual(values=c('#1b9e77','#e7298a')) +
    theme(legend.position = "none")
  ggsave("./figures/cluster0vs1.svg", width = 5, height = 4)
  
  df <- data.frame(cluster = data_subset_0vs1$seurat_clusters, Infection_load = data_subset_0vs1$Infection_load)
  ggplot(df, aes(x = cluster, y= Infection_load, color = cluster)) +
    geom_jitter(size = 2) +
    #geom_violin() +
    scale_y_continuous(trans='log10') +
    theme_classic() +
    ylab("log10 (infection load)") +
    scale_color_manual(values=c('#1b9e77','#66a61e')) +
    theme(legend.position = "none")
  ggsave("./figures/infection_load_0vs_1.svg", width = 5, height = 4)
  
  
  
  Idents(data_subset_0vs1) <- data_subset_0vs1$seurat_clusters
  cluster1.markers <- FindMarkers(data_subset_0vs1, ident.1 = 1, ident.2 = 0, min.pct = 0.25)
  write.csv(cluster1.markers, file = "./tables/cluster0vs1.markers.csv")
  
  
  cluster1_list <- cluster1.markers[cluster1.markers$avg_logFC > 0,]
  cluster0_list <- cluster1.markers[cluster1.markers$avg_logFC < 0,]
  
  
  hm_gene <- c(rownames(cluster1_list)[1:40], rownames(cluster0_list)[1:40])
  scaled <- data_subset_0vs1@assays$RNA@scale.data[hm_gene,]
  cell <- c(rownames(data_subset_0vs1@meta.data[data_subset_0vs1@meta.data$seurat_clusters == 0,]), rownames(data_subset_0vs1@meta.data[data_subset_0vs1@meta.data$seurat_clusters == 1,]))
  scaled <- scaled[,cell]
  
  
  
  svg("./figures/heatmap_0vs1.svg", width = 8, height = 12)
  p2 <- Heatmap(scaled, cluster_columns = FALSE, cluster_rows = TRUE, col = colorRamp2(seq(-2,2,0.1), viridis(41)), show_column_names = FALSE, show_row_names = TRUE, show_heatmap_legend = FALSE)
  print(p2)
  dev.off()
  
  
  
  
  clu1_gene <- rownames(cluster1_list)
  gene.df <- bitr(clu1_gene, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
  ego_clu1 <- enrichGO(gene         = gene.df$ENSEMBL,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05)
  
  
  clu0_gene <- rownames(cluster0_list)
  gene.df <- bitr(clu0_gene, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
  ego_clu0 <- enrichGO(gene         = gene.df$ENSEMBL,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENSEMBL',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05)
  
  
  
  pos <- data.frame(ego_clu1[1:20,], wave = rep("positive",20))
  neg <- data.frame(ego_clu0[1:20,], wave = rep("negative",20))
  #neg$Count <- neg$Count * -1
  neg$p.adjust <- log10(neg$p.adjust)
  pos$p.adjust <- log10(pos$p.adjust) * -1
  df <- rbind(pos, neg)
  
  df <- df[complete.cases(df),]
  df$row <- make.unique(paste(df$ID, ":",df$Description))
  df$row <- factor(df$row, levels = rev(df$row))
  
  ggplot(df, aes(x = row, y = p.adjust, fill = wave)) +
    geom_bar(stat="identity") +
    # scale_y_continuous(trans="log10") +
    scale_x_discrete(position = "top") +
    #geom_hline(yintercept = 0.05, linetype="dashed") + 
    coord_flip() +
    theme_classic() +
    xlab("") +
    ylab("log10(p.adjust)") +
    scale_fill_manual(values=c('#e7298a','#1b9e77')) +
    theme(legend.position = "none") #+
  #ggtitle("clu 0 vs clu 1")
  
  ggsave("./figures/GSE_clu0_clu1.svg", width = 10, height =6)  
  
  
  ggplot(df, aes(x = row, y = p.adjust, color = wave)) + 
    geom_point(aes(size = Count), alpha = 1) +
    scale_size(range = c(0.5, 12)) + # Adjust the range of points size
    theme_classic() +
    coord_flip()
  
  
  
  pos <- data.frame(ego_positive[1:25,], wave = rep("positive",25))
  neg <- data.frame(ego_negative[1:25,], wave = rep("negative",25))
  #neg$Count <- neg$Count * -1
  neg$p.adjust <- log10(neg$p.adjust)
  pos$p.adjust <- log10(pos$p.adjust) * -1
  df <- rbind(pos, neg)
  
  
  
  ########################################
  ########################################
  ########################################
  ### Correlation analysis on other infected clusters >>> these are supplementary figures with correlation analysi
  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  ### Focusing on the analysis of cluster one (1)
  cells <- rownames(data@meta.data[data$seurat_clusters == 1,])
  data_subset <- subset(data, cells = cells)
  
  ### Finding the genes that are correlated with viral load
  correlation <- data.frame(cor = apply(data_subset@assays$RNA@data, 1, cor , y = data_subset$Infection_load), ident = rep("ident",21624), label = rownames(data_subset@assays$RNA@data))
  correlation <- correlation[complete.cases(correlation),]
  
  color <- c(1:20640)
  for (i in 1:20640) {
    if (correlation$cor[i] > 0.4 | correlation$cor[i] < -0.4) {
      color[i] <- "strong"
    } else {
      color[i] <- "weak"
    }
  }
  
  correlation <- data.frame(correlation, condition = color)
  
  ggplot(correlation, aes(x=ident, y=cor)) +
    geom_jitter(aes(x=ident, y=cor, color = condition),size = 1, position = position_jitter(seed = 1)) +
    geom_boxplot(outlier.shape = NA) +
    theme_classic() +
    geom_hline(yintercept = 0.4,  linetype="dashed") +
    geom_hline(yintercept = -0.4,  linetype="dashed") +
    xlab("") +
    scale_color_manual(values=c('#F8766D', '#999999')) +
    ylab("coefiicient of correlation") +
    coord_flip() +
    theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  ggsave("./figures/correlation_boxplot_cluster1.svg", width = 11, height = 2)

  
  
  ### creating a heatmap for highly correlated genes
  cell_order <- names(sort(data_subset$Infection_load, decreasing = FALSE))
  Infection_load_order <- as.numeric(sort(data_subset$Infection_load, decreasing = FALSE))
  time <- data.frame(data$HTO_classification)[cell_order, ]
  
  hm_gene <- correlation[correlation$cor > 0.4 | correlation$cor < -0.4, ]
  norm <- data_subset@assays$RNA@data[rownames(hm_gene), cell_order]
  col <- colnames(norm)
  scaled <- t(apply(norm, 1, scale))
  colnames(scaled) <- col
  
  col_fun = colorRamp2(c(min(data_subset$Infection_load), max(data_subset$Infection_load)), c("white", "red"))
  ha = HeatmapAnnotation(Infection_load = Infection_load_order ,simple_anno_size = unit(1, "cm"), show_legend = FALSE, col = list(Infection_load = col_fun, time = c("Hashtag1-TotalA" = "#e41a1c", "Hashtag2-TotalA" = "#377eb8", "Hashtag3-TotalA" = "#4daf4a", "Hashtag4-TotalA" = "#984ea3")))
  svg("./figures/correlated_gene_heatmap_cluster1.svg", width = 11, height = 20)
  p <- Heatmap(scaled, cluster_columns = FALSE, cluster_rows = TRUE,top_annotation = ha, col = colorRamp2(seq(-2,2,0.1), viridis(41)), show_column_names = FALSE, show_heatmap_legend = FALSE)
  print(p)
  dev.off()
  
  ### GO term analysis of positive and negative genes
  positive <- rownames(correlation[correlation$cor > 0.3,])
  negative <- rownames(correlation[correlation$cor < -0.3,])
  gene.df <- bitr(positive, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
  ego_positive <- enrichGO(gene         = gene.df$ENSEMBL,
                           OrgDb         = org.Hs.eg.db,
                           keyType       = 'ENSEMBL',
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.05)
  
  write.csv(ego_positive, file = "./tables/GSE_positive.csv")
  
  
  gene.df <- bitr(negative, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
  ego_negative <- enrichGO(gene         = gene.df$ENSEMBL,
                           OrgDb         = org.Hs.eg.db,
                           keyType       = 'ENSEMBL',
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.05)
  
  write.csv(ego_negative, file = "./tables/GSE_negative.csv")
  
  pos <- data.frame(ego_positive[1:25,], wave = rep("positive",25))
  neg <- data.frame(ego_negative[1:25,], wave = rep("negative",25))
  #neg$Count <- neg$Count * -1
  neg$p.adjust <- log10(neg$p.adjust)
  pos$p.adjust <- log10(pos$p.adjust) * -1
  df <- rbind(pos, neg)
  
  df <- df[complete.cases(df),]
  df$row <- make.unique(paste(df$Description, ":",df$ID))
  df$row <- factor(df$row, levels = rev(df$row))
  
  ggplot(df, aes(x = row, y = p.adjust, fill = wave)) +
    geom_bar(stat="identity") +
    #scale_y_continuous(trans="log10") +
    #scale_x_discrete(position = "top") +
    #geom_hline(yintercept = 0.05, linetype="dashed") + 
    coord_flip() +
    theme_classic() +
    geom_vline(xintercept = 0.5) + 
    xlab("") +
    ylab("log10(p.adjust)") +
    scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
    theme(legend.position = "none")
  
  ggsave("./figures/GSE_pos_neg_cluster1.svg", width = 10, height = 7) 
  
  
  #######################################################################################################################################
  #######################################################################################################################################
  #######################################################################################################################################
  ### Focusing on the analysis of cluster three (3)
  cells <- rownames(data@meta.data[data$seurat_clusters == 3,])
  data_subset <- subset(data, cells = cells)
  
  ### Finding the genes that are correlated with viral load
  correlation <- data.frame(cor = apply(data_subset@assays$RNA@data, 1, cor , y = data_subset$Infection_load), ident = rep("ident",21624), label = rownames(data_subset@assays$RNA@data))
  correlation <- correlation[complete.cases(correlation),]
  
  color <- c(1:18590)
  for (i in 1:18590) {
    if (correlation$cor[i] > 0.4 | correlation$cor[i] < -0.4) {
      color[i] <- "strong"
    } else {
      color[i] <- "weak"
    }
  }
  
  correlation <- data.frame(correlation, condition = color)
  
  ggplot(correlation, aes(x=ident, y=cor)) +
    geom_jitter(aes(x=ident, y=cor, color = condition),size = 1, position = position_jitter(seed = 1)) +
    geom_boxplot(outlier.shape = NA) +
    theme_classic() +
    geom_hline(yintercept = 0.4,  linetype="dashed") +
    geom_hline(yintercept = -0.4,  linetype="dashed") +
    xlab("") +
    scale_color_manual(values=c('#F8766D', '#999999')) +
    ylab("coefiicient of correlation") +
    coord_flip() +
    theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  ggsave("./figures/correlation_boxplot_cluster3.svg", width = 11, height = 2)
  
  
  
  ### creating a heatmap for highly correlated genes
  cell_order <- names(sort(data_subset$Infection_load, decreasing = FALSE))
  Infection_load_order <- as.numeric(sort(data_subset$Infection_load, decreasing = FALSE))
  time <- data.frame(data$HTO_classification)[cell_order, ]
  
  hm_gene <- correlation[correlation$cor > 0.4 | correlation$cor < -0.4, ]
  norm <- data_subset@assays$RNA@data[rownames(hm_gene), cell_order]
  col <- colnames(norm)
  scaled <- t(apply(norm, 1, scale))
  colnames(scaled) <- col
  
  col_fun = colorRamp2(c(min(data_subset$Infection_load), max(data_subset$Infection_load)), c("white", "red"))
  ha = HeatmapAnnotation(Infection_load = Infection_load_order ,simple_anno_size = unit(1, "cm"), show_legend = FALSE, col = list(Infection_load = col_fun, time = c("Hashtag1-TotalA" = "#e41a1c", "Hashtag2-TotalA" = "#377eb8", "Hashtag3-TotalA" = "#4daf4a", "Hashtag4-TotalA" = "#984ea3")))
  svg("./figures/correlated_gene_heatmap_cluster3.svg", width = 11, height = 20)
  p <- Heatmap(scaled, cluster_columns = FALSE, cluster_rows = TRUE,top_annotation = ha, col = colorRamp2(seq(-2,2,0.1), viridis(41)), show_column_names = FALSE, show_heatmap_legend = FALSE)
  print(p)
  dev.off()
  
  ### GO term analysis of positive and negative genes
  positive <- rownames(correlation[correlation$cor > 0.3,])
  negative <- rownames(correlation[correlation$cor < -0.3,])
  gene.df <- bitr(positive, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
  ego_positive <- enrichGO(gene         = gene.df$ENSEMBL,
                           OrgDb         = org.Hs.eg.db,
                           keyType       = 'ENSEMBL',
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.05)
  
  write.csv(ego_positive, file = "./tables/GSE_positive.csv")
  
  
  gene.df <- bitr(negative, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
  ego_negative <- enrichGO(gene         = gene.df$ENSEMBL,
                           OrgDb         = org.Hs.eg.db,
                           keyType       = 'ENSEMBL',
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.05)
  
  write.csv(ego_negative, file = "./tables/GSE_negative.csv")
  
  pos <- data.frame(ego_positive[1:25,], wave = rep("positive",25))
  neg <- data.frame(ego_negative[1:25,], wave = rep("negative",25))
  #neg$Count <- neg$Count * -1
  neg$p.adjust <- log10(neg$p.adjust)
  pos$p.adjust <- log10(pos$p.adjust) * -1
  df <- rbind(pos, neg)
  
  df <- df[complete.cases(df),]
  df$row <- make.unique(paste(df$Description, ":",df$ID))
  df$row <- factor(df$row, levels = rev(df$row))
  
  ggplot(df, aes(x = row, y = p.adjust, fill = wave)) +
    geom_bar(stat="identity") +
    #scale_y_continuous(trans="log10") +
    #scale_x_discrete(position = "top") +
    #geom_hline(yintercept = 0.05, linetype="dashed") + 
    coord_flip() +
    theme_classic() +
    geom_vline(xintercept = 0.5) + 
    xlab("") +
    ylab("log10(p.adjust)") +
    scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
    theme(legend.position = "none")
  
  ggsave("./figures/GSE_pos_neg_cluster3.svg", width = 10, height = 7)
    
