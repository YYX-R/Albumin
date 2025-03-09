##install R packages
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(devtools)
library(harmony)
library(reticulate)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readr)
library(stringr)
library(ggpubr)
library(patchwork) 
library(ggplot2)


###get multiple samples data
path <- getwd()
full_path <- paste(path, "/single_cell_yyx/", sep = "")
print(full_path)
dir_name=list.files("single_cell_yyx")
scRNAlist <- list()
for(i in 1:length(dir_name)){
  counts <- Read10X(data.dir = paste(full_path, dir_name[i], sep = ""))
  scRNAlist[[i]] <- CreateSeuratObject(counts, project = dir_name[i],min.cells = 3, min.features = 300)
}


####Batch calculation of mitochondrial and red blood cell ratios####
for(i in 1:length(scRNAlist)){
  sc <- scRNAlist[[i]]
  sc[["mt_percent"]] <- PercentageFeatureSet(sc, pattern = "^mt-")
  HB_genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") 
  HB_m <- match(HB_genes, rownames(sc@assays$RNA))
  HB_genes <- rownames(sc@assays$RNA)[HB_m]  
  HB_genes <- HB_genes[!is.na(HB_genes)] 
  sc[["HB_percent"]] <- PercentageFeatureSet(sc, features=HB_genes)   
  scRNAlist[[i]] <- sc
  rm(sc)
}




####Batch drawing of violin plots before quality control####
violin_before <- list()
for(i in 1:length(scRNAlist)){
  violin_before[[i]] <- VlnPlot(scRNAlist[[i]],
                                features = c("nFeature_RNA", "nCount_RNA", "mt_percent","HB_percent"), 
                                pt.size = 0.01, 
                                ncol = 4) 
}
violin_before_merge <- CombinePlots(plots = violin_before,nrow=length(scRNAlist),legend='none')
violin_before_merge




####Batch filtering of cells, MT, HB genes####
scRNAlist <- lapply(X = scRNAlist, FUN = function(x){
  x <- subset(x, 
              subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & 
                mt_percent < 10 & 
                HB_percent < 3 & 
                nCount_RNA < quantile(nCount_RNA,0.95) & 
                nCount_RNA > 1000)})
View(scRNAlist[[1]]@meta.data)


####merge samples####
scRNAlist <- merge(x=scRNAlist[[1]],y=scRNAlist[-1])
## count the number of cells
table(scRNAlist[[]]$orig.ident)

colnames(scRNAlist@meta.data)
Idents(scRNAlist) <- scRNAlist@meta.data$orig.ident
Idents(scRNAlist) <- factor(Idents(scRNAlist), levels = c("WT", "Cre"))


###Draw a violin plot after quality control###
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
WT_color <- "#A3CFF7"  
CRE_color <- "#F7A3A3"  
scRNAlist$orig.ident <- factor(scRNAlist$orig.ident, levels = c("WT", "Cre"))
violin_nFeature_RNA <- VlnPlot(scRNAlist, 
                               features = "nFeature_RNA", 
                               pt.size = 0,  
                               group.by = "orig.ident",  
                               fill = "orig.ident") +  
  scale_fill_manual(values = c("WT" = WT_color, "Cre" = CRE_color)) + 
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Arial", face = "bold"),
    axis.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA)  
  )
violin_nCount_RNA <- VlnPlot(scRNAlist, 
                             features = "nCount_RNA", 
                             pt.size = 0,  
                             group.by = "orig.ident",  
                             fill = "orig.ident") +  
  scale_fill_manual(values = c("WT" = WT_color, "Cre" = CRE_color)) + 
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Arial", face = "bold"),
    axis.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA)  
  )
violin_mt_percent <- VlnPlot(scRNAlist, 
                             features = "mt_percent", 
                             pt.size = 0,  
                             group.by = "orig.ident", 
                             fill = "orig.ident") + 
  scale_fill_manual(values = c("WT" = WT_color, "Cre" = CRE_color)) +  
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Arial", face = "bold"),
    axis.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA)  
  )
combined_violin_plot <- violin_nFeature_RNA | violin_nCount_RNA | violin_mt_percent
combined_violin_plot


ggsave(file = "combined_violin_plot.png",plot =combined_violin_plot,he = 5,wi = 10,dpi = 600 )






####Data normalization, screening of highly variable genes and PCA dimensionality reduction####
scRNAlist <- NormalizeData(scRNAlist) %>%
  FindVariableFeatures(selection.method = "vst",nfeatures = 5000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 30, verbose = T)
a=DimPlot(scRNAlist,reduction = "pca",group.by = "orig.ident")
print(a)
top15 <- head(VariableFeatures(scRNAlist), 15) 
plot1 <- VariableFeaturePlot(scRNAlist) 
plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE, size=3) 
feat_15 <- CombinePlots(plots = list(plot1,plot2),legend = "bottom")
feat_15
ggsave(file = "feat_15.png",plot = feat_15,he = 10,wi = 15 )

scRNAlist <- JoinLayers(scRNAlist)#connect the count data of layers



####cell cycle score####
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(scRNAlist))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(scRNAlist))
scRNAlist <- CellCycleScoring(object=scRNAlist,  g2m.features=g2m_genes,  s.features=s_genes)
scRNAlist=CellCycleScoring(object = scRNAlist, 
                           s.features = s_genes, 
                           g2m.features = g2m_genes, 
                           set.ident = TRUE)
scRNAlist <- CellCycleScoring(object=scRNAlist,  g2m.features=g2m_genes,  s.features=s_genes)
scRNAlist@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal()

####RunHarmony####
scRNA_harmony <- RunHarmony(scRNAlist, group.by.vars = "orig.ident")
scRNA_harmony@reductions[["harmony"]][[1:5,1:5]]
b=DimPlot(scRNA_harmony,reduction = "harmony",group.by = "orig.ident")
b
pca_harmony_integrated <- CombinePlots(list(a,b),ncol=1)
pca_harmony_integrated

save(scRNA_harmony,scRNAlist,file = "scdata2.Rdata")
load("scdata2.Rdata")


####umap/tsne####
ElbowPlot(scRNA_harmony, ndims=50, reduction="harmony")
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:10) %>% FindClusters(resolution = 0.2)
head(scRNA_harmony@meta.data)
scRNA_harmony_filtered <- subset(scRNA_harmony, seurat_clusters != "3")
scRNA_harmony <-scRNA_harmony_filtered 
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:10)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:10)
umap_integrated1 <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "orig.ident")
umap_integrated2 <- DimPlot(scRNA_harmony, reduction = "umap", label = TRUE)
tsne_integrated1 <- DimPlot(scRNA_harmony, reduction = "tsne", group.by = "orig.ident") 
tsne_integrated2 <- DimPlot(scRNA_harmony, reduction = "tsne", label = TRUE)
umap_tsne_integrated <- CombinePlots(list(tsne_integrated1,tsne_integrated2,umap_integrated1,umap_integrated2),ncol=2)
umap_tsne_integrated
ggsave("umap_tsne_integrated.png",umap_tsne_integrated,wi=25,he=15)

scRNA_harmony@meta.data$seurat_clusters = scRNA_harmony@meta.data[["RNA_snn_res.0.2"]]
table(scRNA_harmony@meta.data$seurat_clusters)
scRNA_harmony <- JoinLayers(scRNA_harmony)











###Draw a bar chart of cell population ratios###
cluster_data$orig.ident <- factor(cluster_data$orig.ident, levels = c("WT", "Cre"))
cluster_counts <- cluster_data %>%
  group_by(orig.ident, seurat_clusters) %>%
  summarise(count = n()) %>%
  ungroup()
total_counts <- cluster_counts %>%
  group_by(orig.ident) %>%
  summarise(total = sum(count))

cluster_ratios <- cluster_counts %>%
  left_join(total_counts, by = "orig.ident") %>%
  mutate(ratio = count / total)  

print(cluster_ratios)

library(ggplot2)
ggplot(cluster_ratios, aes(x = orig.ident, y = ratio, fill = as.factor(seurat_clusters))) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Group", y = "Proportion", title = "Cell Proportion Comparison by Cluster") +
  theme_minimal() +  
  scale_fill_manual(values = c("0" = "#E53935", "1" = "#388E3C", "2" = "#4575B4", "4" = "purple")) + 
  theme(
    text = element_text(size = 12, family = "Arial", face = "bold"),  
    axis.text.x = element_text(angle = 0, hjust = 0.5),  
    legend.title = element_blank()
  )



cluster_colors <- c("0" = "#D73027", 
                    "1" = "#1A9850", 
                    "2" = "#4575B4", 
                    "4" = "#8E44AD")
umap_integrated <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "seurat_clusters", label = FALSE) +
  scale_color_manual(values = cluster_colors) + 
  theme(
    text = element_text(size = 12, family = "Arial", face = "bold"), 
    plot.title = element_blank(),  
    axis.text = element_text(size = 8, face = "bold")  
  )
umap_integrated

ggsave("umap_integrated.png",umap_integrated,wi=3,he=3,dpi = 600)

tsne_integrated <- DimPlot(scRNA_harmony, reduction = "tsne", group.by = "seurat_clusters", label = TRUE) +
  scale_color_manual(values = cluster_colors) +  
  theme(
    text = element_text(size = 12, family = "Arial", face = "bold"),  
    plot.title = element_blank(),  
    axis.text = element_text(size = 8, face = "bold")  
  )

tsne_integrated
umap_integrated + tsne_integrated + plot_layout(ncol = 2)

###Markers 
c=VlnPlot(scRNA_harmony,
          features = c("Cd3e","Cd8a","Cd4","Tcf7","Il7r","Sell","Tox","Havcr2","Pdcd1","Entpd1")
          ,  
          pt.size = 0,
          group.by = "seurat_clusters",
          stack = TRUE,
          same.y.lims = TRUE,
          flip = TRUE) + 
  NoLegend() 
ggsave("c.png",c,wi=3,he=6,dpi = 600)




cluster_colors <- c("0" = "#D73027", 
                    "1" = "#1A9850", 
                    "2" = "#4575B4", 
                    "4" = "#8E44AD")

umap_wt_cluster <- DimPlot(subset(scRNA_harmony, orig.ident == "WT"), 
                           reduction = "umap", 
                           group.by = "seurat_clusters", 
                           label = FALSE) + 
  ggtitle("WT") + 
  scale_color_manual(values = cluster_colors) +  
  theme(
    text = element_text(size = 12, family = "Arial", face = "bold"),
    plot.title = element_text(hjust = 0.5),  
    axis.text = element_text(size = 8, face = "bold")
  )
ggsave("umap_wt_cluster.png",umap_wt_cluster,wi=3,he=3,dpi = 600)


umap_cre_cluster <- DimPlot(subset(scRNA_harmony, orig.ident == "Cre"), 
                            reduction = "umap", 
                            group.by = "seurat_clusters", 
                            label = FALSE) + 
  ggtitle("Cre") + 
  scale_color_manual(values = cluster_colors) + 
  theme(
    text = element_text(size = 12, family = "Arial", face = "bold"),
    plot.title = element_text(hjust = 0.5), 
    axis.text = element_text(size = 8, face = "bold")
  )
ggsave("umap_cre_cluster.png",umap_cre_cluster,wi=3,he=3,dpi = 600)

umap_wt_cluster + umap_cre_cluster






genes_of_interest <- c("Gzmb", "Ifng", "Tnf", "Il2","Tfe3")

scRNA_harmony_sub_cluster2 <- subset(scRNA_harmony, seurat_clusters == 2)
scRNA_harmony_sub_cluster4 <- subset(scRNA_harmony, seurat_clusters == 4)

scRNA_harmony$orig.ident <- factor(scRNA_harmony$orig.ident, levels = c("WT", "Cre"))

dot_plot_cluster2 <- DotPlot(scRNA_harmony_sub_cluster2, features = genes_of_interest, group.by = "orig.ident") +  
  scale_size_continuous(range = c(2, 8)) +  
  labs(x = "Cluster", y = "Gene Expression", title = "Gene Expression in Cluster 2") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Arial", face = "bold"),
    axis.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA)  
  )

dot_plot_cluster4 <- DotPlot(scRNA_harmony_sub_cluster4, features = genes_of_interest, group.by = "orig.ident") +  
  scale_size_continuous(range = c(2, 8)) +  
  labs(x = "Cluster", y = "Gene Expression", title = "Gene Expression in Cluster 4") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, family = "Arial", face = "bold"),
    axis.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA) 
  )


combined_plot <- dot_plot_cluster2 / dot_plot_cluster4
combined_plot
ggsave("combined_plot.png",combined_plot,wi=4.5,he=3.2,dpi = 600)
