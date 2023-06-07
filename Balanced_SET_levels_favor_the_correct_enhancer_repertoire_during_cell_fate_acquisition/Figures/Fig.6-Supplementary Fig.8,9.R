library(Seurat)
library(dplyr)
library(patchwork)
library(stringr)
library(ggplot2)
library(openxlsx)
library(pals)
library(RColorBrewer)
library(openxlsx)
library(gprofiler2)
library(data.table)
library(readxl)
library(ComplexHeatmap)

options(future.globals.maxSize = 8000 * 1024^2)

setwd("./")

#Fig. 6 e UMAP plot of Gene Expression and chromatin accessibility


Multiome <- readRDS("/home/zaghi/Documents/scRNA_Multiome/MultiOme.rds") #load single cell gene expression object from Seurat

scATAC_Multiome <- loadArchRProject("Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR/") #load ArchR object for single cell chromatin accessibility

#transfer single cell ATAC UMAP coordinates in Seurat object to plot UMAP


barcodes_scATAC <- as.data.frame(scATAC_Multiome@embeddings@listData[["UMAP"]]@listData[["df"]]) %>% 
  rownames_to_column() %>% 
  rename(barcode_scATAC=1) 

barcodes_scATAC$sample_number <- str_split_fixed(barcodes_scATAC$barcode_scATAC, pattern = "-", n = 2)[,2]
barcodes_scATAC$sample <- str_split_fixed(barcodes_scATAC$barcode_scATAC, pattern = "#", n = 2)[,1]
barcodes_scATAC$barcode <- str_split_fixed(barcodes_scATAC$barcode_scATAC, pattern = "#", n = 2)[,2]
barcodes_scATAC$barcode <- str_split_fixed(barcodes_scATAC$barcode, pattern = "-", n = 2)[,1]

barcodes_scATAC <- barcodes_scATAC %>% 
  rename(sample_numer_ATAC=2)

barcodes_scATAC$sample[barcodes_scATAC$sample == "embryo_mut"] <- "1"
barcodes_scATAC$sample[barcodes_scATAC$sample == "embryo_ctrl"] <- "2"
barcodes_scATAC$sample[barcodes_scATAC$sample == "P2_ctrl"]<- "3"
barcodes_scATAC$sample[barcodes_scATAC$sample == "P2_mut"]<- "4"

barcodes_scATAC$barcode <- paste(barcodes_scATAC$barcode,sep = "-",barcodes_scATAC$sample)
rownames(barcodes_scATAC) <- barcodes_scATAC[,6]
barcodes_scATAC$barcode <- NULL

barcodes_scATAC <- barcodes_scATAC[c(2:3)]

barcodes_scATAC <- barcodes_scATAC %>% 
  dplyr::rename(UMAP_1=1,
                UMAP_2=2)

barcodes_scATAC <- barcodes_scATAC[order(rownames(barcodes)),]
barcodes_scATAC <- barcodes_scATAC[rownames(Multihome@meta.data),]

Multihome[["ATAC"]] <- CreateDimReducObject(embeddings = Multihome[["umap"]]@cell.embeddings, assay = DefaultAssay(Multihome),
                                             key = "cloupe_")

Multihome[["ATAC"]]@cell.embeddings[,1] <- barcodes_scATAC[,1]
Multihome[["ATAC"]]@cell.embeddings[,2] <- barcodes_scATAC[,2]

#plot UMAP with clusters and samples for ATAC & RNA

#ATAC
colors <- kelly(11)
colors[1] <- "brown" 
p1 = DimPlot(Multihome, reduction = "ATAC", label = F, pt.size = 0.3, cols = colors,group.by="Label_cluster") + 
  theme_void()+
  theme(legend.text = element_text(color = "black", size=12))
p1
ggsave("SETBP1_epigenomics/pipeline/plots/UMAP_no_labels_ATAC.png", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 200, height = 200, units = "mm", dpi = 300, limitsize = TRUE)
  p1+theme(plot.title = element_text(color="black", size=7, face="bold.italic"),
        axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=5, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 5),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=5),
        axis.title.y = element_text(face = "bold", color = "black", size = 5),
        legend.text = element_text(face = "bold", color = "black", size = 12),
        legend.position="top", text = element_text(face = "bold", color = "black", size=5),
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "UMAP 1", y = "UMAP 2")
  
  colors <- c("red","blue","darkgreen","purple")
  p1 = DimPlot(Multihome, reduction = "ATAC", label = F, pt.size = 0.3, cols = colors,group.by="Condition") + 
    theme_void()+
    theme(legend.text = element_text(color = "black", size=12))
  p1
  ggsave("SETBP1_epigenomics/pipeline/plots/UMAP_no_labels_ATAC_samples.png", plot = last_plot(), device = NULL, path = NULL,
         scale = 1, width = 200, height = 200, units = "mm", dpi = 300, limitsize = TRUE)
  p1+theme(plot.title = element_text(color="black", size=7, face="bold.italic"),
           axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=5, hjust =1), 
           axis.title.x = element_text(face = "bold", color = "black", size = 5),
           axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=5),
           axis.title.y = element_text(face = "bold", color = "black", size = 5),
           legend.text = element_text(face = "bold", color = "black", size = 12),
           legend.position="top", text = element_text(face = "bold", color = "black", size=5),
           panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "UMAP 1", y = "UMAP 2")

#RNA
colors <- kelly(11)
colors[1] <- "brown"
  p1 = DimPlot(Multihome, reduction = "umap", label = F, pt.size = 0.3, cols = colors,group.by="Label_cluster") + 
    theme_void()+
    theme(legend.text = element_text(color = "black", size=12))
  p1
  ggsave("SETBP1_epigenomics/pipeline/plots/UMAP_no_labels_RNA.png", plot = last_plot(), device = NULL, path = NULL,
         scale = 1, width = 200, height = 200, units = "mm", dpi = 300, limitsize = TRUE)
  p1+theme(plot.title = element_text(color="black", size=7, face="bold.italic"),
           axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=5, hjust =1), 
           axis.title.x = element_text(face = "bold", color = "black", size = 5),
           axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=5),
           axis.title.y = element_text(face = "bold", color = "black", size = 5),
           legend.text = element_text(face = "bold", color = "black", size = 12),
           legend.position="top", text = element_text(face = "bold", color = "black", size=5),
           panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "UMAP 1", y = "UMAP 2")
  
  
  colors <- c("red","blue","darkgreen","purple")
  p1 = DimPlot(Multihome, reduction = "umap", label = F, pt.size = 0.3, cols = colors,group.by="Condition") + 
    theme_void()+
    theme(legend.text = element_text(color = "black", size=12))
  p1
  ggsave("SETBP1_epigenomics/pipeline/plots/UMAP_no_labels_RNA_samples.png", plot = last_plot(), device = NULL, path = NULL,
         scale = 1, width = 200, height = 200, units = "mm", dpi = 300, limitsize = TRUE)
  p1+theme(plot.title = element_text(color="black", size=7, face="bold.italic"),
           axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=5, hjust =1), 
           axis.title.x = element_text(face = "bold", color = "black", size = 5),
           axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=5),
           axis.title.y = element_text(face = "bold", color = "black", size = 5),
           legend.text = element_text(face = "bold", color = "black", size = 12),
           legend.position="top", text = element_text(face = "bold", color = "black", size=5),
           panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "UMAP 1", y = "UMAP 2")
  

#Fig. 6 g Heatmap of chromatin accessibility along pseudotime 

#Select cells containing only ctrl or mutant cells to perform the analysis by genotype

idxPass <- which(scATAC_Multiome$Sample %in% c("embryo_ctrl","P2_ctrl"))
cellsPass <- scATAC_Multiome$cellNames[idxPass]
scATAC_Multiome_ctrl <- scATAC_Multiome[cellsPass, ]

idxPass <- which(scATAC_Multiome$Sample %in% c("embryo_mut","P2_mut"))
cellsPass <- scATAC_Multiome$cellNames[idxPass]
scATAC_Multiome_mut <- scATAC_Multiome[cellsPass, ]

#add to each genotype specific object the specific peaks associated with the specific clusters analyzed 

peaks1 <- read_tsv("Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/PeakCalls/clusters/da_edo/ctrl/AP_RGC") 

peaks2 <- read_tsv("Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/PeakCalls/clusters/da_edo/ctrl/INP") 

peaks3<- read_tsv("Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/PeakCalls/clusters/da_edo/ctrl/ExN_DL") 

peaks4 <- read_tsv("Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/PeakCalls/clusters/da_edo/mut/AP_RGC") 

peaks5 <- read_tsv("Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/PeakCalls/clusters/da_edo/mut/INP") 

peaks6 <- read_tsv("Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/PeakCalls/clusters/da_edo/mut/ExN_DL") 

peak <- rbind.data.frame(peaks1,peaks2,peaks3,
                         peaks4,peaks5,peaks6) %>% 
  makeGRangesFromDataFrame() %>% 
  IRanges::reduce()

scATAC_Multiome_ctrl <- addPeakSet(
  ArchRProj =  scATAC_Multiome_ctrl,
  peakSet = peak,
  genomeAnnotation = getGenomeAnnotation(scATAC_Multiome),
  force = TRUE
)

scATAC_Multiome_ctrl <- addPeakMatrix(scATAC_Multiome_ctrl, force=T)

peaks1 <- read_tsv("Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/PeakCalls/clusters/da_edo/ctrl/AP_RGC") 

peaks2 <- read_tsv("Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/PeakCalls/clusters/da_edo/ctrl/INP") 

peaks3<- read_tsv("Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/PeakCalls/clusters/da_edo/ctrl/ExN_UL") 

peaks4 <- read_tsv("Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/PeakCalls/clusters/da_edo/mut/AP_RGC") 

peaks5 <- read_tsv("Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/PeakCalls/clusters/da_edo/mut/INP") 

peaks6 <- read_tsv("Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/PeakCalls/clusters/da_edo/mut/ExN_UL") 



peak <- rbind.data.frame(peaks1,peaks2,peaks3,
                         peaks4,peaks5,peaks6) %>% 
  makeGRangesFromDataFrame() %>% 
  IRanges::reduce()

scATAC_Multiome_mut <- addPeakSet(
  ArchRProj =  scATAC_Multiome_mut,
  peakSet = peak,
  genomeAnnotation = getGenomeAnnotation(scATAC_Multiome),
  force = TRUE
)

scATAC_Multiome_mut <- addPeakMatrix(scATAC_Multiome_mut, force=T)

#Calculate pseudotime trajectory based on the the clusters of interest

trajectory_neu <- c("AP_RGC", "INP", "ExN_DL")

scATAC_Multiome_ctrl <- addTrajectory(
  ArchRProj = scATAC_Multiome_ctrl, #add pseudotime trajecotry of interest to ArchR object
  name = "Neurons", 
  groupBy = "Clusters",
  trajectory = trajectory_neu, 
  embedding = "UMAP_2", 
  force = TRUE,
  reducedDims =NULL
)

scATAC_Multiome_mut <- addTrajectory(
  ArchRProj = scATAC_Multiome_mut, 
  name = "Neurons", 
  groupBy = "Clusters",
  trajectory = trajectory_neu, 
  embedding = "UMAP_2", 
  force = TRUE
)

#Compute and plot Heatmaps of ATAC peaks Z-score across pseudotime 


trajPM_ctrl  <- getTrajectory(ArchRProj = scATAC_Multiome_ctrl, name = "Neurons", useMatrix = "PeakMatrix", log2Norm = TRUE) #Compute trajectory based on Peaks matrix of control cells

p_trajPM_ctrl <- plotTrajectoryHeatmap(trajPM_ctrl, pal = paletteContinuous(set = "solarExtra"), returnMatrix = T)
#plot peak matrix based on control peaks
png("SETBP1_epigenomics/pipeline/plots/p_trajPM_ctrl.png",pointsize = 1,res=1200,height = 20,width = 20,
    units = "cm")
p_trajPM_ctrl
dev.off()

p_trajPM_mut <- plotTrajectoryHeatmap(trajPM_mut, pal = paletteContinuous(set = "solarExtra"), returnMatrix = T)
#plot peak matrix based on mutant peaks
png("SETBP1_epigenomics/pipeline/plots/p_trajPM_mut.png",pointsize = 1,res=1200,height = 20,width = 20,
    units = "cm")
p_trajPM_mut
dev.off()




library(monocle3)
########Trjectories with monocle3#######

cds_sort <- subset(Multi.object, subset=Condition %in% c("embryo_ctrl", "P2_ctrl") & Label_cluster %in% c("AP_RGC", "INP", "ExN_DL", "ExN_UL", "Late_Prog"))
cds <- as.cell_data_set(cds_sort, cell_metadata=as.data.frame(cds_sort@meta.data))
cds <- preprocess_cds(cds, num_dim = 10, norm_method = "log")
cds <- align_cds(cds, alignment_group = "Condition")

cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
cds <- learn_graph(cds)
cds <- order_cells(cds, reduction_method = "UMAP")


cds@clusters$UMAP$partition <- cds_sort@meta.data$Label_cluster
names(cds@clusters$UMAP$partition) <- rownames(cds_sort@meta.data)
cds@clusters$UMAP$clusters <- cds_sort@meta.data$Label_cluster
names(cds@clusters$UMAP$clusters) <- rownames(cds_sort@meta.data)

traj_ctrl_Cond <- plot_cells(cds, color_cells_by = "Condition", norm_method = "log", label_cell_groups = FALSE)+
  theme(plot.title = element_text(color="black", size=5, face="bold.italic"),
                       axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=5, hjust =1), 
                       axis.title.x = element_text(face = "bold", color = "black", size = 5),
                       axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=5),
                       axis.title.y = element_text(face = "bold", color = "black", size = 5),
                       legend.text = element_text(face = "bold", color = "black", size = 5),
                       legend.position="top",
                       panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "UMAP 1", y = "UMAP 2")+
  ggtitle("ctrl")
traj_ctrl_Clus <- plot_cells(cds, color_cells_by = "cluster", group_cells_by = "cluster", norm_method = "log", label_cell_groups = FALSE)+
  theme(plot.title = element_text(color="black", size=5, face="bold.italic"),
                       axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=5, hjust =1), 
                       axis.title.x = element_text(face = "bold", color = "black", size = 5),
                       axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=5),
                       axis.title.y = element_text(face = "bold", color = "black", size = 5),
                       legend.text = element_text(face = "bold", color = "black", size = 5),
                       legend.position="top",
                       panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "UMAP 1", y = "UMAP 2")+
  ggtitle("ctrl")

Multi.object <- AddMetaData(object = Multi.object, metadata = cds@principal_graph_aux@listData$UMAP$pseudotime, col.name = "Pseudotime_ctrl" )

cds_sort <- subset(Multi.object, subset=Condition %in% c("embryo_mut", "P2_mut") & Label_cluster %in% c("AP_RGC", "INP", "ExN_DL", "ExN_UL", "Late_Prog"))
cds <- as.cell_data_set(cds_sort, cell_metadata=as.data.frame(cds_sort@meta.data))
cds <- preprocess_cds(cds, num_dim = 10, norm_method = "log")
cds <- align_cds(cds, alignment_group = "Condition")

cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds = cds, reduction_method = "UMAP")
cds <- learn_graph(cds)
cds <- order_cells(cds, reduction_method = "UMAP")


cds@clusters$UMAP$partition <- cds_sort@meta.data$Label_cluster
names(cds@clusters$UMAP$partition) <- rownames(cds_sort@meta.data)
cds@clusters$UMAP$clusters <- cds_sort@meta.data$Label_cluster
names(cds@clusters$UMAP$clusters) <- rownames(cds_sort@meta.data)


traj_ctrl_Cond <- plot_cells(cds, color_cells_by = "Condition", norm_method = "log", label_cell_groups = FALSE)+
  theme(plot.title = element_text(color="black", size=5, face="bold.italic"),
                       axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=5, hjust =1), 
                       axis.title.x = element_text(face = "bold", color = "black", size = 5),
                       axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=5),
                       axis.title.y = element_text(face = "bold", color = "black", size = 5),
                       legend.text = element_text(face = "bold", color = "black", size = 5),
                       legend.position="top",
                       panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "UMAP 1", y = "UMAP 2")+
  ggtitle("mut")
traj_ctrl_Clus <- plot_cells(cds, color_cells_by = "cluster", group_cells_by = "cluster", norm_method = "log", label_cell_groups = FALSE)+
  theme(plot.title = element_text(color="black", size=5, face="bold.italic"),
                       axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=5, hjust =1), 
                       axis.title.x = element_text(face = "bold", color = "black", size = 5),
                       axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=5),
                       axis.title.y = element_text(face = "bold", color = "black", size = 5),
                       legend.text = element_text(face = "bold", color = "black", size = 5),
                       legend.position="top",
                       panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "UMAP 1", y = "UMAP 2")+
  ggtitle("mut")

Multi.object <- AddMetaData(object = Multi.object, metadata = cds@principal_graph_aux@listData$UMAP$pseudotime, col.name = "Pseudotime_mut" )



library(scales)
library(ggridges)
library(viridis)
library(hrbrthemes)
library(easyGgplot2)

ctrl_Time <- Multi.object@metadata %>% filter(!is.na(Pseudotime_ctrl))
mut_Time <- Multi.object@metadata %>% filter(!is.na(Pseudotime_mut))

ctrl_Time$Pseudotime_ctrl_scaled <- rescale(ctrl_Time$Pseudotime_ctrl, to = c(0, 100))
mut_Time$Pseudotime_mut_scaled <- rescale(mut_Time$Pseudotime_mut, to = c(0, 100))
ctrl_Time <- ctrl_Time[c(2, 4)]
ctrl_Time$Condition <- "ctrl"
names(ctrl_Time)[2] <- "Pseudotime"
mut_Time <- mut_Time[c(2, 4)]
mut_Time$Condition <- "mut"
names(mut_Time)[2] <- "Pseudotime"

pseudo <- rbind(ctrl_Time, mut_Time)
pseudo["Label_cluster"][pseudo["Label_cluster"] == "ExN_UL"] <- "ExN"
pseudo["Label_cluster"][pseudo["Label_cluster"] == "ExN_DL"] <- "ExN"
pseudo["Label_cluster"][pseudo["Label_cluster"] == "ExN_L1"] <- "ExN"

a <- ggplot2.density(data=filter(pseudo, pseudo$Label_cluster %in% c("ExN", "INP" ,"AP_RGC")), 
                     xName='Pseudotime', groupName='Label_cluster', legendPosition="top",
                     faceting=TRUE, facetingVarNames="Condition", densityFill = 'Label_cluster', 
                     fillGroupDensity = T, colorGroupDensityLine = T, 
                     groupColors = c("#BE0032","cornflowerblue" ,"#F3C300"))+
  #scale_fill_manual(values=kelly(length(unique(pseudo$Label_cluster))+1)[2:4])+
  theme(plot.title = element_text(color="black", size=20, face="bold.italic", hjust = 0),
        axis.text.x = element_text(angle = 45, face = "bold", color = "black", size=20, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 20),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=20),
        axis.title.y = element_text(face = "bold", color = "black", size = 20, vjust = 0),
        legend.text = element_text(face = "bold", color = "black", size = 20),
        legend.position="top",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "", y = "% of Cell")



b <- ggplot2.density(data=filter(pseudo, pseudo$Label_cluster %in% c("ExN", "INP" ,"AP_RGC")), 
                     xName='Pseudotime', groupName='Condition', legendPosition="top",
                     faceting=TRUE, facetingVarNames="Condition", densityFill = 'Label_cluster', 
                     fillGroupDensity = T, colorGroupDensityLine = T, 
                     groupColors = okabe(8)[6:7])+
  theme(plot.title = element_text(color="black", size=5, face="bold.italic"),
        axis.text.x = element_text(angle = 45, face = "bold", color = "black", size=20, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 20),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=20),
        axis.title.y = element_text(face = "bold", color = "black", size = 20, vjust = 0),
        legend.text = element_text(face = "bold", color = "black", size = 20),
        legend.position="top",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "Pseudotime", y = "% of Cell")



c <- ggplot2.density(data=filter(pseudo, pseudo$Label_cluster %in% c("INP", "ExN", "Late_Prog", "AP_RGC", "IN", "ExN_L1", "Astro")),
                     xName='Pseudotime', groupName='Label_cluster', legendPosition="top",
                     faceting=TRUE, facetingVarNames="Condition", densityFill = 'Label_cluster', 
                     fillGroupDensity = T, colorGroupDensityLine = T, 
                     groupColors = watlington(length(unique(pseudo$Label_cluster))))+
  #scale_fill_manual(values=kelly(length(unique(pseudo$Label_cluster))+1)[2:4])+
  theme(plot.title = element_text(color="black", size=20, face="bold.italic", hjust = 0),
        axis.text.x = element_text(angle = 45, face = "bold", color = "black", size=20, hjust =1), 
        axis.title.x = element_text(face = "bold", color = "black", size = 20),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=20),
        axis.title.y = element_text(face = "bold", color = "black", size = 20, vjust = 0),
        legend.text = element_text(face = "bold", color = "black", size = 20),
        legend.position="top",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  labs(x = "", y = "% of Cell")+
  ggtitle("Pseudotime distribution of lineage ExN -> INP -> AP_RGC")


png("ExN_INP_AP_RGC_in_pseudotime.png", res = 330, height = 20, width = 30, units = "in")
c/a/b
dev.off()

########################Diffrential Enhancer Usage#############

library(plyranges)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(dplyr)
library(stringr)
library(viridis)
library(ggplot2)
library(openxlsx)
library(pheatmap)



files <- list.files("/home/zaghi/Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/PeakCalls/clusters/ctrl/", full.names = T)
WT_peaks <- lapply(files, read.table,  header = T)
names(WT_peaks) <- list.files("/home/zaghi/Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/PeakCalls/clusters/ctrl/", full.names = F)

files <- list.files("/home/zaghi/Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/PeakCalls/clusters/mut/", full.names = T)
MUT_peaks <- lapply(files, read.table,  header = T)
names(MUT_peaks) <- list.files("/home/zaghi/Documents/scATAC_Multiome_Bam/scATAC_multiome_ArchR_correct/PeakCalls/clusters/mut/", full.names = F)


for (i in 1:length(WT_peaks)) {
  WT_peaks[[i]] <- WT_peaks[[i]] %>% makeGRangesFromDataFrame(keep.extra.columns = T)
  MUT_peaks[[i]] <- MUT_peaks[[i]] %>% makeGRangesFromDataFrame(keep.extra.columns = T)
}

common_peaks <- list()

for (i in 1:length(WT_peaks)) {
  common_peaks[[i]] <- subsetByOverlaps(WT_peaks[[i]], MUT_peaks[[i]])
  names(common_peaks)[i] <- names(WT_peaks)[i]
}

for (i in 1:length(WT_peaks)) {
  common_peaks[[i]] <- common_peaks[[i]] %>%
    annotatePeak(tssRegion=c(-10000, 2000), 
                 TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene, 
                 annoDb="org.Mm.eg.db") %>%
    as.data.frame()
  
  WT_peaks[[i]] <- WT_peaks[[i]] %>%
    annotatePeak(tssRegion=c(-10000, 2000), 
                 TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene, 
                 annoDb="org.Mm.eg.db") %>%
    as.data.frame()
  
  MUT_peaks[[i]] <- MUT_peaks[[i]] %>%
    annotatePeak(tssRegion=c(-10000, 2000), 
                 TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene, 
                 annoDb="org.Mm.eg.db") %>% 
    as.data.frame()
  
}

WT_gene_count <- list()

for (i in names(WT_peaks)) {
  WT_gene_count[[i]] <- WT_peaks[[i]] %>% 
    group_by(SYMBOL) %>%
    count() %>%
    as.data.frame()
}

MUT_gene_count <- list()

for (i in names(WT_peaks)) {
  MUT_gene_count[[i]] <- MUT_peaks[[i]] %>% 
    group_by(SYMBOL) %>%
    count() %>%
    as.data.frame()
}

common_gene_count <- list()

for (i in names(WT_peaks)) {
  common_gene_count[[i]] <- common_peaks[[i]] %>% 
    group_by(SYMBOL) %>%
    count() %>%
    as.data.frame()
}


all_counts <- list()


for (i in names(WT_gene_count)) {
  all_counts[[i]] <- full_join(full_join(WT_gene_count[[i]], MUT_gene_count[[i]], by = "SYMBOL"), common_gene_count[[i]], by = "SYMBOL")
  all_counts[[i]][is.na(all_counts[[i]])] <- 0 
  names(all_counts[[i]])[2:4] <- c("WT_count", "MUT_count", "Common_count")
}


for (i in names(all_counts)) {
  all_counts[[i]]$Unique_WT <- all_counts[[i]]$WT_count - all_counts[[i]]$Common_count
  all_counts[[i]]$Unique_MUT <- all_counts[[i]]$MUT_count - all_counts[[i]]$Common_count
  all_counts[[i]] <- filter(all_counts[[i]], all_counts[[i]]$Unique_MUT >= 0)
}


for (i in 1:length(all_counts)) {
  all_counts[[i]]$Perc_of_common <-  (all_counts[[i]]$Common_count/(all_counts[[i]]$Unique_WT + all_counts[[i]]$Unique_MUT + all_counts[[i]]$Common_count))*100
  all_counts[[i]]$Condition <- names(all_counts)[[i]]
}

for (i in 1:length(all_counts)) {
  all_counts[[i]] <-  filter(all_counts[[i]], all_counts[[i]]$Perc_of_common < 100 | all_counts[[i]]$Perc_of_common > 0 )
}




all_plots <- list()

for (i in 1:10) {
  p <- ggplot(all_counts[[i]], aes(x=all_counts[[i]]$Unique_WT, y=all_counts[[i]]$Unique_MUT, color = all_counts[[i]]$Perc_of_common)) +
    geom_point() +
    scale_color_viridis()+
    geom_jitter(position = "jitter", aes(x=all_counts[[i]]$Unique_WT + 0.2, y=all_counts[[i]]$Unique_MUT + 0.2))+
    geom_smooth(method = "loess")+
    theme(plot.title = element_text(color="black", size=20, face="bold.italic", hjust = 0),
          axis.text.x = element_text(angle = 45, face = "bold", color = "black", size=20, hjust =1), 
          axis.title.x = element_text(face = "bold", color = "black", size = 15),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=20),
          axis.title.y = element_text(face = "bold", color = "black", size = 15, vjust = 0),
          legend.text = element_text(face = "bold", color = "black", size = 10),
          legend.position="top",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "Unique WT", y = "Unique MUT")+
    ggtitle(paste(names(all_counts)[i]))+
    xlim(c(0, max(append(all_counts[[i]]$Unique_WT, all_counts[[i]]$Unique_MUT))))+
    ylim(c(0, max(append(all_counts[[i]]$Unique_WT, all_counts[[i]]$Unique_MUT))))
  ggsave(p, filename = paste(names(all_counts)[i], ".png", sep = ""), width = 10, height = 10)
  
}



to_vio <- bind_rows(all_counts)
ggplot(filter(to_vio, to_vio$Common_count > 0), aes(x=Condition, y=Common_count)) + 
  geom_violin()+
  #geom_density()
  ylim(c(0, 10))

for (i in names(all_counts)) {
  density(all_counts[[i]]$MUT_count)
}

library(readxl)
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

DEGS <- read_excel_allsheets("DEGS_intra_cluster.xlsx")


for (i in intersect(names(DEGS), names(all_counts))) {
  all_counts[[i]] <- dplyr::left_join(all_counts[[i]], DEGS[[i]], by="SYMBOL")
}

for (i in intersect(names(DEGS), names(all_counts))) {
  all_counts[[i]] <- all_counts[[i]] %>%
    mutate(DEGS = case_when(
      avg_log2FC > 0 ~ "UP",
      avg_log2FC < 0 ~ "DOWN",
      is.na(avg_log2FC) ~ "NO_DEGS"))
  
}


for (i in intersect(names(DEGS), names(all_counts))) {
  p <- ggplot(all_counts[[i]], aes(x=all_counts[[i]]$Unique_WT, y=all_counts[[i]]$Unique_MUT, color = as.factor(DEGS))) +
    theme_classic()+
    geom_point() +
    geom_jitter(position = "jitter", aes(x=all_counts[[i]]$Unique_WT + 0.2, y=all_counts[[i]]$Unique_MUT + 0.2)) +
    labs(x = "Unique WT", y = "Unique MUT")+
    ggtitle(paste(i))+
    xlim(c(0, max(append(all_counts[[i]]$Unique_WT, all_counts[[i]]$Unique_MUT))))+
    ylim(c(0, max(append(all_counts[[i]]$Unique_WT, all_counts[[i]]$Unique_MUT))))+
    scale_colour_manual(name = "Special Points",
                        values = c("blue", alpha("grey", 0.1), "red"))
  ggsave(p, filename = paste(i, ".png", sep = ""), width = 10, height = 10)
}



ggplot(all_counts[[1]], aes(x=all_counts[[1]]$Unique_WT, y=all_counts[[1]]$Unique_MUT, color = as.factor(DEGS))) +
  theme_classic()+
  geom_point() +
  labs(x = "Unique WT", y = "Unique MUT")+
  geom_jitter(position = "jitter", aes(x=all_counts[[1]]$Unique_WT + 0.2, y=all_counts[[1]]$Unique_MUT + 0.2))+
  xlim(c(0, max(append(all_counts[[i]]$Unique_WT, all_counts[[i]]$Unique_MUT))))+
  ylim(c(0, max(append(all_counts[[i]]$Unique_WT, all_counts[[i]]$Unique_MUT))))+
  scale_colour_manual(name = "Special Points",
                      values = c("blue", alpha("grey", 0.1), "red"))


for (i in  intersect(names(DEGS), names(all_counts))) {
  all_counts[[i]] <- all_counts[[i]][order(all_counts[[i]]$DEGS),]
}

to_write <- list()

for (i in intersect(names(DEGS), names(all_counts))) {
  df <- all_counts[[i]] %>%
    filter(DEGS!="NO_DEGS")
  to_write[[i]] <- df
}

openxlsx::write.xlsx(to_write, "Peaks_common_by_gene.xlsx", overwrite = F)





