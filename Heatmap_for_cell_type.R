setwd("C:/Users/Jacky/Desktop/PDAC/scRNA")

seu_obj <- readRDS("./RDS/Total_Harmony.rds")
library(ClusterGVis)
library(org.Hs.eg.db)
library(Seurat)
library(dplyr)

levels(seu_obj)
Idents(seu_obj) <- seu_obj@meta.data$celltype

#寻找标记基因
seu_obj.markers.all <- Seurat::FindAllMarkers(seu_obj,
                                           only.pos = TRUE,
                                           min.pct = 0.25,
                                           logfc.threshold = 0.25)

#选取top20的标记基因
seu_obj.markers <- seu_obj.markers.all %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 20, wt = avg_log2FC)

head(seu_obj.markers)

st.data <- prepareDataFromscRNA(object = seu_obj,
                                diffData = seu_obj.markers,
                                showAverage = TRUE)

str(st.data)

enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.25,
                        topn = 5,
                        seed = 5201314)

head(enrich)
write.table(enrich, file = "enrich.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

#挑选需要展示的特征基因
markGenes = unique(seu_obj.markers$gene)[sample(1:length(unique(seu_obj.markers$gene)),75,
                                             replace = F)]

markGenes

visCluster(object = st.data,
           plot.type = "line")

pdf('sc1.pdf',height = 10,width = 6,onefile = F)
visCluster(object = st.data,
           plot.type = "heatmap",
           column_names_rot = 75,
           markGenes = markGenes,
           cluster.order = c(1:10)) #cluster.order 参数调整聚类顺序。
dev.off()

pdf('sc2.pdf',height = 10,width = 14,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 75,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           cluster.order = c(1:10),
           go.col = rep(jjAnno::useMyCol("stallion",n = 9),each = 5),
           add.bar = T)
dev.off()

st.data <- prepareDataFromscRNA(object = seu_obj,
                                diffData = seu_obj.markers,
                                showAverage = TRUE,
                                keep.uniqGene = FALSE,
                                sep = "_")

df <- st.data$wide.res #比如 CD3D 这个基因有三个重复,表示在三个亚群都有显著性差异。

# line plot
visCluster(object = st.data,
           plot.type = "line") #我们再来看折线图数量就都是 20 了。
#这些基因在热图上也都会展示:
# heatmap plot
pdf('sc3.pdf',height = 10,width = 6,onefile = F)
visCluster(object = st.data,
           plot.type = "heatmap",
           column_names_rot = 45,
           markGenes = markGenes,
           cluster.order = c(1:10))
dev.off()

seu_obj.markers1 <- seu_obj.markers.all %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 6, wt = avg_log2FC)

# retain duplicate diff gene in multiple clusters
st.data <- prepareDataFromscRNA(object = seu_obj,
                                diffData = seu_obj.markers1,
                                showAverage = FALSE)

# check
str(st.data)

# heatmap plot
pdf('sc4.pdf',height = 10,width = 8,onefile = F)
visCluster(object = st.data,
           plot.type = "heatmap",
           markGenes = unique(seu_obj.markers1$gene),
           column_title_rot = 45,
           cluster.order = 1:10
) #记得设置 show_column_names = F。
dev.off()

pdf('sc5.pdf',height = 10,width = 8,onefile = F)
new.cluster.ids <- levels(seu_obj)

visCluster(object = st.data,
           plot.type = "heatmap",
           markGenes = unique(seu_obj.markers1$gene),
           column_title_rot = 45,
           cluster.order = 1:10,
           #show_column_names = F,
           sample.cell.order = rev(new.cluster.ids),
           sample.col = jjAnno::useMyCol("paired",n = 10))
dev.off()

#添加富集注释:
# add GO annotation
pdf('sc6.pdf',height = 12,width = 16,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_title_rot = 45,
           markGenes = unique(seu_obj.markers1$gene),
           markGenes.side = "left",
           annoTerm.data = enrich,
           #show_column_names = F,
           line.side = "left",
           cluster.order = c(1:10),
           add.bar = T)
dev.off()

dev.on()
library(cols4all)
mycol <- c4a('vivid',12) #ѡȡ??ɫ

p1 <- DimPlot(seu_obj, reduction = "SeuratUMAP2D", label = TRUE, cols = mycol)
p2 = DimPlot(seu_obj, reduction = "umap", group.by='orig.ident', cols = mycol)
plotc <- p1+p2
plotc

dev.new()

devtools::install_github("junjunlab/seu_objtoolVis")
library(seu_objtoolVis)

clusterCornerAxes(object = seu_obj,reduction = 'umap',
                  noSplit = T)

clusterCornerAxes(object = seu_obj,reduction = 'umap',
                  noSplit = T,arrowType = 'open')


clusterCornerAxes(object = seu_obj,reduction = 'umap',
                  noSplit = T,
                  cornerTextSize = 3.5,
                  themebg = 'bwCorner',
                  addCircle = TRUE,
                  cicAlpha = 0.2,
                  nbin = 200)

library(RColorBrewer) 
library(viridis)
library(wesanderson)
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
pie(rep(6,n), col=sample(color, n))
col_vector
col_vector =c(wes_palette("Darjeeling1"), wes_palette("GrandBudapest1"), wes_palette("Cavalcanti1"), wes_palette("GrandBudapest2"), wes_palette("FantasticFox1"))
pal <- wes_palette("Zissou1", 10, type = "continuous")
pal2 <- wes_palette("Zissou1", 5, type = "continuous")
pal[3:10]

# 画dotplot
final.markers <- c('TOX','PDCD1','LAG3','TIGIT','ZNF683','ITGAE','CXCR6','HSPA1A','HSPA1B','HSPH1','ENO1','HMGA1','GTF3A','GZMK','CXCR3','FCGR3A','CX3CR1','FGFBP2','FCER1G','IKZF2','TYROBP','TRDC','TRGC1','TRDV1')

seu_obj$seurat_clusters <- factor(x = seu_obj$seurat_clusters, levels = c('progenitors','circulating','basal','active'))
Idents(seu_obj) <- "integrated_snn_res.0.5"
library(ggplot2)

##  RotatedAxis() scale_colour_gradientn都是比较重要的
DotPlot(seu_obj, features = final.markers , dot.scale = 10) + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1), axis.text.y = element_text(face="bold")) + 
  scale_colour_gradientn(colours = pal)+ theme(legend.position="right")  + labs(title = "cluster markers", y = "", x="")

library(stringr)  
genes_to_check <- c('TOX','PDCD1','LAG3','TIGIT','ZNF683','ITGAE','CXCR6','HSPA1A','HSPA1B','HSPH1','ENO1','HMGA1','GTF3A','GZMK','CXCR3','FCGR3A','CX3CR1','FGFBP2','FCER1G','IKZF2','TYROBP','TRDC','TRGC1','TRDV1')
genes_to_check=str_to_upper(unique(genes_to_check))
genes_to_check

th=theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
p_all_markers <- DotPlot(seu_obj, features = genes_to_check,
                         assay='RNA' ,group.by = 'customclassif' )  + coord_flip()+th
p_all_markers

ggsave(p_all_markers, filename="check_all_marker_by_seurat_celltype.pdf",height = 7,width = 9)

data<-p_all_markers$data

colnames(data)

colnames(data)<-c("AverageExpression_unscaled","Precent Expressed","Features","celltype","Average Expression")

unique(data$`Precent Expressed`)

p = ggplot(data,aes(celltype,Features,size = `Precent Expressed` ))+
  geom_point(shape=21,aes(fill= `Average Expression`),position =position_dodge(0))+
  theme_minimal()+xlab(NULL)+ylab(NULL) +
  scale_size_continuous(range=c(1,10))+theme_bw()+
  scale_fill_gradient(low = "grey", high = "#E54924")+
  theme(legend.position = "right",legend.box = "vertical", #图例位置
        legend.margin=margin(t= 0, unit='cm'),
        legend.spacing = unit(0,"in"),
        axis.text.x  = element_text(color="black",size=16,angle = 45, 
                                    vjust = 0.5, hjust=0.5),#x轴
        axis.text.y  = element_text(color="black",size=12),#y轴
        legend.text = element_text(size =12,color="black"),#图例
        legend.title = element_text(size =12,color="black"),#图例
        axis.title.y=element_text(vjust=1,  
                                  size=16)
  )+labs(x=" ",y = "Features");p

p1 = ggplot(data,aes(celltype,Features,size = `Precent Expressed` ))+
  geom_point(shape=21,aes(fill= `Average Expression`),position =position_dodge(0))+
  theme_minimal()+xlab(NULL)+ylab(NULL) +
  scale_size_continuous(range=c(1,10))+theme_bw()+
  scale_fill_gradient(low = "#E54924", high = "#498EA4")+
  theme(legend.position = "right",legend.box = "vertical",
        legend.margin=margin(t= 0, unit='cm'),
        legend.spacing = unit(0,"in"),
        axis.text.x  = element_text(color="black",size=16,angle = 45, 
                                    vjust = 0.5, hjust=0.5),
        axis.text.y  = element_text(color="black",size=12),
        legend.text = element_text(size =12,color="black"),
        legend.title = element_text(size =12,color="black"),
        axis.title.y=element_text(vjust=1,  
                                  size=16)
  )+labs(x=" ",y = "Features")

p1



























