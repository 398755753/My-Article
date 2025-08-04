getwd()
list.files("/public/home/CD4 Apt-CS/High_CD4")
high_CD4.data <- Read10X(data.dir = "/public/home/CD4 Apt-CS/High_CD4")
dim(high_CD4.data)
high_CD4 <- CreateSeuratObject(counts = high_CD4.data,project = "high_CD4",min.cells = 3,min.features = 200)
high_CD4

list.files("/public/home/CD4 Apt-CS/PBS_Control")
PBS.data <- Read10X(data.dir = "/public/home/CD4 Apt-CS/PBS_Control")
dim(PBS.data)
PBS <- CreateSeuratObject(counts = PBS.data,project = "PBS",min.cells = 3,min.features = 200)
PBS

combined <- merge(high_CD4,y=PBS)

#查看线粒体基因表达
combined[["percent.mt"]]<-PercentageFeatureSet(combined,pattern="^mt-")
VlnPlot(combined,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
plot1<-FeatureScatter(combined,feature1 = "nCount_RNA",feature2 = "percent.mt")
plot2<-FeatureScatter(combined,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
CombinePlots(plots=list(plot1,plot2))
combined <-subset(combined,subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10)
#数据标准化
combined <-NormalizeData(combined,normalization.method = "LogNormalize",scale.factor = 10000)
#高变基因
combined <- FindVariableFeatures(combined, selection.method = "vst",nfeatures = 2000)
top10 <- head(VariableFeatures(combined),10)
plot1<-VariableFeaturePlot(combined)
plot2<-LabelPoints(plot=plot1,points=top10, repel =TRUE)
CombinePlots(plots = list(plot1,plot2))
#归一化，选取的高变基因进行scale
combined <- ScaleData(combined, features = VariableFeatures(combined))

#PCA分析
combined <- RunPCA(combined, verbose = T)
print(combined[["pca"]], dims = 1:5, nfeatures = 5) 

library(harmony)
combined  <- RunHarmony(combined, group.by.vars = "orig.ident")

DimPlot(combined, reduction = "pca")
DimHeatmap(combined, dims = 1:2, cells = 500, balanced = TRUE)
ElbowPlot(combined)


#细胞聚类分析
combined <- FindNeighbors(combined, dims = 1:25)
library(clustree)
combined <- FindClusters(object = combined,resolution = c(seq(.1,1.6,.45)))
clustree(combined@meta.data, prefix = "RNA_snn_res.")

combined <- FindClusters(object = combined,resolution = 0.5)
combined <- RunUMAP(combined,dims = 1:25)
DimPlot(combined, reduction = "umap",group.by = "orig.ident")
DimPlot(combined, reduction = "umap",label = TRUE)
save(combined,file = "combined_V4_clusterunnamed.rData")



#---------------差异marker基因分析---------------------
#combined2 <- JoinLayers(combined)
combined.markers <- FindAllMarkers(combined,only.pos = TRUE,min.pct=0.25,logfc.threshold = 0.25)
top5 <- combined.markers %>% group_by(cluster) %>% top_n(n=5,wt=avg_log2FC)
DoHeatmap(combined, features = top5$gene, label = T) 
write.csv(top10,"top10.csv")



new.cluster.ids <- c("0"="Tumor cells","1"="Tumor cells","2"="Tumor cells","5"="Tumor cells","12"="Tumor cells",
                     "8"="Fibroblasts","11"="Fibroblasts","15"="Fibroblasts",
                     "13"="Endothelial",
                     "10"="T cells",
                     "14"="B&DC cells",
                     "3"="Monocyte","6"="Monocyte",
                     "4"="Macrophage","7"="Macrophage","9"="Macrophage","16"="Macrophage")
combined2_named <- RenameIdents(combined, new.cluster.ids) 
combined2_named$celltype <- combined2_named@active.ident
DimPlot(combined2_named, group.by = "celltype",label = F,pt.size = 0.005)
save(combined2_named,file = "combined2_V4_clusternamed.rData")
combined2_named_TSNE <- RunTSNE(combined2_named,dims = 1:25)
DimPlot(combined2_named_TSNE, reduction = "tsne", group.by = "seurat_clusters",label = F,pt.size = 1)

#把内皮和成纤维合并
new.cluster.ids <- c("0"="Tumor cells","1"="Tumor cells","2"="Tumor cells","5"="Tumor cells","12"="Tumor cells",
                     "8"="Stromal cells","11"="Stromal cells","15"="Stromal cells",
                     "13"="Stromal cells",
                     "10"="T cells",
                     "14"="B&DC cells",
                     "3"="Monocyte","6"="Monocyte",
                     "4"="Macrophage","7"="Macrophage","9"="Macrophage","16"="Macrophage")
combined2_named <- RenameIdents(combined, new.cluster.ids) 
combined2_named$celltype <- combined2_named@active.ident
DimPlot(combined2_named, group.by = "celltype",label = F,pt.size = 0.005)
save(combined2_named,file = "combined2_V4_clusternamed.rData")
combined2_named_TSNE <- RunTSNE(combined2_named,dims = 1:25)
DimPlot(combined2_named_TSNE, reduction = "tsne", group.by = "seurat_clusters",label = F,pt.size = 1)

combined_T_TSNE <- RunTSNE(combined_T,dims = 1:25)
colors <- c("#E9A86D","#FA6666","#D8A7DD","#70B073",
            "grey","grey","grey","grey",
            "grey","grey","grey","grey","grey")
DimPlot(combined_T_TSNE, reduction = "tsne", label = F,pt.size = 1.8, cols = colors)

F8766D
#E98C5D
#53B400
colors <- c("#F8766D","#D0A300","#00C094",
            "#00B6EB","#B39BFF","#FB71DA","grey","grey",
            "grey","grey","grey","grey","grey")
DimPlot(combined2_named_TSNE,reduction = "tsne", group.by = "celltype",pt.size = 1,cols = colors)
#导出图片9*12


#---------------T细胞差异基因------------------------------------

combined_T2 <-combined_T
meta.data <- combined_T2@meta.data
meta.data$orig.ident <- ifelse(meta.data$orig.ident == "high_CD4", "Post",
                               ifelse(meta.data$orig.ident == "PBS", "Pre", meta.data$orig.ident))
# 合并两列并添加下划线
meta.data$celltype <- paste(meta.data$orig.ident, "_", meta.data$celltype, sep="")

# 将更新后的metadata赋回Seurat对象
combined_T2@meta.data <- meta.data

# 确保Seurat对象更新成功
combined_T2

combined_CD4T <- combined_T[,(combined_T$celltype %in% c("CD4 T"))]
combined_CD4T.markers <- FindMarkers(combined_CD4T,only.pos = TRUE,
                                     ident.1 = "PBS", ident.2 = "high_CD4",min.pct=0.25,logfc.threshold = 0.25)
top10 <- combined_CD4T.markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
DoHeatmap(combined_CD4T, features = top10$gene, label = T) +NoLegend()



#


#---------------marker dimplot---------------------------------------------------
Tumor_marker <- c("Epcam", "Krt19")
Fibroblasts_marker <- c("Col1a1","Col1a2","Pdgfrb","Tagln")
Endothelial_marker <- c("Pecam1","Vwf","Plvap")
T_marker <- c("Cd3d","Cd3e")
B_DC_marker <- c("Cd79a","Iglc2","Xcr1","Siglech")
Monocyte_marker <- c("S100a8","S100a9")
Macrophage_marker <- c("Apoe","C1qa","Cd68","Ccr5")
levels(combined2_named) <- c("Tumor cells","Fibroblasts","Endothelial","T cells","B&DC cells","Monocyte","Macrophage")
feature <- list("Tumor_marker"=Tumor_marker,"Fibroblasts_marker"=Fibroblasts_marker,"Endothelial_marker"=Endothelial_marker,
                "T_marker"=T_marker,"B_DC_marker"=B_DC_marker,"Monocyte_marker"=Monocyte_marker,
                "Macrophage_marker"=Macrophage_marker)
p <-DotPlot(object = combined2_named,features = feature)
p1 <- ggplot(p$data,aes(x=features.plot,y=id))+
  geom_point(aes(size=pct.exp,color=avg.exp.scaled))+
  facet_grid(facets = ~feature.groups,switch = "x",scales = "free_x",space = "free_x")+
  theme_classic()+
  scale_color_gradient2(low = "blue",mid = "white",high = "red")+
  theme(axis.text.x = element_text(angle = 90))




levels(combined) <- c("Tumor cells","Fibroblasts","Endothelial","T cells","B&DC cells","Monocyte","Macrophage")

p <-DotPlot(object = combined)
p1 <- ggplot(p$data,aes(x=features.plot,y=id))+
  geom_point(aes(size=pct.exp,color=avg.exp.scaled))+
  facet_grid(facets = ~feature.groups,switch = "x",scales = "free_x",space = "free_x")+
  theme_classic()+
  scale_color_gradient2(low = "blue",mid = "white",high = "red")+
  theme(axis.text.x = element_text(angle = 90))




markers <- c("Cd4","Foxp3") 
DotPlot(combined_T,features = markers)+coord_flip()
VlnPlot(combined_T, features = c("Cd69"),split.by = c("orig.ident"),y.max = 5)

markers <- c("Cd4","Foxp3")
DotPlot(combined2_Immu,group.by = c("orig.ident","celltype"),features = markers)+coord_flip()
DotPlot(combined_T_PBS,features = markers)+coord_flip()
DotPlot(combined_T_high_CD4,features = markers)+coord_flip()

NK <- combined_T[,combined_T@meta.data$celltype == "NK"]
library(ggpubr)
combined_T@meta.data$orig.ident <- factor(combined_T@meta.data$orig.ident, levels = c("PBS","high_CD4"))
comparisons = list(c("high_CD4","PBS" ))
p1 <- VlnPlot(subset(combined_T,Cd69>=0),features = "Cd69",group.by = "orig.ident",y.max = 3) +
  stat_compare_means(comparisons = comparisons) + labs(title = "T_Cd69")+geom_boxplot(width=.2,col="black",fill="white")+ 
  NoLegend()
p1




#---------------平均表达值热图------------------------------------------

markers <- c("Cd80","Cd40","Cd86","Ebi3","Dhrs9","Ccr7","Ido1","Lamp3","Ccl19","Pdl1")
DotPlot(combined_Immu_Avg,features = markers)+coord_flip()


combined_Immu_nogama <- combined_Immu[, combined_Immu@meta.data$celltype != "Other T"]
#8.13修改至全细胞
keep <-c("B&DC cells","Monocyte","Macrophage")
combined_APC_Avg <- combined2_Immu[,combined2_Immu$celltype %in% keep]


dat=AverageExpression(combined_APC_Avg, group.by = c("orig.ident","celltype"))$RNA
test <- dat[apply(dat, 1, function(x) sd(x)!=0),]
colnames(test)
new_order <- c("PBS_Macrophage","high_CD4_Macrophage",
               "PBS_Monocyte","high_CD4_Monocyte",
               "PBS_B&DC cells","high_CD4_B&DC cells")  # 若使用列索引
# 使用 new_order 重新排列列
test_sorted <- test[, new_order]

#gene.set=unique(c("Tnf","Nos2","Il1b","Il12a","Cd86","Ccl2","Ccl5","Ptgs2",
                  "Cd14","Cxcl2",
                  "Cxcl3","Cd40","Ebi3","Dhrs9",
                  "Cd80","Cxcl5","Cxcl9","Ccr7","Ido1",
                  "Cd69","Gzma",
                  "Cxcr3","Cxcl10","Prf1"))

#gene.set=unique(c("Tnf","Nos2","Il1b","Il12a","Cd86","Ccl2","Ccl5","Ptgs2",
#                  "Cd14","Il1b",
#                 "Cd86","Cd40","Ebi3","Dhrs9",
#                  "Cd80","Cd86","Cd40","Ccr7","Ido1",
#                  "Cd69","Gzma","Prf1",
#                  "Cxcr3","Gzma","Prf1"))


gene.set=unique(c("Cd74","H2-Ab1","H2-Aa","H2-Eb1"))

#test_CD4 <- test_sorted[,c("PBS_CD4 T","high_CD4_CD4 T","high_CD4_NK")]
#gene.set=unique(c("Ctla4","Havcr2","Lag3","Pdcd1","Tigit"))
#gene.set=unique(c("Cd69","Icos","Cd226","Tnfrsf14","Tnfrsf9","Cd28"))

#new_order <- c("PBS_CD8 T","high_CD4_CD8 T","PBS_CD4 T",
#              "high_CD4_CD4 T","PBS_NK","high_CD4_NK",
#               "PBS_Other T","high_CD4_Other T")
#test_sorted <- test[, new_order]
p1<-pheatmap(test_sorted[gene.set,], border_color = "white",
         scale = 'row',
         clustering_method = "ward.D2",
         cluster_rows = F,cluster_cols = F,
         show_rownames = T,
         cellwidth = 20,cellheight = 20)
         #breaks = c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01)),
         #gaps_col = c(2,4,6,8,10),  
         #main="costimulatory molecular genes",
         #legend_breaks=seq(-3,3,1))
p1
p2<-pheatmap(test_CD4[gene.set,], border_color = "white",
             scale = 'row',
             clustering_method = "ward.D2",
             cluster_rows = F,cluster_cols = F,
             show_rownames = T,
             cellwidth = 20,cellheight = 20,
             breaks = c(seq(-1,-0.1,by=0.00001),seq(0,1,by=0.00001)),
             #gaps_col = c(2,4,6,8,10),  
             main="cytotoxicity associated genes",
             legend_breaks=seq(-1,1,0.5))
p2


annotation_col = data.frame( CellType = factor(rep(c("CD4 T", "Macrophage","Monocyte","B&DC cells",
                                                     "CD8 T","NK"), 2)), Time = 1:2 )
rownames(annotation_col) = paste("Test", 1:12, sep = "")
ann_colors = list(CellType = c("CD4 T" = "#1B9E77", "Macrophage" = "#D95F02", "Monocyte"="#7570B3",
                               "B&DC cells"="#E7298A", "CD8 T"= "#66A61E") )
pheatmap(test_sorted, annotation_col = annotation_col,
          annotation_colors = ann_colors, main = "Title")


'HLA-DRB1','HLA-DQB1','HLA-DRA','HLA-DQA1','HLA-DRB5'





grouped_counts_all <- combined2_named@meta.data %>%
  group_by(orig.ident, celltype) %>%
  summarize(count = n())
FeaturePlot(combined2_named,features = c("Fscn1","Ccr7","Ccl22","Tnfrsf9","Sema7a","Stat4","Il12b"))


#---------------CD4 T 热图------------------------------------

keep <-c("CD4 T")
combined_CD4T <- combined_T[,combined_T$celltype %in% keep]


dat=AverageExpression(combined_CD4T, group.by = c("orig.ident","celltype"))$RNA
test <- dat[apply(dat, 1, function(x) sd(x)!=0),]
colnames(test)
new_order <- c("PBS_CD4 T","high_CD4_CD4 T")  
# 若使用列索引
# 使用 new_order 重新排列列
test_sorted <- test[, new_order]

gene.set=unique(c("Foxp3","Cd74","H2-Ab1","H2-Aa","H2-Eb1"))

#gene.set=unique(c("Tnf","Nos2","Il1b","Il12a","Cd86","Ccl2","Ccl5","Ptgs2",
#                  "Cd14","Il1b",
#                 "Cd86","Cd40","Ebi3","Dhrs9",
#                  "Cd80","Cd86","Cd40","Ccr7","Ido1",
#                  "Cd69","Gzma","Prf1",
#                  "Cxcr3","Gzma","Prf1"))


#gene.set=unique(c("Il1a","Il1b","Tnf","Ccl2","Ccl5","Ptgs2","Cd68","Itgam","Cxcl10","Cd274"))

#test_CD4 <- test_sorted[,c("PBS_CD4 T","high_CD4_CD4 T","high_CD4_NK")]
#gene.set=unique(c("Ctla4","Havcr2","Lag3","Pdcd1","Tigit"))
#gene.set=unique(c("Cd69","Icos","Cd226","Tnfrsf14","Tnfrsf9","Cd28"))

#new_order <- c("PBS_CD8 T","high_CD4_CD8 T","PBS_CD4 T",
#              "high_CD4_CD4 T","PBS_NK","high_CD4_NK",
#               "PBS_Other T","high_CD4_Other T")
#test_sorted <- test[, new_order]
p1<-pheatmap(test_sorted[gene.set,], border_color = "white",
             scale = 'row',
             clustering_method = "ward.D2",
             cluster_rows = F,cluster_cols = F,
             show_rownames = T,
             cellwidth = 20,cellheight = 20,
             breaks = c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01)),
             #gaps_col = c(2,4,6,8,10),  
             #main="costimulatory molecular genes",
             legend_breaks=seq(-2,2,1))
p1
p2<-pheatmap(test_CD4[gene.set,], border_color = "white",
             scale = 'row',
             clustering_method = "ward.D2",
             cluster_rows = F,cluster_cols = F,
             show_rownames = T,
             cellwidth = 20,cellheight = 20,
             breaks = c(seq(-1,-0.1,by=0.00001),seq(0,1,by=0.00001)),
             #gaps_col = c(2,4,6,8,10),  
             main="cytotoxicity associated genes",
             legend_breaks=seq(-1,1,0.5))
p2




#---------------计算细胞占比------------------------------------------
#百分比堆叠柱状图
celltypes <- c("CD4 T","CD 8T","NK","Other T","B&DC cells","Monocyte","Macrophage")
counts_high_CD4 <- c(3950, 951, 56,148,67,2110,1492) # 样本A中各类细胞的数量
counts_PBS <- c(11284,953,178,580,68,2920,3874) # 样本B中各类细胞的数量

# 构建数据框
df <- data.frame(
  celltype = rep(celltypes), # 每种细胞类型重复两次，对应两个样本
  sample = c(rep("high_CD4", length(celltypes)), rep("PBS", length(celltypes))), # 样本标识
  count = c(counts_high_CD4, counts_PBS) # 各个样本中细胞的数量
)

# 首先按celltype和sample统计细胞数，并计算每个样本中细胞类型的百分比
df_total <- df %>%
  group_by(celltype) %>%
  summarize(total_count = sum(count))
df_joined <- df_sum %>%
  left_join(df_total, by = "celltype")
df_percent <- df_joined %>%
  mutate_at(vars(matches("count")), list(percent = ~ ./total_count * 100),
            .names = "{col}_percent") %>%
  select(-total_count)
df_plot <- df_percent %>%
  pivot_wider(names_from = "sample", values_from = c(count, count_percent))

df_plot_percent <- df_plot %>%
  pivot_longer(cols = c(count_PBS, count_high_CD4), names_to = "sample", values_to = "count") %>%
  mutate(sample = case_when(
    sample == "count_PBS" ~ "PBS",
    sample == "count_high_CD4" ~ "High CD4"
  ))
ggplot(df_plot_percent, aes(x = celltype, y = count, fill = sample)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Cell Type", y = "Percentage", fill = "Sample") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#同一样本中细胞占比
table(combined_Immu$orig.ident)#查看各组细胞数
prop.table(table(Idents(combined_Immu)))
table(Idents(combined_Immu), combined_Immu$orig.ident)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(combined_Immu), combined_Immu$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)

allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
library(ggplot2)


ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  scale_fill_manual(values = allcolour)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))


initial_ratio <- 0
Cellratio_processed <- Cellratio %>%
  group_by(Var1) %>%
  arrange(Var2) %>%
  mutate(net_change = c(initial_ratio, diff(Freq)))


#---------------计算细胞占比2----------------------------------
status=c('Untreated','Apt-CS')
Tumor_cells=c(56.83,45.02)
Stromal_cell=c(5.69,11.48)
Immune_cell=c(37.48,43.50)
#monocyte=c(0.1471,0.2405)
#macrophage=c(0.1951,0.1700)
#BDC=c(0.0034,0.0076)
#CD4T=c(0.0603,0.1216)
#CD8T=c(0.8983,0.4189)
#OtherT=c(0.0224,0.2230)
#NK=c(0.0190,0.2365)
#创建数据框
data=data.frame(status=status,"Stromal_cell"=Stromal_cell,"Immune_cell"=Immune_cell,"Tumor_cells"=Tumor_cells)
data$Status <- factor(data$status, c('Untreated','Apt-CS'))
#ggplot2画图需要宽数据变成长数据 
melt.data <- melt(data, variable.name = 'Cell', value.name = 'ratio')
#配色
colors=c("#EC8D5A","#ED6F6A","#5E95CF")
#colors=c("#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
 #           "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
  #          "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
   #         "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
#绘制柱状堆叠图的大致轮廓
p=ggplot(melt.data ,aes(x = Status, y = ratio, fill = Cell)) + 
  geom_bar(stat="identity",width = 0.6) + 
  scale_fill_manual(values =colors ) + #添加柱状堆叠图颜色
  theme_bw() + 
  theme(axis.text = element_text(colour = 'black'), 
        axis.text.x = element_text(angle = 0,size=15), 
        axis.title.y = element_text(vjust = 1,size=15),
        axis.ticks.length = unit(1,'mm'),
        #panel.grid = element_blank() , #设置主题背景
        #panel.border = element_blank(),
  plot.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.border = element_blank(),
strip.text.y = element_blank(), 
strip.text = element_text(size=10,face="bold"),
strip.background = element_blank(),
axis.line = element_line(color = 'black'))+
  labs(y="Percentage(%)",x=" ",fill="cell type")
p



#---------------进一步聚类T-----------------------------------------
combined_T <- combined2_named[,combined2_named$celltype %in% ("T cells")]
save(combined_T,file = "combined_T.rData")

#高变基因
combined_T <- FindVariableFeatures(combined_T, selection.method = "vst",nfeatures = 2000)
top10 <- head(VariableFeatures(combined_T),10)
plot1<-VariableFeaturePlot(combined_T)
plot2<-LabelPoints(plot=plot1,points=top10, repel =TRUE)
CombinePlots(plots = list(plot1,plot2))
#归一化，选取的高变基因进行scale
combined_T <- ScaleData(combined_T, features = VariableFeatures(combined_T))

#PCA分析
combined_T <- RunPCA(combined_T, verbose = T)
print(combined_T[["pca"]], dims = 1:5, nfeatures = 5) 

library(harmony)
combined_T  <- RunHarmony(combined_T, group.by.vars = "orig.ident")

DimPlot(combined_T, reduction = "pca")
DimHeatmap(combined_T, dims = 1:2, cells = 500, balanced = TRUE)
ElbowPlot(combined_T)


#细胞聚类分析
combined_T <- FindNeighbors(combined_T, dims = 1:25)
library(clustree)
combined_T <- FindClusters(object = combined_T,resolution = c(seq(.1,1.6,.45)))
clustree(combined_T@meta.data, prefix = "RNA_snn_res.")

combined_T <- FindClusters(object = combined_T,resolution = 1)
combined_T <- RunUMAP(combined_T,dims = 1:25)
DimPlot(combined_T, reduction = "umap",label=TRUE)

markers <- c("Cd4","Cd8a","Cd8b1","Cd3d","Cd3e","Nkg7","Ncr1","Cd69","Ptprc")
DotPlot(combined_T,features = markers)+coord_flip()
VlnPlot(combined_T, features = c("Trgc1","Trdv2"))

new.cluster.ids <- c("0"="CD8 T","1"="CD8 T","2"="CD8 T","3"="CD8 T","4"="CD4 T",
                     "6"="Other T","5"="NK")
combined_T <- RenameIdents(combined_T, new.cluster.ids) 
combined_T$celltype <- combined_T@active.ident
DimPlot(combined_T, group.by = "celltype",label = T)

colors <- c("#F79581","#94BEE4","#F69CC7","#A989C1",
            "grey","grey","grey","grey","grey",
            "grey","grey","grey","grey","grey")

colors <- c("#5CBDC4","#FA6666","#E5AECC","#D9E077",
            "grey","grey","grey","grey",
            "grey","grey","grey","grey","grey")
DimPlot(combined_T, group.by = "celltype",label = F, cols = colors)

combined_T_TSNE <- RunTSNE(combined_T,dims = 1:25)
colors <- c("#E9A86D","#FA6666","#D8A7DD","#70B073",
            "grey","grey","grey","grey",
            "grey","grey","grey","grey","grey")
DimPlot(combined_T_TSNE, reduction = "tsne", label = F,pt.size = 1.8, cols = colors)


save(combined_T,file = "combined_T_clusternamed.rData")

combined_T_high_CD4 <- combined_T_TSNE[,(combined_T_TSNE$orig.ident %in% c("high_CD4"))]
DimPlot(combined_T_high_CD4, group.by = "celltype", cols = colors)+
  ggtitle("high CD4")
combined_T_PBS <- combined_T_TSNE[,(combined_T_TSNE$orig.ident %in% c("PBS"))]
DimPlot(combined_T_PBS, group.by = "celltype", cols = colors)+
  ggtitle("PBS")

#combined_T_high_CD4 <- RunTSNE(combined_T_high_CD4,dims = 1:25)
#colors <- c("#8BC6D8","#FA6666","#D9E077","#E5AECC",
#            "grey","grey","grey","grey",
#            "grey","grey","grey","grey","grey")
#DimPlot(combined_T_PBS, reduction = "tsne", label = F,pt.size = 2, cols = colors)

FeaturePlot(combined_T,features = c("Cd4","Cd8a","Ifng","Gzmb","Tnf","Gata3","Rora","Cd69","Ptprc"))


#---------------CD4CD69占比-------------------------------
FeaturePlot(combined2_named,features = c("Cd69"))

combined_T_high_CD4_TSNE <- RunTSNE(combined_T_high_CD4,dims = 1:25)
combined_T_PBS_TSNE <- RunTSNE(combined_T_PBS,dims = 1:25)

FeaturePlot(combined_T_TSNE,reduction = "tsne",features = c("Cd4","Cd69"))
FeaturePlot(combined_T_high_CD4_TSNE,reduction = "tsne",features = c("Cd4","Cd74","H2-Ab1","H2-Aa","H2-Eb1","H2-Eb2"))
FeaturePlot(combined_T_PBS_TSNE,reduction = "tsne",features = c("Cd4","Cd74","H2-Ab1","H2-Aa","H2-Eb1","H2-Eb2"))
markers <- c("Cd4","Foxp3")
DotPlot(combined_T,group.by = c("orig.ident","celltype"),features = markers)+coord_flip()
DotPlot(combined_T_PBS,features = markers)+coord_flip()
DotPlot(combined_T_high_CD4,features = markers)+coord_flip()


#创建一个新的seurat文件，其中仅包含免疫细胞细胞
keep <-c("B&DC cells","T cells","Monocyte","Macrophage")
combined2_Immu <- combined2_named[,combined2_named$celltype %in% keep]
#combined2_Immu <- merge(combined2_Immu, y = combined_T)
#save(combined2_Immu,file = "combined2_Immu_renamedT.rData")
FeaturePlot(combined2_Immu,features = c("Cd4","Cd69"))

combined2_Immu_high_CD4 <- combined2_Immu[,(combined2_Immu$orig.ident %in% c("high_CD4"))]
combined2_Immu_PBS <- combined2_Immu[,(combined2_Immu$orig.ident %in% c("PBS"))]
FeaturePlot(combined2_Immu_high_CD4,features = c("Cd4","Cd69"))
FeaturePlot(combined2_Immu_PBS,features = c("Cd4","Cd69"))


#---------------计算某个基因表达的几种图-----------------
levels(combined2_named@meta.data$orig.ident) <- c("PBS", "high_CD4")
VlnPlot(combined2_Immu, features = c("Cd4"),
        split.by = c("celltype"),y.max = 4)
combined2_named_PBS <- combined2_named[,combined2_named$orig.ident %in% "PBS"]
combined2_named_high_CD4 <- combined2_named[,combined2_named$orig.ident %in% "high_CD4"]

combined_Tumor_high_CD4 <- combined_Tumor[,(combined_Tumor$orig.ident %in% c("high_CD4"))]
combined_Tumor_PBS <- combined_Tumor[,(combined_Tumor$orig.ident %in% c("PBS"))]
DimPlot(combined_Tumor_high_CD4, group.by = "celltype", cols = colors)+ggtitle("Post")
DimPlot(combined_Tumor_PBS, group.by = "celltype", cols = colors)+ggtitle("Pre")

FeaturePlot(combined_Tumor_high_CD4,features = c("Cd74","H2-Ab1","H2-Aa","H2-Eb1","H2-Eb2"))
FeaturePlot(combined_Tumor_PBS,features = c("Cd74","H2-Ab1","H2-Aa","H2-Eb1","H2-Eb2"))

markers <- c("Cd74","H2-Ab1","H2-Aa","H2-Eb1","H2-Eb2")
DotPlot(combined_Tumor,group.by = "orig.ident",features = markers)+coord_flip()
DotPlot(combined_T_high_CD4,features = markers)+coord_flip()



# #-------------- 确定你感兴趣的聚类编号
cluster_of_interest <- "0"  # 例如聚类 0

# 获取基因表达矩阵
expr_matrix <- GetAssayData(object = combined2_named, layer = "counts")

# 筛选出特定聚类的细胞索引
#cells_in_cluster <- which(combined2_named@meta.data$celltype == c("T cells","Tumor cells","Stromal cells",
                                                                  "Monocyte","Macrophage","B&DC cells"))

# 提取特定基因（例如 "Cd4"）的表达值
expr_Cd69 <- expr_matrix["Cd69",]

# 获取 orig.ident 列
orig_ident_values <- combined2_named@meta.data$orig.ident
# 获取 celltype 列
celltype_values1 <- combined2_named@meta.data$celltype
celltype_values <- combined_T@meta.data$celltype

# 将表达值和 orig.ident 保存为数据框
expr_Cd69_df <- data.frame(
  cell_id = colnames(expr_matrix),
  expression = as.vector(expr_Cd69),
  orig.ident = orig_ident_values,
  celltype1 = celltype_values1

)


# 保存为 CSV 文件
write.csv(expr_Cd69_df, file = "Cd69_expression_in_T_cells.csv")

cluster_of_interest <- "0"  # 例如聚类 0

# 获取基因表达矩阵
expr_matrix <- GetAssayData(object = combined2_named, layer = "counts")

# 筛选出特定聚类的细胞索引
cells_in_cluster <- which(combined2_named@meta.data$celltype == "T cells")

# 提取特定基因（例如 "Cd4"）的表达值
expr_Cd69 <- expr_matrix["Cd69", cells_in_cluster]

# 获取 orig.ident 列
orig_ident_values <- combined2_named@meta.data$orig.ident
# 获取 celltype 列
celltype_values <- combined2_named@meta.data$celltype

# 将表达值和 orig.ident 保存为数据框
expr_Cd69_df <- data.frame(
  cell_id = colnames(expr_matrix)[cells_in_cluster],
  expression = as.vector(expr_Cd69),
  orig.ident = orig_ident_values,
  celltype = celltype_values
)

# 保存为 CSV 文件
write.csv(expr_Cd69_df, file = "Cd69_expression_in_T_cells.csv")


#创建一个新的seurat文件，其中T细胞被详细注释
keep <-c("Tumor cells","B&DC cells","Stromal cells","Monocyte","Macrophage")
combined2_NoT <- combined2_named[,combined2_named$celltype %in% keep]
combined2_renamedT <- merge(combined_T, y = combined2_NoT)

#进行Normalize归一化后
normalized_data <- GetAssayData(object = combined2_renamedT, slot = "data")
gene_expression <- normalized_data["Cd69", ]
orig_ident_values <- combined2_renamedT@meta.data$orig.ident
celltype_values <- combined2_renamedT@meta.data$celltype
gene_expression_df <- data.frame(cell_id = colnames(normalized_data), 
                                 expression = as.vector(gene_expression),
                                 orig.ident = orig_ident_values,
                                 celltype = celltype_values)

VlnPlot(combined2_Immu, features = c("Cd4"),split.by = c("orig.ident"),y.max = 5)
write.csv(gene_expression_df, file = "Cd69_expression_in_all_cells.csv")
#---------------进一步聚类B------------------------------------
combined_BDC <- combined2_named[,combined2_named$celltype %in% ("B&DC cells")]
save(combined_BDC,file = "combined_BDC.rData")

#高变基因
combined_BDC <- FindVariableFeatures(combined_BDC, selection.method = "vst",nfeatures = 2000)
top10 <- head(VariableFeatures(combined_BDC),10)
plot1<-VariableFeaturePlot(combined_BDC)
plot2<-LabelPoints(plot=plot1,points=top10, repel =TRUE)
CombinePlots(plots = list(plot1,plot2))
#归一化，选取的高变基因进行scale
combined_BDC <- ScaleData(combined_BDC, features = VariableFeatures(combined_BDC))

#PCA分析
combined_BDC <- RunPCA(combined_BDC, verbose = T)
print(combined_BDC[["pca"]], dims = 1:5, nfeatures = 5) 

library(harmony)
combined_BDC  <- RunHarmony(combined_BDC, group.by.vars = "orig.ident")

DimPlot(combined_BDC, reduction = "pca")
DimHeatmap(combined_BDC, dims = 1:2, cells = 500, balanced = TRUE)
ElbowPlot(combined_BDC)


#细胞聚类分析
combined_BDC <- FindNeighbors(combined_BDC, dims = 1:25)
library(clustree)
combined_BDC <- FindClusters(object = combined_BDC,resolution = c(seq(.1,1.6,.45)))
clustree(combined_BDC@meta.data, prefix = "RNA_snn_res.")

combined_BDC <- FindClusters(object = combined_BDC,resolution = 1.5)
combined_BDC <- RunUMAP(combined_BDC,dims = 1:25)
DimPlot(combined_BDC, reduction = "umap",label=TRUE)

markers <- c("Cd4","Cd79a","Iglc2","Jchain","Xcr1","Siglech","Mctp2","Fscn1", "Ccr7","Lilrb4", "Itgax")

markers <- c("Cd4")

DotPlot(combined2_renamedT,features = markers)+coord_flip()
VlnPlot(combined2_renamedT, features = c("Cd4"))



new.cluster.ids <- c("2"="CD4+")
combined_BDC <- RenameIdents(combined_BDC, new.cluster.ids) 
combined_BDC$celltype <- combined_BDC@active.ident
DimPlot(combined_BDC, group.by = "celltype",label = T)
save(combined_BDC,file = "combined_BDC_clusternamed.rData")

combined_BDC_high_CD4 <- combined_BDC[,(combined_BDC$orig.ident %in% c("high_CD4"))]
DimPlot(combined_BDC_high_CD4, group.by = "celltype")+
  ggtitle("high CD4")
combined_BDC_PBS <- combined_BDC[,(combined_BDC$orig.ident %in% c("PBS"))]
DimPlot(combined_BDC_PBS, group.by = "celltype")+
  ggtitle("PBS")
#---------------计算比值-----------------------------------------
grouped_counts_all <- combined2_named@meta.data %>%
  group_by(orig.ident, celltype) %>%
  summarize(count = n())

# 提取样本1和样本2中cell1和cell2的数量
high_CD4_CD4_count <- grouped_counts %>%
  filter(orig.ident == "high_CD4" & celltype == "CD4 T") %>%
  pull(count)
high_CD4_CD8_count <- grouped_counts %>%
  filter(orig.ident == "high_CD4" & celltype == "CD8 T") %>%
  pull(count)

PBS_CD4_count <- grouped_counts %>%
  filter(orig.ident == "PBS" & celltype == "CD4 T") %>%
  pull(count)
PBS_CD8_count <- grouped_counts %>%
  filter(orig.ident == "PBS" & celltype == "CD8 T") %>%
  pull(count)

# 计算两个样本中cell1/cell2的比例
ratio_high_CD4 <- high_CD4_CD4_count / high_CD4_CD8_count
ratio_PBS <- PBS_CD4_count / PBS_CD8_count

# 输出比例
cat("样本1中cell1/cell2的比例是：", ratio_high_CD4, "\n")
cat("样本2中cell1/cell2的比例是：", ratio_PBS, "\n")


df_ratios <- data.frame(Sample = c("high_CD4", "PBS"),
                        Ratio = c(ratio_high_CD4, ratio_PBS))
ggplot(df_ratios, aes(x = Sample, y = Ratio)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.5) +
  geom_text(aes(label = round(Ratio, 4)), vjust = -0.5) +
  labs(title = "CD4T/CD8T",
       x = "Sample",
       y = "Ratio")



#---------------计算CD4\CD8 in all--------------------------
# 提取样本1和样本2中cell1和cell2的数量
high_CD4_CD4_count <- grouped_counts %>%
  filter(orig.ident == "high_CD4" & celltype == "CD4 T") %>%
  pull(count)
high_CD4_CD8_count <- grouped_counts %>%
  filter(orig.ident == "high_CD4" & celltype == "CD8 T") %>%
  pull(count)
high_CD4_NK_count <- grouped_counts %>%
  filter(orig.ident == "high_CD4" & celltype == "NK") %>%
  pull(count)
high_CD4_BDC_count <- grouped_counts_all %>%
  filter(orig.ident == "high_CD4" & celltype == "B&DC cells") %>%
  pull(count)
high_CD4_Mac_count <- grouped_counts_all %>%
  filter(orig.ident == "high_CD4" & celltype == "Macrophage") %>%
  pull(count)
high_CD4_Mono_count <- grouped_counts_all %>%
  filter(orig.ident == "high_CD4" & celltype == "Monocyte") %>%
  pull(count)

PBS_CD4_count <- grouped_counts %>%
  filter(orig.ident == "PBS" & celltype == "CD4 T") %>%
  pull(count)
PBS_CD8_count <- grouped_counts %>%
  filter(orig.ident == "PBS" & celltype == "CD8 T") %>%
  pull(count)
PBS_NK_count <- grouped_counts %>%
  filter(orig.ident == "PBS" & celltype == "NK") %>%
  pull(count)
PBS_BDC_count <- grouped_counts_all %>%
  filter(orig.ident == "PBS" & celltype == "B&DC cells") %>%
  pull(count)
PBS_Mac_count <- grouped_counts_all %>%
  filter(orig.ident == "PBS" & celltype == "Macrophage") %>%
  pull(count)
PBS_Mono_count <- grouped_counts_all %>%
  filter(orig.ident == "PBS" & celltype == "Monocyte") %>%
  pull(count)

# 计算两个样本中cell1/CD3T的比例
ratio_high_CD4 <- high_CD4_CD8_count / 113
ratio_PBS <- PBS_CD8_count / 569

# 计算两个样本中cell1/All immu的比例
ratio_high_CD4 <- high_CD4_Mono_count / 3817
ratio_PBS <- PBS_Mono_count / 7442

# 计算两个样本中cell1/All的比例
ratio_high_CD4 <- high_CD4_Mono_count / 8774
ratio_PBS <- PBS_Mono_count / 19857

# 输出比例
cat("样本1中cell1/cell2的比例是：", ratio_high_CD4, "\n")
cat("样本2中cell1/cell2的比例是：", ratio_PBS, "\n")

df_ratios <- data.frame(Sample = c("PBS", "high_CD4"),
                        Ratio = c(ratio_PBS, ratio_high_CD4))
ggplot(df_ratios, aes(x = Sample, y = Ratio)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.5) +
  geom_text(aes(label = round(Ratio, 6)), vjust = -0.5) +
  labs(title = "Mono/ALL",
       x = "Sample",
       y = "Ratio")+
  theme_classic()

Sample=c('Untreated','Apt-CS')
Ratio_CD4CD8=c(6.72,29.03)
Ratio_CD4CD3=c(6.15,15.93)
Ratio_NK=c(0.15,0.92)
Ratio_BDC=c(0.91,1.76)
Ratio_Mono=c(35.24,55.28)
#创建数据框
data=data.frame(Sample=Sample,Ratio_Mono=Ratio_Mono)
data$Sample <- factor(data$Sample, c('Untreated','Apt-CS'))
#ggplot2画图需要宽数据变成长数据 
melt.data <- melt(data, variable.name = 'Cell', value.name = 'ratio')
p<-ggplot(melt.data, aes(x = Sample, y = Ratio_Mono)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.5) +
  geom_text(aes(label = round(Ratio_Mono, 6)), vjust = -0.5) +
  labs(title = "Mono/Immune cell",
       x = "Sample",
       y = "Ratio")+
  theme_classic()
p

#---------------细胞通讯-----------------------------------------------

library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(SeuratData)
options(stringsAsFactors = FALSE)

combined_high_CD4 <- combined2_named[,(combined2_named$orig.ident %in% c("high_CD4"))]
DimPlot(combined_high_CD4, group.by = "celltype")+
  ggtitle("high CD4")
combined_PBS <- combined2_named[,(combined2_named$orig.ident %in% c("PBS"))]
DimPlot(combined_PBS, group.by = "celltype")+
  ggtitle("PBS")

table(combined2_named$celltype)
Idents(combined2_named) <- 'celltype'
combined2 <- subset(combined2_named,idents=c("Tumor cells", "Stromal cells","T cells",
                                             "B&DC cells","Monocyte","Macrophage") )
combined2$celltype <-as.factor(as.character(combined2$celltype))
table(combined2$orig.ident)                    
Idents(combined2) <- 'orig.ident'                    
high_CD4 <- subset(combined2,idents="high_CD4")                    
PBS <- subset(combined2,idents="PBS")                     
                    
cco.high_CD4 <- createCellChat(high_CD4@assays$RNA@data,meta = high_CD4@meta.data,group.by = "celltype")                    
cco.PBS <- createCellChat(PBS@assays$RNA@data,meta = PBS@meta.data,group.by = "celltype")                    
save(cco.high_CD4,cco.PBS,file = "cco.rda")
                    
cellchat <- cco.high_CD4
cellchat@DB <- CellChatDB.mouse
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat,PPI.mouse)
cellchat <- computeCommunProb(cellchat,raw.use = TRUE,population.size = TRUE)
cellchat <- filterCommunication(cellchat,min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat,slot.name = "netP")
cellchat <- computeNetSimilarity(cellchat,type = "functional")
#cellchat <- netEmbedding(cellchat,type = "functional")
#cellchat <- netClustering(cellchat,type = "functional")
cellchat <- computeNetSimilarity(cellchat,type = "structural")                    
#cellchat <- netEmbedding(cellchat,type = "structural")
#cellchat <- netClustering(cellchat,type = "structural")           
cco.high_CD4 <- cellchat
saveRDS(cco.high_CD4,"cco.high_CD4.rds")
#修改high_CD4为PBS重复上述行

cco.list <- list(Pre=cco.PBS,Post=cco.high_CD4)
cellchat <- mergeCellChat(cco.list,add.names = names(cco.list),cell.prefix = TRUE)

#柱状图对比相互作用数量和强度
gg1 <- compareInteractions(cellchat, show.legend = F,group = c(1,2),measure = "count")
gg2 <- compareInteractions(cellchat,show.legend = F,group = c(1,2),measure = "weight")
p <- gg1 + gg2
p
ggsave("Overview_number_strength.pdf",p,width = 6,height = 4)

#两组对比后的数量和强度差异网络图
par(mfrow=c(1,2))
netVisual_diffInteraction(cellchat,weight.scale = T,edge.width.max = 25) #edge.width.max 修改线的粗细
netVisual_diffInteraction(cellchat,weight.scale = T,measure = "weight",edge.width.max = 25)


par(mfrow=c(1,2))
weight.max <- getMaxWeight(cco.list,attribute = c("idents","count"))
for (i in 1:length(cco.list)) {
  netVisual_circle(cco.list[[i]]@net$count,weight.scale = T,label.edge = F,
                   edge.weight.max = weight.max[2],edge.width.max = 12,
                   title.name = paste0("Number of interactions - ",names(cco.list)[i]))
}

#信号通路
gg1 <- rankNet(cellchat,mode = "comparison",stacked = T,do.stat = TRUE)
gg2 <- rankNet(cellchat,mode = "comparison",stacked = F,do.stat = TRUE)
p <- gg1 + gg2
p



cellchat.Immu.Tumor <- identifyOverExpressedGenes(cellchat.Immu.Tumor, group.dataset = "datasets", pos.dataset = "Post", 
                                       features.name = "Post", only.pos = FALSE, thresh.pc = 0.1, 
                                       thresh.fc = 0.1, thresh.p = 1) 
net <- netMappingDEG(cellchat.Immu.Tumor, features.name = "Post")
net.up <- subsetCommunication(cellchat.Immu.Tumor, net = net, datasets = "Post",
                              ligand.logFC = 0.1,receptor.logFC = NULL)
#上调通路的弦图
netVisual_chord_gene(cco.high_CD4_Immu_Tumor, sources.use = 2, targets.use = c(1:7), slot.name = 'netP', 
                     net = net.up,lab.cex = 0.8, small.gap = 1, 
                     title.name = paste0("Up-regulated signaling in Apt-CS")) 
netVisual_chord_cell(cco.high_CD4_Immu_Tumor, sources.use = 2, targets.use = c(1:7), slot.name = 'netP', 
                     net = net.up,lab.cex = 0.8, small.gap = 3.5, 
                     title.name = paste0("Up-regulated signaling in LS"))
#所有通路的弦图
netVisual_chord_gene(cco.high_CD4_Immu_Tumor, sources.use = 2, targets.use = c(1:7), 
                     lab.cex = 0.8,
                     small.gap = 1,
                     slot.name = "netP",title.name = "Chord diagram  2: show pathway") 





#某一具体通路
df.net <- subsetCommunication(cellchat)
df.netP <- subsetCommunication(cellchat,slot.name = "netP")
pathways.show <- unique(df.netP$Post$pathway_name)

pathways.show <- c("PD-L1")
weight.max <- getMaxWeight(cco.list,slot.name = c("netP"),attribute = pathways.show)
par(mfrow=c(1,2),xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_aggregate(cco.list[[i]],signaling = pathways.show,layout = "circle",
                      edge.weight.max = weight.max[1],edge.width.max = 10,
                      signaling.name = paste(pathways.show,names(cco.list)[i]))
}
#某一具体通路的和弦图
par(mfrow=c(1,2),xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_aggregate(cco.list[[i]],signaling = pathways.show,layout = "chord",
                      pt.title = 3, title.space = 0.05, vertex.label.cex = 0.6,
                      signaling.name = paste(pathways.show,names(cco.list)[i]))
}
#
p = plotGeneExpression(cellchat,signaling = "PD-L1")
ggsave("PD-L1_GeneExpression_vln.pdf",p,width = 8,height = 8)
p

netVisual_aggregate(cco.high_CD4_Immu_Tumor, signaling = pathways.show, layout = "chord")

#点图
netVisual_bubble(cellchat.Immu.Tumor, sources.use = 2, targets.use = c(1:7), comparison = c(1, 2), 
                 max.dataset = 2,title.name = "Increased signaling", angle.x = 45,remove.isolate = T )
netVisual_bubble(cellchat.Immu.Tumor, sources.use = 2, targets.use = c(1,5,6), comparison = c(1, 2), 
                 max.dataset = 2,angle.x = 45,remove.isolate = T )


#和弦图
par(mfrow=c(1,2),xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_chord_gene(cco.list[[i]],sources.use = 6,targets.use = c(1:7),
                       lab.cex = 0.6,legend.pos.x = 10,legend.pos.y = 2)
  
}


#---------------单样本--------------------
par(mfrow=c(1,2))
weight.max <- getMaxWeight(cco.list,attribute = c("idents","count"))
for (i in 1:length(cco.list.Immu.Tumor)) {
  netVisual_circle(cco.list.Immu.Tumor[[i]]@net$count,weight.scale = T,label.edge = F,
                   edge.weight.max = weight.max[2],edge.width.max = 12,
                   title.name = paste0("Number of interactions - ",names(cco.list.Immu.Tumor)[i]))
}

groupSize <- as.numeric(table(cco.high_CD4@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cco.high_CD4@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions of high CD4")
netVisual_circle(cco.high_CD4@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength of high CD4")
groupSizePBS <- as.numeric(table(cco.PBS@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cco.PBS@net$count, vertex.weight = groupSizePBS, weight.scale = T, label.edge= F, title.name = "Number of interactions of PBS")
netVisual_circle(cco.PBS@net$weight, vertex.weight = groupSizePBS, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength of PBS")

#去除tumor
par(mfrow=c(1,2))
s.cell <- c("B&DC cells","Endothelial","Fibroblasts","Macrophage","Monocyte","T cells")
count1 <- cco.list[[1]]@net$count[s.cell,s.cell]
count2 <- cco.list[[2]]@net$count[s.cell,s.cell]
weight.max <- max(max(count1),max(count2))
netVisual_circle(count1,weight.scale = T,label.edge = T,edge.weight.max = weight.max,
                 edge.width.max =12,title.name = paste0("Number of interactions-",names(cco.list)[1])) 
netVisual_circle(count2,weight.scale = T,label.edge = T,edge.weight.max = weight.max,
                 edge.width.max =12,title.name = paste0("Number of interactions-",names(cco.list)[2]))

#每种细胞发出的信号
mat <- cco.high_CD4@net$count
par(mfrow=c(3,3),xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0,nrow=nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,]<-mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,arrow.width = 0.2,
                   arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
  
#某一具体通路
pathways.show <- c("CD40")
par(mfrow=c(1,1))
netVisual_aggregate(cco.high_CD4,signaling = pathways.show,layout = "circle")
p = plotGeneExpression(cco.high_CD4,signaling = "PD-L1")
p
markers <- c("Clec2d")
DotPlot(combined_T_NoNK,features = markers)+coord_flip()
  
#尝试单方向箭头

df.net <- subsetCommunication(cco.high_CD4, sources.use = 6, targets.use = c(1:7))
cell.levels <- levels(factor(df.net$cell.type))  # 或者使用实际列名

netVisual_chord_cell(cco.high_CD4_Immu_Tumor, net = cco.high_CD4_Immu_Tumor@netP,
                     sources.use = 6, targets.use = c(1:7)) 

netVisual_chord_gene(cco.high_CD4,
                     sources.use = c(1,7),
                     targets.use = 6,
                     lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cco.high_CD4_Immu_Tumor, sources.use = 2, targets.use = c(1:7), 
                     lab.cex = 0.8,
                     small.gap = 1,
                     slot.name = "netP",title.name = "Chord diagram  2: show pathway") 
pathways.show <- c("IL1","CD137","CD80")
netVisual_chord_gene(cco.high_CD4, signaling = pathways.show)
pathways.show <- c("IL1")
par(mfrow=c(1,1))
netVisual_aggregate(cco.high_CD4, signaling = pathways.show, layout = "chord")

netVisual_chord_gene(cco.high_CD4, sources.use = c(1,6), targets.use = 1, legend.pos.x = 30)
       
netVisual_chord_cell(cco.high_CD4,signaling = "CD80")

#---------------T细胞亚群通讯---------------------------------
combined_Mono <- combined2_named[,combined2_named$celltype %in% ("Monocyte")]
save(combined_Mono,file = "combined_Mono.rData")
combined_Immu <- merge(combined_T, y = c(combined_Mono,combined_Mac,combined_BDC))
combined_Immu


combined_Immu_high_CD4 <- combined_Immu[,(combined_Immu$orig.ident %in% c("high_CD4"))]
combined_Immu_PBS <- combined_Immu[,(combined_Immu$orig.ident %in% c("PBS"))]


table(combined_Immu$celltype)
Idents(combined_Immu) <- 'celltype'
combined8 <- subset(combined_Immu,idents=c("CD4 T","CD8 T","NK","Other T","B&DC cells","Monocyte","Macrophage") )
combined8$celltype <-as.factor(as.character(combined8$celltype))
table(combined8$orig.ident)   
Idents(combined8) <- 'orig.ident'                    
high_CD4_Immu <- subset(combined8,idents="high_CD4")                    
PBS_Immu <- subset(combined8,idents="PBS")                     

cco.high_CD4_Immu <- createCellChat(high_CD4_Immu@assays$RNA@data,meta = high_CD4_Immu@meta.data,group.by = "celltype")                    
cco.PBS.Immu <- createCellChat(PBS_Immu@assays$RNA@data,meta = PBS_Immu@meta.data,group.by = "celltype")                    
save(cco.high_CD4_Immu,cco.PBS.Immu,file = "cco_Immu.rda")

cellchat <- cco.high_CD4_Immu
cellchat@DB <- CellChatDB.mouse
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat,PPI.mouse)
cellchat <- computeCommunProb(cellchat,raw.use = TRUE,population.size = TRUE)
cellchat <- filterCommunication(cellchat,min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat,slot.name = "netP")
cellchat <- computeNetSimilarity(cellchat,type = "functional")
#cellchat <- netEmbedding(cellchat,type = "functional")
#cellchat <- netClustering(cellchat,type = "functional")
cellchat <- computeNetSimilarity(cellchat,type = "structural")                    
#cellchat <- netEmbedding(cellchat,type = "structural")
#cellchat <- netClustering(cellchat,type = "structural")           
cco.high_CD4_Immu <- cellchat
saveRDS(cco.high_CD4_Immu,"cco.high_CD4_Immu.rds")
#修改high_CD4为PBS重复上述行

cco.list.Immu <- list(Pre=cco.PBS.Immu,Post=cco.high_CD4_Immu)
cellchat.Immu <- mergeCellChat(cco.list.Immu,add.names = names(cco.list.Immu),cell.prefix = TRUE)

#柱状图对比相互作用数量和强度
gg1 <- compareInteractions(cellchat.Immu, show.legend = F,group = c(1,2),measure = "count")
gg2 <- compareInteractions(cellchat.Immu,show.legend = F,group = c(1,2),measure = "weight")
p <- gg1 + gg2
p
ggsave("Overview_number_strength_TBDC.pdf",p,width = 6,height = 4)

#两组对比后的数量和强度差异网络图-做不出来
par(mfrow=c(1,2))
netVisual_diffInteraction(cellchat.Immu,weight.scale = T)
netVisual_diffInteraction(cellchat.Immu,weight.scale = T,measure = "weight")
#方法二
par(mfrow=c(1,2))
weight.max <- getMaxWeight(cco.list.Immu,attribute = c("idents","count"))
for (i in 1:length(cco.list.Immu)) {
  netVisual_circle(cco.list.Immu[[i]]@net$count,weight.scale = T,label.edge = F,
                   edge.weight.max = weight.max[2],edge.width.max = 12,
                   title.name = paste0("Number of interactions - ",names(cco.list.Immu)[i]))
}

#信号通路
gg1 <- rankNet(cellchat.Immu,mode = "comparison",stacked = T,do.stat = TRUE)
gg2 <- rankNet(cellchat.Immu,mode = "comparison",stacked = F,do.stat = TRUE)
p <- gg1 + gg2
p

#某一具体通路
pathways.show <- c("MHC-II")
par(mfrow=c(1,1))
netVisual_aggregate(cco.PBS.TT,signaling = pathways.show,layout = "circle")
p = plotGeneExpression(cco.high_CD4.TBDC,signaling = "CD137")
p

Idents(combined_T) <- 'celltype'
combined_T_NoNK <- subset(combined_T,idents=c("CD4 T","CD8 T","Other T") )
markers <- c("Cd28","Ctla4")
DotPlot(combined_T_NoNK,features = markers)+coord_flip()

netVisual_chord_gene(cco.high_CD4_Immu, sources.use = 2, targets.use = c(1:7), lab.cex = 0.5,
                     slot.name = "netP",title.name = "Chord diagram  2: show pathway") 

#---------------T and Macrophage------------------
combined_Mac <- combined2_named[,combined2_named$celltype %in% ("Macrophage")]
save(combined_Mac,file = "combined_Mac.rData")
combined_T_Mac <- merge(combined_T, y = combined_Mac)
combined_T_Mac

combined_T_Mac_high_CD4 <- combined_T_Mac[,(combined_T_Mac$orig.ident %in% c("high_CD4"))]
combined_T_Mac_PBS <- combined_T_Mac[,(combined_T_Mac$orig.ident %in% c("PBS"))]


table(combined_T_Mac$celltype)
Idents(combined_T_Mac) <- 'celltype'
combined5 <- subset(combined_T_Mac,idents=c("CD4 T","CD8 T","NK","Other T","Macrophage") )
combined5$celltype <-as.factor(as.character(combined5$celltype))
table(combined5$orig.ident)   
Idents(combined5) <- 'orig.ident'                    
high_CD4_TMac <- subset(combined5,idents="high_CD4")                    
PBS_TMac <- subset(combined5,idents="PBS")                     

cco.high_CD4.TMac <- createCellChat(high_CD4_TMac@assays$RNA@data,meta = high_CD4_TMac@meta.data,group.by = "celltype")                    
cco.PBS.TMac <- createCellChat(PBS_TMac@assays$RNA@data,meta = PBS_TMac@meta.data,group.by = "celltype")                    
save(cco.high_CD4.TMac,cco.PBS.TMac,file = "cco_TMac.rda")

cellchat <- cco.PBS.TMac
cellchat@DB <- CellChatDB.mouse
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat,PPI.mouse)
cellchat <- computeCommunProb(cellchat,raw.use = TRUE,population.size = TRUE)
cellchat <- filterCommunication(cellchat,min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat,slot.name = "netP")
cellchat <- computeNetSimilarity(cellchat,type = "functional")
#cellchat <- netEmbedding(cellchat,type = "functional")
#cellchat <- netClustering(cellchat,type = "functional")
cellchat <- computeNetSimilarity(cellchat,type = "structural")                    
#cellchat <- netEmbedding(cellchat,type = "structural")
#cellchat <- netClustering(cellchat,type = "structural")           
cco.PBS.TMac <- cellchat
saveRDS(cco.PBS.TMac,"cco.PBS.TMac.rds")
#修改high_CD4为PBS重复上述行

cco.list.TMac <- list(Pre=cco.PBS.TMac,Post=cco.high_CD4.TMac)
cellchat.TMac <- mergeCellChat(cco.list.TMac,add.names = names(cco.list.TMac),cell.prefix = TRUE)

#柱状图对比相互作用数量和强度
gg1 <- compareInteractions(cellchat.TMac, show.legend = F,group = c(1,2),measure = "count")
gg2 <- compareInteractions(cellchat.TMac,show.legend = F,group = c(1,2),measure = "weight")
p <- gg1 + gg2
ggsave("Overview_number_strength_TMac.pdf",p,width = 6,height = 4)

#两组对比后的数量和强度差异网络图-做不出来
par(mfrow=c(1,2))
netVisual_diffInteraction(cellchat.TMac,weight.scale = T)
netVisual_diffInteraction(cellchat.TMac,weight.scale = T,measure = "weight")
#方法二
par(mfrow=c(1,2))
weight.max <- getMaxWeight(cco.list.TMac,attribute = c("idents","count"))
for (i in 1:length(cco.list.TMac)) {
  netVisual_circle(cco.list.TMac[[i]]@net$count,weight.scale = T,label.edge = F,
                   edge.weight.max = weight.max[2],edge.width.max = 12,
                   title.name = paste0("Number of interactions - ",names(cco.list.TMac)[i]))
}

#信号通路
gg1 <- rankNet(cellchat.TMac,mode = "comparison",stacked = T,do.stat = TRUE)
gg2 <- rankNet(cellchat.TMac,mode = "comparison",stacked = F,do.stat = TRUE)
p <- gg1 + gg2
p

#某一具体通路
pathways.show <- c("MHC-II")
par(mfrow=c(1,1))
netVisual_aggregate(cco.high_CD4.TBDC,signaling = pathways.show,layout = "circle")
p = plotGeneExpression(cco.PBS.TT,signaling = "MHC-II")
p
markers <- c("H2-Aa","H2-Ab1","Cd4")
DotPlot(combined_T_Mac_high_CD4,features = markers)+coord_flip()
netVisual_bubble(cellchat.TBDC, sources.use = 1, targets.use = 2, comparison = c(1, 2), 
                 max.dataset = 2,title.name = "Increased signaling", angle.x = 45,remove.isolate = T )

#---------------T and Tumor------------------
combined_CD4T <- combined_T[,combined_T$celltype %in% ("CD4 T")]
#save(combined_CD8T,file = "combined_CD8T.rData")
combined_CD4T_Tumor <- merge(combined_Tumor, y = combined_CD4T)
combined_CD4T_Tumor

combined_CD4T_Tumor_high_CD4 <- combined_CD4T_Tumor[,(combined_CD4T_Tumor$orig.ident %in% c("high_CD4"))]
combined_CD4T_Tumor_Tumor_PBS <- combined_CD4T_Tumor[,(combined_CD4T_Tumor$orig.ident %in% c("PBS"))]


table(combined_CD4T_Tumor$celltype)
Idents(combined_CD4T_Tumor) <- 'celltype'
combined10 <- subset(combined_CD4T_Tumor,idents=c("CD4 T","Tumor cells") )
combined10$celltype <-as.factor(as.character(combined10$celltype))
table(combined10$orig.ident)   
Idents(combined10) <- 'orig.ident'                    
high_CD4_CD4T_Tumor <- subset(combined10,idents="high_CD4")                    
PBS_CD4T_Tumor <- subset(combined10,idents="PBS")                     

cco.high_CD4_CD4T_Tumor <- createCellChat(high_CD4_CD4T_Tumor@assays$RNA@data,
                                          meta = high_CD4_CD4T_Tumor@meta.data,group.by = "celltype")                    
cco.PBS_CD4T_Tumor <- createCellChat(PBS_CD4T_Tumor@assays$RNA@data,meta = PBS_CD4T_Tumor@meta.data,group.by = "celltype")                    
save(cco.high_CD4_CD4T_Tumor,cco.PBS_CD4T_Tumor,file = "cco_CD4T_Tumor.rda")

cellchat <- cco.PBS_CD4T_Tumor
cellchat@DB <- CellChatDB.mouse
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat,PPI.mouse)
cellchat <- computeCommunProb(cellchat,raw.use = TRUE,population.size = TRUE)
cellchat <- filterCommunication(cellchat,min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat,slot.name = "netP")
cellchat <- computeNetSimilarity(cellchat,type = "functional")
cellchat <- computeNetSimilarity(cellchat,type = "structural")                    
cco.PBS_CD4T_Tumor <- cellchat
saveRDS(cco.PBS_CD4T_Tumor,"cco.PBS_CD4T_Tumor.rds")
#修改high_CD4为PBS重复上述行

cco.list.CD4T.Tumor <- list(Pre=cco.PBS_CD4T_Tumor,Post=cco.high_CD4_CD4T_Tumor)
cellchat.CD4T.Tumor <- mergeCellChat(cco.list.CD4T.Tumor,
                                     add.names = names(cco.list.CD4T.Tumor),
                                     cell.prefix = TRUE)

#柱状图对比相互作用数量和强度
gg1 <- compareInteractions(cellchat.CD4T.Tumor, show.legend = F,group = c(1,2),measure = "count")
gg2 <- compareInteractions(cellchat.CD4T.Tumor,show.legend = F,group = c(1,2),measure = "weight")
p <- gg1 + gg2
p
#ggsave("Overview_number_strength_CD4T_BDC.pdf",p,width = 6,height = 4)

#两组对比后的数量和强度差异网络图-做不出来
par(mfrow=c(1,2))
netVisual_diffInteraction(cellchat.CD4T.Tumor,weight.scale = T,edge.width.max = 50)
netVisual_diffInteraction(cellchat.CD4T.Tumor,weight.scale = T,measure = "weight",edge.width.max = 20)
#方法二
par(mfrow=c(1,2))
weight.max <- getMaxWeight(cco.list.Immu.Tumor,attribute = c("idents","count"))
for (i in 1:length(cco.list.Immu.Tumor)) {
  netVisual_circle(cco.list.Immu.Tumor[[i]]@net$count,weight.scale = T,label.edge = F,
                   edge.weight.max = weight.max[2],edge.width.max = 12,
                   title.name = paste0("Number of interactions - ",names(cco.list.Immu.Tumor)[i]))
}

#信号通路
gg1 <- rankNet(cellchat.Immu.Tumor,mode = "comparison",stacked = T,do.stat = TRUE)
gg2 <- rankNet(cellchat.Immu.Tumor,mode = "comparison",stacked = F,do.stat = TRUE)
p <- gg1 + gg2
p

#某一具体通路
pathways.show <- c("FN1")
par(mfrow=c(1,1))
netVisual_aggregate(cco.high_CD4_Immu_Tumor,signaling = pathways.show,layout = "circle")
p = plotGeneExpression(cco.high_CD4_Immu_Tumor,signaling = "LCK")
p
markers <- c("Col1a1","Col1a2")
DotPlot(combined2_named,features = markers)+coord_flip()
netVisual_bubble(cellchat.Immu.Tumor, sources.use = 1, targets.use = 2, comparison = c(1, 2), 
                 max.dataset = 2,title.name = "Increased signaling", angle.x = 45,remove.isolate = T )


#---------------通路富集--------------------
#单细胞只保留免疫细胞，转成RNAseq形式做GSEA、GO通路富集
library(BiocManager)
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)

DefaultAssay(combined2_Immu) <- 'RNA'
Idents(combined2_Immu) <- 'orig.ident'

markers <- FindMarkers(combined2_Immu, ident.1 = "PBS", ident.2 = "high_CD4", 
                       only.pos = T,logfc.threshold = 0.25)
head(markers)
dim(markers)
markers = markers %>% rownames_to_column('gene') %>% filter(p_val_adj < 0.05)
head(markers)
dim(markers)
OrgDb = "org.Mm.eg.db"# 根据物种来指定
gene_convert <- bitr(markers$gene, fromType = 'SYMBOL', 
                     toType = 'ENTREZID', OrgDb = OrgDb)
head(markers)
markers = markers%>%inner_join(gene_convert,by=c("gene"="SYMBOL"))
head(markers)

#GO富集
ont = "BP"# posible value: BP, CC, MF, all
go.results <- enrichGO(markers$ENTREZID, keyType="ENTREZID",ont="BP",
                       OrgD = OrgDb, readable = FALSE)
head(go.results)
go.results <- enrichGO(markers$ENTREZID, keyType="ENTREZID",ont="ALL",
                       OrgD = OrgDb, readable = TRUE)
head(go.results)
dotplot(go.results,showCategory = 40,label_format=10000) #点图
barplot(go.results,showCategory = 40,label_format=10000) #柱状图

temp <- go.results@result
head(temp)
# 保存元数据到CSV文件
write.csv(temp, file = "go.results.csv", row.names = FALSE)

#前22条做柱状图，y轴为-log10(qvalue)
go.results_bp <- go.results[order(go.results$qvalue,decreasing = F),]
go.results_bp$Description <- factor(go.results_bp$Description, 
                                    levels = go.results_bp$Description)
go.results_top22 <- go.results_bp[1 : 22,]
ggplot(data=go.results_top22, aes(x=reorder(Description,-qvalue),y=-log10(qvalue))) + 
  geom_bar(stat="identity", width=0.8,fill='#DF4521') + 
  coord_flip() +  xlab("GO term") + ylab("-log10(qvalue)") + 
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())





#KEGG富集
organism = "mmu"# 对应KEGG数据库里的物种缩写
kegg.results <- enrichKEGG(markers$ENTREZID, organism = organism)
head(kegg.results)
kegg.results <- setReadable(kegg.results, OrgDb = OrgDb, keyType='ENTREZID')
head(kegg.results)
dotplot(kegg.results,label_format=10000) #点图
barplot(kegg.results, showCategory = 10,label_format=10000) #柱状图





#GSEA富集
min.pct = 0.01## 至少多少比例的细胞表达这个基因，过滤一些只在极少数细胞中有表达的基因
logfc.threshold = 0.01## 过滤掉在两组中几乎没有差异的基因
markers.for.gsea <- FindMarkers(combined2_Immu, ident.1 = "PBS", 
                                ident.2 = "high_CD4", 
                                min.pct = 0.01, logfc.threshold=0.01)
# GSEA 要求输入的是一个排好序的列表
Markers_genelist <- markers.for.gsea$avg_log2FC
names(Markers_genelist)= rownames(markers.for.gsea)
head(Markers_genelist)
Markers_genelist <- sort(Markers_genelist, decreasing = T)
head(Markers_genelist)
m_df = msigdbr(species = 'Mus musculus' , category = "C2")
head(m_df)
mf_df = m_df %>% dplyr::select(gs_name,gene_symbol) 
colnames(mf_df)<-c("term","gene")
gsea.results <- GSEA(Markers_genelist, TERM2GENE = mf_df)
head(gsea.results)
gseaplot(gsea.results, gsea.results@result$ID[1])
gseaplot(gsea.results, gsea.results@result$ID[2])
write_csv(gsea.results %>% data.frame, "gsea_results.csv")
write.csv(gsea.results, file = "gsea.results.csv")

dotplot(gsea.results)

m_df2 = msigdbr(species = 'Mus musculus' , category = "C7") #C7是免疫专用基因集，但是好像没什么区别
mf_df2 = m_df2 %>% dplyr::select(gs_name,gene_symbol) 
colnames(mf_df2)<-c("term","gene")
gsea.results2 <- GSEA(Markers_genelist, TERM2GENE = mf_df)
head(gsea.results2)
gseaplot(gsea.results, gsea.results@result$ID[1])

dotplot(gsea.results2)



#另一个方法GSEAGo
df <- markers[ , ]       # 移除第一列
rownames(df) <- df[, 7]  # 使用第一列作为行名
df <- df[ , -2]          # 保留avg_log2FC、p_val_adj、ENTREZID三列，并且ENTREZID作为行名
geneList=df$avg_log2FC
names(geneList)=df$ENTREZID 
geneList=sort(geneList,decreasing = T)
head(geneList)

ego <- gseGO(geneList     = geneList,
             OrgDb        = org.Mm.eg.db,
             ont          = "BP")

go=DOSE::setReadable(ego, OrgDb='org.Mm.eg.db',keyType='ENTREZID')
go=ego@result
pro='test_gsea'

go_gse=go
sortgo<-go_gse[order(go_gse$enrichmentScore, decreasing = T),]
head(sortgo)
dim(sortgo)
#write.table(sortkk,"gsea_output2.txt",sep = "\\t",quote = F,col.names = T,row.names = F)
#可以根据自己想要的通路画出需要的图
library(enrichplot)
plot <- gseaplot2(ego,#数据
          row.names(sortgo)[109],#画那一列的信号通路
          
          base_size = 20,#字体大小
          color = "#DF4521",#线条的颜色
          pvalue_table = TRUE,#加不加p值
          ES_geom="line")#是用线，还是用d点
  
plot

#GseaVis包，可以画出带NES值的图
ego2 <- setReadable(ego,
                   OrgDb = "org.Mm.eg.db",
                   keyType = "ENTREZID")
gseaNb(object = ego2,
       geneSetID = 'GO:0050863',
       subPlot = 3,
       addPval = T)



#write.csv(go,file = 'gse_go.csv') 

g_gseago<- ggplot(ego, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue)) + 
  geom_bar(stat="identity") + 
  scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
  scale_x_discrete(name ="Pathway names") +
  scale_y_continuous(name ="log10P-value") +
  coord_flip() + theme_bw(base_size = 15)+
  theme(plot.title = element_text(hjust = 0.5),  axis.text.y = element_text(size = 15))+
  ggtitle("Pathway Enrichment") 
g_gseago




 

