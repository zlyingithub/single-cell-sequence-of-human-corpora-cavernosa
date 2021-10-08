setwd("F:/学习/数据分析/single cell sequence/测序数据/human/ED海绵体分析/分析结果")
save.image(file = "ED分析.Rdata")
library(cowplot)
library(Seurat)
library(dplyr)
library(monocle)
library(stringr)
library(pheatmap)
library(vegan)  
library(SingleR)
library(RColorBrewer)
library(GENIE3)
library(CellChat)
color_plot<-brewer.pal(12, "Paired") #display.brewer.all()
color_plot<-brewer.pal(8, "Dark2")
color_plot<-brewer.pal(8, "Set2")
color_plot<-brewer.pal(8, "Set1")
color_plot<-c("#e41a1c","#377eb8","#7c02a2") #3组病人配色

color_plot<-c("#f47c4e","#e264b1","#6382c6","#4ab694","#9dd147","#efc922")

LZ037_expr <- Read10X(data.dir = "F:/学习/数据分析/single cell sequence/测序数据/human/ED海绵体分析/P21001/AllCells/1_Expression_Matrix")
LZ037 <- CreateSeuratObject(counts = LZ037_expr, project = "LZ037")
LZ037$sample <- "LZ037"
LZ037$sample1 <- "Normal_1"
LZ037$age<-"58 yo"
LZ037$tech<-"10X"
LZ037$disease<-"Normal"
remove(LZ037_expr)
LZ037[["percent.mt"]] <- PercentageFeatureSet(object = LZ037, pattern = "^MT-")

LZ040_expr <- Read10X(data.dir = "F:/学习/数据分析/single cell sequence/测序数据/human/ED海绵体分析/P9199_P21171/P9199/1_Expression_Matrix/filtered_feature_bc_matrix")
LZ040 <- CreateSeuratObject(counts = LZ040_expr, project = "LZ040")
LZ040$sample <- "LZ040"
LZ040$sample1 <- "Psychical_1"
LZ040$age<-"28 yo"
LZ040$tech<-"10X"
LZ040$disease<-"Psychical"
remove(LZ040_expr)
LZ040[["percent.mt"]] <- PercentageFeatureSet(object = LZ040, pattern = "^MT-")

LZ041_expr <- Read10X(data.dir = "F:/学习/数据分析/single cell sequence/测序数据/human/ED海绵体分析/P9199_P21171/P21171/1_Expression_Matrix/filtered_feature_bc_matrix")
LZ041 <- CreateSeuratObject(counts = LZ041_expr, project = "LZ041")
LZ041$sample <- "LZ041"
LZ041$sample1 <- "Vascular_1"
LZ041$age<-"31 yo"
LZ041$tech<-"10X"
LZ041$disease<-"Vascular"
remove(LZ041_expr)
LZ041[["percent.mt"]] <- PercentageFeatureSet(object = LZ041, pattern = "^MT-")

LZ043_expr <- Read10X(data.dir = "F:/学习/数据分析/single cell sequence/测序数据/human/ED海绵体分析/zlhos_1_zlhos_2/zlhos_1/1_Expression_Matrix/zlhos_1_filtered_feature_bc_matrix")
LZ043 <- CreateSeuratObject(counts = LZ043_expr, project = "LZ043")
LZ043$sample <- "LZ043"
LZ043$sample1 <- "Normal_2"
LZ043$age<-"55 yo"
LZ043$tech<-"10X"
LZ043$disease<-"Normal"
remove(LZ043_expr)
LZ043[["percent.mt"]] <- PercentageFeatureSet(object = LZ043, pattern = "^MT-")

LZ044_expr <- Read10X(data.dir = "F:/学习/数据分析/single cell sequence/测序数据/human/ED海绵体分析/zlhos_1_zlhos_2/zlhos_2/1_Expression_Matrix/zlhos_2_filtered_feature_bc_matrix")
LZ044 <- CreateSeuratObject(counts = LZ044_expr, project = "LZ044")
LZ044$sample <- "LZ044"
LZ044$sample1 <- "Normal_3"
LZ044$age<-"65 yo"
LZ044$tech<-"10X"
LZ044$disease<-"Normal"
remove(LZ044_expr)
LZ044[["percent.mt"]] <- PercentageFeatureSet(object = LZ044, pattern = "^MT-")

LZ045_expr <- Read10X(data.dir = "F:/学习/数据分析/single cell sequence/测序数据/human/ED海绵体分析/P9345/1_Expression_Matrix/P9345_filtered_feature_bc_matrix")
LZ045 <- CreateSeuratObject(counts = LZ045_expr, project = "LZ045")
LZ045$sample <- "LZ045"
LZ045$sample1 <- "Vascular_2"
LZ045$age<-"28 yo"
LZ045$tech<-"10X"
LZ045$disease<-"Vascular"
remove(LZ045_expr)
LZ045[["percent.mt"]] <- PercentageFeatureSet(object = LZ045, pattern = "^MT-")

LZ046_expr <- Read10X(data.dir = "F:/学习/数据分析/single cell sequence/测序数据/human/ED海绵体分析/P21278/1_Expression_Matrix/P21278_filtered_feature_bc_matrix")
LZ046 <- CreateSeuratObject(counts = LZ046_expr, project = "LZ046")
LZ046$sample <- "LZ046"
LZ046$sample1 <- "Diabetes_1"
LZ046$age<-"46 yo"
LZ046$tech<-"10X"
LZ046$disease<-"Diabetes"
remove(LZ046_expr)
LZ046[["percent.mt"]] <- PercentageFeatureSet(object = LZ046, pattern = "^MT-")

allsample_combined_cca<-merge(x=LZ037,y=c(LZ037,LZ040,LZ041,LZ043,LZ044,LZ045))

allsample_combined_cca<-RenameIdents(object = allsample_combined_cca, `LZ003` = "LZ003",`lz013` = "BD",`lz014` = "BD",`lz015` = "BD",`LZ009` = "other_10X",`LZ007` = "other_10X",`LZ008` = "other_10X",`LZ011` = "other_10X",`lz026` = "OLD",`lz033` = "OLD",`lz034` = "QD")
allsample_combined_cca<-RenameIdents(object = allsample_combined_cca, `Normal` = "Normal",`Psychical` = "Vascular",`Vascular` = "Vascular",`Diabetes` = "Diabetes")
allsample_combined_cca<-RenameIdents(object = allsample_combined_cca, `Normal_1` = "Normal_1",`Normal_2` = "Normal_2",`Normal_3` = "Normal_3",
                                     `Psychical_1` = "non-DM_1",`Vascular_1` = "non-DM_2",`Vascular_2` = "non-DM_3",`Diabetes_1` = "Diabetes_1")
allsample_combined_cca$sample2<-Idents(allsample_combined_cca)
allsample_combined_cca$disease2<-Idents(allsample_combined_cca)
allsample_combined_cca$batch1<-Idents(allsample_combined_cca)

allsample_combined_cca$sample1<-factor(allsample_combined_cca$sample1,levels = c("Normal_1","Normal_2","Normal_3","Psychical_1","Vascular_1","Vascular_2","Diabetes_1"),ordered = F)
allsample_combined_cca$seurat_12clusters<-factor(allsample_combined_cca$seurat_12clusters,levels = c("EndothelialC","LEC","ENF","Fibroblast","SMC","MFB1","MFB2","MFB3","SchwannC","M1","M2","T"),ordered = F)
allsample_combined_cca<-RenameIdents(object = allsample_combined_cca, `EndothelialC` = "EndothelialC",`LEC` = "EndothelialC",`ENF` = "Fibroblast",`Fibroblast` = "Fibroblast",
                                     `SMC` = "SMC",`MFB1` = "SMC",`MFB2` = "SMC",`MFB3` = "SMC",`SchwannC` = "SchwannC",`M1` = "MAC",`M2` = "MAC",`T` = "T")
#输出表达矩阵
MARKER<-allsample_combined_cca1@assays$RNA@counts
MARKER<-data.frame(MARKER)
colnames(MARKER)<-paste(allsample_combined_cca1@meta.data[,4],rownames(allsample_combined_cca1@meta.data))

#批次效应处理
allsample_combined_cca_list <- SplitObject(allsample_combined_cca, split.by = "sample1")
for (i in 1:length(allsample_combined_cca_list)) {
  allsample_combined_cca_list[[i]] <- NormalizeData(allsample_combined_cca_list[[i]], verbose = FALSE)
  allsample_combined_cca_list[[i]] <- FindVariableFeatures(allsample_combined_cca_list[[i]], selection.method = "vst", 
                                                           nfeatures = 1000, verbose = FALSE)
}
allsample_combined_cca_anchors <- FindIntegrationAnchors(object.list = allsample_combined_cca_list, dims = 1:30,anchor.features = 1000)
allsample_combined_cca <- IntegrateData(anchorset = allsample_combined_cca_anchors, dims = 1:30)


#质量控制
allsample_combined_cca <- subset(x = allsample_combined_cca, subset = nFeature_RNA > 800 & nFeature_RNA < 7000 & percent.mt < 10 & nCount_RNA <30000)
FeaturePlot(allsample_combined_cca, features = c("nFeature_RNA","nCount_RNA" ,"percent.mt"), cols = c("blue", "yellow"),max.cutoff="q95",min.cutoff = "q10",pt.size = 0.1, reduction = "tsne",ncol = 3)
VlnPlot(allsample_combined_cca,features = c("nFeature_RNA","nCount_RNA" ,"percent.mt"),pt.size = 0.1,group.by = "sample1")
#标准化
#DefaultAssay(allsample_combined_cca) <- "RNA"
#DefaultAssay(allsample_combined_cca) <- "integrated"
#allsample_combined_cca<-NormalizeData(allsample_combined_cca)
#allsample_combined_cca<-FindVariableFeatures(allsample_combined_cca, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
allsample_combined_cca <- ScaleData(allsample_combined_cca)
allsample_combined_cca <- RunPCA(allsample_combined_cca, npcs = 40)
#聚类
allsample_combined_cca <- FindNeighbors(object = allsample_combined_cca, dims = 1:20)
allsample_combined_cca <- FindClusters(object = allsample_combined_cca, resolution = 1)
#降维
allsample_combined_cca <- RunUMAP(object = allsample_combined_cca, dims = 1:30,n.neighbors = 35L,seed=15)
allsample_combined_cca <- RunTSNE(object = allsample_combined_cca, dims = 1:30,n.neighbors = 35L,seed=17)

# Visualization
DimPlot(object = allsample_combined_cca, reduction = "tsne", group.by = "age2",label.size = 8,pt.size = 0.1)+ theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line = element_blank())
DimPlot(object = allsample_combined_cca, reduction = "tsne", group.by = "sample1",label.size = 8,pt.size = 0.1)+ theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line = element_blank())
DimPlot(object = allsample_combined_cca, reduction = "tsne",  label = T,label.size = 6,pt.size = 0.1,cols = color_plot)+ theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line = element_blank())

#SingleR 自动命名 ref_Human_all在注释文件夹中
allsamples_for_SingleR <- GetAssayData(allsample_combined_cca, slot="data")
clusters <- allsample_combined_cca@meta.data$seurat_clusters
pred.hesc <- SingleR(test = allsamples_for_SingleR, ref = ref_Human_all, 
                     labels = ref_Human_all$label.main,
                     #因为样本主要为免疫细胞（而不是全部细胞），因此设置为label.fine
                     
                     #这里我们为上一步分的9个cluster注释celltype
                     assay.type.test = "logcounts", assay.type.ref = "logcounts")

allsample_combined_cca$singleR<-as.character(pred.hesc$labels)
DimPlot(allsample_combined_cca,reduction = "tsne" ,group.by="singleR",label = T,pt.size = 0.1,label.size = 5)
table(pred.hesc$labels)

#看细胞亚型基因表达
Idents(allsample_combined_cca)<-allsample_combined_cca$seurat_41clusters
allsample_combined_cca_subset<-subset(x = allsample_combined_cca, idents = "SC")#亚型
Idents(allsample_combined_cca_subset)<-allsample_combined_cca_subset$age2 #分组
VlnPlot(allsample_combined_cca_subset,features = c("MKI67"),pt.size = 0.1,group.by = "seurat_INC6clusters",cols = color_plot_qualitative)+NoLegend()+geom_smooth(aes(group=1),size=2,method = 'lm',formula = y~ns(x,3),color="black")+theme(axis.title=element_blank(),axis.text = element_text(size = 35),plot.title =element_text(size = 45,face="plain") )
FeaturePlot(object = allsample_combined_cca, reduction = "tsne", features = c("CTGF"), cols = c("#eeeeee", "darkred"), ncol = 1,pt.size = 0.1)+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line = element_blank(),legend.position = c(.85,.15))

VlnPlot(allsample_combined_cca,features = c("gene_SASP"),pt.size = 0.1,group.by = "seurat_11clusters",cols = color_plot)+geom_boxplot(fill = "black",color="gray90",size=1,width=0.05,outlier.size = 0)+NoLegend()+theme(axis.title=element_blank(),axis.text = element_text(size = 35),plot.title =element_text(size = 45,face="plain") )
#带趋势线
VlnPlot(allsample_combined_cca_germcells,features = c("percent.glycolysis_gene"),pt.size = 0.1)+geom_boxplot(fill = "black",color="gray90",size=1,width=0.05,outlier.size = 0)+NoLegend()+theme(axis.title=element_blank(),axis.text = element_text(size = 35),plot.title =element_text(size = 45,face="plain") )+stat_summary(aes(group=1),fun.y=mean, geom="smooth", shape=1,size=2,color="gray")+ scale_y_continuous(limits = c(0,2))

allsample_combined_cca<-subset(x = allsample_combined_cca, idents = c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","15","16","17","18","19","20","21"))
allsample_combined_cca <- RenameIdents(object = allsample_combined_cca, `0` = "Fibroblast", `1` = "EndothelialC",`2` = "EndothelialC",
                                       `3` = "Fibroblast",`4` = "Fibroblast",`5` = "EndothelialC",`6` = "EndothelialC",`7` = "Fibroblast",`8` = "SMC",`9` = "Fibroblast",`10` = "Fibroblast",`11` = "EndothelialC",
                                       `12` = "Fibroblast",`13` = "ENF",`14` = "LEC",`15` = "T",`16` = "M1",`17` = "Fibroblast",`18` = "MFB1",`19` = "MFB2",`20` = "M2",`21` = "MFB3",`22` = "SchwannC",`23`="T",`24`="T")
allsample_combined_cca$seurat_noa9clusters<-Idents(allsample_combined_cca)
allsample_combined_cca$age<-factor(allsample_combined_cca$age,levels = c("OA","iNOA","KS","AZFa_Del"),ordered = F)
allsample_combined_cca$seurat_10clusters<-factor(allsample_combined_cca$seurat_10clusters,levels = c("SPG","SPC","SPT","SC","LC","PTM","EC","VSM","MAC","LYM"),ordered = F)
MARKERS <- FindAllMarkers(object = allsample_combined_cca, only.pos = F, min.pct = 0.1, logfc.threshold = 1)
MARKERS <- FindMarkers(object = allsample_combined_cca, only.pos = F,ident.1 = "M1",ident.2 = "M2", min.pct = 0.1, logfc.threshold = 0.25)
write.csv(MARKERS,"6群DEGs.csv")
#计算基因平均值
MARKERS <- AverageExpression(allsample_combined_cca, return.seurat = T)
#计算metadata里的平均值（表格，按照xx分组list属性，平均值或其他）
MARKERS<-aggregate(allsample_combined_cca_10X_KSandOA@meta.data,list(allsample_combined_cca_10X_KSandOA@meta.data$seurat_11clusters,allsample_combined_cca_10X_KSandOA@meta.data$age),mean)

#基因表达量和占比气泡图
DotPlot(allsample_combined_cca_subset, features = c("MYOC","DEFB1","LEPR","PPP1R14A","KERA","BMP7"), cols = c("lightgrey", "red"))+scale_x_discrete(position = "top")+ coord_flip()

#计算差异基因
Idents(allsample_combined_cca)<-allsample_combined_cca$seurat_11clusters
listfori<-c("Fibroblast","EndothelialC","SMC")
for (i in listfori) {
  allsample_combined_cca_subset<-subset(x = allsample_combined_cca, idents = i)
  Idents(allsample_combined_cca_subset)<-allsample_combined_cca_subset$sample1
  MARKERS <- FindMarkers(object = allsample_combined_cca_subset, only.pos = F,ident.1 = "Normal_1",ident.2 = "Vascular_2", min.pct = 0.1, logfc.threshold = 0.25)
  write.csv(MARKERS,paste("Normal_1vsVascular_2差异基因",i,".csv",sep=""))
}

#计算基因集（多个基因）表达
listfori<-c("gene_SASP","gene_collagen","gene_Interleukins","gene_Chemokines","gene_Otherinflammatorymolecules","gene_Growthfactors","gene_Proteases","gene_L_R","gene_PEGNOS","gene_Insolublefactors")
for (i in listfori) {
  allsample_combined_cca<-AddModuleScore(allsample_combined_cca,features = get(i),name = i)
}


#计算随年龄变化回归分析(只分析成人的样本)
DefaultAssay(allsample_combined_cca) <- "RNA"
Idents(allsample_combined_cca)<-allsample_combined_cca$age1
allsample_combined_cca <- RenameIdents(object = allsample_combined_cca, `2 yo` = 1, `8 yo` = 1, `11 yo` = 1, `20+ yo` = 2, `51 yo` = 3, `67 yo` = 4, `81 yo` = 5)
allsample_combined_cca$age3<-Idents(allsample_combined_cca)
allsample_combined_cca_adult<-subset(x = allsample_combined_cca, idents = c(2,3,4,5))#只分析成年后的样本
#建立表格数据
Idents(allsample_combined_cca_adult)<-allsample_combined_cca_adult$seurat_10clusters
for(j in c("SPG" ,"SPC","SPT","SC","LC","PTM","EC","VSM","MAC","LYM")){
  allsample_combined_cca_subset<-subset(x = allsample_combined_cca_adult, idents = j)
  allsample_combined_cca_subset<-NormalizeData(allsample_combined_cca_subset)
  gene_subset<-row.names(subset(allsample_combined_cca_subset@assays$RNA@meta.features,vst.mean>0.1))
  MARKERS<-t(data.frame(allsample_combined_cca_subset@assays$RNA@data))[,gene_subset]
  MARKER<-cbind(allsample_combined_cca_subset@meta.data[,c(2,3,19:35)],MARKERS) #在回归分析中加入的额外参数，需要根据meta.data修改
  MARKER$age3<-as.numeric(MARKER$age3) #将年龄等级数据转换为数值
  a<-colnames(MARKER)
  #回归分析（每次做前都建立空表格）
  c<-read.csv("空表格.csv",row.names=1) #建立空表格
  for(i in a){
    b<-lm(get(i)~age3,MARKER) #填入相应的Y~X
    b<-summary(b)
    b1<-t(data.frame(b$coefficients[2,])) #导出相关系数
    b1$R_squared<-data.frame(b$r.squared) #导出r值
    b1<-data.frame(b1)
    colnames(b1)<-c("a","Std.Error","t","p","r.squared")
    rownames(b1)<-i
    c<-rbind(c,b1)
  }
  write.csv(c,paste0("所有基因随衰老变化回归分析",j,".csv"))
}
#作图
#回归图
ggplot(data = MARKER, mapping = aes(x = age3, y = RPS26),geom="jitter")+ 
  geom_point(aes(color = nFeature_RNA,size= nCount_RNA),alpha = 1/10)+ #以drv为分组设置点的颜色
  geom_smooth(method = 'lm', formula = y ~ x,colour="ORANGE",size=2)+ scale_y_continuous(limits = c(1,5))
#

#画各个样本比例
ggplot(data = allsample_combined_cca_subset@meta.data, mapping = aes(x = disease2, fill = seurat_SMC4clusters))+geom_bar(stat = 'count', position = 'fill')+scale_fill_manual(values=color_plot) + labs(x = '', y = 'Rate')+theme_classic()

#代谢分析
#读取基因
gene_X<-read.csv("F:/学习/数据分析/single cell sequence/注释文件/human/chr_X_gene.csv",header = F,sep = "\t",,stringsAsFactors=F)[,1]
gene_Y<-read.csv("F:/学习/数据分析/single cell sequence/注释文件/human/chr_Y_gene.csv",header = F,sep = "\t",,stringsAsFactors=F)[,1]

gene_OXPHOS<-list(read.csv("F:/学习/数据分析/single cell sequence/注释文件/human/氧化磷酸化-无线粒体基因 go 基因列表.csv",header = F,sep = "\t",,stringsAsFactors=F)[,1])
gene_glycolysis<-list(read.csv("F:/学习/数据分析/single cell sequence/注释文件/human/经典糖酵解 go 基因列表.csv",header = F,sep = "\t",,stringsAsFactors=F)[,1])
gene_triglyceride_metabolic<-list(read.csv("F:/学习/数据分析/single cell sequence/注释文件/human/甘油三酯代谢.csv",header = F,sep = "\t",,stringsAsFactors=F)[,1])
gene_mitochondria<-list(read.csv("F:/学习/数据分析/single cell sequence/注释文件/human/线粒体基因.csv",header = F,sep = "\t",,stringsAsFactors=F)[,1])

#计算X，Y染色体基因表达，注意这里使用的genelist应该都要在object的RNA@meta.features里面存在，所以应当筛选
gene_ALL<-row.names(allsample_combined_cca@assays$RNA@meta.features)
gene_X<-subset(gene_X,gene_X %in% gene_ALL)
gene_Y<-subset(gene_Y,gene_Y %in% gene_ALL)

#计算比例
allsample_combined_cca[["percent.X_gene"]] <- PercentageFeatureSet(object = allsample_combined_cca, features= gene_X)
allsample_combined_cca[["percent.Y_gene"]] <- PercentageFeatureSet(object = allsample_combined_cca, features= gene_Y)


#基因表达相关性回归分析
MARKER<-allsample_combined_cca_10X_KS_SC@meta.data
MARKERS<-allsample_combined_cca_10X_KS_SC@assays$RNA@data
MARKERS<-data.frame(MARKERS)
MARKERS<-t(MARKERS)
MARKER<-cbind(MARKER,MARKERS)

ggplot(data = MARKER, mapping = aes(x = XIST, y = percent.chr_neXi_gene))+ 
  geom_point(aes(color = percent.chr_X_gene,size= nCount_RNA))+ #以drv为分组设置点的颜色
  geom_smooth(method = 'lm', formula = y ~ x,colour="ORANGE",size=2)+labs(caption ="y = 4.092-0.698x  p-value = 8.185e-14")
MARKERS<-lm(percent.chr_eXi_gene~XIST,MARKER)
summary(MARKERS) 

#简单回归相关性作图 
FeatureScatter(allsample_combined_cca_subset, feature1 = "PC_1", feature2 = "SRD5A2",cols = color_plot,pt.size = 0.1)+
  geom_smooth(aes(group=1),size=3,method = 'lm',formula = y~ns(x,5),color="grey")+ 
  scale_y_continuous(limits = c(0.05,0.2))+ scale_x_continuous(limits = c(-15,18))


#绘制火山图
MARKER<-read.csv("KS SSC 所有基因与all X回归分析.csv",header = T,stringsAsFactors=F)

ggplot(data = MARKER, aes(x = r, y = -log10(p.value),size=mean,color=UPDOWN)) +
  geom_point(alpha=0.6)+scale_size(limits=c(0,4)) + scale_color_manual(values=c("#90bff9", "grey","#f2b77c"))+
  xlim(c(-1.1, 0.7)) +
  ylim(c(0, 4.7)) +
  geom_vline(xintercept=c(-0.1,0.1),lty=4,col="gray",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="gray",lwd=0.8)


#各组平均表达热图
cgene<-c("GFRA1","SYCP3","TNP1","SOX9","STAR","MYH11","VWF","NOTCH3","CD163","TPSAB1","CD3D")
MARKER <- AverageExpression(allsample_combined_cca, return.seurat = T)
MARKERS<-ScaleData(MARKERS,features = top30$gene)

top30 <- MARKERS %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
DoHeatmap(MARKERS, features = c("USP9Y","DDX3Y","UTY","HSFY1","HSFY2","RBMY1A1","RBMY1B","RBMY1C","RBMY1D","RBMY1E","RBMY1J","DAZ1","DAZ2","DAZ3","DAZ4","BPY2","CDY1B","PRY","CSPG4P1Y","DAZL","BOLL","SRY","ZFX","ZFY"), size = 3,slot = "data", disp.min = 0, disp.max = 5)+scale_fill_gradient2(low = "steelblue", mid = "white", high = "darkorange")
DoHeatmap(MARKERS, features = top30$gene, size = 6,draw.lines = F)+scale_fill_gradient2(low = "#02197b", mid = "white", high = "#b10000")

#计算多个基因打分
IMMUNE_MARKER<-list(c("PTPRC","TPSAB1","CD163","CD3D","CD4","CD8A","MS4A1", "CD79A","GNLY","ITGAX"))
gene_androgen_synthesis<-list(c("STAR","CYP11A1","CYP17A1","HSD17B6","SCARB1","MED1","HSD3B2", "HSD17B3","SRD5A1","SRD5A2","SRD5A3"))
allsample_combined_cca<-AddModuleScore(object = allsample_combined_cca,features = gene_paternally_imprinted,name = 'gene_paternally_imprinted')
#由于有些通路没有差异的基因比较多，只取差异基因计算
allsample_combined_cca<-AddModuleScore(object = allsample_combined_cca,features = intersect(gene_paternally_imprinted),name = 'gene_paternally_imprinted')




#计算差异度
#建立表格数据
Idents(allsample_combined_cca)<-allsample_combined_cca$age2
allsample_combined_cca_subset<-subset(x = allsample_combined_cca, idents = c("childbearing","agedness"))
MARKER<-AverageExpression(allsample_combined_cca_subset,add.ident = "seurat_10clusters")
MARKER<-MARKER$integrated
#MARKER<-MARKER$RNA
#MARKER<-MARKER[allsample_combined_cca_8yo@assays$RNA@var.features,]或者取前1000个差异基因
MARKER<-t(MARKER)
MARKER<-vegdist(MARKER,'bray')
MARKER<-as.matrix(MARKER)
write.csv(MARKER,"各年龄各细胞bray差异度.csv") #手动修改成作图的两列
MARKER<-read.csv("细胞差异度10cluster.csv",,header = T)
#标准化各组
MARKER$cell<-factor(MARKER$cell,levels = rev(c("SPG","SPC","SPT","SC","LC","PTM","EC","VSM","MAC","LYM")))
ggplot(data = MARKER, mapping = aes(x=Dis, y = cell))+ 
  geom_point(aes(color = Dissimilarity.Jaccard,size= Dissimilarity.Bary))+scale_color_gradient(low = "gray", high = "red")+
  scale_size_continuous(range=c(4,10))+theme(axis.title=element_blank(),axis.text = element_text(size = 20))

#小提琴图
Idents(allsample_combined_cca)<-allsample_combined_cca$seurat_10clusters
allsample_combined_cca_subset<-subset(x = allsample_combined_cca, idents = "LC")
Idents(allsample_combined_cca_subset)<-allsample_combined_cca_subset$age2
VlnPlot(allsample_combined_cca_subset$seurat_INC6clusters,features = c("gene_paternally_imprinted1"),pt.size = 0.1,group.by = "age2",cols = color_plot_qualitative)+NoLegend()+geom_smooth(aes(group=1),size=2,method = 'lm',formula = y~ns(x,8),color="black")+theme(axis.title=element_blank(),axis.text = element_text(size = 35),plot.title =element_text(size = 45,face="plain") )

allsample_combined_cca_subset <- subset(x = allsample_combined_cca_subset, subset = percent.Y_gene <0.3)

#画带透明度的Featureplot
i<-c("GFRA1")
P1<-FeaturePlot(object = allsample_combined_cca, reduction = "umap", features = c(i), ncol = 1,pt.size = 0.1)+theme(axis.title=element_blank(),axis.text = element_text(size = 0),legend.position = c(.85,.15))
ggplot(P1$data,aes(UMAP_1,UMAP_2,color=get(i)))+geom_point(size=0.1,alpha=0.7)+scale_color_gradient(low = "#f3f7f7",high = "#5c0000")+P1$theme+labs(title=i)+theme(plot.title = element_text(hjust = 0.5))


#构建互动图形，手动选择细胞
plot <- FeaturePlot(pbmc3k.final, features = "MS4A1")
HoverLocator(plot = plot, information = FetchData(pbmc3k.final, vars = c("ident", "PC_1", "nFeature_RNA")))
pbmc3k.final <- RenameIdents(pbmc3k.final, DC = "CD14+ Mono")
plot <- DimPlot(pbmc3k.final, reduction = "umap")
select.cells <- CellSelector(plot = plot)

#转录因子调控关系
exprMatr<-as.matrix(allsample_combined_cca_subset[["RNA"]]@data)[unique(MARKERS$gene),] 
head(exprMatr[1:5,1:5])
regulators <- read.csv("F:/学习/数据分析/single cell sequence/注释文件/human/Homo_sapiens_TF.csv",header = F,stringsAsFactors=F)[,1]
regulators<-intersect(regulators,unique(MARKERS$gene))
weightMat <- GENIE3(exprMatr, regulators=regulators, nCores=4, verbose=TRUE)
linkList <- getLinkList(weightMat, reportMax=5)
linkList <- getLinkList(weightMat, threshold=0.1)

#手动气泡图
MARKER<-read.csv("F:/学习/数据分析/single cell sequence/测序数据/human/ED海绵体分析/FIGURE/FIG 3 ENC/signaling.csv",header = T,stringsAsFactors=F)

MARKER$IPA_Pathways<-factor(MARKER$IPA_Pathways,levels = c("RhoGDI Signaling","Oxidative Phosphorylation","PPAR Signaling","HIPPO signaling",
                                                              "Apelin Adipocyte Signaling Pathway","TGF-β Signaling","Renin-Angiotensin Signaling","Leukocyte Extravasation Signaling","Glycolysis I","Calcium Signaling","Insulin Receptor Signaling","VEGF Signaling","NGF Signaling","Remodeling of Epithelial Adherens Junctions","Integrin Signaling"),ordered = F)
MARKER$disease<-factor(MARKER$disease,levels = c("non-DM","Diabetes"),ordered = F)



ggplot(MARKER,aes(disease,IPA_Pathways,p.value,z.score))+ 
  geom_point(aes(size=p.value,color=z.score))+
  scale_colour_gradient2(low="lightblue",mid="lightgrey",high="red")+
  scale_size_continuous(range=c(6,10))+
  scale_y_discrete(position = "right")+ theme_bw()+theme(panel.grid.major=element_line(colour=NA))

#手动分组型小提琴加箱图(加漂亮的箱图)
i<-c("gene_IGFreceptorsignaling1")#输入基因名称
MARKER<-data.frame(allsample_combined_cca_ENC[["RNA"]]@data)[i,]
MARKER<-t(MARKER)
rownames(MARKER)<-rownames(allsample_combined_cca_ENC@meta.data)
MARKER<-cbind(allsample_combined_cca_ENC@meta.data,MARKER)
p<-VlnPlot(allsample_combined_cca_ENC,features = i, split.by = "disease2",pt.size = 0.1,cols = color_plot)+ scale_fill_manual(values = alpha(color_plot, .8))+NoLegend()+theme(axis.title=element_blank(),axis.text.x=element_text(angle=0,vjust=0.5,hjust = 0.5,colour="black",family="Arial",face="plain",size=40), axis.text.y=element_text(family="Arial",size=30,face="plain"),plot.title =element_text(size = 30,face="plain") )
#修改分组方式
p+geom_boxplot(data = MARKER,aes(x=seurat_ENC4clusters, y=get(i)),width=0.1,size=1,position=position_dodge(0.9), alpha = 0.5,color="#eeeeee")
