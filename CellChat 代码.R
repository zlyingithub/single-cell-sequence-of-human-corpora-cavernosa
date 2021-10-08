library(CellChat)
library(ggplot2)
library(ggalluvial)
library(Seurat)
library(NMF)
library(ggalluvial)
Idents(allsample_combined_cca)<-allsample_combined_cca$seurat_8clusters
allsample_combined_cca_subset<-subset(x = allsample_combined_cca, idents = c("INC","SMC","SWC","MAC","T"))
allsample_combined_cca_subset$seurat_5v4clusters<-Idents(allsample_combined_cca_subset)
allsample_combined_cca_ENC$seurat_5v4clusters<-allsample_combined_cca_ENC$seurat_ENC4clusters
allsample_combined_cca_subset<-merge(allsample_combined_cca_ENC,allsample_combined_cca_subset)
#创建cellchat对象
cellchat <- createCellChat(allsample_combined_cca_subset, meta = allsample_combined_cca_subset@meta.data, group.by = "seurat_5v4clusters", do.sparse = T)
cellchat <- addMeta(cellchat, allsample_combined_cca_subset@meta.data, meta.name = unique(colnames(allsample_combined_cca_subset@meta.data)))
cellchat <- setIdent(cellchat, ident.use = "seurat_5v4clusters") # set "labels" as default cell identity
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

#加载与设定需要的CellChatDB数据库

CellChatDB <- CellChatDB.human
# use Secreted signaling for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
#CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor")
#CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact")
#cellchat@DB <- CellChatDB.use # set the used database in the object
cellchat@DB <- CellChatDB #所有互作模式
#预处理表达数据以进行细胞间相互作用分析
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
cellchat <- identifyOverExpressedGenes(cellchat) 
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
#推断细胞间相互作用网络与分析
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")



#可视化
pathways.show <- c("VEGF") 
vertex.receiver = seq(1,3) # a numeric vector
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = groupSize) # Hierarchy plot
netVisual_aggregate(cellchat, signaling = "TGFb", layout = "circle", vertex.size = groupSize) # Circle plot
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use = color_plot)

#查看配对对整个信号贡献
netAnalysis_contribution(cellchat, signaling = "LIGHT")

#小提琴
plotGeneExpression(cellchat, signaling = "NOTCH")

C2_gene<-c("ACKR1","LIFR","LTBR","IL1R1","CD34","SELP","SELE","gene_MHC21")

#所有通路贡献
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, signaling =unique(cellchat@DB$interaction$pathway_name), pattern = "outgoing",width = 6,height = 32,font.size = 10)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, signaling =unique(cellchat@DB$interaction$pathway_name), pattern = "incoming",width = 6,height = 32,font.size = 10)
ht1 + ht2
#只看感兴趣通路
cellchat_signaling<-c("TGFb","BMP","IL6","ANGPT","KIT","SEMA3","BTLA","VCAM","VEGF","CCL","LIFR","OSM","IL1","LIGHT","NGF","CALCR","CD34","SELPLG","WNT","CX3C","CSF3","IFN-I","TRAIL","CD40","SPP1","VISFATIN","ANGPTL","NPR1","NPR2","HGF","FN1","CD46","CDH5","NOTCH","PECAM1","PTPRM","SN","GDF","MSTN","GDNF","ACTIVIN","EGF","FGF","TWEAK","SELE")
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, signaling =cellchat_signaling, pattern = "outgoing",width = 6,height = 32,font.size = 10)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, signaling =cellchat_signaling, pattern = "incoming",width = 6,height = 32,font.size = 10)
ht1 + ht2

#按输入输出聚类
netAnalysis_signalingRole_scatter(cellchat)
netAnalysis_signalingRole_network(cellchat, signaling = "WNT", width = 8, height = 2.5, font.size = 10)

#相互作用热图
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = "SELE", color.heatmap = "Reds",font.size = 16,font.size.title = 20)
