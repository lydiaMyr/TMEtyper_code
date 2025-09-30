library(Seurat)
library(ggplot2)
load("F:/workspace/TME_project_file/seurat_obj.rda")
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
DimPlot(seurat_obj, reduction = "umap", label = TRUE) + 
  ggtitle("UMAP Clustering") + theme_minimal()
# 定义免疫细胞标记基因列表
immune_markers <- list(
  T_cells = c("CD3D", "CD3E"),
  CD8_T = c("CD8A", "CD8B", "GZMK"),
  CD4_T = c("CD4", "IL7R"),
  Tregs = c("FOXP3", "IL2RA"),
  B_cells = c("CD19", "MS4A1"),
  NK_cells = c("NCAM1", "NKG7"),
  Macrophages = c("CD68", "CD163"),
  DCs = c("CD1C", "CLEC9A")
)

# 可视化标记基因表达（点图）
DotPlot(seurat_obj, features = unlist(immune_markers)) +
  RotatedAxis() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red")


new_cluster_ids <- c(
  "0" = "CD4 T",
  "1" = "CD8 T",
  "2" = "CD8 T",
  "3" = "NK",
  "4" = "B cell",
  "5" = "CD8 T",
  "6" = "NK",
  "7" = "Macrophage",
  "8" = "Macrophage",
  "9" = "CD8 T",
  "10" = "B cell",
  "11" = "Macrophage",
  "12" = "B cell",
  "13" = "B cell"
)
seurat_obj$cell_type <- plyr::revalue(as.character(Idents(seurat_obj)), new_cluster_ids)
pdf("F:/workspace/TME_project_file/单细胞/large_cell_sub.pdf",6.5,6)
DimPlot(seurat_obj, group.by = "cell_type", label = TRUE) + 
  ggtitle("Cell Type Annotation")
dev.off()

marker_genes <- c("PAX8","CDX2","PLA2G2D","DSG3","CSF3R","LRP2","GPM6B","CXCL9")
pdf("F:/workspace/TME_project_file/单细胞/CXCL9_feature_plot.pdf",6.5,6)
FeaturePlot(seurat_obj, features = "CXCL9")
dev.off()
# 绘制热图
DotPlot(seurat_obj,,features = marker_genes) + 
  RotatedAxis() + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red")
DoHeatmap(seurat_obj, features = marker_genes) +
  scale_fill_gradientn(colors = c("blue", "white", "red"))


seurat_obj_new = seurat_obj
Idents(seurat_obj_new) = "cell_type"
cell_sub_id = as.vector(unlist(Idents(seurat_obj_new)))
names(cell_sub_id) = names(Idents(seurat_obj_new))
macro_name = as.vector(unlist(Idents(macro_cell_new)))
names(macro_name) = names(Idents(macro_cell_new))
cell_sub_id[names(macro_name)]=macro_name
Idents(seurat_obj_new) =factor(cell_sub_id)

#细胞通讯分析
library(CellChat)
seurat_obj_new$cell_type=cell_sub_id
CellChatDB <- CellChatDB.human  # 或 CellChatDB.mouse
showDatabaseCategory(CellChatDB)  # 查看数据库分类
cellchat@DB <- subsetDB(CellChatDB, search = "Secreted Signaling")  # 选择分泌信号通路
cellchat <- createCellChat(object = seurat_obj_new, 
                           meta = seurat_obj_new@meta.data, 
                           group.by = "cell_type")  # 关键参数：指定分组列名

cellchat <- subsetData(cellchat)                     # 筛选信号相关基因
cellchat <- identifyOverExpressedGenes(cellchat)     # 鉴定差异表达基因
cellchat <- identifyOverExpressedInteractions(cellchat) # 高表达配体-受体对

cellchat <- computeCommunProb(cellchat, raw.use = TRUE)  # 配体-受体水平概率
cellchat <- filterCommunication(cellchat, min.cells = 10) # 过滤低可信通讯
df.net <- subsetCommunication(cellchat)                  # 提取结果表格


cellchat <- computeCommunProbPathway(cellchat)           # 信号通路水平概率
cellchat <- aggregateNet(cellchat)                       # 聚合网络

pdf("F:/workspace/TME_project_file/单细胞/cell_chat_result.pdf",6,6)

#netVisual_circle(cellchat@net$count, title.name = "交互次数")
#netVisual_circle(cellchat@net$weight, title.name = "交互强度")
netVisual_heatmap(cellchat, measure = "weight")  
dev.off()

df.net <- subsetCommunication(cellchat, 
                              sources.use = "TAM", 
                              targets.use = c("M1","M2","CD4 T","CD8 T","NK","Mreg","B cell"))
genes=unique(c(df.net$ligand,df.net$receptor))

# 绘制TAMs发送的信号通路（如VEGF）
pdf("F:/workspace/TME_project_file/单细胞/cell_chat_ligand_receptor.pdf",4,6)
pairLR = data.frame(interaction_name=df.net$interaction_name)
netVisual_bubble(
  cellchat,
  sources.use = "TAM",    # 发送者
  targets.use = sort(c("M1","M2","CD4 T","CD8 T","NK","Mreg","B cell")),   # 接收者
  pairLR.use = pairLR,    # 指定基因对
  sort.by.target = TRUE,         # 按接收者排序
  title.name = "TAM Interactions"                    # X轴标签旋转角度
)

df.net.m2 <- subsetCommunication(cellchat, 
                              sources.use = "M2", 
                              targets.use = c("M1","TAM","CD4 T","CD8 T","NK","Mreg","B cell"))
# 绘制TAMs发送的信号通路（如VEGF）
pairLR = data.frame(interaction_name=df.net$interaction_name)
netVisual_bubble(
  cellchat,
  sources.use = "M2",    # 发送者
  targets.use = sort(c("M1","TAM","CD4 T","CD8 T","NK","Mreg","B cell")),   # 接收者
  pairLR.use = pairLR,    # 指定基因对
  sort.by.target = TRUE,         # 按接收者排序
  title.name = "M2 Interactions"                    # X轴标签旋转角度
)


df.net.mreg <- subsetCommunication(cellchat, 
                              sources.use = "Mreg", 
                              targets.use = c("M1","TAM","CD4 T","CD8 T","NK","M2","B cell"))
# 绘制TAMs发送的信号通路（如VEGF）
pairLR = data.frame(interaction_name=df.net$interaction_name)
netVisual_bubble(
  cellchat,
  sources.use = "Mreg",    # 发送者
  targets.use = sort(c("M1","TAM","CD4 T","CD8 T","NK","M2","B cell")),   # 接收者
  pairLR.use = pairLR,    # 指定基因对
  sort.by.target = TRUE,         # 按接收者排序
  title.name = "Mreg Interactions"                    # X轴标签旋转角度
)

dev.off()

group_names <- levels(cellchat@idents)  # 获取所有细胞群名称
tam_index <- which(group_names == "TAM")  # 假设TAM的标签为"TAM"
df.net_tam_sender <- subsetCommunication(cellchat,
  sources.use = "TAM",  # TAM作为发送者
  targets.use = group_names[-tam_index]  # 所有非TAM细胞作为接收者
)
df.net_tam_receiver <- subsetCommunication(cellchat,
  sources.use = group_names[-tam_index],  # 所有非TAM细胞作为发送者
  targets.use = "TAM"  # TAM作为接收者
)
par(mfrow = c(1,2))


# 计算网络中心性（必须步骤）
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")


# 绘制全局交互网络
pdf("F:/workspace/TME_project_file/单细胞/TAM_cell_chat_result.pdf",8,8)
netVisual_circle(
  cellchat@net$count,
  sources.use = "TAM",
  targets.use = c("M1","M2","CD4 T","CD8 T","B cell","NK","Mreg"),
  title.name = "TAM Interaction Network"  # 正确参数名
)
dev.off()



# 需要先计算信号通路的层级结构
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# 绘制层级图
pdf("F:/workspace/TME_project_file/单细胞/TAM_cell_chat_gene_result.pdf",6,6)
netVisual_aggregate(
  cellchat,
  signaling = "CCL",
  layout = "hierarchy",       # 层级布局
  vertex.receiver = c(1,2)  # 接收细胞群的索引（根据group_names顺序）
)
dev.off()
#macrophage
macro_cell = subset(seurat_obj, idents = c(7,8,11))
macro_cell <- FindClusters(macro_cell, resolution = 0.1)
macro_cell <- RunUMAP(macro_cell, dims = 1:20)
pdf("F:/workspace/TME_project_file/单细胞/macro_subset_plot.pdf",4,4)
DimPlot(macro_cell, reduction = "umap", label = TRUE) + 
  ggtitle("UMAP Clustering") + theme_minimal()
dev.off()


pdf("F:/workspace/TME_project_file/单细胞/CXCL9_macro_feature_plot.pdf",6.5,6)
FeaturePlot(macro_cell, features = "CXCL9")
dev.off()


markers <- FindAllMarkers(
  object = macro_cell,
  assay = "RNA",           # 使用哪个assay
  slot = "data",           # 使用counts/data/scale.data
  test.use = "wilcox",     # 检验方法：wilcox/MAST/roc等
  latent.vars = "nCount_RNA", # 考虑技术偏差
  pseudocount.use = 1,     # 伪计数
  verbose = TRUE           # 显示进度
)


markers_m2_GZMB <- FindMarkers(macro_cell, 
                      ident.1 = 1, 
                      ident.2 = 2,
                      test.use = "wilcox", # 使用Wilcoxon检验
                      min.pct = 0.1,      # 在任一cluster中至少10%细胞表达的基因
                      logfc.threshold = 0.25) # logFC阈值

# 3. 查看结果
head(markers_m2_GZMB[order(markers_m2_GZMB$p_val_adj), ]) # 按调整p值排序

# 提取基因列表和log2FC值
gene_list <- markers_m2_GZMB$avg_log2FC
names(gene_list) <- rownames(markers_m2_GZMB)
gene_list <- sort(gene_list, decreasing = TRUE) # 按log2FC排序
# 进行GSEA分析
ids <- bitr(names(gene_list), 
            fromType = "SYMBOL", 
            toType = "ENTREZID", 
            OrgDb = org.Hs.eg.db)

# 合并log2FC值并处理重复
gene_list_entrez <- gene_list[ids$SYMBOL]
names(gene_list_entrez) <- ids$ENTREZID
gene_list_entrez <- na.omit(gene_list_entrez)
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)
gsea_result <- gseKEGG(
    geneList = gene_list_entrez,
    organism = "hsa", # 人类，其他物种参见KEGG代码
    keyType = "ncbi-geneid",
    nPerm = 1000,    # 置换次数
    minGSSize = 10,   # 最小基因集大小
    maxGSSize = 500,  # 最大基因集大小
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH", # 多重检验校正方法
    verbose = FALSE
)

# 查看结果
head(gsea_result@result)
tt = gsea_result@result
tt%>%dplyr::filter(Description%in%c("TNF signaling pathway","NOD-like receptor signaling pathway","Toll-like receptor signaling pathway","C-type lectin receptor signaling pathway","Complement and coagulation cascades","Cytokine-cytokine receptor interaction"))->tt_filter
library(enrichplot)
pdf("F:/workspace/TME_project_file/单细胞/GZMB_macropahge_gene_enrich.pdf",6,6)
gseaplot2(gsea_result, geneSetID = tt_filter$ID,title="Anti-tumor pathways")
dev.off()

 # 也可用"dot"
# ego <- enrichGO(gene = sapply(row.names(top5_list[[3]]),function(x) as.vector(unlist(strsplit(x,"[.]")))[1]), 
#                OrgDb = "org.Hs.eg.db", 
#                keyType = "SYMBOL")
markers_m2_GZMB%>%dplyr::filter(avg_log2FC<(-1))%>%dplyr::filter(p_val_adj<0.05)->GZMB_up
gene_ls = row.names(GZMB_up)
gene.df <- bitr(gene_ls, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

ego <- enrichKEGG(gene = gene.df$ENTREZID, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.05)               
dotplot(ego)

ego <- enrichGO(gene = row.names(GZMB_up),
               OrgDb = "org.Hs.eg.db", 
               keyType = "SYMBOL")
dotplot(ego)              


top5_list <- lapply(split(markers, markers$cluster), function(x) {
  head(x[order(x$p_val_adj, -x$avg_log2FC), ], 5)
})
marker_top5=c()
for(i in seq(0,3)){
  marker_top5 = c(marker_top5,top5_list[[i+1]][,"gene"])
}
pdf("F:/workspace/TME_project_file/单细胞/macro_marker.pdf",8,6)
DotPlot(macro_cell_new, 
        features = unique(marker_top5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
new_cluster_ids <- c(
  "0" = "M2",
  "1" = "M1",
  "2" = "Mreg",
  "3" = "TAM"
)
macro_cell$cell_type <- plyr::revalue(as.character(Idents(macro_cell)), new_cluster_ids)
pdf("F:/workspace/TME_project_file/单细胞/macro_subtype_plot.pdf",4.5,4)
DimPlot(macro_cell_new, group.by = "cell_type", label = TRUE) + 
  ggtitle("Cell Type Annotation")

dev.off()
library(ggplot2)
library(clusterProfiler)
library("org.Hs.eg.db")
top_list <- lapply(split(markers, markers$cluster), function(x) {
  head(x[order(x$p_val_adj, -x$avg_log2FC), ], 1000)
})
ego <- enrichGO(gene = sapply(row.names(top5_list[[3]]),function(x) as.vector(unlist(strsplit(x,"[.]")))[1]), 
               OrgDb = "org.Hs.eg.db", 
               keyType = "SYMBOL")
gene_ls = sapply(row.names(top5_list[[3]]),function(x) as.vector(unlist(strsplit(x,"[.]")))[1])
gene.df <- bitr(gene_ls, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

ego <- enrichKEGG(gene = gene.df$ENTREZID, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.05)               
dotplot(ego)
item=data.frame(ego)

item$GeneRatio<-item$Count/4

pdf("F:/workspace/TME_project_file/单细胞/TAM_enrich.pdf",8,5)
ggplot(item[1:10,],aes(x=GeneRatio,y=reorder(Description,GeneRatio)))+
  geom_point(aes(size=Count,color=p.adjust))+
  scale_color_gradient(low = "#f85a40", high = "#7552cc")+
  theme_bw()+theme(axis.text=element_text(size = 16))+
  scale_size_continuous(name="Count",
                        range = c(1,10))+
  labs(y="Terms")
dev.off()

macro_cell_new = macro_cell
Idents(macro_cell_new) = "cell_type"
pdf("F:/workspace/TME_project_file/单细胞/macro_CXCL9_sub_plot.pdf",4.5,4)
FeaturePlot(macro_cell_new,features = "CXCL9")
dev.off()
FeaturePlot(macro_cell_new,features = genes_ls)
library(monocle3)
cds <- new_cell_data_set(
  expression_data = GetAssayData(macro_cell_new, assay = "RNA", slot = "counts"),
  cell_metadata = macro_cell_new@meta.data,
  gene_metadata = data.frame(gene_short_name = rownames(macro_cell_new), 
                            row.names = rownames(macro_cell_new))
)
save.image(file="F:/workspace/TME_project_file/scRNA_analysis_cellchat.rda",compress=TRUE,ascii = FALSE)


load("F:/workspace/TME_project_file/scRNA_analysis_cellchat.rda")

library("Matrix")
library("Seurat")
library("monocle3")
library("ggplot2")
# 预处理
cds <- preprocess_cds(cds, 
                     method = "PCA", 
                     num_dim = 50, 
                     norm_method = "log", 
                     pseudo_count = 1)
cds <- reduce_dimension(cds, 
                       preprocess_method = "PCA", 
                       reduction_method = "UMAP",  # 必须指定为 UMAP
                       umap.min_dist = 0.1)

cds <- cluster_cells(cds,
                    resolution = 1e-5,  # 低分辨率聚焦主轨迹
                    cluster_method = "leiden",
                    random_seed = 42)
# 学习轨迹图结构
cds <- learn_graph(cds, 
                  use_partition = FALSE, 
                  close_loop = FALSE,
                  learn_graph_control = list(ncenter = 100))

#根据umap图选择根节点
# 绘制基础轨迹图（未指定伪时间）
plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           cell_size = 1)
# 交互式选择根节点
cds <- order_cells(cds, reduction_method = "UMAP")

# 操作流程：
# 1. 运行上述代码后，R 会弹出 UMAP 图窗口
# 2. 点击轨迹起点附近的节点（如最早发育状态的细胞群）
# 3. 按 [Esc] 或关闭窗口完成选择
# 选择轨迹根节点（根据M1/M2标记表达）
# 轨迹图（按伪时间着色）
pdf("F:/workspace/TME_project_file/单细胞/pseudotime_plot.pdf",6.5,6)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 3,
           cell_size = 1) +
  scale_color_viridis_c()
dev.off()


# 轨迹图（按伪时间着色）
plot_cells(cds,
           color_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 3,
           cell_size = 1)


pdf("F:/workspace/TME_project_file/单细胞/CXCL9_pseudotime_exp.pdf",5,3)
plot_genes_in_pseudotime(cds[c("CXCL9","CXCL10","CXCL11","CCL17","CCL18","CCL22"), ], 
                         color_cells_by = "pseudotime",
                         ncol = 1)
dev.off()
# 绘制并应用颜色
pdf("F:/workspace/TME_project_file/单细胞/cluster_pseudotime_exp.pdf",6.5,6)

plot_cells(cds, 
           color_cells_by = "cell_type",
           label_cell_groups = TRUE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 3,
           cell_size = 1) 
dev.off()


m1_genes <- c("CD86", "CXCL9","IRF5")
m2_genes <- c("CD163", "CCL18","IRF4")
pdf("F:/workspace/TME_project_file/单细胞/M1_M2_pseudotime_exp.pdf",4,5)
plot_genes_in_pseudotime(cds[m1_genes, ], 
                         color_cells_by = "pseudotime",
                         ncol = 1)
plot_genes_in_pseudotime(cds[m2_genes, ], 
                         color_cells_by = "pseudotime",
                         ncol = 1)  
plot_genes_in_pseudotime(cds[c("IGJ","PTPRS","GZMB")], 
                         color_cells_by = "pseudotime",
                         ncol = 1)
dev.off()
                     
pdf("F:/workspace/TME_project_file/单细胞/cellchat_gene_pseudotime.pdf",10,10)
plot_genes_in_pseudotime(cds[genes[1:5], ], 
                         color_cells_by = "pseudotime",
                         ncol = 1)
plot_genes_in_pseudotime(cds[genes[6:10], ], 
                         color_cells_by = "pseudotime",
                         ncol = 1)
plot_genes_in_pseudotime(cds[genes[11:15], ], 
                         color_cells_by = "pseudotime",
                         ncol = 1)
plot_genes_in_pseudotime(cds[genes[16:20], ], 
                         color_cells_by = "pseudotime",
                         ncol = 1)
plot_genes_in_pseudotime(cds[genes[21:28], ], 
                         color_cells_by = "pseudotime",
                         ncol = 1)
dev.off()
#患者信息补充
cell_info = read.table("F:/workspace/TME_project_file/单细胞/single_cell_info",header=T,sep="\t")
sample_cell=table(cell_info[,5])
cell_sub = Idents(macro_cell_new)
library(dplyr)
cell_id_macro = cell_info%>%dplyr::filter(title%in%names(cell_sub))
sample_cell_count = table(cell_info[,5])
sample_response = unique(cell_info[,c(5,6)])
row.names(sample_response) = sample_response[,1]
new_df = data.frame(id=cell_id_macro[,"title"],cell_type=cell_sub[cell_id_macro[,"title"]],source=cell_id_macro[,5])
sub_count = table(new_df$cell_type,new_df$source)
sub_ratio = data.frame(apply(sub_count,1,function(x) x/sample_cell))

sub_ratio$response = sample_response[row.names(sub_ratio),2]
sub_ratio$id = row.names(sub_ratio)
sub_ratio$point = sapply(sub_ratio$id,function(x) as.vector(unlist(strsplit(x,"_")))[[1]])

library(ggpubr)
pdf("F:/workspace/TME_project_file/单细胞/sample_cell_boxplot.pdf",3.5,4)
for(i in seq(1,4)){
    df1 = data.frame(cell=sub_ratio[,i],group=sub_ratio[,5])
    df2 = data.frame(cell=sub_ratio[,i],group=sub_ratio[,7])
    p1 = ggplot(df1,aes(x=group,y=cell))+geom_violin(aes(fill=group,color=group))+geom_boxplot(width=0.1)+scale_color_brewer(palette="Pastel1")+scale_fill_brewer(palette="Pastel1",)+stat_compare_means(method="wilcox.test",label = "p.format",color="black")+theme_classic()+labs(x="",y="",title=colnames(sub_ratio)[i])+theme(axis.text = element_text(color="black",size=12),legend.position="none")
    p2 = ggplot(df2,aes(x=group,y=cell))+geom_violin(aes(fill=group,color=group))+geom_boxplot(width=0.1)+scale_color_brewer(palette="Pastel1")+scale_fill_brewer(palette="Pastel1",)+stat_compare_means(method="wilcox.test",label = "p.format",color="black")+theme_classic()+labs(x="",y="",title=colnames(sub_ratio)[i])+theme(axis.text = element_text(color="black",size=12),legend.position="none")
    print(p1)
    print(p2)
}
dev.off()