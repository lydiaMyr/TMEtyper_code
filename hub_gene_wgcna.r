load("F:/workspace/nianmomianyi (1)/nianmomianyi/TCGA_cancer_exp_mt_tpm_final.rda")
load("F:/workspace/TME_project_file/cluster_assign_new_tumor.rda")
c1_sam=names(which(cluster_assign_new=="TME1"))
c2_sam=names(which(cluster_assign_new=="TME2"))
c3_sam=names(which(cluster_assign_new=="TME3"))
c4_sam=names(which(cluster_assign_new=="TME4"))
c5_sam=names(which(cluster_assign_new=="TME5"))
c6_sam=names(which(cluster_assign_new=="TME6"))
c7_sam=names(which(cluster_assign_new=="TME7"))
  #protein  coding 

# 安装并加载包

library(org.Hs.eg.db)

all_genes <- keys(org.Hs.eg.db, keytype = "SYMBOL")
gene_info <- select(org.Hs.eg.db, 
                   keys = all_genes,
                   keytype = "SYMBOL",
                   columns = c("ENSEMBL", "GENETYPE","SYMBOL"))

# 筛选蛋白质编码基因
protein_coding_genes <- gene_info[gene_info$GENETYPE == "protein-coding", ]
head(protein_coding_genes)
gene_ls = protein_coding_genes$SYMBOL
TME1_exp = c()
TME2_exp = c()
TME3_exp = c()
TME4_exp = c()
TME5_exp = c()
TME6_exp = c()
TME7_exp = c()

for(ca in names(exp_tpm_cancer_ls_geneID)){
    exp = exp_tpm_cancer_ls_geneID[[ca]]
    row.names(exp) = exp[,1]
    exp = exp[,-1]
    exp1 = t(apply(exp,1,function(x) as.numeric(x)))
    colnames(exp1) = colnames(exp)
    exp_new = exp1[intersect(row.names(exp1),gene_ls),]
    s1 = intersect(c1_sam,colnames(exp_new))
    s2 = intersect(c2_sam,colnames(exp_new))
    s3 = intersect(c3_sam,colnames(exp_new))
    s4 = intersect(c4_sam,colnames(exp_new))
    s5 = intersect(c5_sam,colnames(exp_new))
    s6 = intersect(c6_sam,colnames(exp_new))
    s7 = intersect(c7_sam,colnames(exp_new))
    if(length(s1)>0){
        TME1_exp = cbind(TME1_exp,exp_new[,s1])
    }
    #  if(length(s2)>0){
    #     TME2_exp = cbind(TME2_exp,exp_new[,s2])
    # }
    #  if(length(s3)>0){
    #     TME3_exp = cbind(TME3_exp,exp_new[,s3])
    # }
    #  if(length(s4)>0){
    #     TME4_exp = cbind(TME4_exp,exp_new[,s4])
    # }
    #  if(length(s1)>0){
    #     TME5_exp = cbind(TME5_exp,exp_new[,s5])
    # }
    #  if(length(s1)>0){
    #     TME6_exp = cbind(TME6_exp,exp_new[,s6])
    # }
    #  if(length(s1)>0){
    #     TME7_exp = cbind(TME7_exp,exp_new[,s7])
    # }
}

TME_sub_exp = list(TME1 = TME1_exp, TME2 = TME2_exp,TME3 = TME3_exp,TME4 = TME4_exp,TME5 = TME5_exp,TME6 = TME6_exp,TME7 = TME7_exp )
save(TME_sub_exp,file="F:/workspace/TME_project_file/文章修改/TME_sub_exp.rda")


library(WGCNA)
# 过滤低表达基因
TME1_exp = TME1_exp[,-which(colnames(TME1_exp)=="")]
exp = log2(TME1_exp+1)
gsg <- WGCNA::goodSamplesGenes(exp, verbose = 3)
if (!gsg$allOK) {
  exp <- exp[gsg$goodGenes, gsg$goodSamples]
}

# 可选：去除异常样本
sampleTree <- hclust(dist(t(exp)), method = "average")
plot(sampleTree, main = "Sample clustering", sub="", xlab="")

# 数据转换(建议使用log2转换后的数据)
# 如果数据未标准化，可以执行：
# exp <- log2(exp + 1)
# 选择软阈值
powers <- c(1:20)
sft <- pickSoftThreshold(t(exp), powerVector = powers, verbose = 5)

# 绘制结果
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity")

# 根据上图选择合适的power值
softPower <- sft$powerEstimate  # 或手动指定如softPower <- 6

# 构建网络
net <- blockwiseModules(t(exp), power = softPower,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendrogram = FALSE,
                       saveTOMs = TRUE, saveTOMFileBase = "hub_TOM",
                       verbose = 3)

# 查看模块数量
table(net$colors)


# 转换为颜色标签
moduleColors <- labels2colors(net$colors)

# 绘制模块树
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)



# 计算模块特征基因(MEs)
MEs0 <- moduleEigengenes(t(exp), moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
hub_g="PAX8"
hub_expr <- exp[hub_g, ]  # 或STAT5B，根据实际基因名调整
traitData <- data.frame(hub = as.numeric(hub_expr))
rownames(traitData) <- colnames(exp)
# 计算模块与STAT5表达量的相关性
moduleTraitCor <- cor(MEs, traitData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, ncol(exp))

# 显示结果
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

# 绘制相关性热图
pdf("F:/workspace/TME_project_file/文章修改/TME1_reg_module.pdf")
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
              xLabels = "STAT5 expression",
              yLabels = names(MEs),
              ySymbols = names(MEs),
              colorLabels = FALSE,
              colors = blueWhiteRed(50),
              textMatrix = textMatrix,
              setStdMargins = FALSE,
              cex.text = 0.5,
              zlim = c(-1,1),
              main = "Module-trait relationships")
dev.off()

# 找出与STAT5最相关的模块(选择相关性最高且p值<0.05的模块)
tme_cor_modules <- which(moduleTraitPvalue < 0.05)  # 显著相关的模块
tme_cor_values <- moduleTraitCor[tme_cor_modules, 1]
id = row.names(moduleTraitPvalue)[tme_cor_modules]
# 选择相关性最强的模块
top_module <- id[which.max(abs(tme_cor_values))]
top_module_color <- gsub("ME", "", top_module)

# 提取该模块所有基因
module_genes <- rownames(exp)[moduleColors == top_module_color]

# 查看模块基因数量
length(module_genes)

# 保存结果
TME="TME1"
write.table(module_genes, paste("F:/workspace/TME_project_file/文章修改/",TME,"_associated_module_genes.txt",sep=""),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
