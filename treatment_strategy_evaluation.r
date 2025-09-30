load("F:/workspace/nianmomianyi (1)/nianmomianyi/TCGA_cancer_exp_mt_tpm_final.rda")
exp_mt_new = c()
com_gene = read.table("F:/workspace/TME_project_file/文章修改/estimate_gene_ls.txt",header=T)

for(i in names(exp_tpm_cancer_ls_geneID)){
    exp = exp_tpm_cancer_ls_geneID[[i]]
    row.names(exp) = exp[,1]
    exp = exp[,-1]
    genes = intersect(row.names(exp),com_gene$Gene)
    exp_mt_new = cbind(exp_mt_new,exp[genes,])
}
exp_mt_new1 = t(apply(exp_mt_new,1,function(x) as.numeric(x)))
colnames(exp_mt_new1) = colnames(exp_mt_new)

input_file <- "F:/workspace/TME_project_file/input_expression.gct"
write.table(exp_mt_new1, file=input_file, sep="\t", quote=FALSE, col.names=NA)


# 运行ESTIMATE评分
output_score <- "F:/workspace/TME_project_file/estimate_scores1.gct"
estimateScore(input.ds = input_file,
              output.ds = output_score,
              platform = "illumina")  # 或 "illumina"

# 读取结果
scores <- read.table(output_score, skip=2, header=TRUE, row.names=1,check.names = F)
scores <- t(scores)  # 转置为样本x评分
#row.names(scores) = colnames(exp_mt_new1)


load("F:/workspace/TME_project_file/cluster_assign_new_tumor.rda")

# in.file <- system.file("extdata", "sample_input.gct", package="estimate")
# out.file <- tempfile(pattern="estimate", fileext=".gct")
# estimateScore(in.file, out.file)
#PAX8 CDX2 PLA2G2D DSG3 CSF3R LRP2 GPM6B
c1_sam=names(which(cluster_assign_new=="TME1"))
c2_sam=names(which(cluster_assign_new=="TME2"))
c3_sam=names(which(cluster_assign_new=="TME3"))
c4_sam=names(which(cluster_assign_new=="TME4"))
c5_sam=names(which(cluster_assign_new=="TME5"))
c6_sam=names(which(cluster_assign_new=="TME6"))
c7_sam=names(which(cluster_assign_new=="TME7"))


exp_mt_gene = c()
com_gene = c("PAX8","CDX2","PLA2G2D","DSG3","CSF3R","LRP2","GPM6B","CD274","CTLA4","IFNG")
for(i in names(exp_tpm_cancer_ls_geneID)){
    exp = exp_tpm_cancer_ls_geneID[[i]]
    row.names(exp) = exp[,1]
    exp = exp[,-1]
    genes = intersect(row.names(exp),com_gene)
    exp_mt_gene = cbind(exp_mt_gene,exp[genes,])
}

cluster_hub_ls = list(PAX8=c1_sam,CDX2=c2_sam,PLA2G2D=c3_sam,
DSG3=c4_sam,CSF3R=c5_sam,LRP2=c6_sam,GPM6B=c7_sam)
immune_score=as.numeric(scores[,2][-1])
names(immune_score) = row.names(scores)[-1]
sam_id = names(immune_score)
sam_id_new = sapply(sam_id,function(x) gsub("[.]","-",x))
names(immune_score) = sam_id_new

ICP_score = rbind(as.numeric(exp_mt_gene["CD274",]),as.numeric(exp_mt_gene["CTLA4",]),as.numeric(exp_mt_gene["IFNG",]))
colnames(ICP_score) = colnames(exp_mt_gene)
row.names(ICP_score) = c("CD274","CTLA4","IFNG")

load("F:/workspace/TME_project_file/文章修改/sample_TMB.rda")
#drug target score
load("F:/workspace/TME_project_file/target vs drug/durg_gene_mut_fre_all.rda")
target_mut_fre = gene_fre_all[-1,]

score_mt = c()
for(hub_gene in names(cluster_hub_ls)){
    sam = cluster_hub_ls[[hub_gene]]
    hub_exp = exp_mt_gene[hub_gene,sam]
    im_score = immune_score[sam]
    check_score1 = ICP_score[1,sam]
    check_score2 = ICP_score[2,sam]
    check_score3 = ICP_score[3,sam]
    target_score = target_score_fun(sam)
    score_mt = rbind(score_mt,c(median(im_score),median(check_score1),median(check_score2),median(check_score3),target_score$mut_fre,target_score$tmb))
}
colnames(score_mt) = c("Immune Score","CD274","CTLA4","INFG","Drug target","TMB")
row.names(score_mt) = c("PAX8","CDX2","PLA2G2D","DSG3","CSF3R","LRP2","GPM6B")
score_mt_norm = scale(score_mt)
library(scales)
rescaled_data <- rescale(score_mt_norm, to = c(0, 1))
#绘制雷达图
plot_data = data.frame(rescaled_data)
library(ggradar)
plot_data = cbind(rep(0,7),plot_data)
plot_data[,1]=c("PAX8","CDX2","PLA2G2D","DSG3","CSF3R","LRP2","GPM6B")
colnames(plot_data)[1]="Group"
pdf("F:/workspace/TME_project_file/文章修改/雷达图.pdf",4,4)
ggradar(plot.data = plot_data[1,],group.line.width=0.5,group.point.size=2)
ggradar(plot.data = plot_data[2,],group.line.width=0.5,group.point.size=3)
ggradar(plot.data = plot_data[3,],group.line.width=0.5,group.point.size=3)
ggradar(plot.data = plot_data[4,],group.line.width=0.5,group.point.size=3)
ggradar(plot.data = plot_data[5,],group.line.width=0.5,group.point.size=3)
ggradar(plot.data = plot_data[6,],group.line.width=0.5,group.point.size=3)
ggradar(plot.data = plot_data[7,],group.line.width=0.5,group.point.size=3)
dev.off()




sam_mut = row.names(target_mut_fre)
sam_mut_part = sapply(sam_mut,function(x) paste(as.vector(unlist(strsplit(x,"-")))[1:3],collapse="-"))
names(sam_mut) = sam_mut_part

sam_mut1 = names(tmb_ls)
sam_mut_part2 = sapply(sam_mut1,function(x) paste(as.vector(unlist(strsplit(x,"-")))[1:3],collapse="-"))
names(sam_mut1) = sam_mut_part2
target_score_fun = function(sam_id){
    sam_id_part1 = sapply(sam_id,function(x) paste(as.vector(unlist(strsplit(x,"-")))[1:3],collapse="-"))
    sam_id = intersect(sam_id_part1,sam_mut_part)
    mt = target_mut_fre[sam_mut[sam_id],]
    sam_id = intersect(sam_id_part1,sam_mut_part2)
    tmb = median(tmb_ls[sam_mut1[sam_id]])
    mt_fre = apply(mt,1,function(x) sum(x))
    mt_fre_all = length(which(mt_fre>0))/length(sam_id)
    return(list(mut_fre=mt_fre_all,tmb=tmb))
}

library(ggradar)
my_radar<- data.frame(matrix(runif(33),
                             nrow =3,
                             dimnames =list(
                               rnames=c("模式1","模式2","模式3"),#设定行名
                               cnames=c("模式","蔬菜类","水果类","坚果类","豆制品","肉类及其制品",#设定列名
                                 "水产品","奶制品","禽蛋类","茶及饮料","主食类"))))
my_radar[,1]<-c("模式A","模式B","模式C")
my_radar
plot_data = cbind(rep(0,7),plot_data)
plot_data[,1]=c("PAX8","CDX2","PLA2G2D","DSG3","CSF3R","LRP2","GPM6B")
colnames(plot_data)[1]="Group"
ggradar(plot.data = plot_data[1,])
ggradar(plot.data = plot_data[1,])
