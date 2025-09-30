#jacard index hclust
load("E:/project/remodeling tumor microenvironment/NMF_test/NMF_procedure.rda")
literature_sig<-read.xlsx("E:/project/remodeling tumor microenvironment/NMF_test/TME_related_pathway_done.xlsx",sheet=1)
literature_geneSet<-list()
for(i in seq(1,nrow(literature_sig))){
  tmp<-as.vector(unlist(literature_sig[i,]))
  title<-tmp[1]
  genes<-unique(tmp[-1])
  if(length(grep(TRUE,is.na(genes)))>0){
    genes<-genes[-(which(is.na(genes)==TRUE))]
  }
  literature_geneSet[[title]]<-genes
}

load("E:/project/remodeling tumor microenvironment/pathway/leading_edge_filtered_pathway.rda")#final_path_list
#ls()
load("E:/project/remodeling tumor microenvironment/pathway/literature_keyword_pathway.rda")
pathway_nes_ensembl<-c()
path_ensembl_name<-c()
for(key in names(keyword_path_list)){
    path<-intersect(keyword_path_list[[key]],row.names(pahtway_nes))
    if(length(path)==1){
      path_ensembl_name<-c(path_ensembl_name,key)
      pathway_nes_ensembl<-rbind(pathway_nes_ensembl,pahtway_nes[path,])
    }
    if(length(path)>1){
      path_ensembl_name<-c(path_ensembl_name,key)
      tmp<-apply(pahtway_nes[path,],2,function(x) (median(x)))
      pathway_nes_ensembl<-rbind(pathway_nes_ensembl,tmp)
    }
}
row.names(pathway_nes_ensembl)<-path_ensembl_name

pathway_item_liter<-read.table("E:/project/remodeling tumor microenvironment/NMF_test/literature_item_list.txt",header=F,sep="\t")

literature_geneSet_ense<-list()
for(i in seq(1,nrow(pathway_item_liter))){
  tmp<-as.vector(unlist(pathway_item_liter[i,]))
  title<-tmp[1]
  genes<-unique(tmp[-1])
  if(length(genes)<4){
    genes<-genes[-(length(genes))]
  }
  if(length(genes)==4){
    if(nchar(genes[4])==0)
      genes<-genes[-(length(genes))]
  }
  literature_geneSet_ense[[title]]<-genes
}

#pathway binary
keyword_path_list_all<-c(keyword_path_list,literature_geneSet_ense)
#write.table(names(keyword_path_list),file="E:/project/remodeling tumor microenvironment/NMF_test/GSVA_path_item.txt",sep="\t",quote=F)
keyword_path_list_new<-keyword_path_list_all
keyword_path_list_new[["chemokine"]]<-unique(c(keyword_path_list_new[["chemokine"]],keyword_path_list_new[["chemokines"]]))
keyword_path_list_new[["cytokine"]]<-unique(c(keyword_path_list_new[["cytokine"]],keyword_path_list_new[["cytokines"]]))
keyword_path_list_new[["Epithelial"]]<-unique(c(keyword_path_list_new[["Epithelial"]],keyword_path_list_new[["epithelial"]]))
keyword_path_list_new[["Mesenchymal"]]<-unique(c(keyword_path_list_new[["Mesenchymal"]],keyword_path_list_new[["mesenchymal"]]))
keyword_path_list_new[["mTOR"]]<-unique(c(keyword_path_list_new[["mTOR"]],keyword_path_list_new[["MTOR"]]))
keyword_path_list_new[["vessel"]]<-unique(c(keyword_path_list_new[["vessel"]],keyword_path_list_new[["vessels"]]))
keyword_path_list_new1<-list()
for(key in names(keyword_path_list_new)){
    if(!(key%in%c("chemokines","cytokines","epithelial","mesenchymal","MTOR","vessels")))
    keyword_path_list_new1[[key]]<-keyword_path_list_new[[key]]
}
#length(keyword_path_list_new1)
keyword_path_list_all<-keyword_path_list_new1

#keyword_path_list_all<-literature_geneSet_ense
pathway_nes_binary<-c()
path_ensembl_name<-c()
pahtway_nes_all<-rbind(pahtway_nes,pahtway_nes_liter)
for(key in names(keyword_path_list_all)){
    path<-intersect(keyword_path_list_all[[key]],row.names(pahtway_nes_all))
    if(length(path)>2){
     # path_ensembl_name<-c(path_ensembl_name,key)
      tmp<-apply(pahtway_nes_all[path,],2,function(x) {x[which.max(x)]=1;x[which(x!=1)]=0;return(x)})
      row.names(tmp)<-paste(key,row.names(tmp),sep="|")
      pathway_nes_binary<-rbind(pathway_nes_binary,tmp)
    }
}
apply(pathway_nes_binary,1,function(x) sum(x))->pathway_result_static


pathway_nes_binary_filter<-pathway_nes_binary[names(which(pathway_result_static>10)),]
#pheatmap::pheatmap(t(pathway_nes_binary_filter))
cell_state_fl<-list.files("E:/project/remodeling tumor microenvironment/NMF_test/ecotyper_output_skcm/Carcinoma_Cell_States")
cell_state_fl<-cell_state_fl[grep("txt",cell_state_fl)]
cell_state_result<-c()
cell_state_binary<-c()
for(fl in cell_state_fl){
  ct<-read.table(paste("E:/project/remodeling tumor microenvironment/NMF_test/ecotyper_output_skcm/Carcinoma_Cell_States/",fl,sep=""),header=T,sep="\t",row.names = 1)
  cell_name<-as.vector(unlist(strsplit(fl,"_")))[1]
  colnames(ct)<-paste(cell_name,colnames(ct),sep="_")
  cell_state_result<-rbind(cell_state_result,t(ct))

  tmp<-apply(ct,1,function(x) {x[which.max(x)]=1;x[which(x!=1)]=0;return(x)})
  cell_state_binary<-rbind(cell_state_binary,tmp)

}

sam_name<-colnames(cell_state_binary)
sapply(sam_name,function(x) gsub("\\.","-",x))->sam_name_new
colnames(cell_state_binary)<-sam_name_new
binary_feature<-rbind(pathway_nes_binary_filter,cell_state_binary[,colnames(pathway_nes_binary_filter)])

library(pheatmap)
pheatmap(cell_state_binary)

# library(Matrix)
# jaccard <- function(m) {
#     ## common values:
#     A = tcrossprod(m)
#     ## indexes for non-zero common values
#     im = which(A > 0, arr.ind=TRUE)
#     ## counts for each row
#     b = rowSums(m)

#     ## only non-zero values of common
#     Aim = A[im]

#     ## Jacard formula: #common / (#i + #j - #common)
#     J = sparseMatrix(
#           i = im[,1],
#           j = im[,2],
#           x = Aim / (b[im[,1]] + b[im[,2]] - Aim),
#           dims = dim(A)
#     )

#     return( J )
# }
# #install.packages("ade4"
# jaccard(binary_feature)->re_test1
# library(ade4)
# dist.binary(binary_feature, method = 1, diag = FALSE, upper = FALSE)->re_test

library(jaccard)
#jaccard(binary_feature, y, center = FALSE, px = NULL, py = NULL)
#calculation process
Jaccard_mt<-c()
for(i in seq(1,nrow(binary_feature))){
  apply(binary_feature,1,function(x) jaccard::jaccard(as.numeric(x),as.numeric(binary_feature[i,])))->re
  Jaccard_mt<-rbind(Jaccard_mt,re)
  #Jaccard_mt_update[i,which(re>0.01)]=0
}
row.names(Jaccard_mt)<-row.names(binary_feature)
#apply(binary_feature,1,function(x) jaccard_similarity(x,binary_feature[1,]))->re

#apply(binary_feature,1,function(x) phyper_fun(x,binary_feature[1,]))->re




phyper_fun<-function(x,y){
  stat<-x+y
  all<-length(x)
  a<-sum(x)
  b<-sum(y)
  c<-all-max(a,b)
  d<-length(which(stat==2))
  pva<-1-phyper(d-1,max(a,b),c,min(a,b))
  return(pva)
}

Jaccard_mt_update<-Jaccard_mt
for(i in seq(1,nrow(binary_feature))){
  apply(binary_feature,1,function(x) phyper_fun(x,binary_feature[i,]))->re
  Jaccard_mt_update[i,which(re>0.05)]=0
  #Jaccard_mt_update[which(re>0.01,i)]=0
}

length(which(Jaccard_mt_update>0))
length(which(Jaccard_mt>0))
#length(Jaccard_mt)
Jaccard_mt_dist<-dist(Jaccard_mt_update)
hclust(Jaccard_mt_dist,method = "average")->cluster_re

#plot(cluster_re)
pdf("E:/project/remodeling tumor microenvironment/NMF_test/hclust_plot.pdf",16,12)
plot(cluster_re)
#确定最佳聚类个数为7
cut_avg <- cutree(cluster_re, k = 4)
rect.hclust(cluster_re , k = 4)
abline(h = -1, col = 'red')
dev.off()


table(cluster_re$order)
library(factoextra)
library(NbClust)

write.table(Jaccard_mt_update,file="E:/project/remodeling tumor microenvironment/NMF_test/Jaccard_mt_update.txt",sep="\t",quote=F)
write.table(Jaccard_mt,file="E:/project/remodeling tumor microenvironment/NMF_test/Jaccard_mt.txt",sep="\t",quote=F)
pdf("E:/project/remodeling tumor microenvironment/NMF_test/subtitle_cluster.pdf",4,4)
fviz_nbclust(Jaccard_mt_update[row.names(cell_state_binary),row.names(cell_state_binary)], kmeans, method = "silhouette",k.max=15)+
  labs(subtitle = "Silhouette method")
dev.off()
test<-Jaccard_mt_update
diag(test)<-0
pheatmap(test[row.names(cell_state_binary),row.names(cell_state_binary)],clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean", clustering_method = "average",color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
fviz_nbclust(Jaccard_mt_update, hcut, method = "wss",k.max=15) +
    geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")
set.seed(123)
fviz_nbclust(Jaccard_mt_update, hcut, nstart = 25,  method = "gap_stat", nboot = 50, k.max=30)+
  labs(subtitle = "Gap statistic method")

#pheatmap::pheatmap(t(pathway_nes_binary_filter))


#LUAD 测试
ecotyper_re<-read.table("E:/project/remodeling tumor microenvironment/NMF_test/TCGA_cell_state_assignments.tsv",header=T,sep="\t")
