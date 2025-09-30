#TME Clust netword plot
# random graph
library(GGally)
library(ggnet)
library(network)
library(sna)
library(ggplot2)
load("E:/project/remodeling tumor microenvironment/TME_feature_network.rda")
net = rgraph(10, mode = "graph", tprob = 0.5)
net = network(net, directed = FALSE)

# vertex names
network.vertex.names(net) = letters[1:10]
net %v% "phono" = ifelse(letters[1:10] %in% c("a", "e", "i"), "vowel", "consonant")
net %v% "color" = ifelse(net %v% "phono" == "vowel", "steelblue", "tomato")
ggnet2(net, color = "color")

ggnet2(net, color = "phono", palette = c("vowel" = "steelblue", "consonant" = "tomato"))
#ggnet2(net, color = ifelse(net %v% "phono" == "vowel", "steelblue", "tomato"))

jacard_mt_fun(pathway_nes_binary)->Jaccard_mt_update
hc = hclust(as.dist(1-Jaccard_mt_update), method = "complete")
n_clust = choose_clusters(Jaccard_mt_update, "jaccard", range = 2:(nrow(Jaccard_mt_update) - 1))
  #print(n_clust)
  #n_clust = 5
clust =  hclusCut(Jaccard_mt_update, n_clust)$cluster
new_ft<-c()
  #n_index<-intersect(names(which(table(clust)>2)),names(which(table(clust)<11)))
n_index<-names(which(table(clust)>2))
for(i in n_index){
    new_ft<-c(new_ft,(names(which(clust==i))))
}
  #Jaccard_mt_update<-Jaccard_mt_update[row.names(pathway_nes_binary_filter),row.names(pathway_nes_binary_filter)]
ft<-clust_item
color_ls = as.vector(unlist(sapply(names(ft),function(x) as.vector(unlist(strsplit(x,"\\|")))[1])))
cols = c(brewer.pal(8, "Set2"), brewer.pal(12, "Paired"), brewer.pal(8, "Pastel1"), brewer.pal(8, "Set3"))
names(cols)=color_ls
Jaccard_mt_update<-Jaccard_mt_update[new_ft,new_ft]
table(final_clust)
plot_ls<-list()
plot_index=0
for(i in c(2,4,5,6,8,10,12)){
    plot_index<-plot_index+1
    item<-names(clust_item)[which(clust_item==i)]
    test1<-Jaccard_mt_update[item,item]
    net = network(test1, directed = FALSE)
    # vertex names
    #network.vertex.names(net) = letters[1:10]
    feature_path<-as.vector(unlist(as.vector(unlist(sapply(item,function(x) as.vector(unlist(strsplit(x,"\\|")))[1])))))
    # net %v% "color" = as.vector(unlist(cols[feature_path]))
    # net %v% "label" = as.vector(unlist(sapply(item,function(x) as.vector(unlist(strsplit(x,"\\|")))[2])))

    # ggnet2(net,label = "label", color = "color")

    test1_weighted = network(test1,
                matrix.type = "adjacency",
                ignore.eval = FALSE,
                names.eval = "weights")
    test1_weighted %v% "color" = as.vector(unlist(cols[feature_path]))
    test1_weighted %v% "label" = item#as.vector(unlist(sapply(item,function(x) as.vector(unlist(strsplit(x,"\\|")))[2])))
    p<-ggnet2(test1_weighted, label = "label", color = "color",edge.size = "weights",palette = as.vector(unlist(cols[feature_path])))
    plot_ls[[plot_index]]<-p
}
library(easyGgplot2)
#ggnet2(test1_weighted, label = TRUE)
pdf("E:/project/remodeling tumor microenvironment/Cluster_clinical_info_analysis/TME_feature_network_re_color.pdf",12,6)
ggplot2.multiplot(plotlist=plot_ls, cols=4)
#ggnet2(test1_weighted,label = "label", color = "color",palette = "Set2")
dev.off()
#self maed guide
guide_df<-data.frame(Lable=color_ls,value=rep(1,length(color_ls)),color=cols)
ggplot(guide_df,aes(x=Lable,y=value))+geom_point(aes(color=Lable))+scale_color_manual(values=cols)


save.image(file="E:/project/remodeling tumor microenvironment/TME_feature_network.rda",ascii=FALSE,compress=TRUE )

#cluster signature gene
file_ls<-list.files("E:/project/remodeling tumor microenvironment/TCGA_TPM")
exp1<-read.table("E:/project/remodeling tumor microenvironment/TCGA_TPM/TCGA-ACC_exp_tpm.txt",header=T,sep="\t",check.names=FALSE)
cancer_exp_mt<-exp1[,intersect(colnames(exp1),names(TME_clust_assign_final))]
for(fl in file_ls[-1]){
  exp<-read.table(paste("E:/project/remodeling tumor microenvironment/TCGA_TPM/",fl,sep=""),header=T,sep="\t",check.names=FALSE)
  gene<-intersect(row.names(cancer_exp_mt),row.names(exp))
  sam<-intersect(colnames(exp),names(TME_clust_assign_final))
  cancer_exp_mt<-cbind(cancer_exp_mt,exp[gene,sam])
}
cancer_exp_mt_filter<-cancer_exp_mt[,included_sample]
diff_gene_ls<-list()
diff_gene_ls_sig<-list()
for(i in c(2,4,5,6,8,10,12)){
  wilcox_p_va<-c()
  clust_id<-paste("TME_clust_",i,sep="")
  sam<-intersect(included_sample,names(which(TME_clust_assign_final==clust_id)))
  sam_other<-setdiff(included_sample,sam)
  for(g in row.names(cancer_exp_mt_filter)){
    a<-as.vector(unlist(cancer_exp_mt_filter[g,sam]))
    b<-as.vector(unlist(cancer_exp_mt_filter[g,sam_other]))
    re<-wilcox.test(a,b)
    wilcox_p_va<-c(wilcox_p_va,round(re$p.value,4))
  }
  names(wilcox_p_va)<-row.names(cancer_exp_mt_filter)
  diff_gene_ls[[clust_id]]<-wilcox_p_va
  diff_gene_ls_sig[[clust_id]]<-names(which(wilcox_p_va<0.05))
}
save.image(file="E:/project/remodeling tumor microenvironment/TME_feature_network0805.rda",ascii=FALSE,compress=TRUE )