#pancancer test
#Immune cell fraction
library(openxlsx)
library(GSVA)
load("panCancer_immune_eco.rda")
load("literature_geneSet.rda")
load("pathway_leading_edge_result_update.rda")
#cancer_immune_pathway
load("cancer_related_pathway_update.rda")



library(dplyr)
filter_pathway_leading_edge<-list()
all_sig_pathway<-c()
pathway_rank<-rep(0,length(Cancer_immune_pathway))
names(pathway_rank)<-names(Cancer_immune_pathway)
for(ca in names(leading_edge_result)){
    data.frame(leading_edge_result[[ca]]%>%dplyr::filter(padj<0.05 ))->re
    filter_pathway_leading_edge[[ca]]<-re
    pathway_rank[re[,1]]<-pathway_rank[re[,1]]+1
    all_sig_pathway<-c(all_sig_pathway,re$pathway)
}
all_sig_pathway<-intersect(all_sig_pathway,names(immune_path_sig))
names(which(pathway_rank[all_sig_pathway]>2))->all_sig_pathway_filter

filter_pathway_leading_edge<-list()
all_sig_pathway<-c()
pathway_rank<-rep(0,length(immune_path_sig))
names(pathway_rank)<-names(immune_path_sig)
for(ca in names(leading_edge_result)){
    data.frame(leading_edge_result[[ca]]%>%dplyr::filter(padj<0.05))->re
    filter_pathway_leading_edge[[ca]]<-re
    pathway_rank[re[,1]]<-pathway_rank[re[,1]]+1
    all_sig_pathway<-c(all_sig_pathway,re$pathway)
}
all_sig_pathway<-intersect(all_sig_pathway,names(immune_path_sig))
names(which(pathway_rank[all_sig_pathway]>2))->all_sig_pathway_filter

keyword_path_list_filter<-list()
item_count<-c()
for(item in names(keyword_path_list_all)){
  tt<-intersect(keyword_path_list_all[[item]],all_sig_pathway_filter)
  ct<-length(tt)
  item_count<-c(item_count,ct)
  if(ct>2){
    keyword_path_list_filter[[item]]<-tt
  }
}
names(item_count)<-names(keyword_path_list_all)

path_gene_rank<-list()
for(pa in keyword_path_list_filter){
    gene_ls<-immune_path_sig[[pa]]
    gene_rank<-rep(0,length(gene_ls))
    names(gene_rank)<-gene_ls
    for(ca in names(leading_edge_result)){
        leading_info<-data.frame(leading_edge_result[[ca]])
        if(nrow(leading_info)>0){
            row.names(leading_info)<-as.vector(unlist(leading_info[,"pathway"]))
            if(pa%in%row.names(leading_info)){
                leading_gene<-as.vector(unlist(strsplit(as.vector(unlist(leading_info[pa,"leadingEdge"])),",")))
                gene_rank[leading_gene]<-gene_rank[leading_gene]+1
            }

        }
    }
    path_gene_rank[[pa]]<-gene_rank
}

#path_gene_rank
path_gene_filter<-list()
pa_filter_count<-c()
for(pa in names(path_gene_rank)){
    gene_rank<-path_gene_rank[[pa]]
    path_gene_filter[[pa]]<-names(sort(gene_rank,decreasing=T)[1:20])
    #pa_filter_count<-c(pa_filter_count,length(path_gene_filter[[pa]]))
}
final_path_list<-c(path_gene_filter,eco_immune_sig)


setdiff(as.vector(unlist(keyword_path_pos_neg_list_filter)),names(final_path_list))->miss_item
final_path_list<-c(final_path_list,literature_geneSet)

score_pos_neg_matrix<-c()
#keyword_path_list_pos_neg<-list()
for(item in names(keyword_path_list_filter)){
  path_item<-keyword_path_list_filter[[item]]
  sign_ls<-cbind(rep(item,length(path_item)),path_item)
  for(ca in names(leading_edge_result)){
    score<-data.frame(leading_edge_result[[ca]])
    row.names(score)<-score$pathway
    nes<-score[path_item,"ES"]
    sign_ls<-cbind(sign_ls,nes)
    #keyword_path_list_pos_neg[[paste(item,"pos",sep="_")]]<-unique(c(keyword_path_list_pos_neg[[paste(item,"pos",sep="_")]],keyword_path_list_all[[item]][which(nes>0)]))
    #keyword_path_list_pos_neg[[paste(item,"neg",sep="_")]]<-unique(c(keyword_path_list_pos_neg[[paste(item,"neg",sep="_")]],keyword_path_list_all[[item]][which(nes<0)]))
  }
  score_pos_neg_matrix<-rbind(score_pos_neg_matrix,sign_ls)
}
colnames(score_pos_neg_matrix)<-c("Keyword","path",names(leading_edge_result))
score_pos_neg_matrix_pheatmap<-score_pos_neg_matrix
score_pos_neg_matrix_pheatmap[which(score_pos_neg_matrix<0)]=-1
score_pos_neg_matrix_pheatmap[which(score_pos_neg_matrix>0)]=1
apply(score_pos_neg_matrix_pheatmap,1,function(x) length(which(abs(as.numeric(x))==1)))->tag_count
score_pos_neg_matrix_pheatmap_filter<-score_pos_neg_matrix[-which(tag_count==2),]
score_pos_neg_matrix_pheatmap_new<-t(apply(score_pos_neg_matrix_pheatmap_filter,1,function(x) as.numeric(x)))
row.names(score_pos_neg_matrix_pheatmap_new)<-score_pos_neg_matrix_pheatmap_filter[,2]
colnames(score_pos_neg_matrix_pheatmap_new)<-colnames(score_pos_neg_matrix_pheatmap_filter)
score_pos_neg_matrix_pheatmap_new[which(score_pos_neg_matrix_pheatmap_new<0)]=-1
score_pos_neg_matrix_pheatmap_new[which(score_pos_neg_matrix_pheatmap_new>0)]=1

for(item in names(keyword_path_list_filter)){
  path_item<-intersect(row.names(score_pos_neg_matrix_pheatmap_new),keyword_path_list_filter[[item]])
  if(length(path_item)>3){
     pheatmap(t(score_pos_neg_matrix_pheatmap_new[path_item,-c(1,2)]),cluster_rows=FALSE, cluster_cols=FALSE)
    #path_item<-path_item[-1]
  }
 }

keyword_path_list_pos_neg<-list()
for(item in names(keyword_path_list_filter)){
  path_item<-intersect(keyword_path_list_filter[[item]],row.names(score_pos_neg_matrix_pheatmap_new))
  score<-score_pos_neg_matrix_pheatmap_new[path_item,-c(1,2)]
  for(pa in path_item){
    es<-score[pa,]
    pos_count<-length(which(es>0))
    neg_count<-length(which(es<0))
    if(length(abs(es)>0)>0){
      if(pos_count>neg_count){
          keyword_path_list_pos_neg[[paste(item,"pos",sep="_")]]<-unique(c(keyword_path_list_pos_neg[[paste(item,"pos",sep="_")]],pa))
      }else{
          keyword_path_list_pos_neg[[paste(item,"neg",sep="_")]]<-unique(c(keyword_path_list_pos_neg[[paste(item,"neg",sep="_")]],pa))
      }
          
    }
  }

}

keyword_path_pos_neg_list_filter<-list()
item_count<-c()
for(item in names(keyword_path_list_pos_neg)){
  tt<-keyword_path_list_pos_neg[[item]]
  ct<-length(tt)
  item_count<-c(item_count,ct)
  if(ct>2){
    keyword_path_pos_neg_list_filter[[item]]<-tt
  }
}


#pathway_ls<-unique(as.vector(unlist(keyword_path_list)))
path_gene_rank<-list()
for(pa in all_sig_pathway){
    gene_ls<-Cancer_immune_pathway[[pa]]
    gene_rank<-rep(0,length(gene_ls))
    names(gene_rank)<-gene_ls
    for(ca in names(leading_edge_result)){
        leading_info<-data.frame(leading_edge_result[[ca]])
        if(nrow(leading_info)>0){
            row.names(leading_info)<-as.vector(unlist(leading_info[,"pathway"]))
            if(pa%in%row.names(leading_info)){
                leading_gene<-as.vector(unlist(strsplit(as.vector(unlist(leading_info[pa,"leadingEdge"])),",")))
                gene_rank[leading_gene]<-gene_rank[leading_gene]+1
            }

        }
    }
    path_gene_rank[[pa]]<-gene_rank
}
#path_gene_rank
path_gene_filter<-list()
pa_filter_count<-c()
for(pa in names(path_gene_rank)){
    gene_rank<-path_gene_rank[[pa]]
    path_gene_filter[[pa]]<-names(sort(gene_rank,decreasing=T)[1:20])
    #pa_filter_count<-c(pa_filter_count,length(path_gene_filter[[pa]]))
}
final_path_list<-path_gene_filter


eco_keyword_ls<-list()
eco_immune_sig<-list()
for(i in seq(1,nrow(ecotype_cell_signature))){
  ct<-ecotype_cell_signature[i,]
  cell<-as.vector(unlist(ct[1]))
  cell_state<-paste(cell,ct[3],sep="_")
  eco_keyword_ls[[cell]]<-unique(c(eco_keyword_ls[[cell]],cell_state))
  eco_immune_sig[[cell_state]]<-c(eco_immune_sig[[cell_state]],as.vector(unlist(ct[2])))
}



fl_ls<-list.files("E:/project/remodeling tumor microenvironment/FPKM_exp/")
ca_immune_fra<-list()
for(fl in fl_ls){
    ca<-as.vector(unlist(strsplit(as.vector(unlist(strsplit(fl,"_")))[1],"-")))[2]
    exp<-read.table(paste("E:/project/remodeling tumor microenvironment/FPKM_exp/",fl,sep=""),header=T,sep="\t",row.names = 1)
    exp<-as.matrix(exp)
    gsva(exp,pamr_final_marker_add,method="ssgsea")->immune_nes
    ca_immune_fra[[ca]]<-immune_nes
}
save.image(file="E:/project/remodeling tumor microenvironment/panCancer_immune.rda",ascii=FALSE,compress=TRUE )
load("E:/project/remodeling tumor microenvironment/panCancer_immune.rda")

immune_nes_all<-c()
for(ca in names(ca_immune_fra_eco)){
  immune_nes_all<-rbind(immune_nes_all,t(ca_immune_fra_eco[[ca]]))
}
immune_nes<-t(immune_nes_all)

# T_sub2_score<-immune_nes[T_sub2$Cell,]
# T_sub2_score_norm<-apply(T_sub2_score,2,function(x)(x/sum(x)))

immune_nes_norm<-rbind(l1_cell_score_norm,B_sub_score_norm,Mac_sub_score_norm,DC_sub_score_norm,T_sub1_score_norm)
binary_feature_immune_self<-rbind(B_sub_binary,Mac_sub_binary,DC_sub_binary,T_sub_binary)


#nes score combine
sam_name<-colnames(immune_nes)
sapply(sam_name,function(x) gsub("\\.","-",x))->sam_name_new
colnames(immune_nes)<-sam_name_new

nes_score_all<-rbind(immune_nes[,colnames(pahtway_nes_all)],pahtway_nes_all)


ens_ls<-list.files("EnrichmentScore_result/")
ens_feature<-c()
for(fl in ens_ls){
    load(paste("EnrichmentScore_result/",fl,sep=""))
    ens_feature<-rbind(ens_feature,t(pahtway_nes_liter))
}
pahtway_nes_all<-t(ens_feature)
for(key in names(keyword_path_list_all)){
    path<-intersect(keyword_path_list_all[[key]],row.names(pahtway_nes_all))
    if(length(path)>2){
     # path_ensembl_name<-c(path_ensembl_name,key)
      tmp<-pahtway_nes_all[path,]
      row.names(tmp)<-paste(key,row.names(tmp),sep="|")
  
    }#else(print(key))
}




pathway_nes_binary<-c()
path_ensembl_name<-c()
keyword_path_list_all<-c(keyword_path_list,literature_geneSet_ense)
#keyword_path_list_all<-c(literature_geneSet_ense)
keyword_path_list_new1<-list()
for(key in names(keyword_path_list_new)){
    if(!(key%in%c("chemokines","cytokines","epithelial","mesenchymal","MTOR","vessels")))
    keyword_path_list_new1[[key]]<-keyword_path_list_new[[key]]
}
#length(keyword_path_list_new1)
keyword_path_list_all<-keyword_path_list_new1
#keyword_path_list_all<-literature_geneSet_ense#,keyword_path_list$Immune)
#keyword_path_list_all$Immune<-keyword_path_list$Immune
#pahtway_nes_all<-pahtway_nes_liter
for(key in names(keyword_path_list_all)){
    path<-intersect(keyword_path_list_all[[key]],row.names(pahtway_nes_all))
    if(length(path)>2){
     # path_ensembl_name<-c(path_ensembl_name,key)
      tmp<-apply(pahtway_nes_all[path,],2,function(x) {x[which.max(x)]=1;x[which(x!=1)]=0;return(x)})
      row.names(tmp)<-paste(key,row.names(tmp),sep="|")
      apply(tmp,1,function(x) sum(x))->tmp_result_static
       pathway_nes_binary<-rbind(pathway_nes_binary,tmp)
      # print(key)
      # print(max(tmp_result_static))
      # if(max(tmp_result_static)<ncol(pahtway_nes_all)*0.75){
      #   pathway_nes_binary<-rbind(pathway_nes_binary,tmp)
      # } 
    }#else(print(key))
}
colnames(pathway_nes_binary)<-colnames(pahtway_nes_all)
pathway_nes_binary_filter<-pathway_nes_binary


#immune cell feature
#re_ls<-c()
run_ls<-setdiff(ca_ls,unique(ecotyper_re$Histology))
cell_state_result_all<-c()
cell_state_binary_all<-c()
for(ca in run_ls[-2]){
    re_ls<-c()
    cell_state_fl<-list.files(paste("Ecotype_run_result/",ca,"_ecotype_output/TCGA-",ca,"_exp_tpm",sep=""))
    for(cell in cell_state_fl){
        fl_name<-paste("Ecotype_run_result/",ca,"_ecotype_output/TCGA-",ca,"_exp_tpm/",cell,"/state_abundances.txt",sep="")
        re_ls<-c(re_ls,fl_name)
    }
    re_ls<-re_ls[-grep("Ecotypes",re_ls)]
    #print(length(re_ls))
    cell_state_result<-c()
    cell_state_binary<-c()
    for(fl in re_ls){
            ct<-t(read.table(paste("",fl,sep=""),header=T,sep="\t",row.names = 1))
            cell_name<-as.vector(unlist(strsplit(fl,"\\/")))[7]
            if(cell_name=="CD4.T.cells"){
                print(ca)
                print(dim(ct))
            }
            colnames(ct)<-paste(cell_name,colnames(ct),sep="_")
            cell_state_result<-rbind(cell_state_result,t(ct))
            tmp<-apply(ct,1,function(x) {x[which.max(x)]=1;x[which(x!=1)]=0;return(x)})
            cell_state_binary<-rbind(cell_state_binary,tmp)
    }    
    cell_state_result_all<-rbind(cell_state_result_all,t(cell_state_result))
    cell_state_binary_all<-rbind(cell_state_binary_all,t(cell_state_binary))
}

ecotyper_re<-read.table("TCGA_cell_state_assignments.tsv",header=T,sep="\t")
ecotyper_re<-t(immune_nes)
cell_state_binary_left<-c()
cell_name_all<-unique(sapply(colnames(ecotyper_re),function(x) {as.vector(unlist(strsplit(x,"_"))[1])}))
for(fl in cell_name_all){
    cell_name<-fl
    ct<-ecotyper_re[,grep(cell_name,colnames(ecotyper_re))[-1]]#read.table(paste("E:/project/remodeling tumor microenvironment/NMF_test/ecotyper_output_skcm/Carcinoma_Cell_States/",fl,sep=""),header=T,sep="\t",row.names = 1)
    tmp<-apply(ct,1,function(x) {x[which.max(x)]=1;x[which(x!=1)]=0;return(x)})
    cell_state_binary_left<-rbind(cell_state_binary_left,tmp)
}
colnames(cell_state_binary_left)<-ecotyper_re$TCGA_ID
cell_comm<-intersect(row.names(cell_state_binary_left),colnames(cell_state_binary_all))
immune_cell_state_binary<-rbind(cell_state_binary_all[,cell_comm],t(cell_state_binary_left[cell_comm,]))
cell_state_result_all_final<-rbind(cell_state_result_all[,cell_comm],ecotyper_re[,cell_comm])

sam_name<-row.names(immune_cell_state_binary)
sapply(sam_name,function(x) gsub("\\.","-",x))->sam_name_new
row.names(immune_cell_state_binary)<-sam_name_new
com_sam<-intersect(colnames(pathway_nes_binary_filter),row.names(immune_cell_state_binary))
binary_feature<-rbind(pathway_nes_binary_filter[,com_sam],t(immune_cell_state_binary[com_sam,]))

#免疫细胞状态特征提取
binary_feature_immune<-binary_feature[colnames(immune_cell_state_binary),]
#keyword_path_list,literature_geneSet_ense#
#通路特征提取
binary_feature_pathway<-binary_feature[row.names(pathway_nes_binary_filter),]
rename_ct<-row.names(binary_feature_pathway)
sapply(rename_ct,function(x) as.vector(unlist(strsplit(x,"\\|")))[2])->pathway_feature
row.names(binary_feature_pathway)<-pathway_feature
#肿瘤免疫反应过程
binary_feature_process<-binary_feature_pathway[intersect(row.names(binary_feature_pathway),as.vector(unlist(literature_geneSet_ense))),]
#TME相关生物学通路
binary_feature_biology<-binary_feature_pathway[intersect(row.names(binary_feature_pathway),as.vector(unlist(keyword_path_list))),]

clust_function<-function(binary_feature,fig_name){
  jacard_mt_fun(binary_feature)->Jaccard_mt_update
  hc = hclust(as.dist(1-Jaccard_mt_update), method = "average")
  n_clust = choose_clusters(Jaccard_mt_update, "jaccard", range = 2:(nrow(Jaccard_mt_update) - 1))
  #print(n_clust)
  #n_clust = 5
  clust =  hclusCut(Jaccard_mt_update, n_clust)$cluster
  new_ft<-c()
  n_index<-names(which(table(clust)>2))
  for(i in n_index){
    new_ft<-c(new_ft,(names(which(clust==i))))
  }
  #Jaccard_mt_update<-Jaccard_mt_update[row.names(pathway_nes_binary_filter),row.names(pathway_nes_binary_filter)]
  Jaccard_mt_update<-Jaccard_mt_update[new_ft,new_ft]
  hc = hclust(as.dist(1-Jaccard_mt_update), method = "average")
  n_clust = choose_clusters(Jaccard_mt_update, fig_name, range = 2:(nrow(Jaccard_mt_update) - 1))
  clust =  hclusCut(Jaccard_mt_update, n_clust)$cluster
  return(list(clust_detail=clust,n_clust=n_clust))
}

sam_name<-colnames(binary_feature_immune_self)
sapply(sam_name,function(x) gsub("\\.","-",x))->sam_name_new
colnames(binary_feature_immune_self)<-as.vector(unlist(sam_name_new))

final_sample<-intersect(colnames(binary_feature_immune_self),colnames(binary_feature_process))

#feature all 
sam_name<-colnames(cell_state_binary_new)
sapply(sam_name,function(x) gsub("\\.","-",x))->sam_name_new
colnames(cell_state_binary_new)<-as.vector(unlist(sam_name_new))
com_sam<-intersect(colnames(pathway_nes_binary_filter),colnames(cell_state_binary_new))
binary_feature<-rbind(pathway_nes_binary_filter[,com_sam],cell_state_binary_new[,com_sam])
clust_function(binary_feature,"All_feature")->All_feature_clust
#sample clust assign
process_clust_item<-All_feature_clust$clust
process_clust_abun<-c()
for(i in seq(1,max(process_clust_item))){
  ft<-names(which(process_clust_item==i))
  ft_path_index<-grep("\\|",ft)
  ft_cell<-ft[-ft_path_index]
  ft_path<-as.vector(unlist(sapply(ft[ft_path_index],function(x) as.vector(unlist(strsplit(x,"\\|")))[2])))
  tmp<-apply(t(nes_score_all)[,c(ft_cell,ft_path)],1,function(x) {mean(x)})
  process_clust_abun<-rbind(process_clust_abun,tmp)
}
row.names(process_clust_abun)<-paste("All_feature_clust_",seq(1,max(process_clust_item)),sep="")
process_assign<-apply(process_clust_abun,2,function(x) names(which.max(x)))



cell_state_binary_new<-cell_state_binary_left
clust_function(cell_state_binary_new,"immune")->immune_clust
#clust_function(binary_feature_immune_self[,c],"immune_self")->immune_self_clust
clust_function(binary_feature_process,"process")->process_clust
clust_function(binary_feature_biology,"biology")->biology_clust

self_immune_pathway_combine_binary<-rbind(binary_feature_immune_self[,final_sample],binary_feature_process[,final_sample],binary_feature_biology[,final_sample])
clust_function(self_immune_pathway_combine_binary,"self_combine")->self_combine_clust
#TME_clust_assignment
immune_fra_normalization<-c()
cell_state_binary_new<-c()
for(fl in cell_state_fl){
    cell_name<-as.vector(unlist(strsplit(fl,"_")))[1]
    ct<-cell_state_result_all_final[,grep(cell_name,colnames(cell_state_result_all_final))]#read.table(paste("E:/project/remodeling tumor microenvironment/NMF_test/ecotyper_output_skcm/Carcinoma_Cell_States/",fl,sep=""),header=T,sep="\t",row.names = 1)
    tmp<-apply(ct,1,function(x) {x/sum(x)})
    immune_fra_normalization<-rbind(immune_fra_normalization,tmp)
    tmp<-apply(ct,1,function(x) {x[which.max(x)]=1;x[which(x!=1)]=0;return(x)})
    cell_state_binary_new<-rbind(cell_state_binary_new,tmp)
}
immune_fra_normalization<-t(immune_fra_normalization)
#colnames(cell_state_binary_left)<-ecotyper_re$TCGA_ID
immune_clust_item<-immune_clust$clust_detail
immune_clust_abun<-c()
for(i in seq(1,max(immune_clust_item))){
  ft<-names(which(immune_clust_item==i))
  tmp<-apply(immune_fra_normalization[,ft],1,function(x) {mean(x)})
  immune_clust_abun<-rbind(immune_clust_abun,tmp)
}
row.names(immune_clust_abun)<-paste("Immune_clust_",seq(1,max(immune_clust_item)),sep="")
immune_assign<-apply(immune_clust_abun,2,function(x) names(which.max(x)))


#process_assign
pathway_nes_score<-t(pahtway_nes_all)
biology_clust_item<-biology_clust$clust_detail
biology_clust_abun<-c()
for(i in seq(1,max(biology_clust_item))){
  ft<-names(which(biology_clust_item==i))
  tmp<-apply(pathway_nes_score[,ft],1,function(x) {median(x)})
  biology_clust_abun<-rbind(biology_clust_abun,tmp)
}
row.names(biology_clust_abun)<-paste("Biology_clust_",seq(1,max(biology_clust_item)),sep="")
Biology_assign<-apply(biology_clust_abun,2,function(x) names(which.max(x)))

#biology assign
process_clust_item<-process_clust$clust_detail
process_clust_abun<-c()
for(i in seq(1,max(process_clust_item))){
  ft<-names(which(process_clust_item==i))
  tmp<-apply(pathway_nes_score[,ft],1,function(x) {mean(x)})
  process_clust_abun<-rbind(process_clust_abun,tmp)
}
row.names(process_clust_abun)<-paste("Process_clust_",seq(1,max(process_clust_item)),sep="")
process_assign<-apply(process_clust_abun,2,function(x) names(which.max(x)))

sam_name<-names(immune_assign)
sapply(sam_name,function(x) gsub("\\.","-",x))->sam_name_new
names(immune_assign)<-as.vector(unlist(sam_name_new))
final_sample<-intersect(names(immune_assign),names(process_assign))
sapply(final_sample,function(x) {paste(immune_assign[x],process_assign[x],Biology_assign[x],sep="|")})->TME_clust_final
names(which(table(TME_clust_final)>100))->tag_item

library(jaccard)
#jaccard(binary_feature, y, center = FALSE, px = NULL, py = NULL)
#calculation process
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

jacard_mt_fun<-function(binary_feature){
  Jaccard_mt<-c()
  for(i in seq(1,nrow(binary_feature))){
    apply(binary_feature,1,function(x) jaccard::jaccard(as.numeric(x),as.numeric(binary_feature[i,])))->re
    Jaccard_mt<-rbind(Jaccard_mt,re)
    #Jaccard_mt_update[i,which(re>0.01)]=0
  }
  row.names(Jaccard_mt)<-row.names(binary_feature)
  Jaccard_mt_update<-Jaccard_mt
  for(i in seq(1,nrow(binary_feature))){
    apply(binary_feature,1,function(x) phyper_fun(x,binary_feature[i,]))->re
    Jaccard_mt_update[i,which(re>0.01)]=0
    #Jaccard_mt_update[which(re>0.01,i)]=0
  }
  return(Jaccard_mt_update)
}






hclusCut <- function(x, k, ...) list(cluster = cutree(hclust(as.dist(1-x), method = "average", ...), k=k))
choose_clusters <- function(data, name, range = 2:10)
{
        #gap = clusGap(jaccard, hclusCut, 50, B = 100)
        silh <- data.frame(K = range, Silhouette = sapply(range, function(k){
                sil <<- silhouette(hclusCut(data, k)$cluster, as.dist(1-data))
                tmp <<- summary(sil)
                tmp$avg.width
        }))

        g2 <- ggplot(silh, aes(x = K, y = Silhouette)) +
                geom_point() +
                #ylab("Sparseness (coef)") +
                geom_line() +
                geom_vline(xintercept = silh[which.max(silh$Silhouette), 1], lty = 2, colour = "red") +
                theme_bw() +
                theme(panel.grid = element_blank()) +
                theme(aspect.ratio = 1)
                #facet_wrap(CellType, ncol = 4) +
}
#Jaccard_mt_update<-Jaccard_mt_update[row.names(pathway_nes_binary_filter),row.names(pathway_nes_binary_filter)]
#Jaccard_mt_update<-Jaccard_mt_update[new_ft,new_ft]
#print(n_clust)
clust =  hclusCut(Jaccard_mt_update, n_clust)$cluster
sil = silhouette(clust, as.dist(1-Jaccard_mt_update))
avg_silhouette = summary(sil)

top_ann = as.data.frame(t(sapply(rownames(Jaccard_mt_update), function(x) {
	s= strsplit(x, "_")[[1]]
	c(paste0(s[-length(s)],collapse = "_"), s[length(s)])
})))

colnames(top_ann) = c("CellType","State")
top_ann$InitialEcotype = as.factor(sprintf("IE%02d", clust))

min_states=3
tb = table(top_ann$InitialEcotype)
tb = tb[tb >= min_states]

top_ann = top_ann[top_ann$InitialEcotype %in% names(tb),]

top_ann$ID = rownames(top_ann)
nm = unique(top_ann$InitialEcotype)
mapping = sprintf("E%d", 1:length(nm))
names(mapping) = nm

top_ann$Ecotype = mapping[as.character(top_ann$InitialEcotype)]
top_ann$Ecotype = ecotype_to_factor(top_ann$Ecotype)
top_ann = top_ann[order(top_ann$Ecotype),]
jaccard<-Jaccard_mt_update
jaccard = jaccard[match(top_ann$ID, rownames(jaccard)), match(top_ann$ID, rownames(jaccard))]

sil <- silhouette(as.numeric(as.character(gsub("E", "", as.character(top_ann$Ecotype)))), as.dist(1-jaccard))
avg_silhouette <<- summary(sil)
write.table(avg_silhouette$avg.width, file.path(output_dir, "silhouette.txt"), sep = "\t", row.names = F)

top_ann$"Cell type" = top_ann$CellType

library(jaccard)
library(ComplexHeatmap)
library(RColorBrewer)
library("viridis")
h <- heatmap_simple(jaccard, name = "ht1", top_annotation = top_ann, top_columns = c("Ecotype", "Cell type", "State"), 
 legend_name = "Jaccard index", width = unit(2, "in"), height = unit(2, "in"),
color_palette = c("white", viridis(4)), 
color_range = c(0, 0.1, 0.2, 0.3))
draw(h, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", 
	adjust_annotation_extension = T, merge_legends = T)	

ord = top_ann$Ecotype
dup = (which(!duplicated(ord)) - 1)
fract = dup / nrow(top_ann)
width =  c(fract[-1], 1) - fract
decorate_heatmap_body("ht1", {
    grid.rect(unit(fract, "native"), unit(1-fract, "native"), unit(width, "native"), unit(width, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", lty = 1, lwd = 2))
})


png(file.path(output_dir, "jaccard_matrix.png"), width = 8, height = 7, res = 300, units = "in", family = "Helvetica")
draw(h, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", 
	adjust_annotation_extension = T, merge_legends = T)	
decorate_heatmap_body("ht1", {
    grid.rect(unit(fract, "native"), unit(1-fract, "native"), unit(width, "native"), unit(width, "native"), hjust = 0, vjust = 1, gp = gpar(col = "white", fill = NA, lty = 1, lwd = 2))
})
tmp = dev.off()


#sample to cluster
#pahtway_nes_all
path_raw<-c()
path_new<-c()
for(key in names(keyword_path_list_all)){
    path<-intersect(keyword_path_list_all[[key]],row.names(pahtway_nes_all))
    if(length(path)>2){
      paste(key,path,sep="|")->name_new
      path_raw<-c(path_raw,path)
      path_new<-c(path_new,name_new)
    }#else(print(key))
}
pahtway_nes_all_update<-pahtway_nes_all[path_raw,]
row.names(pahtway_nes_all_update)<-path_new
pahtway_nes_all_update<-pahtway_nes_all_update[intersect(names(clust),row.names(pahtway_nes_all_update)),colnames(binary_feature)]
pathway_nes_binary_update<-c()
for(key in names(keyword_path_list_all)){
    path<-intersect(keyword_path_list_all[[key]],row.names(pahtway_nes_all[path_raw,]))
    if(length(path)>2){
     # path_ensembl_name<-c(path_ensembl_name,key)
      tmp<-apply(pahtway_nes_all[path,],2,function(x) {x[which.max(x)]=1;x[which(x!=1)]=0;return(x)})
      row.names(tmp)<-paste(key,row.names(tmp),sep="|")
       pathway_nes_binary_update<-rbind(pathway_nes_binary_update,tmp)
      # print(key)
      # print(max(tmp_result_static))
      # if(max(tmp_result_static)<ncol(pahtway_nes_all)*0.75){
      #   pathway_nes_binary<-rbind(pathway_nes_binary,tmp)
      # } 
    }#else(print(key))
}
colnames(pathway_nes_binary_update)<-colnames(pahtway_nes_all)
pathway_nes_binary_update<-pathway_nes_binary_update[intersect(names(clust),row.names(pathway_nes_binary_update)),colnames(binary_feature)]


comm_cell<-intersect(colnames(cell_state_result_all),colnames(ecotyper_re))
immune_left<-ecotyper_re[,comm_cell]
row.names(immune_left)<-ecotyper_re$TCGA_ID
immune_result<-rbind(cell_state_result_all[,comm_cell],immune_left)
sam_name<-row.names(immune_result)
sapply(sam_name,function(x) gsub("\\.","-",x))->sam_name_new
row.names(immune_result)<-sam_name_new
immune_result_new<-immune_result[,intersect(colnames(immune_result),names(clust))]
cell_state_binary_update<-c()
for(fl in cell_state_fl){
    cell_name<-as.vector(unlist(strsplit(fl,"_")))[1]
    ct<-immune_result[,grep(cell_name,colnames(immune_result))]#read.table(paste("E:/project/remodeling tumor microenvironment/NMF_test/ecotyper_output_skcm/Carcinoma_Cell_States/",fl,sep=""),header=T,sep="\t",row.names = 1)
    tmp<-apply(ct,1,function(x) {x[which.max(x)]=1;x[which(x!=1)]=0;return(x)})
    cell_state_binary_update<-rbind(cell_state_binary_update,tmp)
}
 
sample_feature_update<-rbind(cell_state_binary_update[,colnames(pathway_nes_binary_update)],pathway_nes_binary_update)
sample_feature_update<-t(sample_feature_update[names(clust),])


hamming <- function(a, b) {
  return(sum(a != b))
}
clust_mt<-matrix(rep(0,16*138),ncol=138,nrow=16)
row.names(clust_mt)<-paste("TME",seq(1,16),sep="")
colnames(clust_mt)<-names(clust)
for(i in seq(1,16)){
  ft<-names(which(clust==i))
  clust_mt[i,ft]=1
}

hamming_matrix<-c()
for(i in seq(1,nrow(sample_feature_update))){
  apply(clust_mt,1,function(x) hamming(as.numeric(sample_feature_update[i,]),as.numeric(x)))->hamming_re
  hamming_matrix<-rbind(hamming_matrix,hamming_re)
}
row.names(hamming_matrix)<-row.names(sample_feature_update)
apply(hamming_matrix,1,function(x) names(which.min(x)))->sample_hamming_clust
