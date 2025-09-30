#load("TCGA_tsne.RData")

CD45_ME<-list()
Stromal_ME<-list()
Purity_ME<-list()
sig_gene_cancer_ls<-cancer_ls[-c(14,31)]
for(cancer in sig_gene_cancer_ls){
    print(cancer)
    Cor<-readRDS(paste("WGCNA_result/sig_module_gene/",cancer,"_Cor.RDS",sep=""))
    Pva<-readRDS(paste("WGCNA_result/sig_module_gene/",cancer,"_pvalue.RDS",sep=""))
    SigG<-readRDS(paste("WGCNA_result/sig_module_gene/",cancer,"_sig_gene.RDS",sep=""))
    cd45Me<-names(which(Cor[,1]>0.6))
    stromalMe<-names(which(Cor[,2]>0.6))
    purityMe<-names(which(Cor[,3]>0.5))
    
#     cd45Me<-names(which.max(Cor[,1]))
#     cafMe<-names(which.max(Cor[,2]))
#    # mscMe<-names(which.max(Cor[,3]))
#     endoMe<-names(which.max(Cor[,3]))
#     periMe<-names(which.max(Cor[,4]))
#     purityMe<-names(which.max(Cor[,5]))
    for(me in cd45Me){
        if(Cor[cd45Me,1]>0.6){
            CD45_ME[[paste(cancer,"_CD45_",me,sep="")]]=as.vector(unlist(SigG[which(SigG$module==me),1]))
        }    
    }
    for(me in stromalMe){
        if(Cor[stromalMe,2]>0.6){
            Stromal_ME[[paste(cancer,"_stromal_",me,sep="")]]=as.vector(unlist(SigG[which(SigG$module==me),1]))
        }
        
    }
    for(me in purityMe){
        if(Cor[purityMe,3]>0.5){
            Purity_ME[[paste(cancer,"_purity_",me,sep="")]]=as.vector(unlist(SigG[which(SigG$module==me),1]))
        }
        
    }
    #sig_gene<-paste("WGCNA_result/sig_module_gene/",cancer,"_sig_gene.RDS",sep=""))
}

library(clusterProfiler)
library(org.Hs.eg.db)
CD45_enrich_list<-list()
for(i in seq(1,length(CD45_ME))){
    gtmp<-as.vector(unlist(sapply(CD45_ME[[i]],function(x) as.vector(unlist(strsplit(x,"\\|")))[1])))
    trans_file<-read.table("entrez_symbol_trans",header=T,sep="\t")
    trans_file_unique<-na.omit(trans_file[!duplicated(trans_file$SYMBOL),])
    row.names(trans_file_unique)<-trans_file_unique$SYMBOL
    gg<-as.vector(unlist(na.omit(trans_file_unique[gtmp,"ENTREZID"])))
    kk <- enrichKEGG(gene         = gg,
                    organism     = 'hsa',
                    pvalueCutoff = 0.05)
    CD45_enrich_list[[names(CD45_ME)[i]]]<-data.frame(kk)
}
CD45_enrichment_df<-data.frame()
pa_list<-c()
for(me in names(CD45_enrich_list)){
    tmp<-CD45_enrich_list[[me]]
    pa_list<-rbind(pa_list,cbind(tmp[,c("Description","GeneRatio","qvalue")],rep(me,nrow(tmp))))
}
colnames(pa_list)[4]="WGCNA_model"
ratio1<-as.numeric(as.vector(unlist(sapply(pa_list$GeneRatio,function(x) {as.vector(unlist(strsplit(x,"\\/")))[1]}))))
ratio2<-as.numeric(as.vector(unlist(sapply(pa_list$GeneRatio,function(x) {as.vector(unlist(strsplit(x,"\\/")))[2]}))))
pa_list$Ratio<-ratio1/ratio2
ggplot(pa_list,aes(x=Description,y=WGCNA_model))+geom_point(aes(size=Ratio,color=-1*log10(qvalue)))
pa_list$Description<-factor(pa_list$Description,levels=names(sort(table(pa_list$Description),decreasing=TRUE)),ordered=TRUE)
pa_list$WGCNA_model<-factor(pa_list$WGCNA_model,levels=names(sort(table(pa_list$WGCNA_model),decreasing=TRUE)),ordered=TRUE)
pdf("test.pdf",15,8)
ggplot(pa_list,aes(x=Description,y=WGCNA_model))+geom_point(aes(size=Ratio,color=-1*log10(qvalue)))+theme(axis.text.x = element_text(angle=90, hjust=0.7, vjust=0.7))
dev.off()

trans_file<-read.table("entrez_symbol_trans",header=T,sep="\t")
row.names(trans_file)<-trans_file[,1]
CD45_top5_pathway<-names(sort(table(pa_list$Description),decreasing=TRUE))[1:10]
top_path_gene<-c()
CD45_top_path_cytoscape_df<-c()
for(me in names(CD45_enrich_list)){
    tmp<-CD45_enrich_list[[me]]
    tmp%>%dplyr::filter(Description%in%CD45_top5_pathway)%>%
        dplyr::select(Description,geneID)->tmp_select
    for(i in seq(1,nrow(tmp_select))){
        item<-tmp_select[i,]
        g<-as.vector(unlist(sapply(as.vector(unlist(item[,2])),function(x) {strsplit(x,"\\/")})))
        CD45_top_path_cytoscape_df<-rbind(CD45_top_path_cytoscape_df,cbind(rep(item[,1],length(g)),trans_file[g,"SYMBOL"]))
        top_path_gene<-c(top_path_gene,g)
    }
}
CD45_top_path_cytoscape_df_unique<-CD45_top_path_cytoscape_df[!duplicated(interaction(CD45_top_path_cytoscape_df[,1],CD45_top_path_cytoscape_df[,2])),]
write.table(CD45_top_path_cytoscape_df_unique,file="ResultForPaper/CD45_KEGG_path",sep="\t",quote=F)

top_path_gene<-unique(top_path_gene)
top_path_gene_symbol<-trans_file[top_path_gene,]
all_gene<-unique(as.vector(unlist(CD45_ME)))
CD45_fre<-rep(0,length(all_gene))
names(CD45_fre)<-all_gene
for(me in names(CD45_ME)){
    CD45_fre[CD45_ME[[me]]]=CD45_fre[CD45_ME[[me]]]+1
}
# top_path_gene_fre<-CD45_fre[top_path_gene_symbol[,2]]
# geme_model_fre_df<-data.frame(Gene=names(top_path_gene_fre),Count=as.numeric(top_path_gene_fre))
# ggplot(gene_model_fre_df,aes(x=Count,y=Gene))

Stromal_ME<-c(CAF_ME,Endo_ME,MSC_ME,pericyte_ME)
Stromal_enrich_list<-list()
for(i in seq(1,length(Stromal_ME))){
    gtmp<-as.vector(unlist(sapply(Stromal_ME[[i]],function(x) as.vector(unlist(strsplit(x,"\\|")))[1])))
    trans_file<-read.table("entrez_symbol_trans",header=T,sep="\t")
    trans_file_unique<-na.omit(trans_file[!duplicated(trans_file$SYMBOL),])
    row.names(trans_file_unique)<-trans_file_unique$SYMBOL
    gg<-as.vector(unlist(na.omit(trans_file_unique[gtmp,"ENTREZID"])))
    kk <- enrichKEGG(gene         = gg,
                    organism     = 'hsa',
                    pvalueCutoff = 0.05)
    Stromal_enrich_list[[names(Stromal_ME)[i]]]<-data.frame(kk)
}
#Stromal_enrich_df<-data.frame()
pa_list<-c()
for(me in names(Stromal_enrich_list)){
    tmp<-Stromal_enrich_list[[me]]
    pa_list<-rbind(pa_list,cbind(tmp[,c("Description","GeneRatio","qvalue")],rep(me,nrow(tmp))))
}
colnames(pa_list)[4]="WGCNA_model"
ratio1<-as.numeric(as.vector(unlist(sapply(pa_list$GeneRatio,function(x) {as.vector(unlist(strsplit(x,"\\/")))[1]}))))
ratio2<-as.numeric(as.vector(unlist(sapply(pa_list$GeneRatio,function(x) {as.vector(unlist(strsplit(x,"\\/")))[2]}))))
pa_list$Ratio<-ratio1/ratio2
#ggplot(pa_list,aes(x=Description,y=WGCNA_model))+geom_point(aes(size=Ratio,color=-1*log10(qvalue)))
pa_list$Description<-factor(pa_list$Description,levels=names(sort(table(pa_list$Description),decreasing=TRUE)),ordered=TRUE)
pa_list$WGCNA_model<-factor(pa_list$WGCNA_model,levels=names(sort(table(pa_list$WGCNA_model),decreasing=TRUE)),ordered=TRUE)
pdf("ResultForPaper/Stromal_enrich.pdf",15,10)
ggplot(pa_list,aes(x=Description,y=WGCNA_model))+geom_point(aes(size=Ratio,color=-1*log10(qvalue)))+theme(axis.text.x = element_text(angle=90, hjust=0.7, vjust=0.7))
dev.off()
#save.image(file="ES_analysis0105.rda",ascii=FALSE,compress=TRUE)
trans_file<-read.table("entrez_symbol_trans",header=T,sep="\t")
row.names(trans_file)<-trans_file[,1]
Stromal_top5_pathway<-names(sort(table(pa_list$Description),decreasing=TRUE))[1:10]
top_path_gene_stromal<-c()
Stromal_top_path_cytoscape_df<-c()
for(me in names(Stromal_enrich_list)){
    tmp<-Stromal_enrich_list[[me]]
    # tmp%>%dplyr::filter(Description%in%Stromal_top5_pathway)%>%
    #     dplyr::select(geneID)->tmp_select
    # tmp_gene<-as.vector(unlist(sapply(as.vector(unlist(tmp_select)),function(x) {strsplit(x,"\\/")})))
    # top_path_gene_stromal<-c(top_path_gene_stromal,tmp_gene)
    tmp%>%dplyr::filter(Description%in%Stromal_top5_pathway)%>%
        dplyr::select(Description,geneID)->tmp_select
    for(i in seq(1,nrow(tmp_select))){
        item<-tmp_select[i,]
        g<-as.vector(unlist(sapply(as.vector(unlist(item[,2])),function(x) {strsplit(x,"\\/")})))
        Stromal_top_path_cytoscape_df<-rbind(Stromal_top_path_cytoscape_df,cbind(rep(item[,1],length(g)),trans_file[g,"SYMBOL"]))
        top_path_gene_stromal<-c(top_path_gene_stromal,g)
    }
  #  pa_list<-rbind(pa_list,cbind(tmp[,c("Description","GeneRatio","qvalue")],rep(me,nrow(tmp))))
}
Stromal_top_path_cytoscape_df_unique<-Stromal_top_path_cytoscape_df[!duplicated(interaction(Stromal_top_path_cytoscape_df[,1],Stromal_top_path_cytoscape_df[,2])),]
write.table(Stromal_top_path_cytoscape_df_unique,file="ResultForPaper/Stromal_KEGG_path",sep="\t",quote=F)



top_path_gene_stromal<-unique(top_path_gene_stromal)
top_path_gene_stromal_symbol<-trans_file[top_path_gene_stromal,]


purity_enrich_list<-list()
for(i in seq(1,length(Purity_ME))){
    gtmp<-as.vector(unlist(sapply(Purity_ME[[i]],function(x) as.vector(unlist(strsplit(x,"\\|")))[1])))
    trans_file<-read.table("entrez_symbol_trans",header=T,sep="\t")
    trans_file_unique<-na.omit(trans_file[!duplicated(trans_file$SYMBOL),])
    row.names(trans_file_unique)<-trans_file_unique$SYMBOL
    gg<-as.vector(unlist(na.omit(trans_file_unique[gtmp,"ENTREZID"])))
    kk <- enrichKEGG(gene         = gg,
                    organism     = 'hsa',
                    pvalueCutoff = 0.05)
    purity_enrich_list[[names(Purity_ME)[i]]]<-data.frame(kk)
}
#Stromal_enrich_df<-data.frame()
pa_list<-c()
for(me in names(purity_enrich_list)){
    tmp<-purity_enrich_list[[me]]
    pa_list<-rbind(pa_list,cbind(tmp[,c("Description","GeneRatio","qvalue")],rep(me,nrow(tmp))))
}
colnames(pa_list)[4]="WGCNA_model"
ratio1<-as.numeric(as.vector(unlist(sapply(pa_list$GeneRatio,function(x) {as.vector(unlist(strsplit(x,"\\/")))[1]}))))
ratio2<-as.numeric(as.vector(unlist(sapply(pa_list$GeneRatio,function(x) {as.vector(unlist(strsplit(x,"\\/")))[2]}))))
pa_list$Ratio<-ratio1/ratio2
#ggplot(pa_list,aes(x=Description,y=WGCNA_model))+geom_point(aes(size=Ratio,color=-1*log10(qvalue)))
pa_list$Description<-factor(pa_list$Description,levels=names(sort(table(pa_list$Description),decreasing=TRUE)),ordered=TRUE)
pa_list$WGCNA_model<-factor(pa_list$WGCNA_model,levels=names(sort(table(pa_list$WGCNA_model),decreasing=TRUE)),ordered=TRUE)
pdf("ResultForPaper/Purity_enrich.pdf",20,6)
ggplot(pa_list,aes(x=Description,y=WGCNA_model))+geom_point(aes(size=Ratio,color=-1*log10(qvalue)))+theme(axis.text.x = element_text(angle=90, hjust=0.7, vjust=0.7))
dev.off()

trans_file<-read.table("entrez_symbol_trans",header=T,sep="\t")
row.names(trans_file)<-trans_file[,1]
purity_top5_pathway<-names(sort(table(pa_list$Description),decreasing=TRUE))[1:10]
top_path_gene_purity<-c()
Purity_top_path_cytoscape_df<-c()
for(me in names(purity_enrich_list)){
    tmp<-purity_enrich_list[[me]]
    # tmp%>%dplyr::filter(Description%in%Stromal_top5_pathway)%>%
    #     dplyr::select(geneID)->tmp_select
    # tmp_gene<-as.vector(unlist(sapply(as.vector(unlist(tmp_select)),function(x) {strsplit(x,"\\/")})))
    # top_path_gene_stromal<-c(top_path_gene_stromal,tmp_gene)
    tmp%>%dplyr::filter(Description%in%purity_top5_pathway)%>%
        dplyr::select(Description,geneID)->tmp_select
    for(i in seq(1,nrow(tmp_select))){
        item<-tmp_select[i,]
        g<-as.vector(unlist(sapply(as.vector(unlist(item[,2])),function(x) {strsplit(x,"\\/")})))
        Purity_top_path_cytoscape_df<-rbind(Purity_top_path_cytoscape_df,cbind(rep(item[,1],length(g)),trans_file[g,"SYMBOL"]))
        top_path_gene_purity<-c(top_path_gene_purity,g)
    }
  #  pa_list<-rbind(pa_list,cbind(tmp[,c("Description","GeneRatio","qvalue")],rep(me,nrow(tmp))))
}
Purity_top_path_cytoscape_df_unique<-Purity_top_path_cytoscape_df[!duplicated(interaction(Purity_top_path_cytoscape_df[,1],Purity_top_path_cytoscape_df[,2])),]
write.table(Purity_top_path_cytoscape_df_unique,file="ResultForPaper/Purity_KEGG_path",sep="\t",quote=F)

TME_pathway_ls<-rbind(cbind(rep("CD45",length(CD45_top5_pathway)),CD45_top5_pathway),
                      cbind(rep("Stromal",length(Stromal_top5_pathway)),Stromal_top5_pathway),
                      cbind(rep("Tumor",length(purity_top5_pathway)),purity_top5_pathway))
colnames(TME_pathway_ls)<-c("Source","Traget")
write.table(TME_pathway_ls,file="ResultForPaper/TME_path_cytoscape",sep="\t",quote=F)


all_genes<-unique(c(Purity_top_path_cytoscape_df_unique[,2],Stromal_top_path_cytoscape_df_unique[,2],CD45_top_path_cytoscape_df_unique[,2]))
gene_fre_static<-rep(0,length(all_genes))
names(gene_fre_static)<-all_genes
for(g in all_genes){
    if(g %in% Purity_top_path_cytoscape_df_unique[,2]){
        gene_fre_static[g]<-gene_fre_static[g]+1
    }
    if(g %in% Stromal_top_path_cytoscape_df_unique[,2]){
        gene_fre_static[g]<-gene_fre_static[g]+1
    }
    if(g %in% CD45_top_path_cytoscape_df_unique[,2]){
        gene_fre_static[g]<-gene_fre_static[g]+1
    }
}
purity_new<-c()
stromal_new<-c()
cd45_new<-c()
for(g in names(which(gene_fre_static>2))){
    purity_new<-rbind(purity_new,Purity_top_path_cytoscape_df_unique[which(Purity_top_path_cytoscape_df_unique[,2]==g),])
    stromal_new<-rbind(stromal_new,Stromal_top_path_cytoscape_df_unique[which(Stromal_top_path_cytoscape_df_unique[,2]==g),])
    cd45_new<-rbind(cd45_new,CD45_top_path_cytoscape_df_unique[which(CD45_top_path_cytoscape_df_unique[,2]==g),])
}
colnames(purity_new)<-c("Source","Target")
colnames(stromal_new)<-c("Source","Target")
colnames(cd45_new)<-c("Source","Target")
write.table(purity_new,file="ResultForPaper/Purity_sig_gene")
write.table(stromal_new,file="ResultForPaper/Stromal_sig_gene")
write.table(cd45_new,file="ResultForPaper/CD45_sig_gene")


all_gene<-unique(as.vector(unlist(CD45_ME)))
CD45_fre<-rep(0,length(all_gene))
names(CD45_fre)<-all_gene
for(me in names(CD45_ME)){
    CD45_fre[CD45_ME[[me]]]=CD45_fre[CD45_ME[[me]]]+1
}
all_gene<-unique(as.vector(unlist(Stromal_ME)))
Stromal_fre<-rep(0,length(all_gene))
names(Stromal_fre)<-all_gene
for(me in names(Stromal_ME)){
    Stromal_fre[Stromal_ME[[me]]]=Stromal_fre[Stromal_ME[[me]]]+1
}
all_gene<-unique(as.vector(unlist(Purity_ME)))
Purity_fre<-rep(0,length(all_gene))
names(Purity_fre)<-all_gene
for(me in names(Purity_ME)){
    Purity_fre[Purity_ME[[me]]]=Purity_fre[Purity_ME[[me]]]+1
}
names(CD45_fre)<-as.vector(unlist(sapply(names(CD45_fre),function(x) {as.vector(unlist(strsplit(x,"\\|"))[1])})))
names(Stromal_fre)<-as.vector(unlist(sapply(names(Stromal_fre),function(x) {as.vector(unlist(strsplit(x,"\\|"))[1])})))
names(Purity_fre)<-as.vector(unlist(sapply(names(Purity_fre),function(x) {as.vector(unlist(strsplit(x,"\\|"))[1])})))

CD45_high_fre<-names(which(CD45_fre>9))
Stromal_high_fre<-names(which(Stromal_fre>5))
Purity_high_fre<-names(which(Purity_fre>1))


#gene list enrichment a
TMB_gene_ls<-c(CD45_ME,Stromal_ME,Purity_ME)

Pancancer_expression_new<-Pancancer_expression[,-c(which(anno_color=="LAML"),which(anno_color=="UCEC"))]
sample_id<-colnames(Pancancer_expression_new)
sample_id_sub<-sapply(sample_id,function(x) {paste(as.vector(unlist(strsplit(x,"\\."))[1:4]),collapse="-")})
sample_id_sub1<-sapply(sample_id_sub,functioStromal_fren(x) {substring(x, 1, str_length(x)-1)})
rm_l<-length(which(sample_id_sub1%in%names(which(table(sample_id_sub1)>1))))
tumor_sam<-intersect(sample_id_sub1[-which(sample_id_sub1%in%names(which(table(sample_id_sub1)>1)))],row.names(Tumor_purity))
Pancan_expr<-Pancancer_expression_new[,names(sample_id_sub1)[which(sample_id_sub1%in%tumor_sam)]]
#Pancan_expr
gsva(as.matrix(na.omit(Pancan_expr)),TMB_gene_ls,method="ssgsea",ssgsea.norm=TRUE)->TME_ES

# readr::write_rds(TME_ES,
#TME_ES analysis
lymph_ls<-names(CD45_ME)
stromal_ls<-names(Stromal_ME)
purity_ls<-names(Purity_ME)
lymph_es<-apply(TME_ES[lymph_ls,],2,median)
stromal_es<-apply(TME_ES[stromal_ls,],2,median)
purity_es<-apply(TME_ES[purity_ls,],2,median)
sample_TME_type<-rep(0,length(lymph_es))
names(sample_TME_type)<-names(lymph_es)
# for(i in seq(1,length(lymph_es))){
  
# }
SLP_mt<-rbind(stromal_es,lymph_es,purity_es)
SLP_max_tag<-apply(SLP_mt,2,function(x) {which.max(x)})
SL_tag<-stromal_es-lymph_es
LP_tag<-lymph_es-purity_es
SP_tag<-stromal_es-purity_es

T1=intersect(names(which(SLP_max_tag==1)),names(which(LP_tag>0)))
T2=intersect(names(which(SLP_max_tag==1)),names(which(LP_tag<0)))
T3=intersect(names(which(SLP_max_tag==2)),names(which(SP_tag>0)))
T4=intersect(names(which(SLP_max_tag==2)),names(which(SP_tag<0)))
T5=intersect(names(which(SLP_max_tag==3)),names(which(SL_tag>0)))
T6=intersect(names(which(SLP_max_tag==3)),names(which(SL_tag<0)))
sample_TME_type[T1]="G1"
sample_TME_type[T2]="G2"
sample_TME_type[T3]="G3"
sample_TME_type[T4]="G4"
sample_TME_type[T5]="G5"
sample_TME_type[T6]="G6"


Exp_tsne<-Pancan_expr
gene_id_symbol<-row.names(Exp_tsne)
gene_symbol<-as.vector(unlist(sapply(gene_id_symbol,function(x) {as.vector(unlist(strsplit(x,"\\|")))[1]})))
gene_index<-which(gene_symbol%in%c(top_path_gene_symbol[,2],top_path_gene_stromal_symbol[,2]))

library(Rtsne)
library(ggplot2)
tt<-Pancan_expr[gene_index,]
tsne_out<-Rtsne(t(na.omit(tt)),perplexity=100,step=5000)
tsne.data<-data.frame(X=tsne_out$Y[,1],Y=tsne_out$Y[,2])
#pdf("test.pdf")
ggplot(tsne.data,aes(x=X,y=Y))+geom_point()


library(fastcluster)
hc <- hclust(dist(t(tt)), "ave")
plot(hc)
plot(hc, hang = -1)
## Do the same with centroid clustering and squared Euclidean distance,
## cut the tree into ten clusters and reconstruct the upper part of the
## tree from the cluster centers.
hc <- hclust.vector(t(tt), "cen")
# squared Euclidean distances
hc$height <- hc$height^2
memb <- cutree(hc, k = 7)


#dev.off()
# names(anno_color)<-names(sample_id_sub)
# #plot(tsne_out$Y,col=as.factor(group_ls),pch=16)
# group_ls<-as.factor(anno_color[colnames(Pancan_expr)])
# label_index<-table(anno_color[colnames(Pancan_expr)])
# label<-rep(NA,length(group_ls))
# label[label_index]<-names(label_index)
# tsne.data<-data.frame(X=tsne_out$Y[,1],Y=tsne_out$Y[,2],Group=group_ls,Type=label)
# label.df <- data.frame(Group=tsne.data$Group,label=tsne.data$Type)
# label.df_2 <- tsne.data %>% 
#   group_by(Group) %>% 
#   summarize(X =median(X), Y = median(Y)) %>% 
#   left_join(label.df)

tt<-Pancan_expr[gene_index,]
a<-na.omit(scale(tt))

library(HiClimR)
#print(sys.time())
print(Sys.time())
xcor10 <- fastCor(a, upperTri = FALSE)
print(Sys.time())
library(pheatmap2)
bk <- c(seq(-1,-0.1,by=0.01),seq(0,1,by=0.01))
pheatmap2(xcor10[1:2000,1:2000],
         scale = "none",show_rownames=F,show_colnames=F,
         color = c(colorRampPalette(colors = c("blue","white"))(122),colorRampPalette(colors = c("white","red"))(70)),
         breaks=bk)
library(pheatmap)
pdf("test.pdf",18,6)
pheatmap2(xcor10,
         scale = "none",show_rownames=F,show_colnames=F,
         color = c(colorRampPalette(colors = c("blue","white"))(122),colorRampPalette(colors = c("white","red"))(70)),
         breaks=bk)->gene_exp_cor
dev.off()

CD45_ME_id<-names(CD45_ME)
Stromal_ME_id<-names(c(CAF_ME,Endo_ME,pericyte_ME))
Purity_ME_id<-names(Purity_ME)

cancer_name_CD45<-sapply(c(CD45_ME_id),function(x) {as.vector(unlist(strsplit(x,"_")))[1]})
cancer_name_stromal<-sapply(c(Stromal_ME_id),function(x) {as.vector(unlist(strsplit(x,"_")))[1]})
cancer_name_purity<-sapply(c(Purity_ME_id),function(x) {as.vector(unlist(strsplit(x,"_")))[1]})

cancer_gene_model_static<-data.frame(Type=c(rep("CD45",length(unique(cancer_name_CD45))),rep("Stromal",length(unique(cancer_name_stromal))),rep("Tumor",length(unique(cancer_name_purity)))),
                                     Count=c(as.numeric(table(cancer_name_CD45)),as.numeric(table(cancer_name_stromal)),as.numeric(table(cancer_name_purity))),
                                     Cancer=c(names(table(cancer_name_CD45)),names(table(cancer_name_stromal)),names(table(cancer_name_purity))))
pdf("test.pdf",8,6)
ggplot(cancer_gene_model_static,aes(x=Cancer,y=Count,fill=Type))+geom_bar(stat="identity")+
theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1))
dev.off()
# cancer_model_df<-data.frame(Count)


    # ridge_df<-data.frame(Fra=c(stromal_es,lymph_es,purity_es),Sample=rep("Stromal",length(stromal_es)),Tag=c(rep("Lymph",length(Purity)),rep("Leukocyte",length(Leukocyte))))
    # ggplot(ridge_df) +theme(legend.position='none')+
    # geom_density_ridges(
    #   aes(x = Fra, y=Sample ,fill = Tag),
    #   alpha = .8, color = "white"
    # ) + scale_fill_manual(values=c("#e95f5c","#40b3ff"))+theme(axis.text.x=element_text(angle=50,hjust = 0.5,vjust=0.5))+
    #     theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black",size = 0.2),plot.title = element_text(hjust = 0.5,size = 8))+
    # labs(
    #   x = "Fra",
    #   title = cl
    #  # subtitle = "Analysis unit: municipalities (n = 949)",
    #  # caption = "Marc Belzunces (@marcbeldata) | Source: Idescat"
    # )#->
#tt<-Pancan_expr[gene_index,]
library(stringr)
a<-na.omit(t(tt))
sample_id<-colnames(a)
sample_id_sub<-sapply(sample_id,function(x) {paste(as.vector(unlist(strsplit(x,"\\."))[1:4]),collapse="-")})
sample_id_sub1<-sapply(sample_id_sub,function(x) {substring(x, 1, str_length(x)-1)})
Purity<-Tumor_purity[sample_id_sub1,"purity"]
dt<-na.omit(rbind(a,Purity))
#kc<-kmeans(dt,5) 
xcor10 <- fastCor(t(dt), upperTri = FALSE)
kc<-kmeans(xcor10,3)
CD38_index<-grep("CD38",row.names(Pancan_expr))
ADCY2_index<-grep("ADCY2",row.names(Pancan_expr))
ITGA1_index<-grep("ITGA1",row.names(Pancan_expr))[3]

mydata <- t(Pancan_expr[c(CD38_index,ADCY2_index,ITGA1_index),])
#mydata <- d
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,stats::var))
  for (i in 2:15) wss[i] <- sum(kmeans(mydata,
                                       centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
library(fpc)
pamk.best <- pamk(mydata)
cat("number of clusters estimated by optimum average silhouette width:", pamk.best$nc, "\n")
plot(pam(d, pamk.best$nc))


kc<-kmeans(mydata,5)
clusters<-kc$cluster
library(Rtsne)
tsne_out<-Rtsne(t(na.omit(a)),perplexity=100,step=5000)
tsne_out1<-Rtsne(t(na.omit(a)),perplexity=100,step=1000)
tsne_out2<-Rtsne(xcor10,perplexity=30,step=5000)
tsne_out3<-Rtsne(xcor10,perplexity=30,step=1000)
names(anno_color)<-names(sample_id_sub)
#plot(tsne_out$Y,col=as.factor(group_ls),pch=16)
group_ls<-as.factor(clusters)
# label_index<-table(clusters)
# label<-rep(NA,length(group_ls))
# label[label_index]<-names(label_index)
tsne.data<-data.frame(X=tsne_out3$Y[,1],Y=tsne_out3$Y[,2],Group=group_ls)
label.df <- data.frame(Group=tsne.data$Group)
label.df_2 <- tsne.data %>% 
  group_by(Group) %>% 
  summarize(X =median(X), Y = median(Y)) %>% 
  left_join(label.df)
#pdf("test.pdf")

ggplot(tsne.data,aes(x=X,y=Y,color=Group))+geom_point()+theme_minimal()+
  ggrepel::geom_label_repel(data = label.df_2, aes(label = Group))
#dev.off()

cluster_sample_ls<-list(G1=names(which(clusters==1)),G2=names(which(clusters==2)),G3=names(which(clusters==3)),
 G4=names(which(clusters==4)),G5=names(which(clusters==5)))#,G6=names(which(clusters==6)),G7=names(which(clusters==7)))
 #G8=names(which(clusters==8)), G8=names(which(clusters==9)))#,
# G7=names(which(clusters==7)),G8=names(which(clusters==8)),G9=names(which(clusters==9)))
library(survival)
library(survminer)
Time<-c()
Group<-c()
Status<-c()
for(cl in names(cluster_sample_ls)){
    sample_id<-cluster_sample_ls[[cl]]
    sample_id_sub<-sapply(sample_id,function(x) {paste(as.vector(unlist(strsplit(x,"\\."))[1:3]),collapse="-")})
    surv_sam<-intersect(row.names(Survival_info),sample_id_sub)
   # print(table(Sample_cancer[surv_sam,"TCGA.Study"]))
    Time<-c(Time,Survival_info[surv_sam,"DSS.time.cr"])
    Status<-c(Status,Survival_info[surv_sam,"DSS_cr"])
    Group<-c(Group,rep(cl,length(surv_sam)))
}
test.low.high <- list(time = Time, status = Status, group = as.factor(Group))
model1 <- survdiff(Surv(time, status) ~ group, data= test.low.high, na.action=na.exclude)
p_value<-1-pchisq(model1$chisq, df=length(levels(factor(Group)))-1)
fit <- survfit(Surv(time, status) ~ group, data=test.low.high, na.action=na.exclude)
#pdf("G1_G6_survival.pdf",6,6)
p<-ggsurvplot(fit,pval=round(p_value,4),font.x = 14,font.y=14,font.main=14)
print(p)


save.image(file="TCGA_TIL_cluster_20220115.rda",ascii=FALSE,compress=TRUE)
#dev.off()
df1<-Pancancer_expression[,cluster_sample_ls[[1]]]
df2<-Pancancer_expression[,cluster_sample_ls[[2]]]
df3<-Pancancer_expression[,cluster_sample_ls[[3]]]
df4<-Pancancer_expression[,cluster_sample_ls[[4]]]
df5<-Pancancer_expression[,cluster_sample_ls[[5]]]
df6<-Pancancer_expression[,cluster_sample_ls[[6]]]
df7<-Pancancer_expression[,cluster_sample_ls[[7]]]
df_list<-list(df1=df1,df2=df2,df3=df3,df4=df4,df5=df5,df6=df6,df7=df7)
diff_result<-list()
top10<-c()
for(fl in list.files("Cluster_diff_gene")){
    tmp<-as.vector(unlist(strsplit(fl,"\\.")))[1]
    res<-readRDS(paste("Cluster_diff_gene/",tmp,".RDS",sep=""))
    diff_result[[tmp]]<- res
    up_df<-res$up
    up_df<-up_df[which(up_df$FDR!="Inf"),]
    down_df<-res$down
    down_df<-down_df[which(down_df$FDR!="Inf"),]
    up10<-up_df[order(up_df$FDR,decreasing=T),][1:20,"gene"]
    down10<-down_df[order(down_df$FDR,decreasing=T),][1:20,"gene"]
    top10<-c(top10,up10,down10)#,down10)
}
top10_unique<-unique(top10)
top10_unique_symbol<-as.vector(unlist(sapply(top10_unique,function(x) {as.vector(unlist(strsplit(x,"\\|")))[1]})))

all_diff_gene<-unique(as.vector(unlist(diff_result)))
Immune_related_gene<-read.table("immune_related_genes/GeneList.txt",header=T,sep="\t")
Immport_gene_list<-list()
for(pa in unique(Immune_related_gene$Category)){
    Immport_gene_list[[pa]]<-Immune_related_gene[which(Immune_related_gene$Category==pa),"Symbol"]
}

expr<-Pancancer_expression[,as.vector(unlist(cluster_sample_ls))]
gene_ls<-row.names(expr)
gene_ls_symbol<-as.vector(unlist(sapply(gene_ls,function(x) {as.vector(unlist(strsplit(x,"\\|")))[1]})))
expr<-cbind(gene_ls_symbol,expr)
expr<-expr[!duplicated(expr$gene_ls_symbol),]
row.names(expr)<-expr[,1]
expr_new<-expr[,-1]
library(GSVA)
gsva(as.matrix(expr_new),Immport_gene_list,method="ssgsea",ssgsea.norm=TRUE)->Immport_ls_ES
pheatmap2(Immport_ls_ES,scale = "row",cluster_col=F,show_rownames=T,show_colnames=F,annotation_col=annotation_row,annotation_colors=annotation_color,#c("#037ef3","#f85a40","#00c16e","#7552cc","#f48924","#52565e","#caccd1"),#cluster_col=F,cluster_row=F,
         color = c(colorRampPalette(colors = c("blue","white"))(92),colorRampPalette(colors = c("white","red"))(100)),
         breaks=bk)#->ES_cor_pheatmap_SKCM

save.image(file="TCGA_TIL_Immport.rda",ascii=FALSE,compress=TRUE)

top10_imm<-top10_unique_symbol[which(top10_unique_symbol%in%Immune_related_gene$Symbol)]
top10_unique_imm<-top10_unique[which(top10_unique_symbol%in%Immune_related_gene$Symbol)]

gene_fre<-rep(0,length(all_diff_gene))
names(gene_fre)<-all_diff_gene
for(var in names(diff_result)){
    g_ls<-diff_result[[var]]
    gtmp<-c(g_ls$up,g_ls$down)
    gene_fre[gtmp] = gene_fre[gtmp]+1

}
gn<-names(gene_fre)[which(gene_fre>15)]
gn_symbol<-as.vector(unlist(sapply(gn,function(x) as.vector(unlist(strsplit(x,"\\|")))[1])))
legand_receptor<-read.table("immune_legand_receptor/Homo_sapiens.txt",header=T,sep="\t")
Immune_inter<-legand_receptor[intersect(which(legand_receptor$Ligand_symbol%in%gn_symbol),which(legand_receptor$Receptor_symbol%in%gn_symbol)),]
write.table(Immune_inter,file="Immune_inter.txt",sep="\t",quote=F)

annotation_row = data.frame(Cluster=factor(paste("C",clusters,sep="")))
row.names(annotation_row)<-names(clusters)
annotation_color = list(Cluster=c(C1="#037ef3",C2="#f85a40",C3="#00c16e",C4="#7552cc",C5="#f48924",C6="#52565e",C7="#caccd1"))
library(pheatmap2)
mt<-as.matrix(Pancancer_expression[top10_unique,as.vector(unlist(cluster_sample_ls))])
pdf("test1.pdf",12,12)
pheatmap2(mt,scale = "row",show_rownames=T,cluster_col=F,show_colnames=F,annotation_col=annotation_row,annotation_colors=annotation_color,#c("#037ef3","#f85a40","#00c16e","#7552cc","#f48924","#52565e","#caccd1"),#cluster_col=F,cluster_row=F,
         color = c(colorRampPalette(colors = c("blue","white"))(92),colorRampPalette(colors = c("white","red"))(100)),
         breaks=bk)#->ES_cor_pheatmap_SKCM
dev.off()

#cancer cluster distribution
library(dplyr)
library(stringr)
library(stringi)
all_sam<-as.vector(unlist(cluster_sample_ls))
sample_id_sub<-sapply(all_sam,function(x) {paste(as.vector(unlist(strsplit(x,"\\."))[1:3]),collapse="-")})
Sam_cancer<-data.frame(Sam=row.names(Sample_cancer),Cancer=Sample_cancer$TCGA.Study)
Group_sam<-data.frame(Sam=sample_id_sub,Group=c(rep("G1",length(cluster_sample_ls[[1]])),rep("G2",length(cluster_sample_ls[[2]])),
rep("G3",length(cluster_sample_ls[[3]])),rep("G4",length(cluster_sample_ls[[4]])),rep("G5",length(cluster_sample_ls[[5]])),
rep("G6",length(cluster_sample_ls[[6]])),rep("G7",length(cluster_sample_ls[[7]]))))
left_join(Group_sam,Sam_cancer)->Group_cancer_distribution
cancer_sam_sta<-table(Group_cancer_distribution$Cancer,Group_cancer_distribution$Group)
static_mt<-as.vector(unlist(cancer_sam_sta))
Group_cancer_distribution_df<-data.frame(Cancer=rep(row.names(cancer_sam_sta),times=7),Count=static_mt,Group=rep(colnames(cancer_sam_sta),each=nrow(cancer_sam_sta)))
Group_cancer_distribution_df$Count<-as.numeric(Group_cancer_distribution_df$Count)
pdf("Cancer_groupo_distribution_new_2022_0113.pdf",8,4)
p=ggplot(Group_cancer_distribution_df,mapping= aes(x=Cancer,y=Count,fill=Group))+theme(axis.text.x=element_text(angle=50,hjust = 1,vjust=1))+geom_bar(position = 'fill',stat='identity')+
        theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),axis.line = element_line(colour = "black",size = 0.2),plot.title = element_text(hjust = 0.5,size = 8))+
        labs(title =row.names(TME_ES)[i],x=NULL,y=NULL)+theme(axis.text.x =element_text(size=14), axis.text.y=element_text(size=14),plot.title = element_text(size=24))
print(p)
dev.off()




gtmp<-as.vector(unlist(sapply(gn,function(x) as.vector(unlist(strsplit(x,"\\|")))[2])))
# trans_file<-read.table("entrez_symbol_trans",header=T,sep="\t")
# trans_file_unique<-na.omit(trans_file[!duplicated(trans_file$SYMBOL),])
# row.names(trans_file_unique)<-trans_file_unique$SYMBOL
# gg<-as.vector(unlist(na.omit(trans_file_unique[gtmp,"ENTREZID"])))
kk <- enrichKEGG(gene         = gg,
                  organism     = 'hsa',
                    pvalueCutoff = 0.05)

library(foreach)
library(doParallel)
library(parallel)
cl <- makeCluster(6)
registerDoParallel(cl)
dots <- df_list  # 动态参数list
result<-foreach(i=seq(1,6)) %dopar% 
diff_func(dots[[1]],dots[[1+i]],1,1+i)  # 数据与参数组成list传入函数
stopCluster(cl)

cl <- makeCluster(5)
registerDoParallel(cl)
dots <- df_list  # 动态参数list
result<-foreach(i=seq(1,5)) %dopar% 
diff_func(dots[[2]],dots[[2+i]],2,2+i)  # 数据与参数组成list传入函数
stopCluster(cl)

cl <- makeCluster(4)
registerDoParallel(cl)
dots <- df_list  # 动态参数list
result<-foreach(i=seq(1,4)) %dopar% 
diff_func(dots[[3]],dots[[3+i]],3,3+i)  # 数据与参数组成list传入函数
stopCluster(cl)

cl <- makeCluster(3)
registerDoParallel(cl)
dots <- df_list  # 动态参数list
result<-foreach(i=seq(1,3)) %dopar% 
diff_func(dots[[4]],dots[[4+i]],4,4+i)  # 数据与参数组成list传入函数
stopCluster(cl)


cl <- makeCluster(2)
registerDoParallel(cl)
dots <- df_list  # 动态参数list
result<-foreach(i=seq(1,2)) %dopar% 
diff_func(dots[[5]],dots[[5+i]],5,5+i)  # 数据与参数组成list传入函数
stopCluster(cl)

diff_func(dots[[6]],dots[[7]],6,7)

for(dfx in seq(1,length(df_list))){
    print(dfx)
    for(dfy in seq(1,length(df_list))){
        if((dfx!=dfy)&(!(c(paste(dfx,dfy,sep="_"),paste(dfy,dfx,sep="_"))%in%names(diff_result)))){
            diff_re<-diff_func(df_list[[dfx]],df_list[[dfy]])
            diff_result[[paste(dfx,dfy,sep="_")]]=diff_re
        }
    }
}

diff_func<-function(df1,df2,ind1,ind2){
    median1<-apply(df1,1,median)
    median2<-apply(df2,1,median)
    df1_up<-names(which((median1/median2)>2))
    df1_down<-names(which((median1/median2)<0.5))
    mt_all<-cbind(df1,df2)[c(df1_up,df1_down),]
    #mt_down<-cbind(df1,df2)[df1_down,]
    a=ncol(df1)
    ct<-a+ncol(df2)
    pva_ls<-apply(mt_all,1,function(x) {wilcox_function(x[1:a],x[a+1:ct])})
    up_g<-intersect(df1_up,names(which(pva_ls<0.05)))
    down_g<-intersect(df1_down,names(which(pva_ls<0.05)))
    #return(list(up=up_g,down=down_g))
    up_df<-data.frame(gene=up_g,FDR=(median1/median2)[up_g],pva=pva_ls[up_g])
    down_df<-data.frame(gene=down_g,FDR=(median1/median2)[down_g],pva=pva_ls[down_g])
    result<-list(up=up_df,down=down_df)
    readr::write_rds(result,
      file = file.path(paste("Cluster_diff_gene/diff_",paste(ind1,ind2,sep="_"),".RDS",sep="")))
}

wilcox_function<-function(x,y){
    tmp<-wilcox.test(x,y)
    pva<-tmp$p.value
    return(pva)
}



#ADCY2 CD38 ITGA1
ADCY2<-as.numeric(as.vector(unlist(expr_new[307,])))
ADCY2_mutation<-read.table("mutation/ADCY2_mutation",sep="\t",fill=NA,quote="")
ADCY2_mut_sample<-ADCY2_mutation[,16]
ADCY2_id<-sapply(ADCY2_mut_sample,function(x) {paste(as.vector(unlist(strsplit(x,"\\-"))[1:4]),collapse="-")})

sample_id<-colnames(expr_new)
sample_id_sub<-sapply(sample_id,function(x) {paste(as.vector(unlist(strsplit(x,"\\."))[1:4]),collapse="-")})
#sample_id_sub1<-sapply(sample_id_sub,function(x) {substring(x, 1, str_length(x)-1)})
mutation_index<-which(sample_id_sub%in%ADCY2_id)
normal_index<-which(!sample_id_sub%in%ADCY2_id)
G1=colnames(expr_new)[mutation_index]
wilcox.test(ADCY2[normal_index],ADCY2[mutation_index])
G23=ifelse(ADCY2[normal_index] <= median(ADCY2[normal_index]),"Low","High")
names(G23)<-colnames(expr_new)[normal_index]
G2<-colnames(expr_new)[which(ADCY2[normal_index]<10)]
G3<-colnames(expr_new)[intersect(which(ADCY2[normal_index]>26),which(ADCY2[normal_index]<208))]
#cluster_sample_ls<-list(G2=names(which(G23=="Low")),G3=names(which(G23=="High")))#,G6=names(which(clusters==6)),G7=names(which(clusters==7)))
cluster_sample_ls<-list(G1=G1,G2=G2,G3=G3)#,G6=names(which(clusters==6)),G7=names(which(clusters==7)))

CD38<-as.numeric(as.vector(unlist(expr_new[3319,])))
ITGA1<-as.numeric(as.vector(unlist(expr_new[8598,])))
CD38_ITGA1<-CD38+ITGA1+ADCY2
cutgroup <- ifelse(ADCY2 <= median(ADCY2),"Low","High")
names(cutgroup)<-colnames(expr_new)
cluster_sample_ls<-list(G1=names(which(cutgroup=="Low")),G2=names(which(cutgroup=="High")))#,G6=names(which(clusters==6)),G7=names(which(clusters==7)))
 #G8=names(which(clusters==8)), G8=names(which(clusters==9)))#,
# G7=names(which(clusters==7)),G8=names(which(clusters==8)),G9=names(which(clusters==9)))
library(survival)
library(survminer)
Time<-c()
Group<-c()
Status<-c()
for(cl in names(cluster_sample_ls)){
    sample_id<-cluster_sample_ls[[cl]]
    sample_id_sub<-sapply(sample_id,function(x) {paste(as.vector(unlist(strsplit(x,"\\."))[1:3]),collapse="-")})
    surv_sam<-intersect(row.names(Survival_info),sample_id_sub)
   # print(table(Sample_cancer[surv_sam,"TCGA.Study"]))
    Time<-c(Time,Survival_info[surv_sam,"DSS.time.cr"])
    Status<-c(Status,Survival_info[surv_sam,"DSS_cr"])
    Group<-c(Group,rep(cl,length(surv_sam)))
}
test.low.high <- list(time = Time, status = Status, group = as.factor(Group))
model1 <- survdiff(Surv(time, status) ~ group, data= test.low.high, na.action=na.exclude)
p_value<-1-pchisq(model1$chisq, df=length(levels(factor(Group)))-1)
fit <- survfit(Surv(time, status) ~ group, data=test.low.high, na.action=na.exclude)
#pdf("G1_G6_survival.pdf",6,6)
p<-ggsurvplot(fit,pval=round(p_value,4),font.x = 14,font.y=14,font.main=14)
print(p)


#基因表达在不同癌症中对预后的影响
ADCY2<-as.numeric(as.vector(unlist(expr_new[307,])))

pdf("CD38_cancer_survival.pdf",4,4)
for(cancer in cancer_ls[-c(14,31)]){
    Sample_cancer%>%dplyr::filter(TCGA.Study==cancer)->spe_cancer_sample
    sample_id<-colnames(expr_new)
    sample_id_sub<-sapply(sample_id,function(x) {paste(as.vector(unlist(strsplit(x,"\\."))[1:3]),collapse="-")})
    spe_sam<-which(sample_id_sub%in%row.names(spe_cancer_sample))
    ITGA1<-as.numeric(as.vector(unlist(expr_new[3319,spe_sam])))
    if(length(spe_sam>0)){
        cutgroup <- ifelse(ITGA1 <= median(ITGA1),"Low","High")
        names(cutgroup)<-colnames(expr_new)[spe_sam]
        if(length(unique(as.factor(cutgroup))>1)){
            cluster_sample_ls<-list(G1=names(which(cutgroup=="Low")),G2=names(which(cutgroup=="High")))#,G6=names(which(clusters==6)),G7=names(which(clusters==7)))
            Time<-c()
            Group<-c()
            Status<-c()
            for(cl in names(cluster_sample_ls)){
                sample_id<-cluster_sample_ls[[cl]]
                sample_id_sub<-sapply(sample_id,function(x) {paste(as.vector(unlist(strsplit(x,"\\."))[1:3]),collapse="-")})
                surv_sam<-intersect(row.names(Survival_info),sample_id_sub)
            # print(table(Sample_cancer[surv_sam,"TCGA.Study"]))
                Time<-c(Time,Survival_info[surv_sam,"DSS.time.cr"])
                Status<-c(Status,Survival_info[surv_sam,"DSS_cr"])
                Group<-c(Group,rep(cl,length(surv_sam)))
            }
            if(cancer%in%c("TGCT","PCPG","THYM")){
                Status[which(Status==2)]=1
            }
            if(length(which(Status==1))>1){
                test.low.high <- list(time = Time, status = Status, group = as.factor(Group))
                model1 <- survdiff(Surv(time, status) ~ group, data= test.low.high, na.action=na.exclude)
                p_value<-1-pchisq(model1$chisq, df=length(levels(factor(Group)))-1)
                if(p_value<0.1){
                    print(cancer)
                    print(p_value)
                    fit <- survfit(Surv(time, status) ~ group, data=test.low.high, na.action=na.exclude)
                    #pdf("G1_G6_survival.pdf",6,6)
                    p<-ggsurvplot(fit,pval=round(p_value,4),font.x = 14,font.y=14,font.main=14)
                    print(p)
                }

            }
          

        }
        
       # cluster_sample_ls<-list(G1=names(which(cutgroup=="Low")),G2=names(which(cutgroup=="High")))#,G6=names(which(clusters==6)),G7=names(which(clusters==7)))
    #G8=names(which(clusters==8)), G8=names(which(clusters==9
        

    }
   
    # fit <- survfit(Surv(time, status) ~ group, data=test.low.high, na.action=na.exclude)
    # #pdf("G1_G6_survival.pdf",6,6)
    # p<-ggsurvplot(fit,pval=round(p_value,4),font.x = 14,font.y=14,font.main=14)
    # print(p)
}
dev.off()
