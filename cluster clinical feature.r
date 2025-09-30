library(ggplot2)
#cluster clinical feature
load("E:/project/remodeling tumor microenvironment/Cluster_clinical_info_analysis/TCGA_33cancer_clinical_info.rda")
load("E:/project/remodeling tumor microenvironment/Cluster_clinical_info_analysis/Biology_assign.rda")
load("E:/project/remodeling tumor microenvironment/Cluster_clinical_info_analysis/cancer_sample_df.rda")
load("E:/project/remodeling tumor microenvironment/Cluster_clinical_info_analysis/clinical_info_analysis.rda")
library(openxlsx)
#cancer type distribution
cancer_type_dis_plot_df<-c()
sample_cluster<-TME_clust_assign_final[colnames(clust_abun_heatmap)]
for(subtype in unique(sample_cluster)){
    sample<-intersect(names(which(sample_cluster==subtype)),row.names(cancer_sample_df))
    cancer_type_dis_plot_df<-rbind(cancer_type_dis_plot_df,cbind(sample,cancer_sample_df[sample,1],rep(subtype,length(sample))))
}
colnames(cancer_type_dis_plot_df)<-c("Sample","Cancer","Cluster")
static<-table(cancer_type_dis_plot_df[,"Cancer"],cancer_type_dis_plot_df[,"Cluster"])
df_trans<-as.vector(unlist(static))
plot_df<-cbind(df_trans,rep(row.names(static),ncol(static)),rep(colnames(static),each=nrow(static)))
colnames(plot_df)<-c("Count","Cancer","Cluster")
plot_df<-data.frame(plot_df)
plot_df[,"Count"]<-as.numeric(plot_df[,"Count"])
cancer_group<-read.xlsx("E:/project/remodeling tumor microenvironment/Cluster_clinical_info_analysis/Cancer_group.xlsx",sheet=1)
row.names(cancer_group)<-paste("TCGA-",cancer_group[,1],sep="")
cancer_order<-intersect(row.names(cancer_group),row.names(static))
plot_df$Cancer<-factor(plot_df$Cancer,levels=cancer_order,ordered = TRUE)
Tissue<-cancer_group[cancer_order,"Tissue"]
anno<-table(Tissue)[unique(Tissue)]
anno_new<-anno
for(i in seq(2,length(anno))){
    anno_new[i]<-anno[i]+sum(anno[1:(i-1)])
}
pdf("E:/project/remodeling tumor microenvironment/Cluster_clinical_info_analysis/Cancer_TME_distribution0805.pdf",12,6)
ggplot(plot_df,aes(x=Cancer,y=Count,fill=Cluster))+
  geom_bar(stat = 'identity',position='fill',width = 0.5,colour='black')+
  theme(#axis.title =  element_text(size=12,face = "bold"),
        axis.text.x =   element_text(angle=90,  # 横坐标文字旋转九十度
                                     hjust = 1, # 调整横坐标文字位置
                                     size=15)) +annotate("rect", xmin = 0.5, xmax = anno_new[1]+0.5, ymin = 0, ymax = 1,
             alpha = 0, color= "red")+annotate("rect", xmin = anno_new[1]+0.5, xmax = anno_new[2]+0.5, ymin =  0, ymax = 1,
             alpha = 0, color= "red")+annotate("rect", xmin = anno_new[2]+0.5, xmax = anno_new[3]+0.5, ymin =  0, ymax = 1,
             alpha = 0, color= "red")+annotate("rect", xmin = anno_new[3]+0.5, xmax = anno_new[4]+0.5, ymin =  0, ymax = 1,
             alpha = 0, color= "red")+annotate("rect", xmin = anno_new[4]+0.5, xmax = anno_new[5]+0.5, ymin =  0, ymax = 1,
             alpha = 0, color= "red")+annotate("rect", xmin = anno_new[5]+0.5, xmax = anno_new[6]+0.5, ymin =  0, ymax = 1,
             alpha = 0, color= "red")+annotate("rect", xmin = anno_new[6]+0.5, xmax = anno_new[7]+0.5, ymin =  0, ymax = 1,
             alpha = 0, color= "red")+annotate("rect", xmin = anno_new[7]+0.5, xmax = anno_new[8]+0.5, ymin =  0, ymax = 1,
             alpha = 0, color= "red")+annotate("rect", xmin = anno_new[8]+0.5, xmax = anno_new[9]+0.5, ymin =  0, ymax = 1,
             alpha = 0, color= "red")
dev.off()
#sample stage distribution
# ca="SKCM"
# clin_info<-cl[[ca]]
# #stage_info<-clin_info[,"melanoma_clark_level_value"]
# sample_id<-intersect(names(sample_cluster),row.names(cancer_sample_df))
# sample_id_sub<-sapply(sample_id,function(x) {paste(as.vector(unlist(strsplit(x,"-"))[1:3]),collapse="-")})
# names(sample_cluster)<-sample_id_sub
# stage_sam<-intersect(clin_info$"bcr_patient_barcode",sample_id_sub)

# clin_info%>%dplyr::filter(bcr_patient_barcode%in%stage_sam)%>%dplyr::select(bcr_patient_barcode,melanoma_clark_level_value)->sam_stage_info

# cluster<-sample_cluster
# names(cluster)<-sam_id

# sam_stage_info$cluster<-as.vector(unlist(sample_cluster[stage_sam]))
# sam_stage_info%>%dplyr::filter(melanoma_clark_level_value%in%c("I","II","III","IV","V"))->sam_stage_info

pdf("E:/project/remodeling tumor microenvironment/Cluster_clinical_info_analysis/purity_boxplot.pdf",6,8)
ggplot(sample_purity, aes(x=cluster, y=purity, fill=cluster)) + 
  geom_boxplot()+
  labs(title="Tumor purity",x="Cluster", y = "Purity")+scale_fill_brewer(palette="Dark2") + theme_minimal()+stat_compare_means(comparisons = list(c("C1","C3"),c("C1","C4"),c("C1","C5"),c("C1","C6"),c("C3","C4"),c("C3","C5"),c("C3","C6"),c("C4","C5"),c("C4","C6")))
dev.off()
#overall survival 
Sample_survival_info<-c()
for(ca in names(cl)){
    clin_info<-cl[[ca]]
    patient <- clin_info %>%
    #dplyr::filter(vital_status == 'Alive') %>%
    dplyr::select(c(bcr_patient_barcode,gender,days_to_birth,vital_status,days_to_last_followup,race_list,days_to_death))%>%
    reshape::rename(c(bcr_patient_barcode = 'Barcode',
                        gender = 'Gender',
                        vital_status = 'OS',
                        days_to_last_followup='OS.Time',
                        race_list = 'Race',
                        days_to_birth = 'Age' )) %>%
    mutate(OS=ifelse(OS=='Dead',1,0))%>%
    mutate(OS.Time=as.numeric(OS.Time)/365)
    Sample_survival_info<-rbind(Sample_survival_info,patient)
}
#survival plot
library(dplyr)
library(survminer)
library(survival)
surv_plot_df<-c()
for(subtype in unique(sample_cluster)){
    sample_id<-intersect(names(which(sample_cluster==subtype)),row.names(cancer_sample_df))
    sample_id_sub<-sapply(sample_id,function(x) {paste(as.vector(unlist(strsplit(x,"-"))[1:3]),collapse="-")})
    surv_sam<-intersect(Sample_survival_info$Barcode,sample_id_sub)
    Sample_survival_info%>%dplyr::filter(Barcode%in%surv_sam)%>%dplyr::select(OS,OS.Time)->sam_surv_info
    sam_surv_info<-cbind(sam_surv_info,rep(subtype,length(surv_sam)))
    surv_plot_df<-rbind(surv_plot_df,sam_surv_info)
}
colnames(surv_plot_df)<-c("Status","Time","Group")
surv_plot_df$Group<-as.vector(unlist(sapply(surv_plot_df$Group,function(x) as.vector(unlist(strsplit(x,"_")))[3])))
test.low.high <- list(time = as.numeric(surv_plot_df$Time), status = surv_plot_df$Status, group = as.factor(surv_plot_df$Group))
model1 <- survdiff(Surv(time, status) ~ group, data= test.low.high, na.action=na.exclude)
p_value<-1-pchisq(model1$chisq, df=length(levels(factor(test.low.high$group)))-1)
fit <- survfit(Surv(time, status) ~ group, data=test.low.high, na.action=na.exclude)
pdf("E:/project/remodeling tumor microenvironment/Cluster_clinical_info_analysis/Group_survival0805.pdf",6,6)
ggsurvplot(fit,pval=round(p_value,4),font.x = 14,font.y=14,font.main=14)
dev.off()

#tumor purity
library(ggpubr)
purity<-read.table("E:/project/remodeling tumor microenvironment/NMF_test/absolute_purity",header=T,sep="\t",row.names=1)
sam_id<-sapply(names(sample_cluster),function(x) paste(as.vector(unlist(strsplit(x,"-")))[1:4],collapse="-"))
cluster<-sample_cluster
names(cluster)<-sam_id
purity_sam_id<-sapply(as.vector(unlist(purity$sample)),function(x) paste(as.vector(unlist(strsplit(x,"-")))[1:4],collapse="-"))
purity$id_trans<-as.vector(unlist(purity_sam_id))
purity%>%dplyr::filter(id_trans%in%sam_id)%>%dplyr::select(id_trans,purity)->sample_purity
sample_purity$cluster<-as.vector(unlist(cluster[sample_purity$id_trans]))
sample_purity$cluster<-as.vector(unlist(sapply(sample_purity$cluster,function(x) paste("C",as.vector(unlist(strsplit(x,"_")))[3],sep=""))))

pdf("E:/project/remodeling tumor microenvironment/Cluster_clinical_info_analysis/purity_boxplot0905.pdf",6,4)
ggplot(sample_purity, aes(x=cluster, y=purity, fill=cluster)) + 
  geom_boxplot(width=0.5)+
  labs(title="Tumor purity",x="Cluster", y = "Purity")+scale_fill_brewer(palette="Dark2") + theme_minimal()+stat_compare_means(comparisons = list(c("C1","C3"),c("C1","C4"),c("C1","C5"),c("C1","C6"),c("C3","C4"),c("C3","C5"),c("C3","C6"),c("C4","C5"),c("C4","C6")))
dev.off()

save.image(file="E:/project/remodeling tumor microenvironment/Cluster_clinical_info_analysis/clinical_info_analysis.rda",ascii=FALSE,compress=TRUE)
