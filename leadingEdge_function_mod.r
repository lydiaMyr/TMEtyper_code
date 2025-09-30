library(DESeq2)
library(fgsea)
library(dplyr)
#library(BiocParallel)
library(data.table)
library(GSEA)
load("E:/project/remodeling tumor microenvironment/exp_count_cancer_ls_geneID.rda")
 #load("E:/project/remodeling tumor microenvironment/immune_path_sig_top20.rda")
preprocess_fun<-function(ca){
    exp<-exp_count_cancer_ls_geneID[[ca]]
    row.names(exp)<-exp[,1]
    exp<-exp[,-1]
    exp_num<-t(apply(exp,1,function(x) as.numeric(x)))
    colnames(exp_num)<-colnames(exp)
    #colnames(exp_count_cancer_ls_geneID[["TCGA-BRCA"]])
    sample_info<-as.vector(unlist(sapply(colnames(exp_num),function(x) as.vector(unlist(strsplit(x,"-")))[[4]])))
    ifelse(as.vector(sample_info) <"11A","Tumor","Normal")->condition
    colData <- data.frame(row.names=colnames(exp_num), condition)


   
    dds <- DESeqDataSetFromMatrix(countData = exp_num,
                                      colData = colData,
                                      design = ~ condition)



    dds$condition <- relevel(dds$condition, ref = "Normal")  
    keep <- rowSums(counts(dds)) >= 1.5*ncol(exp_num) 
    dds <- dds[keep,]
    dds <- DESeq(dds,quiet = F)
    res <- results(dds,contrast=c("condition", "Tumor", "Normal")) 
    resOrdered <- res[order(res$padj),]  
    tempDEG <- as.data.frame(resOrdered)
    DEG_DEseq2 <- na.omit(tempDEG)
    DEG_DEseq2%>%filter(padj<0.05)->sig_deg
    gene_rank<-sig_deg$`log2FoldChange`
    names(gene_rank)<-row.names(sig_deg)
    return(gene_rank)
    #leading_edge_result[[ca]]<<-fgseaRes
}

gene_rank<-preprocess_fun(ca)
fgseaRes <- fgsea(pathways = immune_path_sig,
                      stats    = gene_rank,
                      minSize  = 15,
                      maxSize  = 500)



# leading_edge_result<<-list()
# for(ca in names(exp_count_cancer_ls_geneID)){
#     print(ca)
#     #if(!(ca %in% c("TCGA-LGG","TCGA-OV")))
#     try(leading_edge_fun(ca))
# }


# save(leading_edge_result,file="pathway_leading_edge_result_final_ls.rda")
