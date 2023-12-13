library(sva)
library(bladderbatch)
library(pamr)
library(limma)
library(EnhancedVolcano)
library(ggplot2)
library(gmodels)
library(ggpubr)
library(reshape2)
library(dplyr)
library(survival)
library(survminer)
library(edgeR)
library(tibble)
library(clinfun)
library(survival)
library(survminer)
library(forestplot)


#install.packages("memoise")
#install.packages("cachem")

########################### pre-process of TCGA LUAD RNA-seq data ########################

# read TCGA LUAD phenotype data
phenotype<-read.csv("phenotype.tsv",sep="\t")
colnames(phenotype)
phenotype$Sample_ID<-phenotype$submitter_id.samples

# read TCGA LUAD RNA-seq data
rnaseq<-read.csv("rna_seq.tsv",sep="\t")
dim(rnaseq)
gsub("\\.","-",colnames(rnaseq)[-1]) %in% phenotype$submitter_id.samples
colnames(rnaseq)[-1]<-gsub("\\.","-",colnames(rnaseq)[-1])
rownames(rnaseq)<-rnaseq$Ensembl_ID

# read Ensembl_ID to gene name annotation table
id_to_gene<-read.csv("id_to_gene.csv")
head(id_to_gene)
rnaseq<-merge(rnaseq,id_to_gene[,c(1,2)],by.x="Ensembl_ID",by.y="id")

# process duplicate genes
duplicate_gene<-names(table(rnaseq$gene)[table(rnaseq$gene)>1])
rnaseq_unique<-rnaseq[!rnaseq$gene %in% duplicate_gene,]
dim(rnaseq_unique)
rnaseq_duplicate<-rnaseq[rnaseq$gene %in% duplicate_gene,]
rnaseq_duplicate_id<-rnaseq_duplicate[,c("Ensembl_ID","gene")]
# change from log2transform to raw count
rnaseq_duplicate_matrix<-subset(rnaseq_duplicate,select=-c(Ensembl_ID,gene))
rnaseq_duplicate_matrix<-(2^rnaseq_duplicate_matrix)-1
rnaseq_duplicate_matrix$gene<-rnaseq_duplicate_id$gene

unique_gene<-c()
dedup<-list()
for (g in unique(rnaseq_duplicate_id$gene)){
  dup_data<-rnaseq_duplicate_matrix[rnaseq_duplicate_matrix$gene==g,]
  dedup_data<-apply(dup_data[,-c(586)],2,sum)
  dedup[[g]]<-dedup_data
}
dedup_matrix<-do.call(rbind,dedup)
View(dedup_matrix)

rnaseq_unique<-rnaseq[!rnaseq$gene %in% duplicate_gene,]
rnaseq_unique_id<-rnaseq_unique[,c("Ensembl_ID","gene")]
rnaseq_unique_matrix<-subset(rnaseq_unique,select=-c(Ensembl_ID,gene))
# change from log2transform to raw count
rnaseq_unique_matrix<-(2^rnaseq_unique_matrix)-1
sum(rownames(rnaseq_unique_id)==rownames(rnaseq_unique_matrix))
nrow(rnaseq_unique_matrix)
rownames(rnaseq_unique_matrix)<-rnaseq_unique_id$gene

# generate RNA-seq raw count table
rnaseq_matrix_count<-rbind(rnaseq_unique_matrix,dedup_matrix)

# only include primary tumor samples
phenotype<-phenotype[phenotype$sample_type.samples=="Primary Tumor",]

# classify TNM stage into 4 groups
phenotype$TNM<-NA
phenotype$TNM[phenotype$pathologic_M %in% c("M1","M1a","M1b")]<-"M"
phenotype$TNM[phenotype$pathologic_N %in% c("N2","N3") & is.na(phenotype$TNM)]<-"N2_N3"
phenotype$TNM[phenotype$pathologic_N %in% c("N1") & is.na(phenotype$TNM)]<-"N1"
phenotype$TNM[phenotype$pathologic_N %in% c("N0") & phenotype$pathologic_M %in% c("M0") & is.na(phenotype$TNM)]<-"N0_M0"
phenotype<-phenotype[!is.na(phenotype$TNM),]
dim(phenotype)

########################## Differential analysis #############################

#' Identify differentially-expressed genes by limma (controlled by tumor purity or other covariates)and generate file for GSEA analysis. If tumor purity is available, it can be treat as a covariate, if not, purity can be assigned as 1
#'
#' @param count_table,RNA-seq raw count table, rows are gene, columns are samples
#' @param feature_table,feature table containing target feature for example tumor, normal
#' @param cutoff, for one gene, at least how many samples with cpm >1
#' @param response,features to be compared, for example, stage I vs.stage IV
#' @param included_samples
#' @param dge, whether differential analysis is performed
#' @param class_level,comparion level
#' @param gsea, whether generating files for subsequent GSEA analysis
#'
#' @return a list containing all results
#' @export
#'
#' @examples
RNA_differental_analysis_limma<-function(count_table=rnaseq_matrix_count,feature_table=phenotype,cutoff=5,response="TNM",included_samples=phenotype$Sample_ID,dge="yes",class_level=c("N0_M0","M"),gsea="yes"){
  count_table<-count_table[,colnames(count_table) %in% included_samples]
  feature_table<-feature_table[feature_table$Sample_ID %in% included_samples,]
  matrix_list<-list()
  group_list<-c()
  purity_list<-c()
  for (c in class_level){
    class_matrix<-count_table[,colnames(count_table) %in% feature_table$Sample_ID[feature_table[,response]==c]]
    matrix_list[[c]]<-class_matrix
    group_list<-c(group_list,rep(c,ncol(class_matrix)))
    purity_list<-c(purity_list,feature_table$tumor_purity[match(colnames(class_matrix),feature_table$Sample_ID)])
  }
  order_matrix<-do.call(cbind,matrix_list)
  sample_info.edger<-factor(group_list,levels=class_level)
  edgeR.DGElist <- DGEList(counts = order_matrix, group=group_list)
  keep<-rowSums(cpm(edgeR.DGElist)>=1)>=cutoff
  edgeR.DGElist.keep<-edgeR.DGElist[keep,]
  edgeR.DGElist.keep$samples$lib.size<-colSums (edgeR.DGElist.keep$counts)
  #only calculate norm.factors, count matrix is not normalized yet
  edgeR.DGElist.keep<-calcNormFactors(edgeR.DGElist.keep,method ="TMM")
  logCPM<-cpm(edgeR.DGElist.keep, log = TRUE, normalized.lib.sizes=TRUE)
  sample_information_edger<-edgeR.DGElist.keep$samples
  result<-list()
  result[["logCPM_matrix"]]<-logCPM
  result[["Sample_class"]]<-sample_information_edger
  result[["group_list"]]<-group_list
  if(dge=="yes"){
    design <- model.matrix(~0+group_list+purity_list)
    real_level<-c(gsub("group_list","",colnames(design)[1:2]),"purity")
    colnames(design)<-c("level1","level2","purity")
    real_level<-data.frame(real_level=real_level,design_level=colnames(design))
    contr.matrix <- makeContrasts(
      level1_vs_level2 =level1-level2,
      levels = colnames(design))
    #design<-model.matrix(~group_list)
    v <- voom(edgeR.DGElist.keep, design, plot=TRUE)
    vfit <- lmFit(v, design)
    vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
    efit <- eBayes(vfit)
    plotSA(efit, main="Final model: Mean-variance trend")
    summary(decideTests(efit))
    DEGenes<-topTable(efit,coef=1, n=Inf)
    #tfit <- treat(vfit, lfc=1)
    #dt <- decideTests(tfit)
    #summary(dt)
    result[["differential_gene_purity_adjusted"]]<-DEGenes
    result[["real_level"]]<-real_level

    design <- model.matrix(~0+group_list)
    real_level<-c(gsub("group_list","",colnames(design)[1:2]))
    colnames(design)<-c("level1","level2")
    real_level_purity_not_adjusted<-data.frame(real_level=real_level,design_level=colnames(design))
    result[["real_level_purity_not_adjusted"]]<-real_level_purity_not_adjusted

    contr.matrix <- makeContrasts(
      level1_vs_level2 =level1-level2,
      levels = colnames(design))
    #design<-model.matrix(~group_list)
    v <- voom(edgeR.DGElist.keep, design, plot=TRUE)
    vfit <- lmFit(v, design)
    vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
    efit <- eBayes(vfit)
    plotSA(efit, main="Final model: Mean-variance trend")
    summary(decideTests(efit))
    DEGenes_purity_not_adjusted<-topTable(efit,coef=1, n=Inf)
    result[["differential_gene_purity_not_adjusted"]]<-DEGenes_purity_not_adjusted
  }
  result[["logCPM"]]<-v$E
  if(gsea=="yes"){
    line1<-c(length(group_list),length(unique(group_list)),1)
    line2<-c("#",unique(group_list))
    line3<-group_list
    max_n<-max(length(line1),length(line2),length(line3))
    length(line1)<-max_n
    length(line2)<-max_n
    length(line3)<-max_n
    cls_file<-rbind(line1,line2,line3)
    cls_file[is.na(cls_file)]<-""
    expression_txt<-v$E
    expression_txt<-as.data.frame(expression_txt)
    expression_txt<-add_column(expression_txt, NAME =rownames(expression_txt), .before = 1)
    expression_txt<-add_column(expression_txt, Description =rep(1,nrow(expression_txt)), .before = 2)
    result[["cls_file"]]<-cls_file
    result[["expression_txt_for_gsea"]]<-expression_txt
  }
  return(result)
}


DE_N0_M0_vs_M<-RNA_differental_analysis_limma()


################# single sample GSEA analysis #############################
immune_list<- split(as.matrix(gene_set)[,1], gene_set[,2])
gsva_matrix_N_M<- gsva(logCPM, immune_list,method='gsva', abs.ranking=TRUE)
median(gsva_matrix_N_M[rownames(gsva_matrix_N_M)== "Type 2 T helper cell",1:225])
median(gsva_matrix_N_M[rownames(gsva_matrix_N_M)== "Type 2 T helper cell",226:314])
median(gsva_matrix_N_M[rownames(gsva_matrix_N_M)== "Type 2 T helper cell",315:385])
median(gsva_matrix_N_M[rownames(gsva_matrix_N_M)== "Type 2 T helper cell",386:410])
wilcox.test(gsva_matrix_N_M[rownames(gsva_matrix_N_M)== "Type 2 T helper cell",386:410],gsva_matrix_N_M[rownames(gsva_matrix_N_M)== "Type 2 T helper cell",315:385])
cor.test(gsva_matrix_N_M[rownames(gsva_matrix_N_M)=="Type 1 T helper cell"],gsva_matrix_N_M[rownames(gsva_matrix_N_M)=="Type 2 T helper cell"])

######################## Analysis of PD-L1 expression in different stage samples ###############
CD274_logCPM<-rbind(data.frame(TNM="N0_M0",CD274=logCPM[rownames(logCPM)== "CD274",1:225]),data.frame(TNM="N1",CD274=logCPM[rownames(logCPM)== "CD274",226:314]),data.frame(TNM="N2_N3",CD274=logCPM[rownames(logCPM)== "CD274",315:385]),data.frame(TNM="M",CD274=logCPM[rownames(logCPM)== "CD274",386:410]))
wilcox.test(logCPM[rownames(logCPM)== "CD274",1:314],logCPM[rownames(logCPM)== "CD274",315:410])
logCPM[rownames(logCPM)== "CD274",colnames(logCPM) %in% c("M.TCGA-93-A4JN-01A","M.TCGA-97-8171-01A" )]
logCPM[rownames(logCPM)== "CD274",colnames(logCPM) %in% c("M.TCGA-93-A4JP-01A","M.TCGA-55-8620-01A","M.TCGA-L9-A5IP-01A","M.TCGA-55-8512-01A","M.TCGA-55-8094-01A" )]

# read TCGA LUAD xCell data
xcell_TCGA<-read.csv("xCell_TCGA_RSEM.txt",sep="\t")
dim(xcell_TCGA)

# read TCGA LUAD patient information table
patient_TCGA<-read.csv("TCGA_patient.csv")

# annotate xCell data
rownames(xcell_TCGA)<-xcell_TCGA$X
xcell_TCGA<-xcell_TCGA[,c(-1)]
colnames(xcell_TCGA)<-gsub("\\.01","",colnames(xcell_TCGA))
colnames(xcell_TCGA)<-gsub("\\.","-",colnames(xcell_TCGA))
patient_TCGA$PATIENT_ID[!patient_TCGA$PATIENT_ID %in% colnames(xcell_TCGA)]
xcell_TCGA<-t(xcell_TCGA)
xcell_TCGA<-as.data.frame(xcell_TCGA)
xcell_TCGA$PATIENT_ID<-rownames(xcell_TCGA)
xcell_TCGA_annotation<-merge(xcell_TCGA,patient_TCGA,by="PATIENT_ID")
dim(xcell_TCGA_annotation)
table(xcell_TCGA_annotation$AJCC_METASTASIS_PATHOLOGIC_PM)

xcell_TCGA_annotation$TNM<-NA
xcell_TCGA_annotation$TNM[xcell_TCGA_annotation$AJCC_METASTASIS_PATHOLOGIC_PM %in% c("M1","M1a","M1b")]<-"M"
xcell_TCGA_annotation$TNM[xcell_TCGA_annotation$AJCC_NODES_PATHOLOGIC_PN %in% c("N2","N3") & is.na(xcell_TCGA_annotation$TNM)]<-"N2_N3"
xcell_TCGA_annotation$TNM[xcell_TCGA_annotation$AJCC_NODES_PATHOLOGIC_PN %in% c("N1") & is.na(xcell_TCGA_annotation$TNM)]<-"N1"
xcell_TCGA_annotation$TNM[xcell_TCGA_annotation$AJCC_NODES_PATHOLOGIC_PN %in% c("N0") & xcell_TCGA_annotation$AJCC_METASTASIS_PATHOLOGIC_PM %in% c("M0") & is.na(xcell_TCGA_annotation$TNM)]<-"N0_M0"
table(xcell_TCGA_annotation$TNM)
xcell_TCGA_annotation_na_rm<-xcell_TCGA_annotation[!is.na(xcell_TCGA_annotation$TNM),]
dim(xcell_TCGA_annotation_na_rm)
xcell_TCGA_annotation_na_rm$TNM_digit<-1
xcell_TCGA_annotation_na_rm$TNM_digit[xcell_TCGA_annotation_na_rm$TNM=="N1"]<-2
xcell_TCGA_annotation_na_rm$TNM_digit[xcell_TCGA_annotation_na_rm$TNM=="N2_N3"]<-3
xcell_TCGA_annotation_na_rm$TNM_digit[xcell_TCGA_annotation_na_rm$TNM=="M"]<-4
# trend test
jonckheere.test(xcell_TCGA_annotation_na_rm$`CD8+ Tcm`, xcell_TCGA_annotation_na_rm$TNM_digit,alternative="decreasing",nperm=500)
xcell_TCGA_annotation_na_rm$TNM_binary<-0
xcell_TCGA_annotation_na_rm$TNM_binary[xcell_TCGA_annotation_na_rm$TNM_digit==4]<-1

# calculate ImmuneScore, StromaScore and MicroenvironmentScore based on the method of xCell paper
microenvironmentScores <- function(adjustedScores) {
  ImmuneScore = apply(adjustedScores[,c('B-cells','CD4+ T-cells','CD8+ T-cells','DC','Eosinophils','Macrophages','Monocytes','Mast cells','Neutrophils','NK cells')],1,sum)/1.5
  StromaScore = apply(adjustedScores[,c('Adipocytes','Endothelial cells','Fibroblasts')],1,sum)/2
  MicroenvironmentScore = ImmuneScore+StromaScore
  adjustedScores = cbind(adjustedScores,ImmuneScore,StromaScore,MicroenvironmentScore)
  return(adjustedScores)
}
xcell_TCGA_annotation_na_rm_immune_score<-microenvironmentScores(xcell_TCGA_annotation_na_rm)
colnames(xcell_TCGA_annotation_na_rm_immune_score)

# stage N0 is assigned as 0, others are assigned as 1
xcell_TCGA_annotation_na_rm_immune_score$TNM_binary_2<-0
xcell_TCGA_annotation_na_rm_immune_score$TNM_binary_2[xcell_TCGA_annotation_na_rm_immune_score$TNM_digit %in% c(2,3,4)]<-1


#' perform wilcoxon, trend or kruskal test for a numerical feature against a categorical feature
#'
#' @param feature_table,containing at least a numerical feature and a categorical feature
#' @param numeric_features,one or more numerical features to be tested
#' @param category_features, one or more categorical features
#' @param test, "wilcox" for Wilcoxon test, "trend" for jonckheere test, "multi" for kruskal test
#'
#' @return median values of numerical features against categorical features as well as test p-values
#' @export
#'
#' @examples

feature_wilcox_trend_test<-function(feature_table=xcell_TCGA_annotation_na_rm,numeric_features,category_features,test="wilcox"){
  result<-list()
  for (n_f in numeric_features){
    for (c_f in category_features){
      if (n_f %in% colnames(feature_table) & c_f %in% colnames(feature_table)){
        n_c_table<-feature_table[,c(n_f,c_f)]
        colnames(n_c_table)<-c("numeric_feature","category_feature")
        group_mean<-aggregate(. ~ category_feature, data = n_c_table, FUN = mean)
        group_mean<-as.data.frame(t(group_mean))
        colnames(group_mean)<-as.character(group_mean[1,])
        group_mean<-group_mean[c(-1),]
        if(test=="wilcox"){
          p_value<-wilcox.test(numeric_feature~category_feature,data=n_c_table)$p.value
          n_c_result<-data.frame(numeric_feature=n_f,category_feature=c_f,pvalue=p_value,test=test)
          nc_result<-cbind(n_c_result,group_mean)
          result[[paste(n_f,c_f,test,sep="_")]]<-nc_result
        }
        if(test=="trend"){
          p_value_increasing<-jonckheere.test(n_c_table$numeric_feature, n_c_table$category_feature,alternative="increasing",nperm=500)$p.value
          p_value_decreasing<-jonckheere.test(n_c_table$numeric_feature, n_c_table$category_feature,alternative="decreasing",nperm=500)$p.value
          p_value_two_side<-jonckheere.test(n_c_table$numeric_feature, n_c_table$category_feature,nperm=500)$p.value
          n_c_result<-data.frame(numeric_feature=n_f,category_feature=c_f,pvalue_decreasing=p_value_decreasing,pvalue_increasing=p_value_increasing,pvalue_two_side=p_value_two_side,test=test)
          nc_result<-cbind(n_c_result,group_mean)
          result[[paste(n_f,c_f,test,sep="_")]]<-nc_result
        }
        if(test=="multi"){
          pvalue<-kruskal.test(n_c_table$numeric_feature, n_c_table$category_feature,data=n_c_table)$p.value
          n_c_result<-data.frame(numeric_feature=n_f,category_feature=c_f,pvalue=pvalue,test=test)
          nc_result<-cbind(n_c_result,group_mean)
          result[[paste(n_f,c_f,test,sep="_")]]<-nc_result
        }

      }
    }
  }
  final_result<-do.call(rbind,result)
  return(final_result)
}

# test the trend of 64 cell subsets and 3 scores against TNM stages
continuous_feature<-colnames(xcell_TCGA_annotation_na_rm_immune_score)[c(2:65,147:150)]
TNM_digit_trend_test<-feature_wilcox_trend_test(feature_table=xcell_TCGA_annotation_na_rm_immune_score,numeric_features=continuous_feature,category_features="TNM_digit",test="trend")
View(TNM_digit_trend_test)

# test the differential expression of 64 cell subsets in M0_N0 (non-metastasis)and others (metastasis)
continuous_feature<-colnames(xcell_TCGA_annotation_na_rm_immune_score)[2:65]
TNM_binary_2_wilcox_test_mean<-feature_wilcox_trend_test(feature_table=xcell_TCGA_annotation_na_rm_immune_score,numeric_features=continuous_feature,category_features="TNM_binary_2",test="wilcox")
TNM_binary_2_wilcox_test_mean$q.value<-p.adjust(TNM_binary_2_wilcox_test_mean$pvalue,method="BH")
TNM_binary_2_wilcox_test_mean$logFC<-log2(TNM_binary_2_wilcox_test_mean$`1`/TNM_binary_2_wilcox_test_mean$`0`)
TNM_binary_2_wilcox_test_mean$log10qvalue<-(-log10(TNM_binary_2_wilcox_test_mean$q.value))

# volcano plot to display differential subsets
pdf("M0_N0_vs_other_differential_volcano_q_0_0_5.pdf",width=8,height=8)
EnhancedVolcano(TNM_binary_2_wilcox_test_mean,lab = as.character(TNM_binary_2_wilcox_test_mean$numeric_feature), x="logFC",y="q.value",
                pCutoff = 0.05, FCcutoff = 0.20,ylim=c(0,2.5),xlim=c(-2.4,1.0))
dev.off()


#' generate box plot with multiple continuous features against a categorical feature
#'
#' @param feature_table,a data frame containing continuous features and categorical features
#' @param include_type, included continuous features
#' @param category_variable, a categorical feature
#'
#' @return a list, p1 ggboxplot style; p2 ggplot style; dataframe original melt tabale
#' @export
#'
#' @examples
box_plot_tcga_validation<-function(feature_table=xcell_TCGA_annotation_na_rm_immune_score,include_type=c("Th1 cells","Th2 cells","Th2_Th1_ratio"),category_variable="TNM_digit"){
  plot_data<-as.data.frame(feature_table[,c(include_type,category_variable)])
  plot_data<-melt(plot_data,id.vars=category_variable)
  colnames(plot_data)<-c("Group","Cell_type","Enrichment_score")
  plot_data$Group<-as.factor(plot_data$Group)
  if(length(unique(plot_data$Cell_type))>=2){
    p1=ggboxplot(plot_data, x ="Group", y = "Enrichment_score",outlier.shape = NA,
                 fill = "Group", palette = "Dark2",add = "jitter",add.params=list(size=0.3,alpha=0.4),legend="none")+facet_wrap(~Cell_type)
    p2=ggplot(plot_data, aes(x=Group, y=Enrichment_score, fill=Group)) +
      geom_boxplot(alpha=0.7) +
      theme(legend.position="none")+facet_wrap(~Cell_type, nrow=1)+theme_bw()
  }
  else{
    p1=ggboxplot(plot_data, x ="Group", y = "Enrichment_score",outlier.shape = NA,
                 fill = "Group", palette = "Dark2",add = "jitter",add.params=list(size=0.3,alpha=0.4),legend="none")
    p2=ggplot(plot_data, aes(x=Group, y=Enrichment_score, fill=Group)) +
      geom_boxplot(alpha=0.7) +
      theme(legend.position="none")+theme_bw()
  }

  result<-list()
  result[["p1"]]<-p1
  result[["p2"]]<-p2
  result[["dataframe"]]<-plot_data
  return(result)
}

############## plot box plots of various cell subsets against Stage #########################

xcell_cell_type<-read.csv("xcell_cell_type.csv")
xcell_cell_type$Cell.types[!xcell_cell_type$Cell.types %in% colnames(xcell_TCGA_annotation_na_rm_immune_score)]
xcell_TCGA_annotation_na_rm_immune_score$Stage<-factor(xcell_TCGA_annotation_na_rm_immune_score$TNM,levels=c("N0_M0","N1","N2_N3","M"))

CD4_cells<-c(as.character(xcell_cell_type$Cell.types[xcell_cell_type$Parent.Child=="CD4+ T-cells"]),"CD4+ T-cells")
CD4_type_box<-box_plot_tcga_validation(include_type=CD4_cells[-7],category_variable="Stage")
pdf("CD4_TNM_box.pdf",width=8,height=6)
CD4_type_box$p1
dev.off()

CD8_cells<-c(as.character(xcell_cell_type$Cell.types[xcell_cell_type$Parent.Child=="CD8+ T-cells"]),"CD8+ T-cells")
CD8_type_box<-box_plot_tcga_validation(include_type=CD8_cells,category_variable="Stage")
pdf("CD8_TNM_box.pdf",width=7,height=5)
CD8_type_box$p1
dev.off()

macrophage_cells<-c(as.character(xcell_cell_type$Cell.types[xcell_cell_type$Parent.Child=="Macrophages"]))
macrophage_cells_box<-box_plot_tcga_validation(include_type=macrophage_cells,category_variable="Stage")
pdf("macrophage_cells_box.pdf",width=6,height=4)
macrophage_cells_box$p1
dev.off()

Epithelial_cells<-c(as.character(xcell_cell_type$Cell.types[xcell_cell_type$Parent.Child=="Epithelial cells"]),"Epithelial cells")
Epithelial_cells_box<-box_plot_tcga_validation(include_type=Epithelial_cells,category_variable="Stage")
pdf("Epithelial_cells_box.pdf",width=7,height=4.5)
Epithelial_cells_box$p1
dev.off()

B_cells<-c(as.character(xcell_cell_type$Cell.types[xcell_cell_type$Parent.Child=="B-cells"]),"B-cells")
B_cells_type_box<-box_plot_tcga_validation(include_type=B_cells,category_variable="Stage")
pdf("B_cells_TNM_box.pdf",width=8,height=6)
B_cells_type_box$p1
dev.off()

DC_cells<-c(as.character(xcell_cell_type$Cell.types[xcell_cell_type$Parent.Child=="DC"]),"DC")
DC_cells_type_box<-box_plot_tcga_validation(include_type=DC_cells,category_variable="Stage")
pdf("DC_cells_TNM_box.pdf",width=7,height=5)
DC_cells_type_box$p1
dev.off()

###################### Overall survival analysis of TCGA LUAD data ###########################
# whether Th2 and Th1 enrichment score can be independent predictor for patient overall survival

# pre-process survival information of the TCGA LUAD data
os_status<-as.character(xcell_TCGA_annotation_na_rm_immune_score$OS_STATUS)
os_s<-c()
for (df in os_status){
  if (df=="[Not Available]"){
    os_s<-c(os_s,NA)
  }
  else{
    df_split<-strsplit(df,":")[[1]][[1]]
    os_s<-c(os_s,df_split)
  }
}

xcell_TCGA_annotation_na_rm_immune_score$os_status<-os_s
xcell_TCGA_annotation_na_rm_immune_score$Th2<-xcell_TCGA_annotation_na_rm_immune_score$`Th2 cells`
xcell_TCGA_annotation_na_rm_immune_score$Th1<-xcell_TCGA_annotation_na_rm_immune_score$`Th1 cells`
xcell_TCGA_annotation_na_rm_immune_score$CD8<-xcell_TCGA_annotation_na_rm_immune_score$`CD8+ T-cells`

xcell_TCGA_annotation_na_rm_immune_score$os_status<-as.numeric(xcell_TCGA_annotation_na_rm_immune_score$os_status)
xcell_TCGA_annotation_na_rm_immune_score$os_months<-xcell_TCGA_annotation_na_rm_immune_score$OS_MONTHS
xcell_TCGA_annotation_na_rm_immune_score$os_months[xcell_TCGA_annotation_na_rm_immune_score$os_months=="[Not Available]"]<-NA
xcell_TCGA_annotation_na_rm_immune_score$os_months<-as.numeric(as.character(xcell_TCGA_annotation_na_rm_immune_score$os_months))

xcell_TCGA_annotation_na_rm_immune_score$SEX<-as.character(xcell_TCGA_annotation_na_rm_immune_score$SEX)
xcell_TCGA_annotation_na_rm_immune_score$SEX[xcell_TCGA_annotation_na_rm_immune_score$SEX=="MALE"]<-"Male"
xcell_TCGA_annotation_na_rm_immune_score$age<-as.numeric(as.character(xcell_TCGA_annotation_na_rm_immune_score$AGE))

xcell_TCGA_annotation_na_rm_immune_score$pStage<-NA
xcell_TCGA_annotation_na_rm_immune_score$pStage[xcell_TCGA_annotation_na_rm_immune_score$AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("Stage I","Stage IA","Stage IB")]<-"I"
xcell_TCGA_annotation_na_rm_immune_score$pStage[xcell_TCGA_annotation_na_rm_immune_score$AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("Stage II","Stage IIA","Stage IIB")]<-"II"
xcell_TCGA_annotation_na_rm_immune_score$pStage[xcell_TCGA_annotation_na_rm_immune_score$AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("Stage IIIA","Stage IIIB")]<-"III"
xcell_TCGA_annotation_na_rm_immune_score$pStage[xcell_TCGA_annotation_na_rm_immune_score$AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("Stage IV")]<-"IIII"

xcell_TCGA_annotation_na_rm_immune_score[,c("OS_STATUS","os_status")]

# multivariable Cox model including Th2,Th1,CD8,stage,age and sex
multi_os<-coxph(Surv(os_months,os_status) ~Th2+Th1+CD8+pStage+age+SEX,data=xcell_TCGA_annotation_na_rm_immune_score,na.action=na.exclude)
as.matrix(coef(summary(multi_os)))
OR<-exp(cbind(OR = coef(multi), confint(multi)))
OR<-as.data.frame(OR)

# forest plot for multivariable Cox model
format(OR$`97.5 %` ,scientific = FALSE)
mean  = c(NA, 22.9, 14.6, 0.002,1.92, 2.82,2.72,1.007,0.924)
lower = c(NA, 1.34, 0.32, 0.000008, 1.27,1.86,1.489,0.99,0.67)
upper = c(NA, 391.8, 657.6, 0.49, 2.90,4.27,5.00,1.02,1.28)

tabletext<-cbind(
  c("Features","Th2", "Th1", "CD8", "StageII", "StageIII","StageIIII","age","SEX"),
  c("P value", "0.03", "0.17", "0.03", "0.002","8.77e-07","0.001","0.40","0.64"),
  c("Hazard Ratio", "22.9", "14.6", "0.002","1.92", "2.82","2.73","1.01","0.92")
)

pdf("TCGA_survival_forest_2.pdf",  width=8, height=4)
forestplot(tabletext, new_page = TRUE,
           mean, lower, upper,
           zero = 1,
           clip=c(0, 25),
           xlog= FALSE, line.margin=0.5,boxsize=0.25,
           col=fpColors(box="royalblue",line="darkblue"))
dev.off()

#compare TMB in primary tumors with distant metastasis and with local metastasis in TCGA dataset
tcga_mutation_burden<-read.csv("tcga_mutation_burden.txt",sep="\t")
xcell_TCGA_annotation_na_rm_immune_score$PATIENT_ID %in% tcga_mutation_burden$Patient_ID
tcga_mutation_burden<-tcga_mutation_burden[grepl("-01",tcga_mutation_burden$Tumor_Sample_ID),]
xcell_TCGA_annotation_na_rm_immune_score_TMB<-merge(xcell_TCGA_annotation_na_rm_immune_score,tcga_mutation_burden,by.x="PATIENT_ID",by.y="Patient_ID")

T1<-xcell_TCGA_annotation_na_rm_immune_score_TMB[xcell_TCGA_annotation_na_rm_immune_score_TMB$AJCC_TUMOR_PATHOLOGIC_PT %in% c("T1","T1a","T1b"),]
T2<-xcell_TCGA_annotation_na_rm_immune_score_TMB[xcell_TCGA_annotation_na_rm_immune_score_TMB$AJCC_TUMOR_PATHOLOGIC_PT %in% c("T2","T2a","T2b"),]
T3<-xcell_TCGA_annotation_na_rm_immune_score_TMB[xcell_TCGA_annotation_na_rm_immune_score_TMB$AJCC_TUMOR_PATHOLOGIC_PT %in% c("T3","T4"),]

tmb_tcga_box<-box_plot_tcga_validation(feature_table=xcell_TCGA_annotation_na_rm_immune_score_TMB,include_type="Non.silent.per.Mb",category_variable="Stage")

pdf("tmb_tcga_box.pdf",width=6,height=4)
tmb_tcga_box$p1
dev.off()
jonckheere.test(xcell_TCGA_annotation_na_rm_immune_score_TMB$Non.silent.per.Mb, xcell_TCGA_annotation_na_rm_immune_score_TMB$TNM_digit,alternative="increasing",nperm=500)



#########################MSK data analysis ####################################
#compare TMB in primary tumors with distant metastasis and with local metastasis

MSK_clinical<-read.csv("MSK_clinical.txt",sep="\t")
dim(MSK_clinical)
MSK_clinical_primary<-MSK_clinical[MSK_clinical$Sample.Type=="Primary",]
dim(MSK_clinical_primary)
table(MSK_clinical_primary$Metastatic.Site)
MSK_clinical_two_sample<-MSK_clinical[MSK_clinical$Number.of.Samples.Per.Patient>1,]
dim(MSK_clinical_two_sample)
MSK_clinical_two_sample$Sample.Type<-as.character(MSK_clinical_two_sample$Sample.Type)
MSK_clinical_two_sample$Metastatic.Site<-as.character(MSK_clinical_two_sample$Metastatic.Site)
pleura_id<-unique(MSK_clinical_two_sample$Patient.ID[MSK_clinical_two_sample$Metastatic.Site %in% c("Pleura","Pleural Fluid","Chest Wall","Pericardium")])
lymph_id<-unique(MSK_clinical_two_sample$Patient.ID[MSK_clinical_two_sample$Metastatic.Site %in% c("Lymph Node")])
distant_id<-unique(MSK_clinical_two_sample$Patient.ID[MSK_clinical_two_sample$Metastatic.Site %in% c("Brain","Bone","Adrenal Gland","Liver","Abdomen","Breast","Duodenum","Skin","Soft Tissue","Subcutaneous Tissue","Rectum")])
median(MSK_clinical_two_sample$TMB..nonsynonymous.[MSK_clinical_two_sample$Sample.Type=="Primary" & MSK_clinical_two_sample$Patient.ID %in% distant_id])
wilcox.test(MSK_clinical_two_sample$TMB..nonsynonymous.[MSK_clinical_two_sample$Sample.Type=="Primary" & MSK_clinical_two_sample$Patient.ID %in% pleura_id],MSK_clinical_two_sample$TMB..nonsynonymous.[MSK_clinical_two_sample$Sample.Type=="Primary" & MSK_clinical_two_sample$Patient.ID %in% distant_id])
pleura_msk<-MSK_clinical_two_sample[MSK_clinical_two_sample$Sample.Type=="Primary" & MSK_clinical_two_sample$Patient.ID %in% pleura_id,]
pleura_msk$metastasis<-"intrathoracic"
distant_msk<-MSK_clinical_two_sample[MSK_clinical_two_sample$Sample.Type=="Primary" & MSK_clinical_two_sample$Patient.ID %in% distant_id,]
distant_msk$metastasis<-"ldistant"
pleura_distant_msk<-rbind(pleura_msk,distant_msk)

# primary tumors with distant metastasis tend to have a higher TMB than those with local metastasis
pleura_distant_msk_box<-box_plot_tcga_validation(feature_table=pleura_distant_msk,include_type="TMB..nonsynonymous.",category_variable="metastasis")

pdf("pleura_distant_msk_boxplot.pdf",width=4,heigh=4)
pleura_distant_msk_box$p1
dev.off()
wilcox.test(pleura_distant_msk$TMB..nonsynonymous.[pleura_distant_msk$metastasis=="intrathoracic"],pleura_distant_msk$TMB..nonsynonymous.[pleura_distant_msk$metastasis=="ldistant"])



