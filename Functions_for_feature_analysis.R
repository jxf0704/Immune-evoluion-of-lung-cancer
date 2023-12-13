library(gmodels)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(dplyr)
library(survival)
library(survminer)
library(survMisc)
library(clinfun)
library(DescTools)
library(ggplot2)
library(easyGgplot2)
library(ggpubr)
library(GSA)
library(jcolors)
library(plyr)
library(pROC)
library(Hmisc)

colnames(target_feature_table)

clinical_features<- c("site","tumor.size.cm.","PD_L1.TPS","PD_L1.CPS","CD4","CD8","organ","metastasis","TPS_status","CPS_status","CD4_status","CD8_status", "Smoking.pack.year.","PFS_status_new", "PFS_new" ,"OS_status_new","OS_new","Smoking","CIS","total_TMB","sense_TMB","impact_TMB","ploidy","polyploidy","oncogene_num", "TSG_num","cancer_gene_num","oncogenic_num","oncogenic_2_num","sig_exome_A","sig_exome_B", "sig_exome_C","sig_exome_D","sig_SBS1","sig_SBS2","sig_SBS4","sig_SBS5","sig_SBS7a","sig_SBS7b","sig_SBS10b","sig_SBS13",  "sig_SBS15","sig_SBS23","Genetic_distance","local","histology","TPS_and_CPS","TPS_or_CPS","CD4_and_CD8","CD4_or_CD8","subtype","no_reseg_no_BAE_FCS", "no_reseg_no_BAE_BCS","no_reseg_no_BAE_GCS", "no_reseg_LOH_FCS","no_reseg_LOH_BCS","no_reseg_LOH_GCS","clonal_count subclonal_count","subclonal_percentage")

numeric_feature<-c("tumor.size.cm.","PD_L1.TPS","PD_L1.CPS","CD4","CD8", "Smoking.pack.year.","CIS","total_TMB","sense_TMB","impact_TMB","ploidy","oncogene_num", "TSG_num","cancer_gene_num","oncogenic_num","oncogenic_2_num","sig_exome_A","sig_exome_B", "sig_exome_C","sig_exome_D","sig_SBS1","sig_SBS2","sig_SBS4","sig_SBS5","sig_SBS7a","sig_SBS7b","sig_SBS10b","sig_SBS13",  "sig_SBS15","sig_SBS23","no_reseg_no_BAE_FCS", "no_reseg_no_BAE_BCS","no_reseg_no_BAE_GCS", "no_reseg_LOH_FCS","no_reseg_LOH_BCS","no_reseg_LOH_GCS","clonal_count","subclonal_count","subclonal_percentage")

category_features<- c("site","organ","metastasis","TPS_status","CPS_status","CD4_status","CD8_status", "PFS_status_new", "OS_status_new","Smoking","polyploidy","local","histology","TPS_and_CPS","TPS_or_CPS","CD4_and_CD8","CD4_or_CD8","subtype")

IHC_mumeric_feature_matrix<-target_feature_table[,numeric_feature]

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  result<-data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
  return(result)
}


IHC_mumeric_feature_correlation<-rcorr(as.matrix(IHC_mumeric_feature_matrix),type="spearman")
write.csv(flattenCorrMatrix(IHC_mumeric_feature_correlation$r, IHC_mumeric_feature_correlation$P),"numeric_feature_correlation.csv")

target_feature_table$metastasis


#' linear or logistic regression among various features in a feature table
#'
#' @param feature_table, a data frame containing categorical and numerical features such as sex, age, CD8, CD4, primary, metastasis, etc
#' @param y_field,dependent variable, a vector containing continuous features (linear)or binary categorical features (logistic, for example, CD8 positive vs. CD8 negative)
#' @param x_field_fixed, controlled features in a multivariable regression, for example,sex,age,histology,etc.
#' @param x_field_supplement, independent variables, i.e. features to be investigated
#' @param family, linear for linear regrssion, logistic for logistic regression
#'
#' @return regression result
#' @export
#'
#' @examples, numeric_category_correlation_univariate<-feature_regression(y_field=numeric_feature,x_field_fixed = c(),x_field_supplement=category_features)

feature_regression<-function(feature_table=target_feature_table,y_field=c("CD8","CD4"),x_field_fixed=c("metastasis","histology","sex","age"),x_field_supplement,family="linear"){
  feature_result_total<-list()
  for (y_feature in y_field){
    for (x_feature in x_field_supplement){
      all_table<-feature_table[,c(y_feature,x_field_fixed,x_feature)]
      colnames(all_table)[[1]]<-"response"
      if(family=="linear"){
        lm_gene<-lm(response~.,data=all_table)
        lm_summary<-summary(lm_gene)
        feature_result<-lm_summary$coefficients
        feature_result<-data.frame(feature_result)
        feature_result<-feature_result[!grepl("Intercept",rownames(feature_result)),]
        feature_result$feature<-as.character(rownames(feature_result))
        feature_result$y_feature<-y_feature
        feature_result$x_feature<-x_feature
        feature_result$family<-"linear"
        feature_p_value<-feature_result[nrow(gene_result),]
        feature_result_total[[paste(y_feature,x_feature,sep="_")]]<-feature_result
        }
      if(family=="logistic"){
        glm_gene<-glm(response~.,data=all_table,family="binomial")
        glm_summary<-summary(glm_gene)
        feature_result<-glm_summary$coefficients
        feature_result<-data.frame(feature_result)
        feature_result<-feature_result[!grepl("Intercept",rownames(feature_result)),]
        feature_result$feature<-as.character(rownames(feature_result))
        feature_result$y_feature<-y_feature
        feature_result$x_feature<-x_feature
        feature_result$family<-"logistic"
        feature_p_value<-feature_result[nrow(gene_result),]
        feature_result_total[[paste(y_feature,x_feature,sep="_")]]<-feature_result
        }
    }
  }
  feature_result_total<-do.call(rbind,feature_result_total)
  return(feature_result_total)
}


#' Correlation between clinical features and gene-level alterations (mutation or CNV)
#'
#' @param feature_table
#' @param gene_binary_table, a dataframe containing binary status of gene alteration,yes or no, as well as a column named Sample_ID
#' @param y_field,binary categorical features or continuous features
#' @param x_field,controlled features in multivariable regression, for example, sex,age,Smoking,etc.
#' @param family, linear or logistic
#' @param gene_num_cutoff, cutoff of the frequency of gene alteration in the cohort, to exclude genes with very low frequency
#'
#' @return
#' @export
#'
#' @examples, gene_cnv_cd8_regression<-gene_regression(gene_binary_table=cnv_level2_matrix$cnv_binary_table,y_field="CD8",family="linear",gene_num_cutoff=5)

gene_regression<-function(feature_table=target_feature_table,gene_binary_table=cancer_gene_cnv_binary_matrix,y_field="CD8",x_field=c("metastasis","histology","sex","age","Smoking"),family="linear",gene_num_cutoff=10){
  gene_name<-colnames(gene_binary_table)[-1]
  all_table<-merge(feature_table,gene_binary_table,by="Sample_ID")
  gene_result_total<-list()
  gene_pvalue_total<-list()
  for (gene in gene_name){
    if (length(x_field)==0){
      gene_table<-all_table[,c(y_field,gene)]
    }
    else{
      gene_table<-all_table[,c(y_field,x_field,gene)]
    }
    colnames(gene_table)[[1]]<-"response"
    mutation_freq<-sum(gene_table[,gene])
    if(family=="linear" & mutation_freq>=gene_num_cutoff){
      lm_gene<-lm(response~.,data=gene_table)
      lm_summary<-summary(lm_gene)
      gene_result<-lm_summary$coefficients
      gene_result<-gene_result[!grepl("Intercept",rownames(gene_result)),]
      gene_result<-as.data.frame(gene_result)
      gene_result$feature<-as.character(rownames(gene_result))
      gene_result$gene<-gene
      gene_result$family<-"linear"
      gene_p_value<-gene_result[nrow(gene_result),]
      gene_p_value$mutation_freq<-mutation_freq
      gene_result_total[[gene]]<-gene_result
      gene_pvalue_total[[gene]]<-gene_p_value
    }
    if(family=="logistic"& mutation_freq>=gene_num_cutoff){
      glm_gene<-glm(response~.,data=gene_table,family="binomial")
      glm_summary<-summary(glm_gene)
      gene_result<-glm_summary$coefficients
      gene_result<-gene_result[!grepl("Intercept",rownames(gene_result)),]
      gene_result<-as.data.frame(gene_result)
      gene_result$feature<-as.character(rownames(gene_result))
      gene_result$gene<-gene
      gene_result$family<-"logistic"
      gene_p_value<-gene_result[nrow(gene_result),]
      gene_p_value$mutation_freq<-mutation_freq
      gene_result_total[[gene]]<-gene_result
      gene_pvalue_total[[gene]]<-gene_p_value
    }
}
  gene_result_total<-do.call(rbind,gene_result_total)
  gene_pvalue_total<-do.call(rbind,gene_pvalue_total)
  colnames(gene_pvalue_total)[[4]]<-"p_value"
  gene_pvalue_total$q_value<-p.adjust(gene_pvalue_total$p_value, method = "BH")
  gene_pvalue_total<-gene_pvalue_total[order(gene_pvalue_total$p_value,decreasing=F),]
  gene_pvalue_total_cutoff<-gene_pvalue_total[gene_pvalue_total$mutation_freq>=gene_num_cutoff,]
  gene_pvalue_total_cutoff$q_value<-p.adjust(gene_pvalue_total_cutoff$p_value, method = "BH")
  final_result_list<-list()
  final_result_list[["gene_result_total"]]<-gene_result_total
  final_result_list[["gene_pvalue_total"]]<-gene_pvalue_total
  final_result_list[["gene_pvalue_total_cutoff"]]<-gene_pvalue_total_cutoff
  return(final_result_list)
}


#' Performing Wilcoxon or Fisher test for phenotype-genetype correlation analysis
#'
#' @param feature_table
#' @param gene_binary_table
#' @param y_field
#' @param test, "wilcox" for Wilcoxon rank sum test if y_field is a continuous feature and "fisher" for Fisher exact test if y_field is a binary categorical feature
#' @param gene_num_cutoff
#'
#' @return
#' @export
#'
#' @examples
gene_wilcox_fisher<-function(feature_table=target_feature_table,gene_binary_table=cancer_gene_cnv_binary_matrix,y_field="CD8",test="wilcox",gene_num_cutoff=5){
  gene_name<-colnames(gene_binary_table)[-1]
  all_table<-merge(feature_table,gene_binary_table,by="Sample_ID")
  #gene_result_total<-list()
  gene_test_total<-list()
  total_gene_table<-list()
  for (gene in gene_name){
    gene_table<-all_table[,c(y_field,gene)]
    colnames(gene_table)[[1]]<-"response"
    mutation_freq<-sum(gene_table[,gene])
    if (mutation_freq>=gene_num_cutoff){
      if(test=="wilcox"){
        mutation_median<-median(gene_table$response[gene_table[,gene]==1])
        non_mutation_median<-median(gene_table$response[gene_table[,gene]==0])
        gene_p<-wilcox.test(gene_table$response[gene_table[,gene]==1],gene_table$response[gene_table[,gene]==0])$p.value
        gene_result<-data.frame(gene=gene,mutation_freq=mutation_freq,mutation_median=mutation_median,non_mutation_median=non_mutation_median,pvalue=gene_p,response=y_field)
        gene_test_total[[gene]]<-gene_result
        colnames(gene_table)[[1]]<-y_field
        total_gene_table[[gene]]<-gene_table
      }
      if(test=="fisher"){
        gene_table_fisher<-table(gene_table$response,gene_table[,gene])
        response_1_freq<-sum(gene_table[,gene][gene_table$response==1])/sum(gene_table$response==1)
        response_0_freq<-sum(gene_table[,gene][gene_table$response==0])/sum(gene_table$response==0)
        gene_p<-fisher.test(gene_table_fisher)$p.value
        gene_result<-data.frame(gene=gene,mutation_freq=mutation_freq,response_1_freq=response_1_freq,response_0_freq=response_0_freq,pvalue=gene_p,response=y_field)
        gene_test_total[[gene]]<-gene_result
        colnames(gene_table)[[1]]<-y_field
        total_gene_table[[gene]]<-gene_table
      }
    }
  }
  gene_test_total<-do.call(rbind,gene_test_total)
  gene_test_total$q_value<-p.adjust(gene_test_total$pvalue, method = "BH")
  gene_test_total<-gene_test_total[order(gene_test_total$pvalue,decreasing=F),]
  final_result_list<-list()
  final_result_list[["gene_test_total"]]<-gene_test_total
  final_result_list[["total_gene_table"]]<-total_gene_table
  return(final_result_list)
}


# generate box plot and pairwise Wilcoxon test among groups
box_plot_function_single_group<-function(feature_table=target_feature_table,y_feature="no_reseg_no_BAE_BCS",x_feature="local",shape_feature="histology",outlie_quantile=0.95){
  target_table<-feature_table[,c(y_feature,x_feature,shape_feature)]
  target_table[,x_feature]<-as.factor(target_table[,x_feature])
  pair_list<-list()
  x_feature_type<-as.character(unique(target_table[,x_feature]))
  for(i in 1:(length(x_feature_type)-1)){
    for (j in (i+1):length(x_feature_type)){
      pair_list[[paste(x_feature_type[[i]],x_feature_type[[j]],sep="_")]]<-c(x_feature_type[[i]],x_feature_type[[j]])
    }
  }
  p1=ggboxplot(target_table, x =x_feature, y = y_feature,outlier.shape = NA,
               fill = x_feature, palette = "Dark2",add = "jitter",legend="right",add.params = list(shape = shape_feature))
  p1<-eval(parse(text=paste("p1+stat_compare_means(comparisons=pair_list)",sep = "")))
  outlie_value<-quantile(target_table[,y_feature],outlie_quantile)[[1]]
  p2=ggboxplot(target_table, x =x_feature, y = y_feature,outlier.shape = NA,
               fill = x_feature, palette = "Dark2",add = "jitter",legend="right",add.params = list(shape = shape_feature))+ylim(c(0,outlie_value))
  p2<-eval(parse(text=paste("p2+stat_compare_means(comparisons=pair_list)",sep = "")))
  result<-list()
  result[["outlie_include"]]<-p1
  result[["outlie_remove"]]<-p2
  return(result)
}

# generate multi-group bar plot with facet_wrap feature, and output summarize table to a txt file
bar_plot_multi_group<-function(feature_table=two_cytoband_IHC,fill_feature="chr11_q14",x_feature="TPS_status",facet_feature="metastasis",shape_feature="histology",outlie_quantile=0.95){
  fill_table<-feature_table[,c(fill_feature,x_feature,facet_feature)]
  fill_table_divid<-split(fill_table,fill_table[,facet_feature])
  cross_all_table<-list()
  for (i in 1:length(fill_table_divid)){
    t<-fill_table_divid[[i]]
    cross_table<-table(t[,fill_feature],t[,x_feature])
    cross_Table<-as.data.frame(cross_table)
    colnames(cross_Table)<-c("fill_feature","x_feature","Freq")
    cross_Table[,"facet_feature"]<-names(fill_table_divid)[[i]]
    cross_all_table[[i]]<-cross_Table
  }
  cross_Table<-do.call(rbind,cross_all_table)
  lapply(fill_table_divid,function(x) capture.output(CrossTable(x[,fill_feature], x[,x_feature], digits=2, prop.r=TRUE, prop.c=TRUE, prop.t=FALSE, prop.chisq=FALSE, fisher=TRUE,sresid=FALSE, format=c("SPSS"), dnn = c(fill_feature,x_feature)), file = paste(fill_feature,x_feature,unique(x[,facet_feature])[[1]],".txt",sep="_")))
  p=ggplot(cross_Table, aes(fill=fill_feature, y=Freq, x=x_feature)) + geom_bar(position="fill", stat="identity")+scale_fill_brewer(palette="Dark2")+theme_bw()+facet_wrap(~facet_feature)
  return(p)
}

# generate single-group bar plot and output summarize table to a txt file
bar_plot_single_group<-function(feature_table=target_feature_table,fill_feature,x_feature="metastasis",color_feature,shape_feature="histology",outlie_quantile=0.95){
  fill_table<-feature_table[,c(fill_feature,x_feature)]
  cross_table<-table(fill_table[,fill_feature],fill_table[,x_feature])
  summary_table<-CrossTable(fill_table[,fill_feature], fill_table[,x_feature], digits=2, prop.r=TRUE, prop.c=TRUE, prop.t=FALSE, prop.chisq=FALSE, fisher=TRUE,sresid=FALSE, format=c("SPSS"), dnn = c(fill_feature,x_feature))
  cross_Table<-as.data.frame(cross_table)
  colnames(cross_Table)<-c("fill_feature","x_feature","Freq")
  capture.output(CrossTable(fill_table[,fill_feature], fill_table[,x_feature], digits=2, prop.r=TRUE, prop.c=TRUE, prop.t=FALSE, prop.chisq=FALSE, fisher=TRUE,sresid=FALSE, format=c("SPSS"), dnn = c(fill_feature,x_feature)), file = paste(fill_feature,x_feature,".txt",sep="_"))
  p=ggplot(cross_Table, aes(fill=fill_feature, y=Freq, x=x_feature)) +
    geom_bar(position="fill", stat="identity")+scale_fill_brewer(palette="Dark2")+theme_bw()
  return(p)

}


bar_plot_single_group_gene<-function(feature_table=target_feature_table,fill_feature,x_feature="metastasis",color_feature,shape_feature="histology",outlie_quantile=0.95,gene=c(),gene_table=cancer_gene_mutation_matrix_73){
  gene_status<-gene_table[,c("Sample_ID",gene)]
  gene_status<-merge(feature_table,gene_status,by="Sample_ID",all.x=T)
  fill_table<-gene_status[,c(fill_feature,x_feature)]
  cross_table<-table(fill_table[,fill_feature],fill_table[,x_feature])
  summary_table<-CrossTable(fill_table[,fill_feature], fill_table[,x_feature], digits=2, prop.r=TRUE, prop.c=TRUE, prop.t=FALSE, prop.chisq=FALSE, fisher=TRUE,sresid=FALSE, format=c("SPSS"), dnn = c(fill_feature,x_feature))
  cross_Table<-as.data.frame(cross_table)
  colnames(cross_Table)<-c("fill_feature","x_feature","Freq")
  capture.output(CrossTable(fill_table[,fill_feature], fill_table[,x_feature], digits=2, prop.r=TRUE, prop.c=TRUE, prop.t=FALSE, prop.chisq=FALSE, fisher=TRUE,sresid=FALSE, format=c("SPSS"), dnn = c(fill_feature,x_feature)), file = paste(fill_feature,x_feature,".txt",sep="_"))
  p=ggplot(cross_Table, aes(fill=fill_feature, y=Freq, x=x_feature)) +
    geom_bar(position="fill", stat="identity")+scale_fill_brewer(palette="Dark2")+theme_bw()
  return(p)
}


#' A function to screen genes whose expression/methylation is significantly affected by certain gene mutation by using limma
#'
#' @param rnaseqdata,expression matrix, for example, rna-seq or methylation matrix
#' @param property_data, a data frame including patient information for example, stage,grade,or gene mutation
#' @param study_field, a binary feature, for example, gene mutation status(0 or 1), sex (0 or 1), etc.
#' @param exclude_field
#'
#' @return
#' @export
#'
#' @examples
genemuation_gene_effect<-function(rnaseqdata,property_data,study_field,exclude_field){
  study_id<-which(colnames(property_data)==study_field)
  tumor_id<-which(colnames(property_data)=="tumor")
  if (!is.na(study_id)){
    study_field_property<-property_data[,study_id]
  }
  exclude_field_properties<-list()
  for (i in 1:length(exclude_field)){
    exclude_id<-which(colnames(property_data)==exclude_field[[i]])
    if (!is.na(exclude_id)){
      exclude_field_property<-property_data[,exclude_id]
      exclude_field_properties[[i]]<-exclude_field_property
    }
  }
  exclude_field_properties[[length(exclude_field)+1]]<-as.character(property_data[,tumor_id])
  exclude_field_properties[[length(exclude_field)+2]]<-study_field_property
  property<-do.call(cbind,exclude_field_properties)
  property_colnames<-c(exclude_field,"tumor",study_field)
  colnames(property)<-property_colnames
  if (nrow(property)==144 & ncol(property)==(length(exclude_field)+2)){
    remove_id<-c()
    for (i in 1:nrow(property)){
      if ("1" %in% property[i,1:(ncol(property)-2)] | NA %in% property[i,1:(ncol(property)-2)] | is.na (property[i,ncol(property)])){
        remove_id<-c(remove_id,i)
      }
    }
  }
  if (length(remove_id)>0){
    property<-as.data.frame(property[-remove_id,])}
  else {
    property<-as.data.frame(property)
  }
  if (table(property[,ncol(property)])[[2]]>=3){
    if(length(remove_id)>0){
      property_data<-rnaseqdata[,-remove_id]
    }
    else{
      property_data<-rnaseqdata
    }
    TS <- paste(property$tumor, property[,ncol(property)], sep=".")
    TS <- factor(TS, levels=c("N.0","T.0","N.1","T.1"))
    design <- model.matrix(~0+TS)
    colnames(design) <- levels(TS)
    fit <- lmFit(property_data, design)
    cont.matrix <- makeContrasts(wt_effect=T.0-N.0,mut_effect=T.1-N.1,interaction=(T.1-N.1)-(T.0-N.0),levels=design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2)
    property_interaction<-topTable(fit2,coef=3,number=nrow(property_data))
    property_interaction$genename<-property_interaction$ID
    significant.record<-property_interaction[property_interaction$adj.P.Val<=0.05,]
    significant.record.pvalue0.5<-property_interaction[property_interaction$P.Value<=0.05,]
    rank<-data.frame(genename=property_interaction$genename,score=property_interaction$logFC)
    return(list(property,property_interaction,significant.record,significant.record.pvalue0.5,rank))
    }
}



#' A function to make summary statistic
#'
#' @param conditionvector, a vector, the elements are the condition name (binary variable or continuous variable) for example,age, gene mutation sex,p53 mutation, etc)
#' @param groupvariable, a binary variable, for example, treatment and control, mutation or non_mutation
#' @param conditiontable, a data frame including the above features
#' @param count_testtype, fisher or chi
#'
#' @return
#' @export
#'
#' @examples
condition_summary<-function(conditionvector,groupvariable,conditiontable,count_testtype){
  count_results<-list()
  continuous_results<-list()
  j=1
  s=1
  final<-list()
  for (i in 1:length(conditionvector)){
    condition_ind<-which(colnames(conditiontable)==conditionvector[[i]])
    group_ind<-which(colnames(conditiontable)==groupvariable)
    testtable<-data.frame(condition=conditiontable[,condition_ind],group=conditiontable[,group_ind])
    testtable<-subset(testtable,!is.na(testtable$condition) & !is.na(testtable$group))
    testtable$condition_name<-conditionvector[[i]]
    conditionvalue<-levels(as.factor(testtable[,1]))
    if (length(conditionvalue)<=2){
      condition_count<-dcast(testtable,group~condition,fun.aggregate=length)
      if (nrow(condition_count)==2 & ncol(condition_count)==3){
        condition_count<-as.matrix(condition_count[,c(2,3)])
        if (count_testtype=="fisher") {
          fisher<-fisher.test(condition_count)
          pvalue<-fisher$p.value
        } else if (count_testtype=="chi"){
          chisq<-chisq.test(condition_count)
          pvalue<-chisq$p.value
        }
        condition_positive<-subset(testtable,condition==1)
        condition_p<-dcast(condition_positive,condition_name~group,fun.aggregate=length)
        if (dim(condition_p)[[2]]==2 & colnames(condition_p)[[2]]=="0"){
          condition_p$"1"=0
        } else if (dim(condition_p)[[2]]==2 & colnames(condition_p)[[2]]=="1"){
          condition_p$"0"=0
          condition_p<-condition_p[,c(1,3,2)]
        }
        condition_t<-dcast(testtable,condition_name~group,fun.aggregate=length)
        percentage<-round(condition_p[,-1]/condition_t[,-1]*100,digits=2)
        total_percentage<-round(sum(condition_p[1,-1])/sum(condition_t[1,-1])*100,digits=2)
        percentage_combine<-paste(condition_p[,-1],"(",percentage,")",sep="")
        percentage_combine<-t(as.data.frame(percentage_combine))
        positive_total<-sum(condition_p[1,-1])
        patient_total<-sum(condition_t[1,-1])
        total_percentage_combine<-paste(positive_total,"(",total_percentage,")",sep="")
        result<-cbind(condition_p,percentage,percentage_combine,positive_total,patient_total,total_percentage,total_percentage_combine,condition_t)
        result$pvalue<-pvalue
        count_results[[j]]<-result
        j=j+1
      }
    }
    else if (length(conditionvalue)>2 & (class(testtable$condition)=="numeric" |class(testtable$condition)=="integer" )){
      testtable<-subset(testtable,!is.na(condition) & !is.na(group))
      k<-kruskal.test(condition ~ group, data = testtable)
      pvalue<-k$p.value
      iqrvalue<-function(x){
        percentile<-quantile(x[,1],c(0.25,0.5,0.75),na.rm=TRUE)
        median<-percentile[[2]]
        iqr<-percentile[[3]]-percentile[[1]]
        final<-data.frame(number=length(x[,1]),medianvalue=median,IQRvalue=iqr)
        return(final)
      }
      a<-ddply(testtable,.(group),.fun=iqrvalue)
      colname<-a$group
      a<-a[,-1]
      a<-as.data.frame(t(a))
      colnames(a)<-colname
      totalnumber<-as.data.frame(a[1,])
      medianvalue<-as.data.frame(a[2,])
      iqrvalue<-as.data.frame(a[3,])
      combine<-paste(totalnumber,"(",medianvalue,",", iqrvalue,")",sep="")
      combine<-t(as.data.frame(combine))
      iqrtotal<-quantile(testtable[,1],c(0.25,0.5,0.75),na.rm=TRUE)
      total.median<-iqrtotal[[2]]
      total.iqr<-iqrtotal[[3]]-iqrtotal[[1]]
      totalcount<-sum(totalnumber[1,])
      outcome<-cbind(totalnumber,medianvalue,iqrvalue,combine)
      outcome$totalpatient<-totalcount
      outcome$total_median<-total.median
      outcome$total_iqr<-total.iqr
      outcome$totalcombine<-paste(totalcount,"(",total.median,",",total.iqr,")",sep="")
      outcome$pvalue<-pvalue
      outcome$conditionname<-conditionvector[[i]]
      continuous_results[[s]]<-outcome
      s=s+1
    }
  }
  final[[1]]<-count_results
  final[[2]]<-continuous_results
  return(final)
}


#' A function for univariable and multivariable CoxPH analysis with gene expression level as predictor;according to median value of tested gene, patients are divided into two groups. In multivariable cox analysis, age, stage,grade and expression group are included
#'
#' @param testgene, a gene name
#' @param patientsurvival,a data.frame including patient survival data (css, css_days, or os, os_days)
#'
#' @return
#' @export
#'
#' @examples
gene_expression_coxphmodel<-function(testgene,patientsurvival){
  low50<-which(testgene<=quantile(testgene,c(0.5)))
  high50<-which(testgene>quantile(testgene,c(0.5)))
  testgenerange<-c(rep("0",length(testgene)))
  testgenerange[low50]<-"0"
  testgenerange[high50]<-"1"
  testdataframe<-data.frame(patientsurvival,geneexpression=testgene,expressionrange=testgenerange)
  uni<-coxph(Surv(css_days, css) ~expressionrange,data=testdataframe,na.action=na.exclude)
  multi<-coxph(Surv(css_days, css) ~age_surgery+stage_cat+grade_cat+expressionrange,data=testdataframe,na.action=na.exclude)
  uni.summary<-as.matrix(coef(summary(uni)))
  multi.summary<-as.matrix(coef(summary(multi)))
  uni.co<-uni.summary[,1]
  uni.eco<-uni.summary[,2]
  uni.p<-uni.summary[,5]
  uni.para<-c(uni.co,uni.eco,uni.p)
  multi.co<-multi.summary[,1]
  multi.eco<-multi.summary[,2]
  multi.p<-multi.summary[,5]
  multi.para<-c(multi.co,multi.eco,multi.p)
  genecox<-c(uni.para,multi.para)
  return(genecox)
}


