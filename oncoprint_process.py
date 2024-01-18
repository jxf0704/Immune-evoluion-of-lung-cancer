# -*- coding: utf-8 -*-

import os
import pandas as pd
from itertools import chain
import scipy.stats as stats

os.chdir("XXXXX")
os.listdir(".")

mutation_oncoprint=pd.read_csv("mutation_oncoprint_primary.csv",index_col=False)
gene_pathway=pd.read_csv("oncogenic_pathway.csv",index_col=False)
cnv_oncoprint_primary=pd.read_csv("gene_cnv_oncoprint_primary.csv",index_col=False)
oncoprint_table=mutation_oncoprint
pathway_table=gene_pathway

def  join_mutation(a_list):
    a_list=a_list[a_list.notnull()]
    a_list=list(set(a_list))
    a_list=["SNV" for item in a_list]
    a_list=list(set(a_list))
    join_result=(",").join(a_list)
    return (join_result) 

def  join_cnv(a_list):
    a_list=a_list[a_list.notnull()]
    a_list=list(set(a_list))
    a_list=["SCNA" for item in a_list]
    a_list=list(set(a_list))
    join_result=(",").join(a_list)
    return (join_result)


def gene_pathway_match_2(oncoprint_table, pathway_table,mutation_cnv="mutation",oncoprint_type="oncogenic",keep_oncogenic=10,keep_pathway=5):
    return_result={}
    oncoprint_table=pd.merge(oncoprint_table,pathway_table,on="Gene",how="left")
    oncoprint_table["Pathway"]=oncoprint_table["Pathway"].fillna("Others")
    return_result["original_table"]=oncoprint_table
    if oncoprint_type=="oncogenic":
        keep_gene=list(oncoprint_table[(oncoprint_table["count"]>=keep_oncogenic) & (oncoprint_table["Pathway"] !="Others")]["Gene"])
    
    if oncoprint_type=="total":
        keep_gene=list(oncoprint_table[oncoprint_table["count"]>=keep_oncogenic]["Gene"])
        oncoprint_table_keep_gene=oncoprint_table[oncoprint_table["Gene"].isin(keep_gene)]
    pathway_oncoprint_table=oncoprint_table
    function_dict={}
    if mutation_cnv=="mutation":
        for patient in pathway_oncoprint_table.columns[2:12]:
            function_dict[patient]=join_mutation
    if mutation_cnv=="cnv":
        for patient in pathway_oncoprint_table.columns[2:12]:
            function_dict[patient]=join_cnv        
    pathway_oncoprint_table_grouped = pathway_oncoprint_table.groupby("Pathway")     
    pathway_oncoprint_table_grouped = pathway_oncoprint_table_grouped.aggregate(function_dict)
    pathway_oncoprint_table_grouped=pathway_oncoprint_table_grouped.reset_index()
    count=[]
    if mutation_cnv=="mutation":
        for i in range(pathway_oncoprint_table_grouped.shape[0]):
            num=(pathway_oncoprint_table_grouped.iloc[i,1:11]=="SNV").sum()
            count.append(num)
    if mutation_cnv=="cnv":
        for i in range(pathway_oncoprint_table_grouped.shape[0]):
            num=(pathway_oncoprint_table_grouped.iloc[i,1:11]=="SCNA").sum()
            count.append(num)        
    pathway_oncoprint_table_grouped["count"]= count   
    pathway_oncoprint_table_grouped =pathway_oncoprint_table_grouped.sort_values(by=['count'], ascending=False)
    return_result["original_pathway_table"]=pathway_oncoprint_table_grouped
    keep_pathway_table=pathway_oncoprint_table_grouped[(pathway_oncoprint_table_grouped["count"]>=keep_pathway) & (pathway_oncoprint_table_grouped["Pathway"]!="Others")]
    keep_pathway_table["Pathway"]=keep_pathway_table["Pathway"]+"_p"
    keep_pathway_table=keep_pathway_table.rename(columns={"Pathway": "Gene"})
    oncoprint_table_keep_gene=oncoprint_table_keep_gene.iloc[:,1:13]
    gene_pathway_table=pd.concat([oncoprint_table_keep_gene,keep_pathway_table])
    return_result["gene_pathway_table"]=gene_pathway_table
    return(return_result)


def fisher_test(primary_table,distance_table,at_least_count=4):
    genes=[]
    p_values=[]
    gene_total=[]
    total_counts=[]
    gene_list=list(set(list(primary_table["Gene"])+list(distance_table["Gene"])))
    for gene in gene_list:
        if gene in list(primary_table["Gene"]):
            gene_primary_count=list(primary_table[primary_table["Gene"]==gene].loc[:,"count"])[0]
        else:
            gene_primary_count=0
        if gene in list(distance_table["Gene"]):
            gene_distance_count=list(distance_table[distance_table["Gene"]==gene].loc[:,"count"])[0]
        else:
            gene_distance_count=0
        total_count=gene_primary_count+gene_distance_count
        total_counts.append(total_count)
        gene_total.append(gene)
        if gene_primary_count+gene_distance_count>=at_least_count:
            pvalue = stats.fisher_exact([[10-gene_primary_count, gene_primary_count], [10-gene_distance_count, gene_distance_count]])[1]
            genes.append(gene)
            p_values.append(pvalue)
            
    result_pvalue={"Gene":genes,"p_value":p_values}
    result_pvalue=pd.DataFrame(result_pvalue)
    result_pvalue =result_pvalue.sort_values(by=['p_value'], ascending=True)
    count_total={"Gene":gene_total,"total_count":total_counts}
    count_total=pd.DataFrame(count_total)
    count_total =count_total.sort_values(by=['total_count'], ascending=False)
    return(count_total,result_pvalue)

            
def primary_distance_combine(primary_table,distance_table,keep_top_gene,keep_top_pathway,sig_level,least_count):
    fisher=fisher_test(primary_table=primary_table,distance_table=distance_table,at_least_count=least_count)
    pathway_count=fisher[0][fisher[0]["Gene"].str.contains(pat='_p')]
    gene_count=fisher[0][~fisher[0]["Gene"].str.contains(pat='_p')]
    keep_pathway=list(pathway_count["Gene"])[0:keep_top_pathway]
    keep_gene=list(gene_count["Gene"])[0:keep_top_gene]
    sig_gene=fisher[1][fisher[1]["p_value"]<=sig_level]
    sig_gene=list(sig_gene["Gene"])
    keep_total=keep_pathway+keep_gene+sig_gene
    keep_primary_table=primary_table[primary_table["Gene"].isin(keep_total)]
    keep_distance_table=distance_table[distance_table["Gene"].isin(keep_total)]
    final=pd.merge(keep_primary_table,keep_distance_table,on="Gene",how="outer")
    return(final)

colorcode=pd.read_csv("colorcode.csv")
def oncoprint_color(variations,background_color='#CCCCCC',color=colorcode,start=0,step=15):
    background="alter_fun = list(background = alter_graphic('rect', fill = '{0}'),".format(background_color)
    for variation in variations:
        variation_color="{0} = alter_graphic('rect', fill = col['{1}']),".format(variation,variation)
        background=background+variation_color
        start=start+step
    return(background)  
    
mutation_oncoprint_primary_total=gene_pathway_match_2(oncoprint_table=mutation_oncoprint, pathway_table=gene_pathway,mutation_cnv="mutation",oncoprint_type="oncogenic",keep_oncogenic=0,keep_pathway=0)  
primary_oncoprint_gene=mutation_oncoprint_primary_total["gene_pathway_table"] 
mutation_oncoprint_distance=pd.read_csv("mutation_oncoprint_distance.csv",index_col=False) 
mutation_oncoprint_distance_total=gene_pathway_match_2(oncoprint_table=mutation_oncoprint_distance, pathway_table=gene_pathway,mutation_cnv="mutation",oncoprint_type="oncogenic",keep_oncogenic=0,keep_pathway=0)  
distance_oncoprint_gene=mutation_oncoprint_distance_total["gene_pathway_table"]

cnv_oncoprint_primary=pd.read_csv("gene_cnv_oncoprint_primary.csv",index_col=False)
cnv_oncoprint_distance=pd.read_csv("gene_cnv_oncoprint_distance.csv",index_col=False)  
cnv_oncoprint_primary_total=gene_pathway_match_2(oncoprint_table=cnv_oncoprint_primary, pathway_table=gene_pathway,mutation_cnv="cnv",oncoprint_type="oncogenic",keep_oncogenic=0,keep_pathway=0)  
cnv_oncoprint_distance_total=gene_pathway_match_2(oncoprint_table=cnv_oncoprint_distance, pathway_table=gene_pathway,mutation_cnv="cnv",oncoprint_type="oncogenic",keep_oncogenic=0,keep_pathway=0)  
primary_oncoprint_cnv=cnv_oncoprint_primary_total["gene_pathway_table"]
distance_oncoprint_cnv=cnv_oncoprint_distance_total["gene_pathway_table"]

primary_table=primary_oncoprint_gene
distance_table=distance_oncoprint_gene
at_least_count=4

mutation_combine=primary_distance_combine()
cnv_combine=primary_distance_combine(primary_table=primary_oncoprint_cnv,distance_table=distance_oncoprint_cnv,keep_top_gene=5,keep_top_pathway=5,sig_level=0.05,least_count=4)
cnv_combine["Gene"]="C_"+cnv_combine["Gene"]
mutation_cnv_combine=pd.concat([mutation_combine,cnv_combine],sort=False)


primary_arm_cnv=pd.read_csv("arm_cnv_oncoprint_primary.csv",index_col=False)
distance_arm_cnv=pd.read_csv("arm_cnv_oncoprint_distance.csv",index_col=False)
arm_combine=primary_distance_combine(primary_table=primary_arm_cnv,distance_table=distance_arm_cnv,keep_top_gene=3,keep_top_pathway=0,sig_level=0.05,least_count=4)

mutation_cnv_arm_combine=pd.concat([mutation_cnv_combine,arm_combine],sort=False)
mutation_cnv_arm_combine.to_csv("mutation_cnv_arm_combine.csv",index=False)
mutation_cnv_arm_combine_oncoprint=pd.read_csv("mutation_cnv_arm_combine_oncoprint.csv")

arm_diff=fisher_test(primary_table=primary_arm_cnv,distance_table=distance_arm_cnv,at_least_count=4)
gene_diff=fisher_test(primary_table=primary_oncoprint_gene,distance_table=distance_oncoprint_gene,at_least_count=4)
cnv_diff=fisher_test(primary_table=primary_oncoprint_cnv,distance_table=distance_oncoprint_cnv,at_least_count=4)

primary_distance_diff={}
primary_distance_diff["arm_diff"]=arm_diff
primary_distance_diff["gene_diff"]=gene_diff
primary_distance_diff["cnv_diff"]=cnv_diff

oncoprint_table_total={}
oncoprint_table_total["mutation_table_primary"]=mutation_oncoprint_primary_total
oncoprint_table_total["mutation_table_distance"]=mutation_oncoprint_distance_total
oncoprint_table_total["cnv_table_primary"]=cnv_oncoprint_primary_total
oncoprint_table_total["cnv_table_distance"]=cnv_oncoprint_distance_total

combine_table={}
combine_table["mutation_table_combine"]=mutation_combine
combine_table["cnv_table_combine"]=cnv_combine
combine_table["arm_table_combine"]=arm_combine

import pickle
with open("oncoprint_process", "wb") as f: 
    pickle.dump((primary_distance_diff,oncoprint_table_total,combine_table), f)
    
with open("oncoprint_process", "rb") as f: #load data next time
    diff, indiviual_table,combine_table = pickle.load(f)    




