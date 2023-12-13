# -*- coding: utf-8 -*-
"""
caculate the pair-wise genetic distance between primary tumor and its metastases, for example, primary-brain metastasis distance, according to the literature of Lucy R Yates, et.al (PMID: 26099045) 

"""

import pandas as pd
import os
import statsmodels.stats.proportion as ssp
import scipy.stats
import statistics


def binomial_upper(mutation_table):
    vaf=[]
    upper_value=[]
    for i in range(mutation_table.shape[0]):
        n=int(mutation_table.loc[i,"total_depth"])
        x=int(mutation_table.loc[i,"var_counts"])
        if n!=0:
            af=float(x/n)
            upper_limit=ssp.proportion_confint(x, n, alpha=0.1, method='beta')[1]
            vaf.append(af)
            upper_value.append(upper_limit)
        else:
            af=9999
            upper_limit=9999
            vaf.append(af)
            upper_value.append(upper_limit)
    mutation_table["vaf"]=vaf
    mutation_table["upper_vaf"]=upper_value
    return mutation_table

filedir="XXXXXX"
input_file=os.listdir(filedir)
pairsample=pd.read_csv("XXXX\IHC_paired_sample.csv",encoding = 'unicode_escape')
patient_id=set(list(pairsample["patient"]))


out_dir="XXXXXX"

for f in input_file:
    sample_name=f.split(".")[0]
    f_table=pd.read_csv(os.path.join(filedir,f),sep="\t")
    f_table=binomial_upper(f_table)
    f_table.to_csv(os.path.join(out_dir,sample_name+".csv"))
    

input_dir="XXXXXXX" 
input_files=os.listdir(input_dir)

def pair_wise_distance(sample1,sample2):
    sample1=sample1.reset_index(drop=True)
    sample2=sample2.reset_index(drop=True)
    sample1=sample1[sample1["vaf"]!=9999]
    sample2=sample2[sample2["vaf"]!=9999]
    sample1_sample2=pd.merge(sample1,sample2,on="mutation_id")
    sample1_sample2=sample1_sample2[(sample1_sample2["vaf_x"]>=0.01) | (sample1_sample2["vaf_y"]>=0.01)]
    vaf_sample1=sample1_sample2["vaf_x"]
    #vaf_sample1=[item for item in vaf_sample1 if item !=0]
    average_vaf1=sum(vaf_sample1)/len(vaf_sample1)
    vaf_sample2=sample1_sample2["vaf_y"]
    #vaf_sample2=[item for item in vaf_sample2 if item !=0]
    average_vaf2=sum(vaf_sample2)/len(vaf_sample2)
    distance=[]
    sample1_sample2=sample1_sample2.reset_index(drop=True)
    for i in range(len(sample1_sample2)):
        sample1_vaf=sample1_sample2.loc[i,"vaf_x"]
        sample1_upper=sample1_sample2.loc[i,"upper_vaf_x"]
        sample1_vaf=sample1_vaf/average_vaf1
        sample1_upper=sample1_upper/average_vaf1
        
        sample2_vaf=sample1_sample2.loc[i,"vaf_y"]
        sample2_upper=sample1_sample2.loc[i,"upper_vaf_y"]
        sample2_vaf=sample2_vaf/average_vaf2
        sample2_upper=sample2_upper/average_vaf2
        result=[abs(sample1_vaf-sample2_vaf),abs(sample1_vaf-sample2_upper),abs(sample1_upper-sample2_vaf)]
        result=min(result)
        distance.append(result)
    average_distance=sum(distance)/len(distance)   
    return (average_distance)

primary_samples=[]
metastasis_samples=[]
genetic_distance=[]
for p in patient_id:
    p_samples=pairsample.loc[pairsample["patient"]==p,"Sample_ID"]
    p_samples=list(p_samples)
    primary=p_samples[0]
    metastasis=p_samples[1]
    sample1=pd.read_csv(os.path.join(input_dir,primary+".csv"))
    sample2=pd.read_csv(os.path.join(input_dir,metastasis+".csv"))
    dis=pair_wise_distance(sample1,sample2)
    primary_samples.append(primary)
    metastasis_samples.append(metastasis)
    genetic_distance.append(dis)


sample_dis={}
sample_dis["primary"]=primary_samples
sample_dis["metastasis"]=metastasis_samples
sample_dis["distance"]=genetic_distance
distance=pd.DataFrame.from_dict(sample_dis)
distance.to_csv("primary_metastasis_distance.csv")    
