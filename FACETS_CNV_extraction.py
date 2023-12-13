# -*- coding: utf-8 -*-
"""
Extracting copy number variation information from FACETS output
 
"""

import os
import pandas as pd
os.chdir("XXXXXX")

oncogenic_gene=pd.read_csv("total_oncogenic_gene.csv")

import os
import pandas as pd
from itertools import compress
import math
import numpy


# extract segment-level cnv data
input_dir="XXXXXX"
files=os.listdir(input_dir)
samplename=[item.split(".")[0] for item in files]
Sample_ID=[item.split("-")[0] for item in samplename]
CNV_data_combine=[]
for file in files:
    samplename=file.split(".")[0]
    sample_id=file.split("-")[0]
    seg_cnv=pd.read_csv(os.path.join(input_dir,file))
    seg_cnv["sample_name"]=samplename
    seg_cnv["Sample_ID"]=sample_id
    CNV_data_combine.append(seg_cnv)
    
CNV_data_combine=pd.concat(CNV_data_combine)
    
# extract chromosomal instability score (CIS)
input_dir="XXXXXX"
files=os.listdir(input_dir)
samplename=[item.split(".")[0] for item in files]
Sample_ID=[item.split("-")[0] for item in samplename]
cis_combine=[]
for file in files:
    cis_cnv=pd.read_csv(os.path.join(input_dir,file))
    cis_cnv=list(cis_cnv.columns)[0]
    cis_combine.append(cis_cnv)
cis_data=[item.replace("%","") for item in cis_combine]
cis_data=[float(item) for item in cis_data]  

cis_data_frame=[samplename,Sample_ID,cis_data]
cis_df = pd.DataFrame (cis_data_frame).transpose()
cis_df.columns = ['sample_name','Sample_ID','CIS'] 

# extract arm-level CNVs
input_dir="XXXXXX"
files=os.listdir(input_dir)
arm_cnv_combine=[]
for file in files:
    samplename=file.split(".")[0]
    sample_id=file.split("-")[0]
    arm_cnv=pd.read_csv(os.path.join(input_dir,file),"\t",names=["arm", "change", "percentage"])
    arm_cnv["sample_name"]=samplename
    arm_cnv["Sample_ID"]=sample_id
    arm_cnv_combine.append(arm_cnv)
arm_data_combine=pd.concat(arm_cnv_combine)


# extract gene-level CNVs
input_dir="XXXXXX"
files=os.listdir(input_dir)
gene_cnv_combine=[]

for file in files:
    if file.endswith("allgene.list"):
        samplename=file.split(".")[0]
        sample_id=file.split("-")[0]
        gene_cnv=pd.read_csv(os.path.join(input_dir,file),"\t",usecols=[0,1,2,3,4])
        genes=[]
        segment_start=[]
        segment_end=[]
        cnv_status=[]
        level=[]
        overlapping=[]
        chromosome=[]
        for i in range(gene_cnv.shape[0]):
            info=gene_cnv.loc[i,"chr;seg_start;seg_end;cnv_stat;level;overlap_length"]
            gene_name=gene_cnv.loc[i,"Gene"]
            if type(info)==str:
                info_split=info.split(";")
                chr=info_split[0]
                seg_start=info_split[1]
                seg_end=info_split[2]
                cnv_stat=info_split[3]
                cnv_level=info_split[4]
                overlap=info_split[5]
                genes.append(gene_name)
                segment_start.append(seg_start)
                segment_end.append(seg_end)
                cnv_status.append(cnv_stat)
                level.append(cnv_level)
                overlapping.append(overlap)
                chromosome.append(chr)
        info_list=[genes,chromosome,segment_start,segment_end,cnv_status,level,overlapping]   
        df = pd.DataFrame (info_list).transpose()
        df.columns =['Gene','chr','seg_start','seg_end','cnv_stat','cnv_level','overlap']
        gene_cnv_table=pd.merge(gene_cnv, df,on='Gene')
        gene_cnv_table["Sample_ID"]=sample_id
        gene_cnv_table["Sample_name"]=samplename
        gene_cnv_combine.append(gene_cnv_table)
        
gene_data_combine=pd.concat(gene_cnv_combine)
gene_data_combine_cnv=gene_data_combine[(gene_data_combine['cnv_stat']=='Amplification') | (gene_data_combine['cnv_stat']=='Deletion')]
cancer_gene=pd.read_csv("cancer_gene.csv")
cancer_gene_cnv=pd.merge(gene_data_combine_cnv,cancer_gene,on="Gene")
cancer_gene_cnv=cancer_gene_cnv.loc[:,["Gene","cnv_stat","cnv_level","Sample_ID","Sample_name"]]



