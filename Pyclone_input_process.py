# -*- coding: utf-8 -*-

# script to prepare files for pyclone analysis
# paired-sample model, for each patient, there are two paired samples,for example, primary and distance; pre-treatment and post-treatment
# input for each patient: MAF file of sample1 (pre),MAF file of sample2(post),FACET segment cnv table, tumor purity table
# output for each patient: Mutation_CNV combine table of sample1,Mutation_CNV combine table of sample2, the mutation_id of the two tables are consistent

import os
import pandas as pd
from itertools import compress

# work directory
work_dir="XXXXX"

patient_dir="XXXXXXXXX"
# a dataframe containing two columns (primary,distance) or any related paired samples (pre-treatment, post-treatment),each cell stores sample name
sample_information=pd.read_csv(os.path.join(patient_dir,"sample_information.csv"))
pre_sample=[item.split("-")[0] for item in sample_information["pre"]]
post_sample=[item.split("-")[0] for item in sample_information["post"]]

silent_list=["5'UTR",'Silent',"3'UTR","3'Flank","5'Flank",'Intron','IGR','RNA'] 

# function used to combine mutation table and cnv table
# after combine, the output combine table should contain the following columns (sample_name,mutation_id,ref_counts,var_counts,normal_cn,minor_cn,major_cn,segment)
def snv_cnv_match(snv_table,cnv_table,major_copy=2):
    snv_table["mutation_id"]=snv_table["mutation_site"]
    snv_table["ref_counts"]=snv_table["t_ref_count"]
    snv_table["var_counts"]=snv_table["t_alt_count"]
    snv_table["normal_cn"]=2
    major_copy_number=[]
    minor_copy_number=[]
    segment=[]
    cnv_tag=[]
    for i in range(snv_table.shape[0]):
        chromosome=int(snv_table["Chromosome"].iloc[i])
        start=int(snv_table["Start_Position"].iloc[i])
        segment_chromosome=cnv_table[cnv_table["chrom"]==chromosome]
        segment_chromosome=segment_chromosome.reset_index()
        findit=False
        for j in range(segment_chromosome.shape[0]):
            seg_start=int(segment_chromosome.loc[j,"start"])
            seg_end=int(segment_chromosome.loc[j,"end"])
            total_cn=segment_chromosome.loc[j,"tcn.em"]
            minor_cn=segment_chromosome.loc[j,"lcn.em"]
            seg=segment_chromosome.loc[j,"seg"]
            if start>=seg_start and start<=seg_end:
                findit=True
                if (not pd.isnull(total_cn)) and (not pd.isnull(minor_cn)):
                    major_cn=int(total_cn)-int(minor_cn)
                    major_copy_number.append(major_cn)
                    minor_copy_number.append(int(minor_cn))
                    segment.append(seg)
                    cnv_tag.append("correct")
                  
                else:
                    major_copy_number.append(total_cn)
                    minor_copy_number.append(0)
                    segment.append(seg)
                    cnv_tag.append("minor_NA")
        if findit==False and major_copy==2:
           major_copy_number.append(2)
           minor_copy_number.append(0)
           segment.append(9999)
           cnv_tag.append("not_cover_major2")
        elif findit==False and major_copy==1: 
           major_copy_number.append(1)
           minor_copy_number.append(1)
           segment.append(9999)
           cnv_tag.append("not_cover_major2")
    snv_table["minor_cn"]=minor_copy_number
    snv_table["major_cn"]=major_copy_number
    snv_table["segment"]=segment
    snv_table["cnv_tag"]=cnv_tag
    if major_copy==2:
        snv_table=snv_table.replace({'major_cn': {0: 2}})
    return(snv_table)   


# function used to genereate pyclone input files for each patient (3 files)
def pyclone_input(output_dir,patient_name,snv_table,cnv_dir,purity_dir,keep_silent=True,silent_list=silent_list,major_copy=2):
    snv_table=snv_table[(snv_table["Chromosome"]!="X") & (snv_table["Chromosome"]!="Y") ]
#gene name + "mutation_id" + "HGVSc", mutation_id: chromosome_startsite, for example, 1_2490459. it is used to label unique mutation site
    snv_table["mutation_site"]=snv_table["Hugo_Symbol"].astype(str)+"_"+snv_table["gene_id"].astype(str)+"_"+snv_table["HGVSc"].astype(str)
    if(keep_silent==False):
        snv_table=snv_table[~snv_table["Variant_Classification"].isin(silent_list)]
    sample_name=snv_table["Tumor_Sample_Barcode"]    
    sample_name=[item.split("-")[0] for item in sample_name]
    snv_table["sample_name"]=sample_name
    unique_sample_name=list(set(sample_name))
    for name in unique_sample_name:
        if name in pre_sample:
            pre_snv_table=snv_table[snv_table["sample_name"]==name]
            pre_snv_table=pre_snv_table.reset_index()
        elif name in post_sample:
            post_snv_table=snv_table[snv_table["sample_name"]==name]
            post_snv_table=post_snv_table.reset_index()
    cnv_files=os.listdir(cnv_dir)
    purity_files=os.listdir(purity_dir)
    for name in unique_sample_name: 
        if name in pre_sample:
            name_index=[name in item for item in cnv_files]
            cnv_file=list(compress(cnv_files, name_index))
            pre_cnv_table=pd.read_csv(os.path.join(cnv_dir,cnv_file[0]))
        elif name in post_sample:
            name_index=[name in item for item in cnv_files]
            cnv_file=list(compress(cnv_files, name_index))
            post_cnv_table=pd.read_csv(os.path.join(cnv_dir,cnv_file[0]))
    pre_input=snv_cnv_match(snv_table=pre_snv_table,cnv_table=pre_cnv_table,major_copy=major_copy) 
    post_input=snv_cnv_match(snv_table=post_snv_table,cnv_table=post_cnv_table,major_copy=major_copy) 
    target_dir=output_dir+patient_name
    pre_filename=patient_name+"_before.tsv"
    pre_input.to_csv(os.path.join(target_dir,pre_filename),sep="\t",index=False)  
    post_filename=patient_name+"_post.tsv"
    post_input.to_csv(os.path.join(target_dir,post_filename),sep="\t",index=False)
    if list(pre_input["mutation_id"])==list(post_input["mutation_id"]):
        print(patient_name+": success")
    for name in unique_sample_name:
        if name in pre_sample:
            name_index=[name in item for item in purity_files]
            purity_file=list(compress(purity_files, name_index))
            pre_purity_table=pd.read_csv(os.path.join(purity_dir,purity_file[0]),sep="\t") 
            pre_purity=float(pre_purity_table.columns[1])
        if name in post_sample:
            name_index=[name in item for item in purity_files]
            purity_file=list(compress(purity_files, name_index))
            post_purity_table=pd.read_csv(os.path.join(purity_dir,purity_file[0]),sep="\t") 
            post_purity=float(post_purity_table.columns[1]) 
    purity_data={"sample":["pre","post"],"purity":[pre_purity,post_purity]}
    purity_table=pd.DataFrame.from_dict(purity_data) 
    purity_name=patient_name+"_"+"purity.tsv"
    purity_table.to_csv(os.path.join(target_dir,purity_name),sep="\t",index=False)  

input_mutation_dir="XXXXXX" 
mutation_files=os.listdir(input_mutation_dir)
patientname=[item.split(".")[0] for item in mutation_files]
patientname=list(set(patientname))
output_dir="XXXXXXX"
cnv_dir="XXXXXX"
purity_dir="XXXXXXXX"


for patient in patientname: 
   p_file=[item for item in mutation_files if patient in item]           
   pcontent=pd.read_csv(os.path.join(input_mutation_dir,p_file[0]),sep="\t",index_col=False) 
   pyclone_input(output_dir=output_dir,patient_name=patient,snv_table=pcontent,cnv_dir=cnv_dir,purity_dir=purity_dir,keep_silent=True,silent_list=silent_list,major_copy=2)     
    
    
    
pyclone_input_dir="E:\\xiaojun\\shiguanai\\pyclone\\pyclone_input_all_mutation"    
