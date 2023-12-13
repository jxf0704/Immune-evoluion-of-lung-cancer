# -*- coding: utf-8 -*-
"""
Extract key information from clonal tree and check the quality of clonal tree
"""
import os 
import pandas as pd
import itertools

input_dir="XXXXXXX"
pair_sample=os.listdir(input_dir)

mutation_CCF=[]
cluster_CCF=[]

d="10"
for d in pair_sample:
    path=os.path.join(input_dir,d)
    path=path+"\\tables\\loci.tsv"
    content=pd.read_csv(path,sep="\t")
    content["patient"]=d
    content["type"]="paired"
    mutation_CCF.append(content)
    
for d in pair_sample:
    path=os.path.join(input_dir,d)
    path=path+"\\tables\\cluster.tsv"
    content=pd.read_csv(path,sep="\t")
    content["patient"]=d
    content["type"]="paired"
    cluster_CCF.append(content)    
    
mutation_CCF_combine=pd.concat(mutation_CCF)    
cluster_CCF_combine=pd.concat(cluster_CCF)

sample_id=list(mutation_CCF_combine["sample_id"])
sample_id=[item.split(".")[0] for item in sample_id]
mutation_CCF_combine["Sample_ID"]=sample_id

mutation_id=list(mutation_CCF_combine["mutation_id"])
chromsome=[item.split("--")[0] for item in mutation_id]
start_position=[item.split("--")[1] for item in mutation_id]
gene=[item.split("--")[2] for item in mutation_id]
AA_change=[item.split("--")[3] for item in mutation_id]

mutation_CCF_combine["Chromosome"]=chromsome
mutation_CCF_combine["Start_Position"]=start_position
mutation_CCF_combine["Gene"]=gene
mutation_CCF_combine["AA_change"]=AA_change

mutation_cluster_CCF_combine=pd.merge(mutation_CCF_combine,cluster_CCF_combine,on=["patient","cluster_id","type"])

mutation_sample_CCF=mutation_CCF_combine.drop('patient', 1)
mutation_sample_CCF=mutation_sample_CCF.drop('cluster_id', 1)
mutation_sample_CCF.drop_duplicates(subset=None, keep='first', inplace=True)

mutation_CCF_combine.to_csv("mutation_CCF_combine_pyclone.csv")
mutation_sample_CCF.to_csv("mutation_sample_CCF_pyclone.csv")

mutation_cluster_CCF_combine.to_csv("mutation_cluster_CCF_combine.csv")
cluster_CCF_combine.to_csv("cluster_CCF_combine.csv")

tree_file="XXXXXX"
tree_file=pd.read_csv(tree_file,sep="\t",header=None)    

tree_dataframe=tree_file
tree_file.columns=["parent","child"]
tree_dataframe=tree_file

def clone_order(tree_dataframe):
    a=tree_dataframe["parent"].values.tolist()
    b=tree_dataframe["child"].values.tolist()
    a.extend(b)
    c=set(a)
    parent_dict={}
    for i in c:
        if not (i in list(tree_dataframe["child"])):
            child=list(tree_dataframe["child"][tree_dataframe["parent"]==i])
            parent_dict[i]=child
    #once a clone is only in parent which means it cannot be a child of another parent clone, then after extract its child clones, this parent clone and its child clones are removed from the clone table (because a child clone can only have one parent clone), to generate new_tree_dataframe.
    new_tree_dataframe=tree_dataframe[~tree_dataframe["parent"].isin(parent_dict.keys())]
    return(parent_dict,new_tree_dataframe) 

# repeat clone_order until the all of the records in tree_dataframe are removed
# the function is to get the clone order from a tree file   
def clone_order_total(tree_dataframe):
    total=[]
    while tree_dataframe.shape[0]>0:
        result=clone_order(tree_dataframe)[0]
        total.append(result)
        tree_dataframe=clone_order(tree_dataframe)[1]
    order_list=[0]
    clone_list=[[int(k) for k in total[0]][0]]
    for i in range(len(total)):
        item=total[i]
        childs=set(list(itertools.chain.from_iterable(total[i].values())))
        for item in childs:
            order_list.append(i+1)
            clone_list.append(item)
    result=pd.DataFrame({"Clone":clone_list,"order":order_list})        
    return(result) 



total_file=[]

input_dir="XXXXXXX"


file_paths=[]
file_names=[]
for root,dirname,filename in os.walk(input_dir):
    for name in filename:
        file_names.append(name)
        file_path=os.path.join(root,name)
        file_paths.append(file_path)

total_file=[item.split(".")[0] for item in file_names]
total_file=list(set(total_file))
total_file=[item for item in total_file if item!="Rplots"]

file_paths=set(list(file_paths))       

ccf_final=[] 
clone_count_total=[]
clone_number_total=[]
patients=[]


# for each mutation site in each sample, assign its parent clone, child clone, mutational CCF, clonal CCF, clone order; for each sample, calculate clone number and mutation number in each clone.  
for patient in total_file:
    
        if not patient.startswith("Rplots"):
            ccf_file=[file for file in file_paths if file.endswith("\\"+patient+".result.csv")]
            
            ccf_file=ccf_file[0]
            ccf_file=pd.read_csv(ccf_file)
            
            ccf_file["mutation_id"]=ccf_file["Gene"]+"_"+ccf_file["Chr"].astype(str)+"_"+ccf_file["Position"].astype(str)
            ccf_file.drop(['Chr','Position','Mutation_type','Gene'],axis=1,inplace=True)
            ccf_long=ccf_file.melt(id_vars=['Clone','mutation_id'])
            ccf_long.columns=['Clone','mutation_id', 'Sample_name','CCF_mutation']
            tree_file=[file for file in file_paths if file.endswith("\\"+patient+".GA.consensusTree.new")]
            tree_file=tree_file[0]
           
            tree_file=pd.read_csv(tree_file,sep="\t",header=None)
            tree_file.columns=["parent","child"]
            tree_file_rename=tree_file.rename(columns={'parent':'Clone'})
            ccf_child=ccf_long.merge(tree_file_rename,how="left")
            tree_file_rename=tree_file.rename(columns={'child':'Clone'})
            ccf_parent=ccf_child.merge(tree_file_rename,how="left")
            order_result=clone_order_total(tree_file)
            ccf_parent_order=ccf_parent.merge(order_result,how="left")
            
            ccf_parent_order["Patient"]=patient
            clone_ccf=[file for file in file_paths if file.endswith("\\"+patient+".cluster.cellularity")]
            clone_ccf=clone_ccf[0]
            clone_ccf=pd.read_csv(clone_ccf,sep="\t")
           
            clone_ccf.columns=['Sample_name','Clone', 'CCF_clone','CCF_clone_sd']
            ccf_parent_order=ccf_parent_order.merge(clone_ccf,how="left",on=["Sample_name","Clone"])
           
            ccf_parent_order["tree_id"]=[patient]*ccf_parent_order.shape[0]
            ccf_final.append(ccf_parent_order)
            clone_number=len(set(ccf_parent_order["Clone"]))
            clone_count = ccf_parent_order.groupby('Clone')
            clone_count=clone_count.aggregate({'mutation_id':len})
            clone_count=clone_count.reset_index()
            clone_count.columns=["Clone","count"]
            clone_count["count"]=clone_count["count"]/2
            clone_count["Patient"]=patient
            clone_count_total.append(clone_count)
            clone_number_total.append(clone_number)
            patients.append(patient)
            

ccf_final_all=pd.concat(ccf_final)
clone_count_total_all_mutation=pd.concat(clone_count_total)
clone_number_mutation_cluster=clone_number_total


# extract key information from the above CCF combine table, including parent clone CCF and the sum CCF of its child clone, to check if pigeonhole principle is satisfied. If being violiated, what is the difference of CCF, as well as the contribution of each clone to the violation. In addition, check the clone order of EGFR mutation, because in NSCLC, EGFR mutation is generally founder clone, with top clone order.  
def information_extract(ccf_table):
    clone_ccf=ccf_table[["Clone","child","parent","CCF_clone"]]
    clone_ccf.drop_duplicates(inplace=True)
    parent_clone=clone_ccf[~pd.isna(clone_ccf["child"])]
    parent_clone=list(set(parent_clone["Clone"]))
    parentclone=[]
    parentccf=[]
    childtotalccf=[]
    for clone in parent_clone:
        parent_ccf=list(clone_ccf[clone_ccf["Clone"]==clone]["CCF_clone"])[0]
        child_clone=list(set(clone_ccf[clone_ccf["Clone"]==clone]["child"]))
        child_total_ccf=0
        for child in child_clone:
            child_ccf=list(clone_ccf[clone_ccf["Clone"]==child]["CCF_clone"])[0]
            child_total_ccf=child_total_ccf+child_ccf
        parentclone.append(clone) 
        parentccf.append(parent_ccf)
        childtotalccf.append(child_total_ccf)
    parent_child_ccf={"parent_clone":parentclone,"parent_ccf":parentccf,"child_total_ccf":childtotalccf}  
    parent_child_ccf=pd.DataFrame(parent_child_ccf)
    parent_child_ccf["parent_child"]=parent_child_ccf["parent_ccf"]-parent_child_ccf["child_total_ccf"]
    parent_child_dif=sum(parent_child_ccf["parent_child"])    
    parent_child_dif_count=sum(parent_child_ccf["parent_child"]<0)
    wrong_ccf=list(parent_child_ccf["parent_child"][parent_child_ccf["parent_child"]<0])
    wrong_ccf=[str(item) for item in wrong_ccf]
    wrong_ccf=",".join(wrong_ccf)
    wrong_ccf_percentage=parent_child_dif_count/parent_child_ccf.shape[0]
    gene=ccf_table["mutation_id"]
    gene=[item.split("_")[0] for item in list(gene)]
    if "EGFR" in gene:
        indices = [i for i, x in enumerate(gene) if x == "EGFR"]
        EGFR_order=min(ccf_table.iloc[indices,]["order"])
    else:
        EGFR_order=9999
    result={"parent_child_dif":[parent_child_dif],"parent_child_dif_count":[parent_child_dif_count],"wrong_ccf":wrong_ccf,"wrong_ccf_percentage":[wrong_ccf_percentage],"EGFR_order":[EGFR_order]}    
    result=pd.DataFrame(result)
    return(parent_child_ccf,result)

 
all_patient= list(set(ccf_final_all["Patient"]))  
all_parameter= list(set(ccf_final_all["parameter"])) 

total_sample_clone_ccf=[]
parent_child_total=[]
for patient in all_patient:
    for parameter in all_parameter:
        target_table=ccf_final_result[(ccf_final_result["parameter"]==parameter) & (ccf_final_result["Patient"]==patient)]
        if target_table.shape[0]>=1:
            for sample in list(set(target_table["Sample_name"])):
                sample_table=target_table[target_table["Sample_name"]==sample]
                sample_result=information_extract(sample_table)
                sample_clone_ccf=sample_result[0]
                parent_child=sample_result[1]
                sample_clone_ccf["patient"]=patient
                sample_clone_ccf["parameter"]=parameter
                sample_clone_ccf["Sample_name"]=sample
                total_sample_clone_ccf.append(sample_clone_ccf)
                parent_child["patient"]=patient
                parent_child["parameter"]=parameter
                parent_child["Sample_name"]=sample
                parent_child_total.append(parent_child)
                
parent_child_total_all=pd.concat(parent_child_total)                

