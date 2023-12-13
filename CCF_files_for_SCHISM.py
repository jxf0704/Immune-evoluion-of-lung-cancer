"""
Prepare input CCF files for SCHISM analysis

"""


import os
import argparse
import ConfigParser
from collections import defaultdict
import pandas as pd

parser = argparse.ArgumentParser(description = "extract cellurity info for SCHISM analysis")
parser.add_argument('--pyclone_dir', help = "the result of pyclone in this directory", required = True)
parser.add_argument('--schism_dir_mode_1', help = " the output file for schism:re-clustering by schism", required = True)
parser.add_argument('--schism_dir_mode_2', help = " the output file for schism:using cluster from pyclone", required = True)

argv = parser.parse_args()
pyclone_dir = argv.pyclone_dir.strip()
schism_dir_mode_1 = argv.schism_dir_mode_1.strip()
schism_dir_mode_2 = argv.schism_dir_mode_2.strip()


pyclone_dir="XXXXXX"

patient_names=os.listdir(pyclone_dir)

for patient in patient_names:
    p_dir=os.path.join(pyclone_dir,str(patient),"tables")
    p_file=os.listdir(p_dir)
    loci=os.path.join(p_dir,"loci.tsv")
    cluster=os.path.join(p_dir,"cluster.tsv")
    loci_table=pd.read_csv(loci,sep="\t")
    cluster_table=pd.read_csv(cluster,sep="\t")
    cluster_id=cluster_table["cluster_id"].unique()
    loci_table=loci_table[loci_table["cluster_id"].isin (cluster_id)]
    sampleID=[x.split('.')[0].split("-")[0] for x in loci_table["sample_id"]]
    loci_table["sampleID"]=sampleID
    loci_output=loci_table[["sampleID","mutation_id","cellular_prevalence","cellular_prevalence_std"]]
    loci_output.columns=["sampleID","mutationID","cellularity","sd"]
    cluster_output=loci_table[["mutation_id","cluster_id"]]
    cluster_output.columns=["mutationID","cluster_id"]
    cluster_output=cluster_output.drop_duplicates()
    os.makedirs(os.path.join(schism_dir_mode_1,str(patient)))
    os.makedirs(os.path.join(schism_dir_mode_2,str(patient)))
    loci_output.to_csv(os.path.join(schism_dir_mode_1,str(patient),str(patient)+".clusterEstimates.tsv"),sep="\t",index=False)
    loci_output.to_csv(os.path.join(schism_dir_mode_2,str(patient),str(patient)+".clusterEstimates.tsv"),sep="\t",index=False)
    cluster_output.to_csv(os.path.join(schism_dir_mode_2,str(patient),str(patient)+".mutation-to-cluster.tsv"),sep="\t",index=False)
    

