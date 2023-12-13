# -*- coding: utf-8 -*-

import argparse
import ConfigParser
import os
import pandas as pd

parser = argparse.ArgumentParser(description = 'run pyclone')
parser.add_argument('--pyclone_input_dir',help = "dir of pyclone input files")
parser.add_argument('--pyclone_output_dir',help='dir of pyclone input files')
parser.add_argument('--max_cluster',help='the max cluster count')
parser.add_argument('--min_cluster_size',help='minmum mutations in a cluster')


argv = parser.parse_args()
pyclone_input_dir=argv.pyclone_input_dir.strip()
pyclone_output_dir=argv.pyclone_output_dir.strip()
max_cluster = argv.max_cluster.strip()
min_cluster_size = argv.min_cluster.strip()


pyclone_input_dir="E:\\xiaojun\\shiguanai\\pyclone\\pyclone_input_all_mutation"
patient_name=os.listdir(pyclone_input_dir)
pyclone_output_dir="E:\\xiaojun\\shiguanai\\pyclone\\pyclone_input_all_mutation\\output"

for patient in patient_name:
    sample_files=os.listdir(os.path.join(pyclone_input_dir,patient))
    pre_sample=os.path.join(pyclone_input_dir,patient,sample_files[0])
    post_sample=os.path.join(pyclone_input_dir,patient,sample_files[1])
    purity_table=os.path.join(pyclone_input_dir,patient,sample_files[2])
    pyclone_output_dir_patient=os.path.join(pyclone_output_dir,patient)
    infiles=pre_sample+" "+post_sample
    purity_content=pd.read_csv(purity_table,sep="\t")
    purities=str(purity_content["purity"][0])+" "+str(purity_content["purity"][1])
    cmd = "PyClone run_analysis_pipeline --in_files " + infiles + " --working_dir "+pyclone_output_dir_patient+" --tumour_contents " + purities + " --num_iters 10000 --density pyclone_beta_binomial --seed 10 --burnin 1000 --thin 5 --max_clusters "+ str(max_cluster) + " --min_cluster_size "+str(min_cluster_size)
    print (cmd)
    os.system(cmd)
