# -*- coding: utf-8 -*-

"""
Generate SCHISM YAML files for each patient

"""

import os
import yaml
import copy

yaml_path="XXXXXX"

with open(yaml_path, 'r') as yaml_file:
    existing_data = yaml.safe_load(yaml_file)
    existing_data_AP=copy.copy(existing_data)
    del existing_data_AP['clustering_method_B']
    del existing_data_AP['clustering_method_C']
    existing_data_AP['clustering_method']=existing_data_AP.pop("clustering_method_A")
    existing_data_DBSCAN=copy.copy(existing_data)
    del existing_data_DBSCAN['clustering_method_A']
    del existing_data_DBSCAN['clustering_method_C']
    existing_data_DBSCAN['clustering_method']=existing_data_DBSCAN.pop("clustering_method_B")
    existing_data_kmeans=copy.copy(existing_data)
    del existing_data_kmeans['clustering_method_A']
    del existing_data_kmeans['clustering_method_B']
    existing_data_kmeans['clustering_method']=existing_data_kmeans.pop("clustering_method_C")
    
mode_1_path="C:\\mode_1"    
patient_name=os.listdir(mode_1_path)
# only update working_dir, mutation_cellularity_input and output_prefix; others are default values
for patient in patient_name:
    work_dir=os.path.join(mode_1_path,str(patient))
    ccf_file=str(patient)+".clusterEstimates.tsv"
    existing_data_AP['working_dir'] = work_dir
    existing_data_AP['mutation_cellularity_input'] = ccf_file
    existing_data_AP['output_prefix']=str(patient)
    existing_data_DBSCAN['working_dir'] = work_dir
    existing_data_DBSCAN['mutation_cellularity_input'] = ccf_file
    existing_data_DBSCAN['output_prefix']=str(patient)
    existing_data_kmeans['working_dir'] = work_dir
    existing_data_kmeans['mutation_cellularity_input'] = ccf_file
    existing_data_kmeans['output_prefix']=str(patient)
    output_AP=os.path.join(work_dir,str(patient)+".AP.yaml")
    output_DBSCAN=os.path.join(work_dir,str(patient)+".DBSCAN.yaml")
    output_kmeans=os.path.join(work_dir,str(patient)+".kmeans.yaml")
    with open(output_AP, 'w') as yaml_file:
      yaml.dump(existing_data_AP, yaml_file, default_flow_style=False)
    with open(output_DBSCAN, 'w') as yaml_file:
      yaml.dump(existing_data_DBSCAN, yaml_file, default_flow_style=False)
    with open(output_kmeans, 'w') as yaml_file:
      yaml.dump(existing_data_kmeans, yaml_file, default_flow_style=False)
        
        
    
    
    
