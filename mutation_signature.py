# -*- coding: utf-8 -*-


#from SigProfilerMatrixGenerator import install as genInstall
#genInstall.install('GRCh37', rsync=False,bash=False)
import os
#os.path.dirname(os.path.abspath(__file__))
#import shutil
#shutil.which("wget") 
#os.environ()
#shutil.which("python")

import pandas as pd
from SigProfilerExtractor import sigpro as sig  


input_dir_primary="I:\\xiaojun\\sigprofiler\\signature_input_primary"
input_dir_distance="I:\\xiaojun\\sigprofiler\\signature_input_distance"
input_dir_distance_specific="I:\\xiaojun\\sigprofiler\\signature_input_distance_specific_pair"

#os.chdir("I:\\xiaojun\\sigprofiler\\")

mutation_table_73=pd.read_csv(os.path.join(input_dir_distance_specific,"mutation_table_exon_73_Distance_specific_pair.csv"))
sample_id=set(list(mutation_table_73["Sample_ID"]))
for s in sample_id:
    s_table=mutation_table_73[mutation_table_73["Sample_ID"]==s]
    s_table=s_table.loc[:,["Chromosome","Start_Position","Sample_ID","Reference_Allele","Tumor_Seq_Allele2"]]
    s_table.to_csv(os.path.join(input_dir_distance_specific,s+".vcf"),sep="\t",index=False,header=False)


data_primary = "I:/xiaojun/sigprofiler/signature_input_primary"
data_distance = "I:/xiaojun/sigprofiler/signature_input_distance"
data_distance_specific = "I:/xiaojun/sigprofiler/signature_input_distance_specific_pair"

sig.sigProfilerExtractor("vcf", "result_73_exome_primary", data_primary, startProcess=1, endProcess=8, genome_build = "GRCh37", refgen="GRCh37",mtype = ["96"],exome=True)

sig.sigProfilerExtractor("vcf", "result_73_exome_distance", data_distance, startProcess=1, endProcess=10, genome_build = "GRCh37", refgen="GRCh37",mtype = ["96"],exome=True)

sig.sigProfilerExtractor("vcf", "result_73_exome_distance_specific_pair", data_distance_specific, startProcess=1, endProcess=8, genome_build = "GRCh37", refgen="GRCh37",mtype = ["96"],exome=True)




