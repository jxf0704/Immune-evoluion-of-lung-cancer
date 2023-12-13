# -*- coding: utf-8 -*-

"""
Generate bash script for running SCHISM for each patient and stored as txt files

"""

import os

mode_1_path="XXXXX"    
patient_name=os.listdir(mode_1_path)
runID=10
# first check yaml file to see if working_dir exists, if existing, using the yaml file; if not, write the value of work_dir to it,using AP clustering method as an example to generate SCHISM bash script for each patient. 
 
for patient in patient_name:
    patient=str(patient)
    work_dir=os.path.join(mode_1_path,str(patient))
    txt1='''#! /bin/bash
    # run schism in step-through mode
    
    runSchism=/jxf/anaconda2/bin/runSchism
    baseDir=%s
    currentPath=$(grep 'working_dir' $baseDir/%s_AP.yaml | sed 's/working_dir://g')
    
    if [ -d $currentPath ]; then
    \techo "The working directory specified in %s_AP.yaml exists."
    \techo "proceeding to analysis"
    # run schism in sequential mode
    \tconf="$baseDir/%s_AP.yaml"
    else
    \techo "The working directory specified in %s_AP.yaml does not exists."
    \techo "Using the path to this shell script as the working directory."
    
    \techo "# adjusted path to working directory" >> $baseDir/%s_AP.path.yaml
    \techo "working_dir: $baseDir" >> $baseDir/%s_AP.path.yaml
    \techo "" >> $baseDir/%s_AP.path.yaml
    \tgrep -v "working_dir" $baseDir/%s_AP.yaml >> $baseDir/%s_AP.path.yaml
    # run schism in sequential mode
    \tconf="$baseDir/%s_AP.path.yaml"
    fi
    
    runSchism prepare_for_hypothesis_test -c $conf  
    
    runSchism hypothesis_test -c $conf   ''' % (work_dir,patient,patient,patient,patient,patient,patient,patient,patient,patient,patient)
    
    txt2 = 'runSchism plot_cpov -c $conf  \n\n'
    for i in range(1, runID + 1):
        txt2 += 'runSchism run_ga --config $conf --mode parallel --runID ' + str(i) + '\n'

    txt3 = '''runSchism run_ga --config $conf --mode serial

    runSchism summarize_ga_results --config $conf

    runSchism consensus_tree -c $conf '''


    f = open(os.path.join(work_dir,patient+'_AP.sh'),'w')
    f.write(txt1 + '\n' + txt2 + '\n' + txt3)
    f.close()
