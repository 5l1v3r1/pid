#!/usr/bin/python

import os
import sys
import time
import glob 

####
# retrieve the temp reads file names to align in the specific cluster directory (in ../templates/) given in input
# the given cluster directory path must be relative to the current directory (src)
####
if(len(sys.argv)==3):
    run_dir = str(sys.argv[1])
    run_dir = run_dir.rstrip('/')+'/'
    read_type = sys.argv[2]

    barcode_dir_list = glob.glob(run_dir+'bc_*_analysis')
    for analysis_dir in barcode_dir_list:
        #collect the list of temp reads files to align
        cmd = 'qsub -cwd -l h_rt=00:59:00 -l h_vmem=10G ./src/p8_detect_mutants_indels.py '+analysis_dir+ ' ' +read_type
        print cmd
        os.system(cmd)

else:
    print sys.argv[1]+': usage: '+sys.argv[0]+' <directory of a run> <read_type>'



