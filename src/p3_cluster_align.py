#!/usr/bin/python

import numpy as np
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
import os
import sys
import time
import glob 

####
# retrieve the temp reads file names to align in the specific cluster directory (in ../templates/) given in input
# the given cluster directory path must be relative to the current directory (src)
####
if(len(sys.argv)==2):
    run_dir = str(sys.argv[1])
    run_dir = run_dir.rstrip('/')+'/'

    barcode_dir_list = glob.glob(run_dir+'bc_*_analysis')
    for analysis_dir in barcode_dir_list:
        #collect the list of temp reads files to align
        list_temp_dirs = glob.glob(analysis_dir+'/temp_*')
        #run one job per file to align on the cluster
        for temp_dir in list_temp_dirs:
            cmd = 'qsub -cwd -l h_rt=12:00:00 -l h_vmem=10G ./src/p3_align.py '+temp_dir
            print cmd
            os.system(cmd)

else:
    print auto_file_name+': usage: '+auto_file_name+' <directory of a run>'



