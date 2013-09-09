#!/usr/bin/python

import numpy as np
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
import os
import sys
import time

auto_file_name = str(sys.argv[0])

####
# retrieve the temp reads file names to align in the specific cluster directory (in ../templates/) given in input
# the given cluster directory path must be relative to the current directory (src)
####
if(len(sys.argv)==2):
    cluster_dir = str(sys.argv[1])
    if cluster_dir[-1]!='/':
        cluster_dir=cluster_dir+'/'
    path_to_dir = "../templates/"
    cluster_dir_base_name = cluster_dir[len(path_to_dir):] # "../templates/" to remove
    
    #collect the list of temp reads files to align
    list_temp_files = os.popen(str('ls '+cluster_dir+'*')).readlines()
    list_temp_files = [file_to_strip.strip() for file_to_strip in list_temp_files]
    print list_temp_files
    #run one job per file to align on the cluster
    for cur_file in list_temp_files:
        cmd = 'qsub -cwd -l h_rt=12:00:00 -l h_vmem=10G ./p3_cluster_align_aux.py ./p3_cluster_align_script.py ' + cur_file
        #cmd = './p3_cluster_align_aux.py p3_cluster_align_script.py ' + cur_file
        print cmd
        os.system(cmd)

else:
    print auto_file_name+': usage: '+auto_file_name+' <specific cluster directory relative path (in ../templates/)>'



