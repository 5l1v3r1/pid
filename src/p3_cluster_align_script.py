#!/usr/bin/python

import numpy as np
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
import os
import sys
import time

auto_file_name = str(sys.argv[0])

####
#align the reads of the temp reads file given in input
####
if len(sys.argv)==2:
    relative_path_to_file_to_align = str(sys.argv[1])
    path_to_templates = "../templates/"
    [cluster_dir_basename,file_to_align_basename] = [relative_path_to_file_to_align.split('/')[i] for i in [2,3]] # "../templates/cluster-<date+id>/<date+id>_temp_<barcode>_<pID>.fasta"
    #print cluster_dir_basename
    [prefix_date_and_id,barcode,pID_to_clean] = [file_to_align_basename.split('_')[i] for i in [0,2,3]] #"<date+id>_temp_<barcode>_<pID>.fasta"
    pID = pID_to_clean.split('.')[0] # remove the file extension

    print '****** job for the file: ' + relative_path_to_file_to_align

    temp_file_name = relative_path_to_file_to_align
    align_file_name = path_to_templates+cluster_dir_basename+'/'+prefix_date_and_id+'_align_'+barcode+'_'+pID+'.fasta'
    cline = MuscleCommandline(input = temp_file_name, out = align_file_name)
    #print cline
    #os.system(str('echo '+pID+' > '+align_file_name))
    time_start = time.time()
    #
    cline()
    #
    time_end = time.time()
    
    print 'time: ' + str(time_end - time_start)
    
else:
    print auto_file_name+': miss the temp reads file to align!'
    
