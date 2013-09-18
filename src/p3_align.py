#!/ebio/ag-neher/share/programs/bin/python2


## script that aligns the temp reads files (one per pID/barcode) in the given directory
## input: directory (for a given barcode) containing the directories named temp_*
## output: aligned reads files (one per pID) contained in each of these directories
import numpy as np
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
import os
import sys
sys.path.append('./src/')
import time
import lib_tools as lt
import glob

######
# DEF FUNCTIONS
######


if (len(sys.argv)==2):
    ####
    # parse the input directory name
    ####
    relative_path_to_temp_directory = str(sys.argv[1])
    relative_path_to_temp_directory = relative_path_to_temp_directory.rstrip('/')+'/'

    list_temp_files = glob.glob(relative_path_to_temp_directory+'*.fasta')

    time_start = time.time()
    for fname in list_temp_files:
        align_file_name= lt.trim_extension(fname)+'_aligned.fasta'
        print 'Alignment for file: ' + fname
        cline = MuscleCommandline(input = fname, out = align_file_name)
        #
        cline()
        #
    time_end = time.time()
    print 'computation time: '+str(time_end-time_start)

else:
    print auto_file_name+': usage: '+auto_file_name+' <temp directory containing temp reads files to align (in ../templates/)>'


