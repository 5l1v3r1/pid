#!/usr/bin/python


# Script that aligns, for a given barcode, all the consensus sequences (each one corresponding to a pID)
# Input: consensus sequences file for a given barcode (in ../templates/dir-<date_and_id>_consensus/)
# Output: aligned consensus sequences file for a given barcode (in ../templates/dir-<date_and_id>_align-global/)


import numpy as np
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
import os
import sys
import time
import lib_tools as lt


auto_file_name = str(sys.argv[0])

######
# DEF FUNCTIONS
######


if(len(sys.argv)>1):
    # parse the input consensus sequences file name
    rundir=str(sys.argv[1]).rstrip('/')+'/'
    if len(sys.argv)>2:
        analysis_type = '_'+sys.argv[2]

    barcode_directories = glob.glob(rundir+'bc_*_analysis*')
    for dname in barcode_directories:
        consensus_file_list = glob.glob(dname+'/consensus*fasta')
        time_start = time.time()
        for fname in consensus_file_list:
            aligned_fname = lt.trim_extension(fname)+'_aligned.fasta'

            cline = MuscleCommandline(input = fname, out = aligned_fname)
            cline()

        time_end = time.time()
        print 'alignment computation time: ' + str(time_end - time_start)

else:
    print auto_file_name + ': usage: '+ auto_file_name + ' <directory containing the directories for each barcode containing the consensus sequence files'
