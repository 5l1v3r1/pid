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


if __name__=='__main__':
    if(len(sys.argv)==2):
        # parse the input consensus sequences file name
        relative_path_to_cons_seq_file=str(sys.argv[1])
        path_to_templates = "../templates/"
        cons_seq_file_basename = relative_path_to_cons_seq_file.split('/')[3]
    
        [prefix_date_and_id,file_type,barcode] = [cons_seq_file_basename.split('_')[i] for i in [0,1,2]]
        barcode = barcode.split('.')[0] # remove file extension ".fasta"

        # create (if necessary) the consensus directory
        lt.check_and_create_directory(str(path_to_templates+'dir-'+prefix_date_and_id+'_align-global'))
    
        # align
        aligned_cons_seq_file_name = str(path_to_templates+'dir-'+prefix_date_and_id+'_align-global/'+prefix_date_and_id+'_align-global_'+barcode+'.fasta') 
        print 'Global consensus sequences alignment for barcode: ' + barcode + ' (from file: ' + relative_path_to_cons_seq_file + ' )'
        
        cline = MuscleCommandline(input = relative_path_to_cons_seq_file, out = aligned_cons_seq_file_name)
        
        time_start = time.time()
        #
        cline()
        #
        time_end = time.time()
        print 'alignment computation time: ' + str(time_end - time_start)

    else:
        print auto_file_name + ': usage: '+ auto_file_name + ' <consensus sequences file (in ../templates/dir-<date_and_id>_consensus)>'
