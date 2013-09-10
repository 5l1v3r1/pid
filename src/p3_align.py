#!/usr/bin/python
#/!\#
#### /ebio/ag-neher/share/pograms/EPD/bins/python
#/!\#

## script that aligns the temp reads files (one per pID/barcode) in the given directory
## input: directory (for a given barcode) containing the temp reads files to align (in ../templates/)
## output: aligned reads files (one per pID) contained in a directory "align" (in ../templates/)

import numpy as np
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
import os
import sys
import time
import lib_tools as lt
import glob

auto_file_name = str(sys.argv[0])


######
# DEF FUNCTIONS
######


if (len(sys.argv)==2):
    ####
    # parse the input directory name
    ####
    relative_path_to_temp_directory = str(sys.argv[1])
    if relative_path_to_temp_directory[-1]!='/':
        relative_path_to_temp_directory+='/'
    path_to_templates='../templates/'
    temp_directory_basename = relative_path_to_temp_directory.split('/')[-2]
    print temp_directory_basename
    [prefix_date_and_id, bc]= [temp_directory_basename.split('_')[i] for i in [0,2]]
    prefix_date_and_id = prefix_date_and_id[len('dir-'):] # "dir-" to remove
    print prefix_date_and_id, bc

    ####
    # collect the temp reads files to align for the given barcode
    ####
    list_temp_files = glob.glob(relative_path_to_temp_directory+'*')

    # list of "pIDs" to generate the pID alignment files 
    list_of_pIDs = [fname.split('_')[-1].split('.')[0] for fname in list_temp_files]
    print 'list_of_pIDs created, length: ' + str(len(list_of_pIDs))
    set_of_pIDs = set(list_of_pIDs) 
    print 'no duplicates in pIDs list: '+str(len(set_of_pIDs)==len(list_of_pIDs))

    ####
    # align the reads (temp files) and creates the corresponding aligned files in the align directory (created automatically)
    ####
    #create the directory for the aligned reads
    dir_align_name_bc = str(path_to_templates+'dir-'+prefix_date_and_id+'_align_'+bc)
    dir_temp_name_bc =  str(path_to_templates+'dir-'+prefix_date_and_id+'_temp_'+bc)

    lt.check_and_create_directory(dir_align_name_bc)
    time_start = time.time()
    for pID in list_of_pIDs:
        temp_file_name = dir_temp_name_bc+ '/' +prefix_date_and_id + '_temp_' + bc + '_' + pID + '.fasta'
        align_file_name=dir_align_name_bc+ '/' +prefix_date_and_id + '_align_'+ bc + '_' + pID + '.fasta'
        print 'Alignment for file: ' + temp_file_name
        cline = MuscleCommandline(input = temp_file_name, out = align_file_name)
        #
        cline()
        #
    time_end = time.time()
    print 'computation time: '+str(time_end-time_start)

else:
    print auto_file_name+': usage: '+auto_file_name+' <temp directory containing temp reads files to align (in ../templates/)>'


