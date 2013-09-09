#!/usr/bin/python


#
# script that gathers all the aligned reads for a given barcode in one file
#

import numpy as np
from Bio import AlignIO
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from collections import Counter
from collections import defaultdict
import os
import sys
import time
import datetime
import lib_tools as lt

auto_file_name = str(sys.argv[0])

######
# DEF FUNCTIONS
######


if (len(sys.argv)==2):
    # parse the input aligned files directory name
    relative_path_to_align_dir=str(sys.argv[1])
    if(relative_path_to_align_dir[-1]!='/'):
        relative_path_to_align_dir+='/'
    path_to_templates = "../templates/"
    align_dir_basename = relative_path_to_align_dir.split('/')[2]
    [prefix_date_and_id,file_type,barcode] = [align_dir_basename.split('dir-')[1].split('_')[i] for i in [0,1,2]]

    # create (if necessary) the all-align-in-one-file directory
    lt.check_and_create_directory(str(path_to_templates+'dir-'+prefix_date_and_id+'_all-align-in-one-file'))
    
    # get the list of files to gather 
    list_files_to_gather = os.popen(str('ls '+relative_path_to_align_dir+'*')).readlines()
    list_files_to_gather = [file_to_gather.strip() for file_to_gather in list_files_to_gather]
    
    
    global_list_alignments = []
    
    for i in range(len(list_files_to_gather)):
        list_seq_cur_file = []
        input_file_name = list_files_to_gather[i]
        for seq_record in SeqIO.parse(input_file_name, 'fasta'):
            list_seq_cur_file.append(seq_record)
        
        cur_alignment = MultipleSeqAlignment(list_seq_cur_file)
        global_list_alignments.append(cur_alignment)
    AlignIO.write(global_list_alignments, str(path_to_templates+'dir-'+prefix_date_and_id+'_all-align-in-one-file/'+prefix_date_and_id+'_all-align-in-one-file_'+barcode+'.fasta'), 'fasta')
       
    
else:
    print auto_file_name+': usage: '+auto_file_name+' <aligned reads files directory (in ../templates/)>'
