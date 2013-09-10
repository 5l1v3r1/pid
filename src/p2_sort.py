#!/usr/bin/python
#/!\#
#### /ebio/ag-neher/share/pograms/EPD/bins/python
#/!\#

## script that gathers the reads per pID (generates temp files) to prepare the alignment process
## input: the filter reads file for a given barcode (in ../data/)
## output: the temp files, one for each pID, are created in a barcode corresponding directory (in ../templates/)


import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import datetime
import time
from Bio import SeqIO
from Bio import AlignIO
from Bio.Alphabet import generic_dna
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MuscleCommandline
from Bio.SeqRecord import SeqRecord
from collections import Counter
from collections import defaultdict
import lib_tools as lt


auto_file_name = str(sys.argv[0])


######
class struct_var_set:
    def __init__(self):
        #self.templates={}
        self.barcode = 'NYD'


######
# MAIN
######
res = struct_var_set()

if (len(sys.argv) == 2):
    relative_path_to_filtered_reads_file = str(sys.argv[1])
    path_to_data_file = "../data/"
    path_to_templates = "../templates/"
    lt.check_and_create_directory(path_to_templates)
    filtered_reads_file_basename = relative_path_to_filtered_reads_file.split('/')[-1]
    dict_pIDs = defaultdict(list)
    
    #import the file containing, for the given barcode, the sequences for each pID
    #generate the pID specific temp files for alignments
    with open(relative_path_to_filtered_reads_file, 'r') as input_file:
        [prefix_date_and_id, res.barcode] = [filtered_reads_file_basename.split('_')[i] for i in [0,1]]
        count=0
        for record in SeqIO.parse(input_file, 'fasta'):
            pID = str(record.id.split('_')[0])
            dict_pIDs[pID].append((record.id,record.seq))

        #create the directory for the given barcode
        dir_name_bc = str(path_to_templates+'dir-'+prefix_date_and_id+'_temp_'+res.barcode)
        lt.check_and_create_directory(dir_name_bc)
        for pID in dict_pIDs.keys():
            # write the temp files for each pID in the corresponding barcode directory
            with open(dir_name_bc+'/'+prefix_date_and_id+'_temp_'+res.barcode+'_'+pID +'.fasta', 'w') as output_pID_file:
                for read in dict_pIDs[pID]:
                    output_pID_file.write(str('>'+read[0]+'\n'))
                    output_pID_file.write(str(read[1]+'\n'))
                    count+=1
                    if(count%500==0):
                        print 'count = ' + str(count)
        print 'total : ' + str(count)
        #else:
        #    print 'no file created'
else: 
    print auto_file_name+': usage: '+auto_file_name+' <filtered reads file (in ../data/)>'
