#!/usr/bin/python
#/!\#
#### /ebio/ag-neher/share/pograms/EPD/bins/python
#/!\#

## script that gathers the reads per pID (generates temp files) to prepare the alignment process
## input: the filter reads file for a given barcode (in ../data/)
## output: the temp files, one for each pID, are created in a barcode corresponding directory (in ../templates/)


import numpy as np
import matplotlib.pyplot as plt
import glob
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

batchsize = 10
######
# MAIN
######
if (len(sys.argv) ==3):
    rundir = str(sys.argv[1]).rstrip('/')+'/'
    readtype = sys.argv[2]

    bc_dir_list = glob.glob(rundir+'bc_*_analysis')
    for bc_dir in bc_dir_list:
        print 'working in: '+bc_dir
        reads_file = bc_dir.rstrip('/')+'/'+readtype+'_reads.fasta'
        print ' -- with reads_file: '+reads_file

        # import the file containing, for the given barcode, the sequences for each pID
        # generate the pID specific temp files for alignments
        count=0
        dict_pIDs = defaultdict(list)
        all_pIDs = []
        with open(reads_file, 'r') as input_file:
            print '-- reading: '+reads_file
            for record in SeqIO.parse(input_file, 'fasta'):
                pID = str(record.id.split('_')[0])
                dict_pIDs[pID].append((record.id,record.seq))
                
        all_pIDs = sorted(dict_pIDs.keys())

        
        batch=0
        for pii, pID in enumerate(all_pIDs):
            if((pii%batchsize)==0):
                batch=pii/batchsize
                temp_pid_files_dir =bc_dir+'/temp_'+readtype+'_'+"{0:04d}".format(batch)
                print 'pii: '+str(pii)+', '+temp_pid_files_dir
                lt.check_and_create_directory(temp_pid_files_dir)
            with open(temp_pid_files_dir+'/'+ pID+'.fasta', 'w') as output_pID_file:
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
    print auto_file_name+': usage: '+auto_file_name+' <run directory> <readtype>'
