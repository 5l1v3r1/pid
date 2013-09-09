#!/usr/bin/python

#
# script that creates consensus sequences files for a given barcode (barcode corresponding to the name of the input aligned files directory (in ../templates/))
# input: aligned files directory (in ../templates/)
# output: consensus sequences files in specific directory (../templates/dir-<date_and_id>_consensus/)

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


######
def make_consensus(file_name, pID):
    #print 'consensus for file: '+file_name
    alignment = AlignIO.read(file_name,"fasta")
    len_align = len(alignment)
    nb_reads_fwd = sum([int(alignment[i].id.split('_')[1]) for i in range(len_align)])
    #print nb_reads_fwd
    nb_reads_rev = sum([int(alignment[i].id.split('_')[2]) for i in range(len_align)])
    #print nb_reads_rev
    nb_reads = nb_reads_fwd + nb_reads_rev
    len_seq = alignment.get_alignment_length()
    if(nb_reads >= 3):
        #compute the consensus sequence taking into account the numbers of fwd/rev reads corresponding to the pID
        
        str_consensus = "".join([Counter("".join([alignment[i,j]*(int(alignment[i].id.split('_')[1])+int(alignment[i].id.split('_')[2])) for i in range(len_align)])).most_common(1)[0][0] for j in range(len_seq)])
        
        #remove gaps
        str_consensus_2 = "".join([i!='-' and i or '' for i in str_consensus])

        return (str_consensus_2, nb_reads_fwd, nb_reads_rev)    
    else:
        print 'file: '+file_name+': not enough reads (#reads < 3) to make consensus'
        return ('nothing',0,0)
#######




if (len(sys.argv)==2):
    # parse the input aligned files directory name
    relative_path_to_align_dir=str(sys.argv[1])
    if(relative_path_to_align_dir[-1]!='/'):
        relative_path_to_align_dir+='/'
    path_to_templates = "../templates/"
    align_dir_basename = relative_path_to_align_dir.split('/')[2]
    [prefix_date_and_id,file_type,barcode] = [align_dir_basename.split('dir-')[1].split('_')[i] for i in [0,1,2]]

    # create (if necessary) the consensus directory
    lt.check_and_create_directory(str(path_to_templates+'dir-'+prefix_date_and_id+'_consensus'))
    
    # get the list of files to consensus (one per pID)
    list_files_to_consensus = os.popen(str('ls '+relative_path_to_align_dir+'*')).readlines()
    dict_consensus_seq = {}
    for cur_file in list_files_to_consensus:
        cur_file = cur_file.strip() # remove '\n'
        cur_file_pID = cur_file.split('_')[-1].split('.')[0] #"<prefix_date_and_id>_align_<barcode>_<pID>.fasta"
        cur_file_barcode = cur_file.split('_')[2]
        #print cur_file
        result_consensus_cur_file = make_consensus(cur_file,cur_file_pID)
        if (result_consensus_cur_file != ('nothing',0,0)): # store the consensus sequence
            dict_consensus_seq[cur_file_pID]=result_consensus_cur_file
        #print result_consensus_cur_file
    # write the consensus sequences of each pID (with #occ >= 3) in the specific consensus seq file for the given barcode
    with open(str(path_to_templates+'dir-'+prefix_date_and_id+'_consensus/'+prefix_date_and_id+'_consensus_'+barcode+'.fasta'),'w') as out_file:
        for pID in sorted(dict_consensus_seq.keys()):
            out_file.write('>'+pID+'_'+str(dict_consensus_seq[pID][1])+'_'+str(dict_consensus_seq[pID][2])+'\n')
            out_file.write(dict_consensus_seq[pID][0]+'\n')
        
    
else:
    print auto_file_name+': usage: '+auto_file_name+' <aligned reads files directory (in ../templates/)>'

