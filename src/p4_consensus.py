#!/usr/bin/python

#
# script that creates consensus sequences files for a given barcode (barcode corresponding to the name of the input aligned files directory (in ../templates/))
# input: aligned files directory (in ../templates/)
# output: consensus sequences files in specific directory (../templates/dir-<date_and_id>_consensus/)

import numpy as np
from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import shutil
import sys
import glob
import time
import lib_tools as lt

auto_file_name = str(sys.argv[0])

######
# DEF FUNCTIONS
######


######
def make_consensus(file_name):
    #print 'consensus for file: '+file_name
    alignment = AlignIO.read(file_name,"fasta")
    nb_reads_fwd = np.asarray([lt.parse_read_label(seq.id)[1] for seq in alignment], int)
    #print nb_reads_fwd
    nb_reads_rev = np.asarray([lt.parse_read_label(seq.id)[2] for seq in alignment], int)
    #print nb_reads_rev
    nb_reads = nb_reads_fwd + nb_reads_rev

    alignment = np.asarray(alignment)
    len_align,len_seq = alignment.shape
    alpha = np.array(list('ATCG-N'))
    nuc_counts = np.zeros((alpha.shape[0], len_seq))
    if(np.sum(nb_reads) >= 3):
        #compute the consensus sequence taking into account the numbers of fwd/rev reads corresponding to the pID
        for ni,nuc in enumerate(alpha):
            nuc_counts[ni] = np.sum((alignment == nuc)*np.repeat(np.array([nb_reads]).T, len_seq, axis=1), axis=0)
        str_consensus = "".join(alpha[np.argmax(nuc_counts, axis=0)])
        
        #remove gaps
        #str_consensus_2 = "".join([nuc for nuc in str_consensus if nuc!='-'])
        str_consensus_2 = lt.remove_gaps(str_consensus)

        return (str_consensus_2, np.sum(nb_reads_fwd), np.sum(nb_reads_rev))    
    else:
        print 'file: '+file_name+': not enough reads (#reads < 3) to make consensus'
        return ('nothing',0,0)
#######



if __name__=='__main__':
    if(len(sys.argv)==3):
        # parse the input consensus sequences file name
        run_dir=str(sys.argv[1]).rstrip('/')+'/'
        readtype = sys.argv[2]
        
        barcode_dir_list = glob.glob(run_dir+'bc_*_analysis')

        for analysis_dir in barcode_dir_list:
            print "analyzing directory:", analysis_dir
            temp_directories = glob.glob(analysis_dir+'/temp_'+readtype+'*')
            consensus_fname = analysis_dir+'/consensus_sequences_'+readtype+'.fasta'
            aligned_reads_fname = analysis_dir+'/aligned_reads_'+readtype+'.fasta'

            with open(consensus_fname, 'w') as consensus_file, \
                    open(aligned_reads_fname, 'w') as aligned_reads_file:
                print temp_directories
                for temp_dir in temp_directories:
                    pID_files = glob.glob(temp_dir+'/*aligned.fasta')
                    print pID_files
                    for pID_file in pID_files:
                        pID = lt.get_last_part_of_path(pID_file).split('_')[0]
                        with open(pID_file, 'r') as infile:
                            tmp_aln = AlignIO.read(infile, 'fasta')
                        AlignIO.write(tmp_aln, aligned_reads_file, 'fasta')
                        consensus_seq = make_consensus(pID_file)
                        if consensus_seq[1]+consensus_seq[2]>2:
                            consensus_file.write('>'+lt.read_label(pID, consensus_seq[1], consensus_seq[2])+'\n')
                            consensus_file.write(consensus_seq[0]+'\n')
                                
                    shutil.rmtree(temp_dir)
    else:
        print auto_file_name+': usage: '+auto_file_name+' <run directory> <read type to work on>'

