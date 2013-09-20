#!/usr/bin/python

####
# Script that detects the possible mutants (pID mutations) with indels and generates a new filtered reads file (in ../data) resolving the pID mutations 
# 
# input: - filtered reads file (in ../data/)
#        - consensus sequences file (in ../templates/dir-<date_and_id>_consensus/)
# 
# output: filtered reads file (in ../data/) resolving the pID mutations, new <date_and_id> = <date_and_id>-mutant-indels
#         new read id: <pID-parent>_<nb-reads-fwd>_<nb-reads-rev>_<pID-origin>  
####

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import Counter
from collections import defaultdict
from itertools import imap
import operator
import sys
import datetime
import time
import lib_tools as lt

auto_file_name = str(sys.argv[0])


# distance max allowed for the consensus seq. alignments
DIST_MAX = 3

# minimum #times a pID must occur to have neighbors 
MIN_OCC_MOST_ABUND_PID = 10

max_pid_alignment_dist = 1.6
# anything <1 does not allow any substitions
# 1 allows only one point mutation
# 1.6 allows one indel and one substitution
# alignment parameters: match=1, unmatch=0, opening gap penalty=-0.6, gap extension penalty=-0.3


########
########



if (len(sys.argv) == 3):
    # parse the input data files names
    barcode_dir = str(sys.argv[1]).rstrip('/')+'/'
    readtype = sys.argv[2]
      
    # Create the dictionary of pIDs (from the output files of p1_trim_and_filter)

    ###################################
    # -- READ THE FILTERED READS FILE #
    ###################################
    # dictionary containing the list of reads (values) associated to a pID (keys)
    dict_all_reads, count_reads_added_to_dict_all_reads, count_reads_with_invalid_pID = lt.parse_readfile(barcode_dir+'aligned_reads_'+readtype+'.fasta')
    # control prints
    ####
    print 'count_reads_with_invalid_pID: '+str(count_reads_with_invalid_pID)
    print 'count_reads_added_to_dict_all_reads: '+str(count_reads_added_to_dict_all_reads)
    print '#pIDs different: '+str(len(dict_all_reads.keys()))


    pID_abundance = {}

    for pID,reads in dict_all_reads.iteritems():
        pID_abundance[pID] = sum([x[0]+x[1] for x in reads])

    print '#pIDs occurring only once: ' + str(sum([x==1 for x in pID_abundance.values()]))
    print '#pIDs occurring only twice: ' + str(sum([x==2 for x in pID_abundance.values()]))
    print '#pIDs occurring at least 3 times: ' + str(sum([x>2 for x in pID_abundance.values()]))
    ####


    ########################################
    # -- READ THE CONSENSUS SEQUENCES FILE #
    ########################################
    # dictionary containing the consensus sequence (values) associated to a pID (keys) if it occurs at least 3 times
    dict_cons_seq, a, b = lt.parse_readfile(barcode_dir+'consensus_sequences_'+readtype+'.fasta')
    print 'From consensus seq. file :'
    print '#pIDs occurring at least 3 times: ' + str(len(dict_cons_seq.keys()))


    ###########################
    # -- SEARCH THE NEIGHBORS #
    ###########################
    # dictionary indicating the "neighbor" state of each pID
    dict_neighbor_state_pIDs = {}
    for pID in dict_all_reads.keys():
        dict_neighbor_state_pIDs[pID] = False

    list_occ_pIDs_dec = sorted(pID_abundance.items(), key = lambda x:x[0], reverse=True)
    list_occ_pIDs_inc = list_occ_pIDs_dec[::-1]

    print 'Check list_occ_pIDs_dec and list_occ_pIDs_inc lengths (should be equal): '        
    print len(list_occ_pIDs_dec), len(list_occ_pIDs_inc)

    # dictionary containing the parent pID (values) of a given pID (keys)
    dict_neighbor_of = defaultdict(list)
    dict_neighbors = defaultdict(list)

    lim_pID_a = [list_occ_pIDs_inc[i][0]>=MIN_OCC_MOST_ABUND_PID for i in range(len(list_occ_pIDs_inc))].count(True)
    print 'lim_pID_a:'+str(lim_pID_a)
    print 'nb occ pID_a at pos lim_pID_a-1 in list_occ_pIDs_dec: '+str(list_occ_pIDs_dec[lim_pID_a-1][0])
    print 'nb occ pID_a at pos lim_pID_a in list_occ_pIDs_dec: '+str(list_occ_pIDs_dec[lim_pID_a][0])

    for ii_a, pID_a in enumerate(list_occ_pIDs_dec[:lim_pID_a]):
        # this condition should always evaluate to true
        if dict_neighbor_state_pIDs[pID_a[1]] == False:
            continue
        # loop over rare pid
        for ii,pID_r in enumerate(list_occ_pIDs_inc[:-(ii_a+1)]):
            # check whether this PID is already a neighbor of some other
             if((dict_neighbor_state_pIDs[pID_r[1]] == False) and (pID_r[0] < pID_a[0])):
                pIDs_align_score = align.localms(pID_a[1], pID_r[1], 1, 0, -0.6, -0.3)
                if (len(pIDs_align_score) != 0): # if pID_a[1] and pID_r[1] can be aligned
                    if (pIDs_align_score[0][2] >= len(pID_a[1]) - max_pid_alignment_dist):
                        print 'pID_a: '+str(pID_a)+', pID_r: '+str(pID_r)
                        print pIDs_align_score[0][0], pIDs_align_score[0][1]
                        print dict_neighbor_state_pIDs[pID_a[1]], dict_neighbor_state_pIDs[pID_r[1]]
                        print "possible neighbors"

                        if pID_r[0]>2: 
                            neighbor_candidates = dict_cons_seq[pID_r[1]]
                        else:
                            neighbor_candidates = dict_all_reads[pID_r[1]]


                        seq_a = dict_cons_seq[pID_a[1]][0][2] 
                        #if (pID_r[0] != 2 or (pID_r[0] == 2 and len(dict_all_reads[pID_r[1]])==1)): # pID_r occurs once or more than twice, or twice with identical reads: only one sequence comparison to do
                            #if(pID_r[0] <= 2): # pID_r occurs once or twice (with identical reads): its sequence is in dict_all_reads
                                # seq_r = dict_all_reads[pID_r[0]][0][2]
                                #seq_r = dict_all_reads[pID_r[1]][0][2] # [nb_fwd,nb_rev,seq]
                            #else: # pID_r occurs at least 3 times: its consensus sequence is in dict_cons_seq
                                #seq_r = dict_cons_seq[pID_r[1]][0][2] # [nb_fwd,nb_rev,seq]
                        #else: # pID_r occurs twice with distinct reads: if one read aligns well: pID_r is considered as a mutant
                            #[seq_r_0,seq_r_1] = [dict_all_reads[pID_r[1]][i][2] for i in [0,1]]
                            
                        #if(pID_)


                        for nf, nr, seq_r in neighbor_candidates:
                            if lt.check_neighbor_plausibility(seq_a, seq_r, DIST_MAX):
                                dict_neighbors[pID_a[1]].append(pID_r[1])
                                dict_neighbor_of[pID_r[1]] = pID_a[1]
                                dict_neighbor_state_pIDs[pID_r[1]] = True
                                print "NEIGHBORS !"
                                break


    print 'neighbors:'
    nb_neighbors_stated = 0
    for pID,state in dict_neighbor_state_pIDs.iteritems():
        if state:
            nb_neighbors_stated += 1
            print pID



    ##############################################
    # -- WRITE THE NEIGHBORS FILTERED READS FILE #
    ##############################################
    count_reads_written = 0
    count_neighbors_reads_written = 0

    # write the new filtered reads file (mutant-indels):
    corrected_aligned_reads_fname = barcode_dir+'corrected_reads.fasta'
    with open(corrected_aligned_reads_fname, 'w') as neighbors_filtered_reads_file:
        for pID in dict_all_reads.keys(): # for all pIDs                                
            if not dict_neighbor_state_pIDs[pID]: # if pID "is not a neighbor"
                print 'write reads for pID: ' + pID

                # write the reads corresponding to pID
                for cur_read in dict_all_reads[pID]: # cur_read is a list [nb_occ_fwd,nb_occ_rev,seq]
                    #pid_orig = pID
                    #mut_flag = ''
                    new_read_id = str(pID+'_'+str(cur_read[0])+'_'+str(cur_read[1])+'_'+pID)
                    neighbors_filtered_reads_file.write('>'+new_read_id+'\n') # write the new id of the read
                    neighbors_filtered_reads_file.write(cur_read[2]+'\n') # write the sequence
                    count_reads_written += (cur_read[0]+cur_read[1])

                # write the neighbors if there are
                if (len(dict_neighbors[pID]) > 0):
                    print '#neighbors pIDs: ' + str(len(dict_neighbors[pID]))
                    for pID_n in dict_neighbors[pID]: # for all the neighbors pID_n of pID

                        # write the neighbors pID_n
                        for read_n in dict_all_reads[pID_n]: # read_n is a list [nb_occ_fwd,nb_occ_rev,seq]
                            pid_orig = pID_n
                            mut_flag = '_m'
                            new_read_id = str(pID+'_'+str(read_n[0])+'_'+str(read_n[1])+'_'+pid_orig+mut_flag)
                            neighbors_filtered_reads_file.write('>'+new_read_id+'\n') # write the new id of the read
                            neighbors_filtered_reads_file.write(read_n[2]+'\n') # write the sequence
                            count_reads_written += (read_n[0]+read_n[1])
                            count_neighbors_reads_written += (read_n[0]+read_n[1])


    print 'nb reads written: '+str(count_reads_written)
    print 'nb neighbored reads written: '+str(count_neighbors_reads_written)
    print 'nb neighbored pIDs: '+str(nb_neighbors_stated)


else:
    print auto_file_name + ': usage: '+ auto_file_name + ' <directory with reads from a barcode> <type of read>'
