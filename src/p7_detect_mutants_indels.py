#!/usr/bin/python

####
# Script that detects the possible mutants (pID mutations) with indels and generates a new filtered reads file (in ../data) resolving the pID mutations 
# 
# input: - filtered reads file (in ../data/)
#        - consensus sequences file (in ../templates/dir-<date_and_id>_consensus/)
# 
# output: filtered reads file (in ../data/) resolving the pID mutations, new <date_and_id> = <date_and_id>-mutant-indels
#         new read id: <pID-parent>_<nb-reads-fwd>_<nb-reads-rev>_<pID-origin>[_m if a mutation occurred]  
####

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.pairwise2 import align
from collections import Counter
from collections import defaultdict
from itertools import imap
import operator
import sys
import datetime
import time
#import lib_tools as lt

auto_file_name = str(sys.argv[0])


# distance max allowed for the consensus seq. alignments
DIST_MAX = 3

# minimum #times a pID must occur to have neighbors 
MIN_OCC_MOST_ABUND_PID = 10




########
########


if __name__=='__main__':
    if (len(sys.argv) == 3):
        # parse the input data files names
        relative_path_to_filtered_reads_file = str(sys.argv[1])
        relative_path_to_cons_seq_file = str(sys.argv[2])
    
        path_to_data = '../data/'
        path_to_templates = '../templates/'
    
        filtered_reads_file_basename = relative_path_to_filtered_reads_file.split('/')[2]
        cons_seq_file_basename = relative_path_to_cons_seq_file.split('/')[3]

        [prefix_date_and_id, barcode] = [filtered_reads_file_basename.split('_')[i] for i in [0,1]]
        [prefix_date_and_id_cs, barcode_cs] = [cons_seq_file_basename.split('_')[i] for i in [0,2]]
        barcode_cs = barcode_cs.split('.')[0] # remove ".fasta"

        # check the barcode similarity
        if ((barcode == barcode_cs) and (prefix_date_and_id == prefix_date_and_id_cs)):
        
            # Create the dictionary of pIDs (from the output files of p1_trim_and_filter)
        
            # dictionary containing the list of reads (values) associated to a pID (keys)
            dict_all_reads = defaultdict(list)

            # dictionary containing the consensus sequence (values) associated to a pID (keys) if it occurs at least 3 times
            dict_cons_seq = defaultdict(list)

            count_reads_with_invalid_pID = 0
            count_reads_added_to_dict_all_reads = 0

            ###################################
            # -- READ THE FILTERED READS FILE #
            ###################################
            with open(relative_path_to_filtered_reads_file, 'r') as filtered_reads_file:
                # for all the valid reads (without Ns), add the read to dict_all_reads
                for record in SeqIO.parse(filtered_reads_file, 'fasta'):
                    [pID,nb_reads_fwd,nb_reads_rev] = [record.id.split('_')[i] for i in [0,1,2]]
                    nb_reads_fwd = int(nb_reads_fwd)
                    nb_reads_rev = int(nb_reads_rev)
                    seq = str(record.seq)
                    if (pID.count('N') == 0): # if the pID is not ambiguous
                        dict_all_reads[pID].append([nb_reads_fwd,nb_reads_rev,seq])
                        count_reads_added_to_dict_all_reads += (nb_reads_fwd+nb_reads_rev)
                    else:
                        # print 'Invalid pID: ' + pID
                        count_reads_with_invalid_pID += (nb_reads_fwd+nb_reads_rev)

            # control prints
            # ###
            print 'count_reads_with_invalid_pID: '+str(count_reads_with_invalid_pID)
            print 'count_reads_added_to_dict_all_reads: '+str(count_reads_added_to_dict_all_reads)
            print '#pIDs distinct: '+str(len(dict_all_reads.keys()))
        
            count_pIDs_occ = Counter()
        
            for pID in dict_all_reads.keys():
                nb_reads_diff = len(dict_all_reads[pID]) # = nb of distinct reads for the given pID
                nb_cur_reads = 0
                for i in range(nb_reads_diff):
                    nb_cur_reads += sum([dict_all_reads[pID][i][j] for j in [0,1]])
                count_pIDs_occ[((nb_cur_reads >= 3) and '3+') or ((nb_cur_reads == 2) and '2') or ('1')] += 1

            print '#pIDs occurring only once: ' + str(count_pIDs_occ['1'])#str(len(list_pIDs_occurring_once))
            print '#pIDs occurring only twice: ' + str(count_pIDs_occ['2'])#str(len(list_pIDs_occurring_twice))
            print '#pIDs occurring at least 3 times: ' + str(count_pIDs_occ['3+'])#str(len(list_pIDs_occurring_at_least_3_times))

        
            ########################################
            # -- READ THE CONSENSUS SEQUENCES FILE #
            ########################################
            # normally, this reading of the consensus file should not modify the list of pIDs (keys of dict_all_reads)
            with open(relative_path_to_cons_seq_file, 'r') as cons_seq_file:
                for record in SeqIO.parse(cons_seq_file, 'fasta'):
                    [pID,nb_occ_fwd,nb_occ_rev] = [record.id.split('_')[i] for i in [0,1,2]]
                    nb_occ_fwd = int(nb_occ_fwd)
                    nb_occ_rev = int(nb_occ_rev)
                    cons_seq = str(record.seq)
                    dict_cons_seq[pID].append([nb_occ_fwd,nb_occ_rev,cons_seq])
            print 'From consensus seq. file :'
            print '#pIDs occurring at least 3 times: ' + str(len(dict_cons_seq.keys()))
        

            ###########################
            # -- SEARCH THE NEIGHBORS #
            ###########################
            # dictionary indicating the "neighbor" state of each pID
            dict_neighbor_state_pIDs = {}
            count_nb_pIDs_occ_sup_eq_3_in_filtered_reads_file = 0
            for pID in dict_all_reads.keys():
                dict_neighbor_state_pIDs[pID] = False
        
                # sorted list of the #occurrences of each pID: couples (#occ, pID)
                list_occ_pIDs_dec = [(sum([dict_cons_seq[pID][i][0]+dict_cons_seq[pID][i][1] for i in range(len(dict_cons_seq[pID]))]), pID) for pID in dict_cons_seq.keys()]
                list_occ_pIDs_inc = [(sum([dict_cons_seq[pID][i][0]+dict_cons_seq[pID][i][1] for i in range(len(dict_cons_seq[pID]))]), pID) for pID in dict_cons_seq.keys()]
        
            # -- add the pIDs with #occurrences < 3 from dict_all_reads to list_occ_pIDs
            for pID in dict_all_reads.keys():
                nb_occ_cur_pID = 0
                for i in range(len(dict_all_reads[pID])):
                    nb_occ_cur_pID += sum([dict_all_reads[pID][i][j] for j in [0,1]]) # add nb_occ_fwd and nb_occ_rev
                if nb_occ_cur_pID < 3:
                    list_occ_pIDs_dec.append((nb_occ_cur_pID,pID))
                    list_occ_pIDs_inc.append((nb_occ_cur_pID,pID))
                else:
                    count_nb_pIDs_occ_sup_eq_3_in_filtered_reads_file += 1
                list_occ_pIDs_dec.sort(key = lambda couple: couple[0], reverse=True)
                list_occ_pIDs_inc.sort(key = lambda couple: couple[0])

            print 'Check list_occ_pIDs_dec and list_occ_pIDs_inc lengths (should be equal): '        
            print len(list_occ_pIDs_dec)
            print len(list_occ_pIDs_inc)
        
            print '#####\n#pIDs with occ >=3 in filtered reads file: '+str(count_nb_pIDs_occ_sup_eq_3_in_filtered_reads_file)
            print '# pIDs (with occ >= 3) in consensus seq. file: '+str(len(dict_cons_seq.keys()))
            

            # dictionary containing the parent pID (values) of a given pID (keys)
            dict_neighbor_of = defaultdict(list)
        
            lim_pID_a = [list_occ_pIDs_inc[i][0]>=MIN_OCC_MOST_ABUND_PID for i in range(len(list_occ_pIDs_inc))].count(True)
            print 'lim_pID_a:'+str(lim_pID_a)
            print 'nb occ pID_a at pos lim_pID_a-1 in list_occ_pIDs_dec: '+str(list_occ_pIDs_dec[lim_pID_a-1][0])
            print 'nb occ pID_a at pos lim_pID_a in list_occ_pIDs_dec: '+str(list_occ_pIDs_dec[lim_pID_a][0])

            for pID_a in list_occ_pIDs_dec[:lim_pID_a]:
                if dict_neighbor_state_pIDs[pID_a[1]] == False:
                    for ii,pID_r in enumerate(list_occ_pIDs_inc[:list_occ_pIDs_inc.index(pID_a)]):
                        if((dict_neighbor_state_pIDs[pID_r[1]] == False) and (pID_r[0] < pID_a[0])):
                            pIDs_align_score = align.localms(pID_a[1], pID_r[1], 1, 0, -0.6, -0.3)
                            if (len(pIDs_align_score) != 0): # if pID_a[1] and pID_r[1] can be aligned
                                if (pIDs_align_score[0][2] >= len(pID_a[1])-1-0.6):
                                    print 'pID_a: '+str(pID_a)+', pID_r: '+str(pID_r)
                                    print pIDs_align_score[0][0], pIDs_align_score[0][1]
                                    print dict_neighbor_state_pIDs[pID_a[1]], dict_neighbor_state_pIDs[pID_r[1]]
                                    print "possible neighbors"
                    
                                    if (pID_a[0] >= MIN_OCC_MOST_ABUND_PID): # if pID_a occurs enough times to have neighbors
                                        if (pID_r[0] != 2 or (pID_r[0]==2 and len(dict_all_reads[pID_r[1]])==1)): # pID_r occurs once or more than twice or twice with the same read: only one sequence comparison to do
                                            seq_a = dict_cons_seq[pID_a[1]][0][2] # [nb_fwd,nb_rev,seq]
                                            if(pID_r[0] <= 2): # pID_r occurs once or twice: its sequence is in dict_all_reads
                                                seq_r = dict_all_reads[pID_r[1]][0][2] # [nb_fwd,nb_rev,seq]
                                            else: # pID_r occurs at least 3 times: its consensus sequence is in dict_cons_seq
                                                seq_r = dict_cons_seq[pID_r[1]][0][2] # [nb_fwd,nb_rev,seq]
                                            score = align.localms(seq_r, seq_a, 1, 0, -1, -1)[0]
                                            print 'Alignment score: '+ str(score[2])
                                            print score[0]
                                            print score[1]
                                            if(score[2] >= len(seq_a) - DIST_MAX):
                                                dict_neighbor_of[pID_r[1]] = pID_a[1]
                                                dict_neighbor_state_pIDs[pID_r[1]] = True
                                                print "NEIGHBORS !"

                                        else: # pID_r occurs twice, occurrences with distinct reads: comparison of sequence of pID_a with its 2 sequences in dict_all_reads
                                            seq_a = dict_cons_seq[pID_a[1]][0][2] # [nb_fwd,nb_rev,seq]
                                            seq_r_0 = dict_all_reads[pID_r[1]][0][2] # [nb_fwd,nb_rev,seq]
                                            seq_r_1 = dict_all_reads[pID_r[1]][1][2] # [nb_fwd,nb_rev,seq]
                                            score_0 = align.localms(seq_r_0, seq_a, 1, 0, -1, -1)[0]
                                            score_1 = align.localms(seq_r_1, seq_a, 1, 0, -1, -1)[0]
                                            print 'Alignment score_0: '+ str(score_0[2])
                                            print score_0[0]
                                            print score_0[1]
                                            print 'Alignment score_1: '+ str(score_1[2])
                                            print score_1[0]
                                            print score_1[1]
                                            if((score_0[2] >= len(seq_a) - DIST_MAX) or (score_1[2] >= len(seq_a) - DIST_MAX)):
                                                dict_neighbor_of[pID_r[1]] = pID_a[1]
                                                dict_neighbor_state_pIDs[pID_r[1]] = True
                                                print "NEIGHBORS !"
                                    else:
                                        print 'PB: pID_a does not have the good number of occurrences !!!'

                
            dict_neighbors = defaultdict(list)
            for pID in dict_neighbor_of.keys():
                dict_neighbors[dict_neighbor_of[pID]].append(pID)
        
            print 'neighbors:'
            nb_neighbors_stated = 0
            for pID in dict_neighbor_state_pIDs.keys():
                if dict_neighbor_state_pIDs[pID]==True:
                    nb_neighbors_stated += 1
                    print pID

        
        
            ##############################################
            # -- WRITE THE NEIGHBORS FILTERED READS FILE #
            ##############################################
            count_reads_written = 0
            count_neighbors_reads_written = 0
        
            # write the new filtered reads file (mutant-indels):
            filtered_reads_file_basename_suffix = filtered_reads_file_basename.split('_')[-1].split('.')[0]
            with open(str(path_to_data+prefix_date_and_id+'-mutant-indels'+'_'+barcode+'_'+filtered_reads_file_basename_suffix+'.fasta'), 'w') as neighbors_filtered_reads_file:
                for pID in dict_all_reads.keys(): # for all pIDs                                
                    if not dict_neighbor_state_pIDs[pID]: # if pID "is not a neighbor"
                        print 'write reads for pID: ' + pID

                        # write the reads corresponding to pID
                        for cur_read in dict_all_reads[pID]: # cur_read is a list [nb_occ_fwd,nb_occ_rev,seq]
                            pid_orig = pID
                            new_read_id = str(pID+'_'+str(cur_read[0])+'_'+str(cur_read[1])+'_'+pid_orig)
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
            print 'PB: "barcode" and "prefix_date_and_id" extracted from input filtered reads file and consensus sequences file are different, Stop here!'
    else:
        print auto_file_name + ': usage: '+ auto_file_name + ' <filtered reads file (in ../data/)> <cons. seq. file (in ../templates/dir-<date_and_id>_consensus/)>'
