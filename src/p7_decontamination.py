#!/usr/bin/python
#/!\#
#### /ebio/ag-neher/share/pograms/EPD/bins/python
#/!\#

# script that decontaminates a filtered reads file for a given barcode and run id
# input: - filtered reads file to decontaminate (in ../data/)
#        - consensus sequences file to check contamination (in ../templates/<date-and-id>_consensus/)
#        - reference sequences file to check contamination (in ../data/)
#        - run id (int) corresponding to the associated filtered reads and consensus sequences files
#
# output:- decontaminated filtered reads file (in ../data/)
#        - filtered reads file containing reads that align badly with all the reference sequences (in ../data/)
 
import numpy as np
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from collections import Counter
from collections import defaultdict
import sys
import datetime
import time
import lib_tools as lt


auto_file_name=str(sys.argv[0])

DIST_MAX = 7

if (len(sys.argv)==5):

    # parse the input file names
    barcode_dir = str(sys.argv[1]).rstrip('/')+'/'
    ref_seq_file_name = str(sys.argv[2])
    readtype = sys.argv[3]
    true_seq_id = sys.argv[4]

    # parse the reference sequences file in dict_ref_seq

    dict_ref_seq ={}
    with open(ref_seq_file_name, 'r') as infile:
        for seq in SeqIO.parse(infile, 'fasta'):
            dict_ref_seq[seq.id]=str(seq.seq)


    #2.parse the consensus sequences files for the given prefix_date_and_id in dict_cons_seq
    dict_cons_seq, a,b = lt.parse_readfile(barcode_dir+'consensus_sequence_'+readtype+'.fasta')

    dict_all_reads, nreads, nbadreads = lt.parse_readfile(barcode_dir+'aligned_reads_'+readtype+'.fasta')

    dict_score_alignments = defaultdict(list)
    for pID in dict_all_reads:
        if pID in dict_cons_seq:
            score_align = lt.align_SW_dist(dict_ref_seq[true_seq_id], dict_cons_seq[pID][2])[2]
            print '----'+str(RUN)+': '+barcode+' - '+pID+': '+pID+', align score: '+str(score_align)
            dict_score_alignments[pID].append(score_align)
        else:
            for nf, nb, seq in dict_all_reads[pID]:
                score_align = lt.align_SW_dist(dict_ref_seq[true_seq_id], dict_cons_seq[pID][2])[2]
                print '----'+str(RUN)+': '+barcode+' - '+pID+': '+pID+', align score: '+str(score_align)
                dict_score_alignments[pID].append(score_align)


    #4.print some statistics
    #for bc in barcodes[prefix_date_and_id[RUN-1]]:
    print '** LIM SCORE ALIGN-'+true_seq_id
    print '---- min = '+str(min([dict_score_alignments[i] for i in dict_score_alignments.keys()]))
    print '---- max = '+str(max([dict_score_alignments[i] for i in dict_score_alignments.keys()]))
    average_scores=np.mean([dict_score_alignments[i] for i in dict_score_alignments.keys()])
    print '---- average score = '+str(average_scores)
    print '---- #pids (total) = '+str(len(dict_all_reads.keys()))


    #5.count good and bad alignments
    # for bc in barcodes[prefix_date_and_id[RUN-1]]:

    count_bad_scores = sum([x<len(ref_seq[true_seq_id])-DIST_MAX for x in reads for n dict_score_alignments])
    nb_bad_scores_reads = sum([x<len(ref_seq[true_seq_id])-DIST_MAX for x in reads for reads in dict_score_alignments])
    count_good_scores = sum([x>=len(ref_seq[true_seq_id])-DIST_MAX for x in reads for reads in dict_score_alignments])
    nb_good_scores_reads = sum([x<len(ref_seq[true_seq_id])-DIST_MAX for x in reads for reads in dict_score_alignments])

    count_bad_scores_once_twice = 0
    nb_bad_scores_reads_once_twice = 0
    count_good_scores_once_twice = 0
    nb_good_scores_reads_once_twice = 0

    for pID in dict_score_alignments.keys():
        if (dict_score_alignments[pID] <  dict_len_ref_seq[(RUN,barcode)]-DIST_MAX):
            count_bad_scores += 1
            nb_bad_scores_reads += dict_pIDs_occ[pID]
        else:
            count_good_scores += 1
            nb_good_scores_reads += dict_pIDs_occ[pID]

    for pID in dict_score_alignments_once_twice.keys():
        if (max(dict_score_alignments_once_twice[pID]) <  dict_len_ref_seq[(RUN,barcode)]-DIST_MAX):
            count_bad_scores_once_twice += 1
            nb_bad_scores_reads_once_twice += dict_pIDs_occ_once_twice[pID]
        else:
            count_good_scores_once_twice += 1
            nb_good_scores_reads_once_twice += dict_pIDs_occ_once_twice[pID]

    print 'PIDS OCC >= 3:'
    print '--RUN '+str(RUN)+': '+barcode+': #bad align scores: '+str(count_bad_scores)
    print '--RUN '+str(RUN)+': '+barcode+': #bad align scores (reads): '+str(nb_bad_scores_reads)
    print '--RUN '+str(RUN)+': '+barcode+': #good align scores: '+str(count_good_scores)
    print '--RUN '+str(RUN)+': '+barcode+': #good align scores (reads): '+str(nb_good_scores_reads)

    print '1 <= PIDS OCC <= 2:'
    print '--RUN '+str(RUN)+': '+barcode+': #bad align scores: '+str(count_bad_scores_once_twice)
    print '--RUN '+str(RUN)+': '+barcode+': #bad align scores (reads): '+str(nb_bad_scores_reads_once_twice)
    print '--RUN '+str(RUN)+': '+barcode+': #good align scores: '+str(count_good_scores_once_twice)
    print '--RUN '+str(RUN)+': '+barcode+': #good align scores (reads): '+str(nb_good_scores_reads_once_twice)
    print '****'


    # 6.check contamination
    print '********'
    print 'Check contamination'
    dict_possible_good_align={}#contains the alignment scores of the pIDs cons. seq. for bc with the ref. seq. of the other barcodes
    dict_reclassif={}
    dict_count_reclassif_pIDs=defaultdict(int)
    dict_count_reclassif_reads=defaultdict(int)
    #for bc in barcodes[prefix_date_and_id[RUN-1]]:
    print 'RUN '+str(RUN)+' - barcode '+barcode+':'
    #dict_possible_good_align={}
    #dict_reclassif={}
    #dict_count_reclassif_pIDs=defaultdict(int)
    #dict_count_reclassif_reads=defaultdict(int)
    index_bc = barcodes[prefix_date_and_id[RUN-1]].index(bc)#index of the barcode bc in the barcodes list
    #barcodes_to_test = barcodes[prefix_date_and_id[RUN-1]][:index_bc]+barcodes[prefix_date_and_id[RUN-1]][index_bc+1:]#list of barcodes - {bc}
    barcodes_to_test = [i for i in dict_ref_seq.keys() if (i!=(RUN,barcode))]

    # loop for pIDs with occ >= 3
    print 'loop for pIDs with occ >= 3'
    for pID in dict_score_alignments.keys():
        if (dict_score_alignments[pID] < dict_len_ref_seq[(RUN,barcode)]-DIST_MAX):
            # if the alignment score of cons. seq. for pID with the reference sequence for bc is lower than the threshold
            print '---found 1 bad alignment: '+pID+', '+str(dict_score_alignments[pID])+' < '+str(dict_len_ref_seq[(RUN,barcode)]-DIST_MAX)
            dict_possible_good_align[pID]=[]
            dict_reclassif[pID]=[]

            for (run2,bc2) in barcodes_to_test: # compare the alignment scores of cons. seq. pID with the ref. seq. for the other barcodes
                score_align_with_bc2 = lt.align_SW_dist(dict_ref_seq[(run2,bc2)], dict_cons_seq[pID])[2]
                dict_possible_good_align[pID].append((run2,bc2,score_align_with_bc2))
                if(score_align_with_bc2 >= dict_len_ref_seq[(run2,bc2)]-DIST_MAX):
                    print '------possible reclassif of pID '+pID+' (#occ: '+str(dict_pIDs_occ[pID])+') in run '+str(run2) +', bc2 '+bc2+', with align score: '+str(score_align_with_bc2)+' >= '+str(dict_len_ref_seq[(run2,bc2)]-DIST_MAX)
                    dict_count_reclassif_pIDs[(run2,bc2)]+=1
                    dict_count_reclassif_reads[(run2,bc2)]+=dict_pIDs_occ[pID]
                    dict_reclassif[pID].append((run2,bc2,score_align_with_bc2))

            #for bc3 in barcodes[prefix_date_and_id[RUN%2]]: # compare the alignment scores of cons. seq. pID with the ref. seq. for the barcodes of the other run
            #    score_align_with_bc3 = align_SW_dist(dict_scsp[((RUN%2)+1,bc3)], dict_cons_seq[(RUN,bc)][pID])[2]
            #    dict_possible_good_align[(RUN,bc)][pID].append(((RUN%2)+1,bc3,score_align_with_bc3))
            #    if(score_align_with_bc3 >= dict_lim_min_score_align[((RUN%2)+1, bc3)]):
            #        print '------possible reclassif of pID '+pID+' (#occ: '+str(dict_pIDs_occ[(RUN,bc)][pID])+') in run '+str((RUN%2)+1)+', bc3 '+bc3+', with align score: '+str(score_align_with_bc3)+' >= '+str(dict_lim_min_score_align[((RUN%2)+1, bc3)])
            #        dict_count_reclassif_pIDs[(RUN,bc)][((RUN%2)+1,bc3)]+=1
            #        dict_count_reclassif_reads[(RUN,bc)][((RUN%2)+1,bc3)]+=dict_pIDs_occ[(RUN,bc)][pID]
            #        dict_reclassif[(RUN,bc)][pID].append(((RUN%2)+1,bc3,score_align_with_bc3))


    # loop for pIDs with 1 <= occ <= 2
    print 'loop for pIDs with 1 <= occ <=2'
    for pID in dict_score_alignments_once_twice.keys():
        if (max(dict_score_alignments_once_twice[pID]) < dict_len_ref_seq[(RUN,barcode)]-DIST_MAX):
            # if the alignment score (occ=1) or the best one (occ=2) for pID with the reference sequence for bc is lower than the threshold
            print '---found 1 bad alignment: '+pID+', '+str(max(dict_score_alignments_once_twice[pID]))+' < '+str(dict_len_ref_seq[(RUN,barcode)]-DIST_MAX)
            dict_possible_good_align[pID]=[]
            dict_reclassif[pID]=[]
            #nb_occ_cur_pID = len(dict_filtered_reads_pIDs_once_twice[(RUN,bc)][pID])      
            nb_occ_cur_pID = sum([dict_pIDs_occ_once_twice[pID][i] for i in [0,1]])

            for bc2 in barcodes_to_test: # compare the alignment scores of the sequence(s) of pID with the ref. seq. for the other barcodes
                score_align_with_bc2 = max([lt.align_SW_dist(dict_ref_seq[(run2,bc2)], dict_filtered_reads_pIDs_once_twice[pID][i])[2] for i in range(nb_occ_cur_pID)])
                dict_possible_good_align[pID].append((run2,bc2,score_align_with_bc2))
                if(score_align_with_bc2 >= dict_len_ref_seq[(run2,bc2)]-DIST_MAX):
                    print '------possible reclassif of pID '+pID+' (#occ: '+str(nb_occ_cur_pID)+') in run '+str(run2) +', bc2 '+bc2+', with align score: '+str(score_align_with_bc2)+' >= '+str(dict_len_ref_seq[(run2,bc2)]-DIST_MAX)
                    dict_count_reclassif_pIDs[(run2,bc2)]+=1
                    dict_count_reclassif_reads[(run2,bc2)]+=nb_occ_cur_pID
                    dict_reclassif[pID].append((run2,bc2,score_align_with_bc2))


            #for bc3 in barcodes[prefix_date_and_id[RUN%2]]: # compare the alignment scores of the sequence(s) of pID with the ref. seq. for the barcodes of the other run
            #    score_align_with_bc3 = max([align_SW_dist(dict_scsp[((RUN%2)+1,bc3)], dict_filtered_reads_pIDs_once_twice[(RUN,bc)][pID][i])[2] for i in range(nb_occ_cur_pID)])
            #    dict_possible_good_align[(RUN,bc)][pID].append(((RUN%2)+1,bc3,score_align_with_bc3))
            #    if(score_align_with_bc3 >= dict_lim_min_score_align[((RUN%2)+1, bc3)]):
            #        print '------possible reclassif of pID '+pID+' (#occ: '+str(nb_occ_cur_pID)+') in run '+str((RUN%2)+1)+', bc3 '+bc3+', with align score: '+str(score_align_with_bc3)+' >= '+str(dict_lim_min_score_align[((RUN%2)+1, bc3)])
            #        dict_count_reclassif_pIDs[(RUN,bc)][((RUN%2)+1,bc3)]+=1
            #        dict_count_reclassif_reads[(RUN,bc)][((RUN%2)+1,bc3)]+=dict_pIDs_occ_once_twice[(RUN,bc)][pID]
            #        dict_reclassif[(RUN,bc)][pID].append(((RUN%2)+1,bc3,score_align_with_bc3))



    #7.print contamination statistics

    # for bc in barcodes[prefix_date_and_id[RUN-1]]:
    count_pIDs_reclassif = 0
    count_pIDs_no_good_alignments = 0
    count_reads_no_good_alignments = 0
    for pID in dict_reclassif.keys():
        if len(dict_reclassif[pID])>0:
            count_pIDs_reclassif+=1
        else:
            count_pIDs_no_good_alignments +=1
            if pID in dict_pIDs_occ.keys():
                count_reads_no_good_alignments += dict_pIDs_occ[pID]
            else:
                count_reads_no_good_alignments += sum([dict_pIDs_occ_once_twice[pID][i] for i in [0,1]]) #dict_pIDs_occ_once_twice[(RUN,bc)][pID]

    print 'RUN '+str(RUN)+' - barcode '+barcode+': #pIDs with better alignments from other barcodes: '+str(count_pIDs_reclassif)
    print '-----> #pIDs with no good alignments: '+ str(count_pIDs_no_good_alignments)
    print '-----> #reads with no good alignments: '+ str(count_reads_no_good_alignments)
    print '##########'
    #for bc in barcodes[prefix_date_and_id[RUN-1]]:
    #print 'RUN '+str(RUN)+' - barcode '+barcode+':'
    #print '--> #pIDs from:  '+str(dict_count_reclassif_pIDs[(RUN,bc)])
    #print '--> #reads from: '+str(dict_count_reclassif_reads[(RUN,bc)])



    #8.decontaminate the filtered reads (neighbors-indels) file
    #dict_count_pIDs_to_bad_aligned_file={}
    #dict_count_pIDs_from_contamination={}
    #dict_count_pIDs_to_decontaminated_file={}

    # for bc in barcodes[prefix_date_and_id[RUN-1]]:
    dict_count_pIDs_to_bad_aligned_file=defaultdict(int)
    dict_count_pIDs_from_contamination=defaultdict(int)
    dict_count_pIDs_to_decontaminated_file=defaultdict(int)
    list_of_reads_to_write_in_bad_aligned_file=[]
    list_of_reads_to_write_in_decontaminated_file=[]
    with open(relative_filtered_reads_file_path, 'r') as infile:
    #with open(str('../data/'+prefix_date_and_id[RUN-1]+'_'+bc+'_sorted_and_filtered_reads_SW_min_occ_same_pid_1_neighbors_indels.fa'), 'r') as infile:
        for record in SeqIO.parse(infile,'fasta'):
            pID = record.id.split('_')[0]
            if (pID in dict_reclassif.keys()): # if the seq. for pID aligns badly with the ref. seq. of its barcode 
                if not (len(dict_reclassif[pID])>0): # and if it does not align well with the ref. seq. of the other barcodes (run 1 and run 2)
                    # don't know what to do with them: keep them in a file
                    list_of_reads_to_write_in_bad_aligned_file.append(record)
                    dict_count_pIDs_to_bad_aligned_file[pID]+=1
                else: # it aligns well with the ref. seq. of an other barcode (run 1 or run 2)
                    # just throw it away
                    dict_count_pIDs_from_contamination[pID]+=1
            else: # the seq. for pID aligns well with the ref. seq. of its barcode
                # to be written in the filtered reads (neighbors_indels) decontaminated file
                list_of_reads_to_write_in_decontaminated_file.append(record)
                dict_count_pIDs_to_decontaminated_file[pID]+=1

    # write the filtered reads (neighbors_indels) decontaminated file
    suffix_filtered_reads_file = filtered_reads_file_basename.split('_')[-1]
    with open(str(path_to_data+prefix_date_and_id+'-decontaminated'+'_'+barcode+suffix_filtered_reads_file), 'w') as outfile:
        print 'write decontaminated file for run '+str(RUN)+', barcode: '+barcode
        for record in list_of_reads_to_write_in_decontaminated_file:
            outfile.write('>'+str(record.id)+'\n')
            outfile.write(str(record.seq)+'\n')
        print '-- #reads written: '+str(len(list_of_reads_to_write_in_decontaminated_file))

    # write the (neighbors_indels) bad aligned reads file
    with open(str(path_to_data+prefix_date_and_id+'-bad-aligned'+'_'+barcode+suffix_filtered_reads_file), 'w') as outfile:
        print 'write bad aligned reads file for run '+str(RUN)+', barcode: '+bc
        for record in list_of_reads_to_write_in_bad_aligned_file:
            outfile.write('>'+str(record.id)+'\n')
            outfile.write(str(record.seq)+'\n')
        print '-- #reads written: '+str(len(list_of_reads_to_write_in_bad_aligned_file))

    # 9. print statistics
    print '#############'
    #for bc in barcodes[prefix_date_and_id[RUN-1]]:
    print '-RUN '+str(RUN)+', barcode '+bc
    print '---> nb non contaminant pIDs: '+str(len(dict_count_pIDs_to_decontaminated_file.keys()))
    print '---> nb non contaminant reads: '+str(sum([dict_count_pIDs_to_decontaminated_file[i] for i in dict_count_pIDs_to_decontaminated_file.keys()]))
    print '--'
    print '---> nb contaminant pIDs: '+str(len(dict_count_pIDs_from_contamination.keys()))
    print '---> nb contaminant reads: '+str(sum([dict_count_pIDs_from_contamination[i] for i in dict_count_pIDs_from_contamination.keys()]))
    print '--'
    print '---> nb bad aligned but non contaminant pIDs : '+str(len(dict_count_pIDs_to_bad_aligned_file.keys()))
    print '---> nb bad aligned but non contaminant reads: '+str(sum([dict_count_pIDs_to_bad_aligned_file[i] for i in dict_count_pIDs_to_bad_aligned_file.keys()]))


else:
    print auto_file_name+': usage: '+auto_file_name+' <filtered reads file (in ../data/)> <consensus seq. file (in ../templates/<date-and-id>_consensus/)> <ref. seq. file (in ../data/)> <run_number>'


 
    
