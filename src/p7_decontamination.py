#!/usr/bin/python
#/!\#
#### /ebio/ag-neher/share/pograms/EPD/bins/python
#/!\#

# script that decontaminates a filtered reads file for a given barcode and run id

 
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
            dict_ref_seq[seq.description]=str(seq.seq)


    #2.parse the consensus sequences files for the given prefix_date_and_id in dict_cons_seq
    dict_cons_seq, a,b = lt.parse_readfile(barcode_dir+'consensus_sequences_'+readtype+'.fasta')
    
    
    #/!\ : aligned reads: need to trim the gaps at one point
    dict_all_reads, nreads, nbadreads = lt.parse_readfile(barcode_dir+'aligned_reads_'+readtype+'.fasta')

    dict_score_alignments = defaultdict(list)
    for pID in dict_all_reads:
        if pID in dict_cons_seq:
            #print 'FLAG:'
            #print dict_cons_seq[pID]
            tmp_cons_seq = dict_cons_seq[pID][0][2]
            score_align = lt.align_SW_dist(dict_ref_seq[true_seq_id], tmp_cons_seq)[2]
            print '----'+str(true_seq_id)+': '+pID+', align score: '+str(score_align)
            dict_score_alignments[pID].append(score_align)
        else:
            for nf, nb, seq in dict_all_reads[pID]:
                score_align = lt.align_SW_dist(dict_ref_seq[true_seq_id], lt.remove_gaps(seq))[2]
                print '----'+str(true_seq_id)+': '+pID+', align score: '+str(score_align)
                dict_score_alignments[pID].append(score_align)


    #4.print some statistics
    print '** LIM SCORE ALIGN-'+true_seq_id
    #print '---- min = '+str(min([dict_score_alignments[i] for i in dict_score_alignments.keys()]))
    print '---- min = '+str(min(np.hstack(dict_score_alignments.values())))
    #print '---- max = '+str(max([dict_score_alignments[i] for i in dict_score_alignments.keys()]))
    print '---- min = '+str(max(np.hstack(dict_score_alignments.values())))
    #average_scores=np.mean([dict_score_alignments[i] for i in dict_score_alignments.keys()])
    average_scores=np.mean(np.hstack(dict_score_alignments.values()))
    print '---- average score = '+str(average_scores)
    print '---- #pids (total) = '+str(len(dict_all_reads.keys()))


    #5.count good and bad alignments
    

    #count_bad_scores = sum([sum([x<len(dict_ref_seq[true_seq_id])-DIST_MAX for x in reads]) for reads in dict_score_alignments])
    #nb_bad_scores_reads = sum([sum([x<len(dict_ref_seq[true_seq_id])-DIST_MAX for x in reads]) for reads in dict_score_alignments])
    #count_good_scores = sum([sum([x>=len(dict_ref_seq[true_seq_id])-DIST_MAX for x in reads]) for reads in dict_score_alignments])
    #nb_good_scores_reads = sum([sum([x<len(dict_ref_seq[true_seq_id])-DIST_MAX for x in reads]) for reads in dict_score_alignments])

                                                                         
    #print '****'
    #print '--RUN-barcode: '+str(true_seq_id)+': #bad align scores: '+str(count_bad_scores)
    #print '--RUN-barcode: '+str(true_seq_id)+': #bad align scores (reads): '+str(nb_bad_scores_reads)
    #print '--RUN-barcode: '+str(true_seq_id)+': #good align scores: '+str(count_good_scores)
    #print '--RUN-barcode: '+str(true_seq_id)+': #good align scores (reads): '+str(nb_good_scores_reads)

 
    # 6.check contamination
    print '********'
    print 'Check contamination'
    dict_possible_good_align={}#contains the alignment scores of the pIDs cons. seq. for bc with the ref. seq. of the other barcodes
    dict_reclassif={}
    dict_count_reclassif_pIDs=defaultdict(int)
    dict_count_reclassif_reads=defaultdict(int)
    print 'RUN - barcode: '+str(true_seq_id)+':'
    #index_bc = barcodes[prefix_date_and_id[RUN-1]].index(bc)#index of the barcode bc in the barcodes list
    barcodes_to_test = [ref_seq_bc for ref_seq_bc in dict_ref_seq.keys() if (ref_seq_bc != true_seq_id)]

    


    for pID,align_scores in dict_score_alignments.iteritems():
        if(max(align_scores) < len(dict_ref_seq[true_seq_id]) - DIST_MAX):
            # if the alignment score of cons. seq. for pID with the reference sequence for bc is lower than the threshold
            print '----found 1 bad alignment: '+pID+', '+str(dict_score_alignments[pID])+' < '+str(len(dict_ref_seq[true_seq_id])-DIST_MAX)
            dict_possible_good_align[pID]=[]
            dict_reclassif[pID]=[]
            is_a_thrice_or_more = pID in dict_cons_seq.keys()
            is_a_twice_with_distinct_reads = len(align_scores)==2
            for ref_seq_id_2 in barcodes_to_test:# compare the alignment scores of cons. seq. with the ref. seq. for the other barcodes
                if is_a_thrice_or_more:
                    list_seq = [dict_cons_seq[pID][0][2]]
                else:
                    list_seq = [seq[2] for seq in dict_all_reads[pID]]
                score_align_with_ref_seq_id_2 = max([lt.align_SW_dist(dict_ref_seq[ref_seq_id_2], lt.remove_gaps(seq))[2] for seq in list_seq])

                dict_possible_good_align[pID].append((ref_seq_id_2,score_align_with_ref_seq_id_2))
                if(score_align_with_ref_seq_id_2 >= len(dict_ref_seq[ref_seq_id_2])-DIST_MAX):
                    print '------possible reclassif of pID '+pID+' in ref_seq_id_2 '+str(ref_seq_id_2)+', with align score: '+str(score_align_with_ref_seq_id_2)+' >= '+str(len(dict_ref_seq[ref_seq_id_2])-DIST_MAX)
                    dict_count_reclassif_pIDs[ref_seq_id_2]+=1
                    dict_count_reclassif_reads[ref_seq_id_2]+=sum([sum(map(int,np.array(dict_all_reads[pID])[:,i])) for i in [0,1]]) #dict_pIDs_occ[pID]
                    dict_reclassif[pID].append((ref_seq_id_2,score_align_with_ref_seq_id_2))

  

    # 7.print contamination statistics
    # 8.decontaminate the filtered reads file
    dict_count_pIDs_to_bad_aligned_file=defaultdict(int)
    dict_count_pIDs_from_contamination=defaultdict(int)
    dict_count_pIDs_to_decontaminated_file=defaultdict(int)
    list_of_reads_to_write_in_bad_aligned_file=[]
    list_of_reads_to_write_in_decontaminated_file=[]
    list_of_reads_from_contamination=[]
    count_reads_no_good_alignments = 0
    count_reads_good_alignments = 0

    with open(barcode_dir+'aligned_reads_'+readtype+'.fasta','r') as aligned_reads_file:
        for record in SeqIO.parse(aligned_reads_file,'fasta'):
            pID,nb_reads_fwd,nb_reads_rev,pID_orig = record.id.split('_')[0:4]
            nb_reads_fwd = int(nb_reads_fwd)
            nb_reads_rev = int(nb_reads_rev)
            read = record.seq
            if pID in dict_reclassif.keys(): # bad alignment with the associated ref. seq.
                if len(dict_reclassif[pID])>0: # contamination from an other run/barcode
                    list_of_reads_from_contamination.append(record)
                    dict_count_pIDs_from_contamination[pID]+=(nb_reads_fwd+nb_reads_rev)
                else: # bad alignment with all the ref. seq.
                    list_of_reads_to_write_in_bad_aligned_file.append(record)
                    dict_count_pIDs_to_bad_aligned_file[pID]+=(nb_reads_fwd+nb_reads_rev)
                    count_reads_no_good_alignments += (nb_reads_fwd+nb_reads_rev)
            else: # good alignment with the associated ref. seq.
                list_of_reads_to_write_in_decontaminated_file.append(record)
                dict_count_pIDs_to_decontaminated_file[pID]+= (nb_reads_fwd+nb_reads_rev)
                count_reads_good_alignments += (nb_reads_fwd+nb_reads_rev)
    
        
   # print 'RUN-barcode '+true_seq_id+': #pIDs with better alignments from other barcodes: '+str(count_pIDs_reclassif)
   # print '-----> #pIDs with no good alignments: '+ str(count_pIDs_no_good_alignments)
    print '-----> #reads with no good alignments: '+ str(count_reads_no_good_alignments)
    print '-----> #reads with good alignments: '+ str(count_reads_good_alignments)
    print '##########'

    # write the decontaminated aligned filtered reads file
    with open(barcode_dir+'aligned_reads_'+readtype+'_decontaminated.fasta','w') as outfile:
        print 'write decontaminated file for run-barcode: '+str(true_seq_id)
        for record in list_of_reads_to_write_in_decontaminated_file:
            outfile.write('>'+str(record.id)+'\n')
            outfile.write(str(record.seq)+'\n')
        print '-- #reads written: '+str(len(list_of_reads_to_write_in_decontaminated_file))
                
    # write the bad aligned filtered reads file
    with open(barcode_dir+'aligned_reads_'+readtype+'_bad_aligned.fasta','w') as outfile:
        print 'write bad aligned reads file for run-barcode: '+str(true_seq_id)
        for record in list_of_reads_to_write_in_bad_aligned_file:
            outfile.write('>'+str(record.id)+'\n')
            outfile.write(str(record.seq)+'\n')
        print '-- #reads written: '+str(len(list_of_reads_to_write_in_bad_aligned_file))

    # write the consensus sequences file
    with open(barcode_dir+'consensus_sequences_'+readtype+'_decontaminated.fasta','w') as outfile:
        print 'write decontaminated file for run-barcode: '+str(true_seq_id)
        for pID,rec in dict_cons_seq.iteritems():
            if pID not in dict_reclassif:
                outfile.write('>'+lt.read_label(pID, rec[0][0], rec[0][1])+'\n')
                outfile.write(rec[0][2]+'\n')


    # 9. print decontamination statistics
    print '#############'
    print '-RUN-barcode '+str(true_seq_id)
    print '---> nb non contaminant pIDs: '+str(len(dict_count_pIDs_to_decontaminated_file.keys()))
    print '---> nb non contaminant reads: '+str(sum(dict_count_pIDs_to_decontaminated_file.values()))
    print '--'
    print '---> nb contaminant pIDs: '+str(len(dict_count_pIDs_from_contamination.keys()))
    print '---> nb contaminant reads: '+str(sum(dict_count_pIDs_from_contamination.values()))
    print '--'
    print '---> nb bad aligned but non contaminant pIDs : '+str(len(dict_count_pIDs_to_bad_aligned_file.keys()))
    print '---> nb bad aligned but non contaminant reads: '+str(sum(dict_count_pIDs_to_bad_aligned_file.values()))


else:
    #print auto_file_name+': usage: '+auto_file_name+' <filtered reads file (in ../data/)> <consensus seq. file (in ../templates/<date-and-id>_consensus/)> <ref. seq. file (in ../data/)> <run_number>'
    print auto_file_name+': usage: '+auto_file_name+' <barcode_dir> <ref_seq_file_name> read_type true_seq_id'


 
    
