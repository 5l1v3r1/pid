#!/usr/bin/python
#/!\#
#### /ebio/ag-neher/share/pograms/EPD/bins/python
#/!\#

## script that filters the raw file (fasta or fastq, Ion Torrent or 454) and gather the reads in separated files (one for each barcode)
## input: raw file containing reads to filter (in ../data/)
## aux: configFile.py (in current directory)
## output: filtered reads files (one per barcode) (in ../data/)


import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.pairwise2 import align
from collections import Counter
from collections import defaultdict
import matplotlib.pyplot as plt
import sys
import datetime
import time
import lib_tools as lt

CONFIG_FILE_NAME='configFile.py'

#####
# DEF CLASSES
#####

#####
class struct_var_set:
    def __init__(self):
        #if(CONFIG_FILE_NAME in sys.argv):
        try:
            execfile(CONFIG_FILE_NAME)
            if ('runid' in cfg.keys() and 'p5_virus_match' in cfg.keys() 
                and 'p3_virus_match' in cfg.keys() and 'input_data_file' in cfg.keys() 
                and 'barcodes' in cfg.keys() and 'reverse' in cfg.keys()):

                # MANDATORY variables in the configuration file:
                self.runid = cfg['runid']
                self.p5_virus_match = cfg['p5_virus_match']
                self.p3_virus_match = cfg['p3_virus_match']
                self.match_seqs = [self.p5_virus_match, self.p3_virus_match]
                self.reverse = cfg['reverse']
                #self.p3_adaptor = cfg['p3_adaptor']
                self.barcodes = cfg['barcodes']
                self.input_data_file = cfg['input_data_file']
                # OPTIONAL variables in the configuration file:
                if('p5_cutoff' in cfg.keys()):
                    self.p5_cutoff = cfg['p5_cutoff']
                else:
                    self.p5_cutoff = len(self.p5_virus_match)-3
                if('p3_cutoff' in cfg.keys()):
                    self.p3_cutoff = cfg['p3_cutoff']
                else:
                    self.p3_cutoff = len(self.p3_virus_match)-3
                if('min_occ_same_pid' in cfg.keys()):
                    self.min_occ_same_pid = cfg['min_occ_same_pid']
                else:
                    self.min_occ_same_pid = 3 # must be > 2 (for consensus sequences)
                if('min_length_pid' in cfg.keys()):
                    self.min_length_pid = cfg['min_length_pid']
                else:
                    self.min_length_pid = 8
                if 'barcode_length_range' in cfg:
                    self.bc_length_range = cfg['barcode_length_range']
                else:
                    self.bc_length_range = set(len, self.barcodes)
                # ADDITIONAL variables for trimming and filtering:
                self.count = 0
                
                self.good_reads = {bc:{} for bc in self.barcodes}                
                self.good_reads_2 = {bc:{} for bc in self.barcodes}
                self.good_read_length = {bc:np.zeros(800) for bc in self.barcodes}
                self.good_fwd_read_length = {bc:np.zeros(800) for bc in self.barcodes}
                self.good_tag_only_read_length = {bc:np.zeros(800) for bc in self.barcodes}
                self.bad_read_length = np.zeros(800)
                self.bad_read_length_454 = {}
                self.good_read_bad_pID_length={bc:np.zeros(800) for bc in self.barcodes}
                #self.templates={}
            else:
                print 'Problem: at least one of the following (mandatory) parameters is missing in ' + CONFIG_FILE_NAME + ':'
                print 'runid, p5_virus_match, p3_virus_match, p3_adaptor, barcodes, input_data_file, reverse.' 
        #else:
        except:
            print 'The configuration file ' + CONFIG_FILE_NAME + ' is missing or mistaken !'



#####
# DEF FUNCTIONS
#####

#####
def retrieve_pID_barcode_and_sequence(res, fwd, scores, seq, seq_id):
    '''
    function that retrieves and stores pID, barcode and sequence for good reads
    res is the structure containing the filtering variables
    fwd is a boolean indicating the direction (True = forward, False = backward)
    scores is the scores obtained for the local alignments of primers
    seq is the sequence (reversed if fwd == True)
    seq_id is the sequence name
    '''

    # retrieve the barcode (3 letters upstream the 5' primer)
    if (res.reverse): # 454
        bc = scores[0][0][scores[0][3]-3:scores[0][3]]
    else: # Ion Torrent
        bc = scores[0][0][:4]
    # print bc
    # retrieve the pID (all the letters downstream the 3' primer)
    pID = scores[1][0][scores[1][4]:]
    # print pID
    L = len(seq)

    if (len(pID) >= res.min_length_pid):
        pID = pID[:res.min_length_pid]
        if (bc in res.barcodes):
            if not (pID in res.good_reads[bc].keys()):
                res.good_reads[bc][pID] = []
            
            len_align_5p = len(scores[0][0])
            pos_end_align_5p_plus_one = scores[0][4]
            pos_start_seq = - (len_align_5p - pos_end_align_5p_plus_one)
            #len_align_3p = len(scores[0][0])
            pos_start_align_3p = scores[1][3]
            pos_end_seq = pos_start_align_3p 

            res.good_reads[bc][pID].append([seq_id, fwd==True and 'fwd' or 'rev', str(seq[pos_start_seq:pos_end_seq])])
            res.good_read_length[bc][L] += 1

    else: # pID too small to be valid
        if(bc in res.barcodes):
            res.good_read_bad_pID_length[bc][L] += 1
            
            
#####
def filter_reads(res):
    time_start = time.time()
    with open(str(res.input_data_file), 'r') as seq_file:
        file_format = res.input_data_file.split('.')[-1]
        for record in SeqIO.parse(seq_file, file_format):
            tmp_seq = str(record.seq)
            L = len(tmp_seq)
            
            if(res.count%500==0):
                print str(res.count) + ' reads analyzed'
            
            # try forward
            prim_scores = lt.check_match_SW(tmp_seq, res.match_seqs)
            if (prim_scores[0][2] >= res.p5_cutoff and   #check 5 prime match
                prim_scores[0][3] in res.bc_length_range and              #check barcode length
                prim_scores[1][2] >= res.p3_cutoff       #check 3 prime match
                ):
                retrieve_pID_barcode_and_sequence(res, True, prim_scores, tmp_seq, record.id)
            else:
                if(res.reverse):
                    # try backward
                    tmp_seq = str(Seq(tmp_seq, generic_dna).reverse_complement())
                    prim_scores = lt.check_match_SW(tmp_seq, res.match_seqs)
                    if (prim_scores[0][2] >= res.p5_cutoff and
                        prim_scores[0][3] in res.bc_length_range and             #check barcode length
                        prim_scores[1][2] >= res.p3_cutoff
                        ):
                        retrieve_pID_barcode_and_sequence(res, False, prim_scores, tmp_seq, record.id)
                    else:# bad read
                        if not (L in res.bad_read_length_454.keys()):
                            res.bad_read_length_454[L] = 0
                        res.bad_read_length_454[L] += 1
                else:# bad read
                    if not (L in res.bad_read_length_454.keys()):
                        res.bad_read_length_454[L] = 0
                    res.bad_read_length_454[L] += 1
                    
            res.count+=1
    time_end = time.time()
    print 'computing time: ' + str(time_end - time_start)
    
    for bc in sorted(res.good_reads.keys()):
        print "barcode " + str(bc) + ":"
        #print res.good_reads[bc]
        print "---#good_reads: " + str(len(res.good_reads[bc]))
        
    print 'Total:'
    print res.count, [len(res.good_reads[bc]) for bc in res.barcodes], [int(np.sum(res.good_read_bad_pID_length[bc])) for bc in res.barcodes]

#####
def logfile_output(res, barcode):
    '''
    function that writes a summary of the data file filtering for the barcode
    res is the structure containing the filtering operation results
    barcode
    '''
    lt.check_and_create_directory(res.runid)
    with  open(res.runid+'/bc_'+barcode+'_summary.txt', 'w') as log_file:
        log_file.write('Barcode: ' + str(barcode) + '\n')
        log_file.write('Total number of reads (all barcodes): ' + str(res.count) + '\n')
        log_file.write('Number of pIDs: ' + str(len(res.good_reads[bc])) + '\n')
        log_file.write('Number of good reads with primer match: '+str(sum([len(res.good_reads[bc][i]) for i in res.good_reads[bc].keys()]))+'\n')
        #log_file.write('Number of reads with exact primer match: ' + str(int(sum([sum(res.good_reads[bc][pid].values()) for pid in res.good_reads[bc].keys()]))) + '\n')
        log_file.write('Number of good reads with invalid pID: ' + str(int(np.sum(res.good_read_bad_pID_length[barcode]))) + '\n')
        log_file.write('Number of un recognized reads (all barcodes): ' + str(int(sum([res.bad_read_length_454[length] for length in res.bad_read_length_454.keys()]))) + '\n')

#####
def export_good_reads_to_fasta_file_for_consensus_compress(res, barcode, mosp):
    '''
    function that exports the good reads occurring with at least res.min_occ_same_pid identical pIDs to make consensus
    res is the structure containing the filtering operation results
    barcode
    mosp min_occ_same_pid
    '''
    ####
    # make directory for this run
    ####
    lt.check_and_create_directory(res.runid)
    with open(res.runid+'/bc_'+barcode+'_min_'+str(mosp)+'.fasta', 'w') as out_file:
        for pID in sorted(res.good_reads[barcode].keys()):
            count_seq = res.good_reads[barcode][pID]
            if len(count_seq) >= mosp: # write the reads corresponding to pID if this one appears at least mosp times
                # count_seq = list of [id, 'fwd'/'rev', seq] corresponding to pID
                dict_seq = defaultdict(lambda:{'fwd':0, 'rev':0})
                for seq_read in count_seq: # count forward and reverse reads
                    dict_seq[seq_read[2]][seq_read[1]]+=1
                for read_to_write in dict_seq.keys():
                    nb_fwd = dict_seq[read_to_write]['fwd']
                    nb_rev = dict_seq[read_to_write]['rev']
                    out_file.write('>' + lt.read_label(pID, nb_fwd, nb_rev) + '\n')
                    out_file.write(read_to_write+'\n')
                

#####
# DEF PLOT FUNCTIONS
#####
def plot_read_length_distribution(res):
    plt.figure()
    ax = plt.subplot(111)
    cols = ['r', 'g', 'b','k']# a modifier pour prendre en compte le nombre de barcodes...
    plt.title('dotted: barcode matched, dash: 5prime end ok, solid: 5p and 3p ok')
    for bi,bc in enumerate(res.barcodes):
        plt.plot(res.good_read_length[bc], label = bc, c = cols[bi])
        #plt.plot(res.good_fwd_read_length[bc], ls = '--', c = cols[bi])
        #plt.plot(res.good_tag_only_read_length[bc], ls = ':', c = cols[bi])
        #plt.plot(res.bad_read_length - np.sum([res.good_tag_only_read_length[bc] for bc in res.barcodes], axis = 0), label = 'unmatched', c = 'k')
        ax.set_yscale('log')
        plt.xlabel('read length')
        plt.ylabel('number of reads')
        plt.legend(loc = 2)
        plt.xlim([0,500])
        plt.savefig(res.runid + '/read_length.pdf')

#####
def plot_number_of_reads_per_pID(res):
    plt.figure()
    ax = plt.subplot(111)
    cols = ['r', 'g', 'b','k']# a modifier pour prendre en compte le nombre de barcodes...
    for bi,bc in enumerate(res.barcodes):
        #plt.plot(np.arange(1, 1 + len(res.templates[bc])), sorted([len(res.templates[bc][t]) for t in res.templates[bc]], reverse = True), label = bc, c = cols[bi])
        plt.plot(np.arange(1, 1 + len(res.good_reads[bc])), sorted([len(res.good_reads[bc][pID]) for pID in res.good_reads[bc].keys()], reverse = True),label = bc, c = cols[bi])
        plt.ylabel('number of reads per PID')
        plt.xlabel('PID rank')
        ax.set_yscale('log')
        ax.set_xscale('log')
        plt.ylim([0.8,1.2e4])
        plt.legend(loc=1)
        plt.savefig(res.runid + '/copy_number.pdf')


        

######
# MAIN
######
if __name__=='__main__':
    try:
        execfile(CONFIG_FILE_NAME)
    except:
        print "messed up config file"         

    try:
        res = struct_var_set() # initialisation of the structure that will contain the filtering operation results
    except:
        print 'Problem for initialisation: there is something wrong with the configuration file ' + CONFIG_FILE_NAME + '!\n'

        
    ####
    # Filter the data file
    ####
    print '--analysis of the input data file: ' + res.input_data_file
    filter_reads(res)

    ####
    # Export the log and fasta files
    ####    
    for bc in res.barcodes:
        print '--generating the log file for the barcode ' + str(bc)
        logfile_output(res, bc)

    #for bc in res.barcodes:
        print '--exporting the good reads (barcode: ' + str(bc) +') occurring with at least ' + str(res.min_occ_same_pid) + ' identical pIDs (with length >= ' + str(res.min_length_pid)  + ') to make consensus'
        # generate a fasta file that contains only good reads with at least min_occ_same_pid identical pIDs to make consensus
        export_good_reads_to_fasta_file_for_consensus_compress(res, bc, res.min_occ_same_pid)
    ####
    # Plots
    ####

    ## read length distribution
    plot_read_length_distribution(res)

    ## number of reads per pID
    plot_number_of_reads_per_pID(res)
