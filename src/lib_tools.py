from Bio.pairwise2 import align
from Bio import SeqIO
from collections import defaultdict
import os

###
def check_match_SW(seq, match_seqs):
    '''
    function that checks a read for correct sequences at the left and right end.
    seq is the read as np.array
    match_seqs is a list of sequences to be matched, again as np.array with the expected start pos
    '''
    scores = []
    for qseq in match_seqs:
        scores.append(align.localms(seq, qseq, 1, 0, -0.6, -0.3)[0])
    return scores

###
def align_SW_dist(seq1,seq2):
    score = align.globalms(seq1, seq2, 1, 0, -0.6, -0.3)[0]
    return score

###
def hamming_dist(s1, s2):
    assert len(s1) == len(s2)
    return sum(imap(operator.ne, s1, s2))

###
def check_and_create_directory(dir_name):
    if  os.path.exists(dir_name):
        print 'Achtung: this directory already exists: '+dir_name
    else:
        os.makedirs(dir_name)

def get_last_part_of_path(fname):
    return fname.split('/')[-1]
def get_base_of_path(fname):
    return '/'.join(fname.split('/')[:-1])

def trim_extension(fname):
    return '.'.join(fname.split('.')[:-1])

def get_extention(fname):
    return fname.split('.')[-1]

def read_label(pID, nfwd, nrev, orig_pID=None):
    if orig_pID:
        return '_'.join(map(str, [pID, nfwd, nrev]))
    else:
        return '_'.join(map(str, [pID, nfwd, nrev, orig_pID]))

def parse_read_label(read_label):
    entries = read_label.split('_')
    if len(entries)==3:
        return entries[0], int(entries[1]), int(entries[2]), entries[0]
    elif len(entries)==4:
        return entries[0], int(entries[1]), int(entries[2]), entries[3]
    else:
        print "parse_read_label: invalid read label:", read_label

def parse_readfile(fname):
    count_reads_added_to_dict_all_reads=0
    count_reads_with_invalid_pID=0
    dict_all_reads = defaultdict(list)
    with open(fname, 'r') as reads_file:
        # for all the valid reads (without Ns), add the read to dict_all_reads
        for record in SeqIO.parse(reads_file, 'fasta'):
            pID,nb_reads_fwd,nb_reads_rev = record.id.split('_')[0:3]
            nb_reads_fwd = int(nb_reads_fwd)
            nb_reads_rev = int(nb_reads_rev)
            seq = str(record.seq)
            if (pID.count('N') == 0): # if the pID is not ambiguous
                dict_all_reads[pID].append([nb_reads_fwd,nb_reads_rev,seq])
                count_reads_added_to_dict_all_reads += (nb_reads_fwd+nb_reads_rev)
            else:
                #print 'Invalid pID: ' + pID
                count_reads_with_invalid_pID += (nb_reads_fwd+nb_reads_rev)
    return dict_all_reads, count_reads_added_to_dict_all_reads, count_reads_with_invalid_pID

def make_file_name(name_parts):
    fields = []
    for f in ['bc', 'pid', 'min']:        
        if f in name_parts:
            fields.extend([f, str(name_parts[f])])

    fields.extend(map(str, name_parts['other']))
    return name_parts['runid']+'/'+'_'.join(fields)+'.'+name_parts['ext']

def parse_file_name(fname):
    entries = fname.split('/')[-1].split('_')
    fields=[]
    remainder=0
    for f in ['bc', 'pid', 'min']:
        if entries[remainder]==f: 
            fields[f]=entries[remainder+1]
            try:
                fields[f]=int(fields[f])
            except:
                pass
            remainder+=2

    fields['other']=entries[remainder:]
    return fields

def make_dir_name(runid, barcode, dir_label):
    return 'dir-'+'_'.join([runid, barcode, dir_label])

def parse_file_name(dir_name):
    return dir_name[4:].split('_')[:2]

def check_neighbor_plausibility(seq1, seq2, distance_cutoff, verbose = False):
    score = align.localms(seq1, seq2, 1, 0, -1, -1)[0]
    if verbose:
        print 'Alignment score: '+ str(score[2])
        print score[0]
        print score[1]
    return (score[2] >= len(seq1) - DIST_MAX)

def remove_gaps(seq):
    return seq.replace('-','')
