
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
