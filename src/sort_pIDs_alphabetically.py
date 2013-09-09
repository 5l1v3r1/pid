#!/usr/bin/python
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import defaultdict
import sys


auto_file_name = str(sys.argv[0])

if(len(sys.argv)==2):
    dict_pIDs = defaultdict(list)
    input_file_name = str(sys.argv[1])
    check_extension = input_file_name.split('.')[-1]
    len_extension = len(check_extension)
    input_file_name_root = input_file_name[:-len_extension-1]

    with open(input_file_name, 'r') as infile:
        for record in SeqIO.parse(infile,'fasta'):
            pID = str(record.id)
            dict_pIDs[pID].append((pID,str(record.seq)))
    
    with open(str(input_file_name_root+'-sorted.'+check_extension), 'w') as outfile:
        for pID in sorted(dict_pIDs.keys()):
            outfile.write(str('>'+dict_pIDs[pID][0][0]+'\n'))
            outfile.write(str(dict_pIDs[pID][0][1]+'\n'))
else:
     print auto_file_name + ': usage: '+ auto_file_name + ' <file to sort>'
