#!/usr/bin/python


#
# script that gathers all the aligned reads for a given barcode in one file
# input: aligned files directory (in ../templates/)
# output: gathered aligned reads files in a specific 'all-align-in-one-file' directory (in ../templates/)

from Bio import AlignIO
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Align import MultipleSeqAlignment
import sys
<<<<<<< HEAD
import time
import datetime
import glob
=======
>>>>>>> fabae8322c3b401072c40dadc2418216e17800d5
import lib_tools as lt
import argparse

### parse command line arguments
parser = argparse.ArgumentParser(description='Gather individual alignment files into one big file')
parser.add_argument('--indir', required=True, type=str, help = 'directory with individual alignments')
parser.add_argument('--outfile', default=None, type=str, help = 'file to store the aligned reads')
args = parser.parse_args()
auto_file_name = parser.prog

<<<<<<< HEAD
if __name__=='__main__':
    if (len(sys.argv)==2):
        # parse the input aligned files directory name
        relative_path_to_align_dir=str(sys.argv[1])
        if(relative_path_to_align_dir[-1]!='/'):
            relative_path_to_align_dir+='/'
        path_to_templates = "../templates/"
        align_dir_basename = relative_path_to_align_dir.split('/')[2]
        [prefix_date_and_id,file_type,barcode] = [align_dir_basename.split('dir-')[1].split('_')[i] for i in [0,1,2]]

        # create (if necessary) the all-align-in-one-file directory
        lt.check_and_create_directory(str(path_to_templates+'dir-'+prefix_date_and_id+'_all-align-in-one-file'))
    
        # get the list of files to gather 
        list_files_to_gather = glob.glob(relative_path_to_align_dir+'*')
    
        global_list_alignments = []
    
        for i in range(len(list_files_to_gather)):
            list_seq_cur_file = []
            input_file_name = list_files_to_gather[i]
            for seq_record in SeqIO.parse(input_file_name, 'fasta'):
                list_seq_cur_file.append(seq_record)
            cur_alignment = MultipleSeqAlignment(list_seq_cur_file)
            global_list_alignments.append(cur_alignment)
        AlignIO.write(global_list_alignments, str(path_to_templates+'dir-'+prefix_date_and_id+'_all-align-in-one-file/'+prefix_date_and_id+'_all-align-in-one-file_'+barcode+'.fasta'), 'fasta')
    else:
        print auto_file_name+': usage: '+auto_file_name+' <aligned reads files directory (in ../templates/)>'
=======

# parse the input aligned files directory name
relative_path_to_align_dir = args.indir.rstrip('/')+'/'

if args.outfile is None:
    path_to_templates = "../templates/"
    align_dir_basename = relative_path_to_align_dir.split('/')[-1]
    [prefix_date_and_id,file_type,barcode] = [align_dir_basename.split('dir-')[1].split('_')[i] for i in [0,1,2]]
    outdir = path_to_templates+'dir-'+prefix_date_and_id+'_all-align-in-one-file'
    # create (if necessary) the all-align-in-one-file directory
    lt.check_and_create_directory(outdir)
    outfname = outdir+'/'+prefix_date_and_id+'_all-align-in-one-file_'+barcode+'.fasta'
else:
    outfname = args.outfile

# get the list of files to gather 
list_files_to_gather = glob.glob(relative_path_to_align_dir+'*.fa*')

global_list_alignments = []    
for fname in list_files_to_gather:
    with open(fname, 'r') as pid_alignment:
        pid = '.'.join(fname.split('.')[:-1]).split('_')[-1]
        global_list_alignments.append((pid, AlignIO.read(pid_alignment, 'fasta')))


global_list_alignments.sort(key = lambda x:x[0])
with open(outfname, 'w') as outfile: 
    for pid, aln in global_list_alignments:
        AlignIO.write(aln, outfname, 'fasta')



>>>>>>> fabae8322c3b401072c40dadc2418216e17800d5
