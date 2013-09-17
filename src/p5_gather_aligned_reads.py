#!/usr/bin/python


#
# script that gathers all the aligned reads for a given barcode in one file
#

from Bio import AlignIO
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Align import MultipleSeqAlignment
import sys
import lib_tools as lt
import argparse

### parse command line arguments
parser = argparse.ArgumentParser(description='Gather individual alignment files into one big file')
parser.add_argument('--indir', required=True, type=str, help = 'directory with individual alignments')
parser.add_argument('--outfile', default=None, type=str, help = 'file to store the aligned reads')
args = parser.parse_args()
auto_file_name = parser.prog


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



