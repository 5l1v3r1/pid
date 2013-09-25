import sys
import glob
import subprocess as sp

if len(sys.argv)==3:
    rundir = sys.argv[1].rstrip('/')+'/'
    readtype = sys.argv[2]

    reads_by_barcode = glob.glob(rundir+'bc*.fasta')
    for bc_dir in reads_by_barcode:
        sp.call(['python', 'src/p2_sort.py',bc_file, readtype])

