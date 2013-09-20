import sys
import glob
import subprocess as sp

if len(sys.argv)==2:
    rundir = sys.argv[1].rstrip('/')+'/'
    
    reads_by_barcode = glob.glob(rundir+'bc*.fasta')
    for bc_file in reads_by_barcode:
        sp.call(['python', 'src/p2_sort.py',bc_file])

