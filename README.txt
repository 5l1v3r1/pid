#### ANALYSIS OF PRIMER ID SEQUENCING DATA ####

## step 1:

COMMAND: python src/p1_trim_and_filter.py configfile

this will check each read for 5' and 3' primer matches (soft matching via Smith Waterman alignment) and split the reads according to their bar codes. filtered_read.fasta files will be deposited in labeled directories. 

## step 2:

COMMAND: python src/p2_sort.py run_directory read_type

this will split the reads according to their pIDs into a largish number of temporary directories. All pIDs within each directory will be aligned by a cluster job. The parameter read_type specifies whether this is to be done one the filtered or corrected reads. In any case, the script looks for a file named  read_type+"_reads.fasta"

## step 3:

COMMAND: python src/p3_cluster_align.py run_directory

Starts a cluster job for each of the temp directory in each of the barcode directories inside the run_directory. 



