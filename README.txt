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

## step 4:

COMMAND: python src/p4_consensus.py run_directory read_type

This script goes over all barcodes in the run directory, gathers the aligned read files in the temporary directory of the desired read type, and builds consensus sequences. it also writes all aligned reads into a single file. 

## step 5:

COMMAND: python src/p5_decontamination.py bar_code_directory ref_seqs read_type true_seq_id

This script takes the aligned reads from one barcode and checks whether the individuals reads or the consensus sequence alignes reasonably well to the reference sequence with the true_seq_id. If a read does not, it is checked against all other reference sequences. All reads that don't align well to their own reference sequence are written into an extra file. 

alternatively to submit batch-jobs to the cluster:

COMMAND: python src/p5_decontamination.py run_directory ref_seqs read_type true_seq_id

## step 6:

COMMAND: python src/p6_detect_mutants_indels.py barcode_dir read_type

Check whether PIDs of low abundance reads are less than a certain edit distance from a high abundance one. Designate a neighbor if reads in addition align well. Produce read files with likely_pIDs and original PIDs. 

After this step, the sorting, alignment and consensus steps (2-4) need to be redone with readtype corrected instead of filtered. 
