###
# Configuration parameters for: p1_trim_and_filter
###

#######
# 454 #
#######
cfg={
    'p5_virus_match':'GTAGCATGACAAAAATCTTAGAGCC',
    'p3_virus_match':'CATTRCTTTGGATGGGTATGAA',
    'barcodes':['ACG','CGT','TAC'],
    'input_data_file':'../data/rawdata_reg2.fsa',
    'p5_cutoff': 21,#=len(p5_virus_match)-4
    'p3_cutoff': 18,#=len(p3_virus_match)-4
    'min_occ_same_pid':1,
    'min_length_pid':10,
    'reverse'=True#,
    #'barcode_length_range':range(3,5)
}

##############
# iontorrent #
##############
#cfg={
#    'runid':'2013-09-23-Run1',
#    'p5_virus_match':'TGGCAGTCTAGCAGAAGAAG',
#    'p3_virus_match':'CCTCAGGAGGGGACCCAG',
#    'barcodes':['TACG','ACGT','CGTA','GTAC'],
#    'input_data_file':'data/subsample_100000_iontorrent.fastq',
#    'reverse':False,
#    'min_occ_same_pid':1,
#    'min_length_pid':8,
#    'barcode_length_range':range(3,5)
#}

