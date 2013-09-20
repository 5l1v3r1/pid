###
# Configuration parameters for: p1_trim_and_filter
###
####
# 454
###
#cfg={
#    'prefix_date_and_id':'2013-09-04-Run1',
#    'p5_virus_match':'TGGCAGTCTAGCAGAAGAAG',
#    'p3_virus_match':'CCTCAGGAGGGGACCCAG',
#    'barcodes':['ACGT', 'TACG', 'CGTA'],
#    'input_data_file':'subsample.fastq',
#    'reverse':True,
#    'min_occ_same_pid':1,
#    'min_length_pid':10,
#    'barcode_length_range':range(3,5)
#}
#
###
# iontorrent
###
cfg={
    'runid':'2013-09-04-Run1',
    'p5_virus_match':'TGGCAGTCTAGCAGAAGAAG',
    'p3_virus_match':'CCTCAGGAGGGGACCCAG',
    'barcodes':['TACG','ACGT','CGTA','GTAC'],
    'input_data_file':'data/subsample_iontorrent.fastq',
    'reverse':False,
    'min_occ_same_pid':1,
    'min_length_pid':8,
    'barcode_length_range':range(3,5)
}

