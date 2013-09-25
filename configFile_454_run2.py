###
# Configuration parameters for: p1_trim_and_filter
###

#######
# 454 #
#######
cfg={
    'runid':'PID_454_run2',
    'p5_virus_match':'GTAGCATGACAAAAATCTTAGAGCC',
    'p3_virus_match':'CATTRCTTTGGATGGGTATGAA',
    'barcodes':['ACG','CGT','TAC'],
#    'input_data_file':'../data/rawdata_reg2.fsa',
    'input_data_file':'/ebio/ag-neher/share/data/PID_454/rawdata_reg2.fsa',
    'p5_cutoff': 21,#=len(p5_virus_match)-4
    'p3_cutoff': 18,#=len(p3_virus_match)-4
    'min_occ_same_pid':1,
    'min_length_pid':10,
    'reverse':True,
    'barcode_length_range':range(3,5)
}

