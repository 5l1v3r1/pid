#!/usr/bin/python

# /!\: script that must be run ONLY after the end of the cluster jobs
# - script that cleans the cluster logs in src directory 
#   and moves the temp and aligned reads files from the cluster specific directory ( in ../templates) given in input
#   to their corresponding temp and align directory in ../templates (barcode specific)

import os
import sys
import time
import glob
from collections import Counter
from collections import defaultdict
import lib_tools as lt

auto_file_name = str(sys.argv[0])

        
######

if __name__=='__main__':

    if (len(sys.argv)==2):

        # 1. clean the cluster logs in src directory (normally the current directory)
        os.system('rm -f p3_cluster_align_aux.py.e* p3_cluster_align_aux.py.o*')
        print 'cluster logs p3_cluster_align_aux_py.e* and p3_cluster_align_aux.py.o* deleted'

        # 2. move the temp and aligned reads files from the cluster directory to their specific temp/align directory
        relative_path_to_cluster_dir = str(sys.argv[1])
        if relative_path_to_cluster_dir[-1]!='/':
            relative_path_to_cluster_dir+='/'

        path_to_templates = "../templates/"
        prefix_date_and_id=relative_path_to_cluster_dir.split('/')[2].split('cluster-')[1]
        print prefix_date_and_id
        
        #list_files = os.popen('ls '+relative_path_to_cluster_dir+'*').readlines()
        list_files = glob.glob(relative_path_to_cluster_dir+'*') 
        total_nb_files_in_cluster_dir = len(list_files)
        print total_nb_files_in_cluster_dir

        count_files = Counter()
        dict_files = defaultdict(list)
        
        for cur_file in list_files:
            cur_file = cur_file.strip()
            cur_file_base_name = cur_file.split('/')[-1]
            cur_file_type, cur_file_bc = [cur_file_base_name.split('_')[i] for i in [1,2]]
            count_files[cur_file_type]+=1
            dict_files[(cur_file_type,cur_file_bc)].append(cur_file_base_name)

        print '#files in cluster directory: '+str(count_files)
        # print dict_files
    
        # move the files to their directory
        for type_and_bc in dict_files.keys():
            # create the aligned files directory if necessary
            new_file_directory = str(path_to_templates+'dir-'+prefix_date_and_id+'_'+type_and_bc[0]+'_'+type_and_bc[1]+'/')
            lt.check_and_create_directory(new_file_directory)
            for cur_file in dict_files[type_and_bc]:            
                print 'move file '+relative_path_to_cluster_dir+cur_file+' to the directory: '+new_file_directory
                os.system('mv '+relative_path_to_cluster_dir+cur_file+' '+new_file_directory)
                count_files[type_and_bc[0]]-=1

        print '#remaining files in cluster directory: '+str(count_files)
    
    else:
        print auto_file_name+': usage: '+auto_file_name+' <cluster directory (in ../templates/)>'
