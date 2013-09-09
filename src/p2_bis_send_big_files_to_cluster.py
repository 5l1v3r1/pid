#!/usr/bin/python
#/!\#
#### /ebio/ag-neher/share/pograms/EPD/bins/python
#/!\#

# script that moves temp reads files with a size greater than the given limit from the given temp directory to a specific cluster directory to prepare alignments
# /!\ : there is only ONE cluster directory per prefix_date_and_id, it is not barcode specific...
# input: temp directory, size min
# output: the temp files are moved into a specific cluster directory (with same prefix_date_and_id as the temp directory)


import os
import sys
import datetime
import time
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import Counter
from collections import defaultdict
import lib_tools as lt

auto_file_name = str(sys.argv[0])


if (len(sys.argv)==3):
    dir_to_check = str(sys.argv[1])
    if (dir_to_check[-1]!='/'):
        dir_to_check+='/'
    file_size_min = int(sys.argv[2])
    ####
    #parse the input directory name
    ####
    path_to_templates = '../templates/'
    #temp_dir_basename = dir_to_check[len(path_to_templates):]
    temp_dir_basename = dir_to_check.split('/')[2]
    [prefix_date_and_id, bc]= [temp_dir_basename.split('_')[i] for i in [0,2]]
    prefix_date_and_id = prefix_date_and_id[len('dir-'):] # "dir-" to remove
    #if bc[-1]=='/':
    #    bc= bc[:-1] # "/" to remove
    print prefix_date_and_id, bc
    ####
    # create the specific cluster directory
    ####
    cluster_dir_name = str(path_to_templates+'cluster-'+prefix_date_and_id)
    lt.check_and_create_directory(cluster_dir_name)

    ####
    # list the temp files and check their size
    ####
    list_temp_files = os.popen('ls '+dir_to_check+'*').readlines()
    list_temp_files = [list_temp_files[i].strip() for i in range(len(list_temp_files))]
    #print list_temp_files
    count_nb_files_moved_to_cluster=0
    # check the size of the temp reads file for each pID
    for cur_file in list_temp_files:
        cur_file_size = os.path.getsize(cur_file)
        #print 'size of file '+cur_file+': '+str(cur_file_size)
        if(cur_file_size >= file_size_min):
            count_nb_files_moved_to_cluster +=1
            print 'file '+cur_file+': '
            print ' --> size = '+str(cur_file_size)+' >= '+str(file_size_min)+' ==>  moved to '+ cluster_dir_name
            os.system(str('mv '+cur_file+' '+cluster_dir_name))
    print str(count_nb_files_moved_to_cluster)+' files with size >= '+str(file_size_min)+' were moved from '+dir_to_check+' to '+cluster_dir_name
else:
    print auto_file_name+': usage: '+auto_file_name+' <temp reads file directory (in ../templates)> <file size min.>'


