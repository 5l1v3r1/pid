#!/usr/bin/python
import os
import sys

cmd = '/ebio/ag-neher/share/programs/EPD/bin/python2.7 ' + '/ebio/ag-neher/share/users/ebenard/PID_Karolinska_2/src/'

for arg in sys.argv[1:]:
    cmd+=arg+' '

print cmd
os.system(cmd)
