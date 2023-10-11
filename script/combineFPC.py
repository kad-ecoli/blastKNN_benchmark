#!/usr/bin/python2
docstring='''
combineFPC.py db.list db.F db.P db.C db.3
    combine 3 aspects of GO term annotations to one file

Input:
    db.list    - list of database entries
    db.{F,P,C} - GO annotations for each aspect

Output:
    db.3       - combined annotations
'''

import sys

if len(sys.argv)<4:
    sys.stderr.write(docstring)
    exit()

listFile     =sys.argv[1]
inputFileList=sys.argv[2:-1]
outputFile   =sys.argv[-1]

fp=open(listFile)
target_list=fp.read().splitlines()
fp.close()

goa_dict=dict()
for target in target_list:
    goa_dict[target]=[]

for filename in inputFileList:
    fp=open(filename)
    for line in fp.read().splitlines():
        target,GOterms=line.split('\t')
        goa_dict[target].append(GOterms)
    fp.close()

txt=''
for target in target_list:
    txt+=target+'\t'+','.join(goa_dict[target])+'\n'
fp=open(outputFile,'w')
fp.write(txt)
fp.close()
