#!/usr/bin/python2
docstring='''
make_target.py goa_old.is_a goa_new.is_a out.is_a
    compare annotations in goa_old.is_a and goa_new.is_a
    output the set of no knowledge and limited knowledge targets to test.is_a
'''
import sys

def make_target(oldfile,newfile,outfile):
    goa_old_dict=dict()
    fp=open(oldfile,'r')
    for line in fp.read().splitlines():
        target,GOterms=line.split()
        goa_old_dict[target]=GOterms
    fp.close()

    txt=''
    fp=open(newfile,'r')
    for line in fp.read().splitlines():
        target,GOterms=line.split()
        if target in goa_old_dict:
            continue
        txt+=line+'\n'
    fp.close()

    fp=open(outfile,'w')
    fp.write(txt)
    fp.close()
    return

if __name__=="__main__":
    if len(sys.argv)!=4:
        sys.stderr.write(docstring)
        exit()

    oldfile=sys.argv[1]
    newfile=sys.argv[2]
    outfile=sys.argv[3]
    make_target(oldfile,newfile,outfile)
