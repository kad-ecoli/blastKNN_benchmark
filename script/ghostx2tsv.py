#!/usr/bin/python2
docstring='''
ghostx2tsv.py query.fasta db.fasta input.ghostx output.m6
'''
import os
import sys
import subprocess

bindir=os.path.dirname(os.path.abspath(__file__))

def fasta2len(filename):
    cmd="%s/fasta2len %s"%(bindir,filename)
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    len_dict=dict()
    for line in stdout.decode().splitlines():
        items=line.split('\t')
        len_dict[items[0]]=items[1]
    return len_dict

def ghostx2tsv(query_filename, db_filename, tblout_filename, output_filename):
    query_len_dict=fasta2len(query_filename)
    db_len_dict   =fasta2len(db_filename)

    txt=''
    fp=open(tblout_filename,'r')
    for line in fp.read().splitlines():
        if line.startswith('#'):
            continue
        items   =line.split()
        qacc    =items[0]
        sacc    =items[1]
        seqid   =items[2]
        alnlen  =items[3]
        evalue  =items[10]
        bitscore=items[11]
    
        qlen=query_len_dict[qacc]
        slen=db_len_dict[sacc]
        nident=float(seqid)*float(alnlen)/100.
        txt+='\t'.join((qacc,qlen,sacc,slen,evalue,bitscore,
            alnlen,"%.0f"%nident))+'\n'
    fp.close()

    fp=open(output_filename,'w')
    fp.write(txt)
    fp.close()
    return

if __name__=="__main__":
    if len(sys.argv)!=5:
        sys.stderr.write(docstring)
        exit()

    query_filename =sys.argv[1]
    db_filename    =sys.argv[2]
    tblout_filename=sys.argv[3]
    output_filename=sys.argv[4]

    ghostx2tsv(query_filename, db_filename, 
        tblout_filename, output_filename)
