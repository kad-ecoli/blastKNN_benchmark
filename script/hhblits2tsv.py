#!/usr/bin/python2
docstring='''
hhblits2tsv.py query.fasta db.fasta input.hhr output.m6
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

def hhblits2tsv(query_filename, db_filename, hhr_filename, output_filename):
    query_len_dict=fasta2len(query_filename)
    db_len_dict   =fasta2len(db_filename)

    nident_dict=dict()
    fp=open(hhr_filename,'r')
    blocks=fp.read().split('\n>')
    fp.close()
    
    qacc=blocks[0].splitlines()[0].split()[1]
    txt=''
    for block in blocks[1:]:
        lines=block.splitlines()
        sacc=lines[0]
        items=lines[1].split()
        if len(items)<5:
            continue
        evalue  =items[1].split('=')[1]
        bitscore=items[2].split('=')[1]
        alnlen  =items[3].split('=')[1]
        pident  =items[4].split('=')[1].rstrip('%')
        nident  =float(pident)*0.01*int(alnlen)
        qlen    =query_len_dict[qacc]
        slen    =db_len_dict[sacc]
        
        txt+='\t'.join((qacc,qlen,sacc,slen,evalue,bitscore,
            alnlen,str(nident)))+'\n'

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
    hhr_filename   =sys.argv[3]
    output_filename=sys.argv[4]

    hhblits2tsv(query_filename, db_filename, hhr_filename, output_filename)
