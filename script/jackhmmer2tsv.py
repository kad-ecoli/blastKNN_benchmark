#!/usr/bin/python2
docstring='''
jackhmmer2tsv.py query.fasta db.fasta input.pfam input.tblout output.m6
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

def jackhmmer2tsv(query_filename, db_filename, pfam_filename,
    tblout_filename, output_filename):
    query_len_dict=fasta2len(query_filename)
    db_len_dict   =fasta2len(db_filename)

    nident_dict=dict()
    fp=open(pfam_filename,'r')
    for block in fp.read().split('//'):
        sacc_list=[]
        sacc_dict=dict()
        for line in block.splitlines():
            if line.startswith('#') or len(line)==0:
                continue
            sacc,sequence=line.split()
            sacc=sacc.split('/')[0]
            sacc_list.append(sacc)
            sacc_dict[sacc]=sequence

        for i in range(len(sacc_list)):
            qacc=sacc_list[i]
            if not qacc in query_len_dict:
                continue
            sequence1=sacc_dict[qacc]
            if not qacc in nident_dict:
                nident_dict[qacc]=dict()
            for j in range(len(sacc_list)):
                sacc=sacc_list[j]
                if not sacc in db_len_dict:
                    continue
                if not sacc in nident_dict:
                    nident_dict[sacc]=dict()
                sequence2=sacc_dict[sacc]
                length=0
                nident=0
                for r in range(len(sequence1)):
                    if sequence1[r]!='.' and sequence1[r]!='-':
                        nident+=(sequence1[r]==sequence2[r])
                        length+=1
                nident_dict[qacc][sacc]=(length,nident)
    fp.close()

    txt=''
    fp=open(tblout_filename,'r')
    for line in fp.read().splitlines():
        if line.startswith('#'):
            continue
        items   =line.split()
        sacc    =items[0]
        qacc    =items[2]
        evalue  =items[4]
        bitscore=items[5]
    
        qlen=query_len_dict[qacc]
        slen=db_len_dict[sacc]
        length=0
        nident=0
        #if qacc in nident_dict and sacc in nident_dict[qacc]:
        if sacc in nident_dict[qacc]:
            length,nident=nident_dict[qacc][sacc]
        txt+='\t'.join((qacc,qlen,sacc,slen,evalue,bitscore,
            str(length),str(nident)))+'\n'
    fp.close()

    fp=open(output_filename,'w')
    fp.write(txt)
    fp.close()
    return

if __name__=="__main__":
    if len(sys.argv)!=6:
        sys.stderr.write(docstring)
        exit()

    query_filename =sys.argv[1]
    db_filename    =sys.argv[2]
    pfam_filename  =sys.argv[3]
    tblout_filename=sys.argv[4]
    output_filename=sys.argv[5]

    jackhmmer2tsv(query_filename, db_filename, pfam_filename,
        tblout_filename, output_filename)
