#!/usr/bin/python2
docstring='''parse_blast.py blast.m6 db.3 submission.tsv score

input:
    blast.m6 - blast output from the following flag:
               -outfmt '6 qacc qlen sacc slen evalue bitscore length nident'
    db.3     - GO annotation
    score    - weight for template
                1: bitscore
                2: bitscore * nident / max(qlen, slen)
                3: nident / qlen
                4: nident / slen
                5: nident / length
                6: nident / max(qlen, slen)
                7: 1
               
               13: nident / qlen
               14: nident / slen
               15: nident / length
               16: nident / max(qlen, slen)

output:
    submission.tsv - GO term prediction
'''
import sys
from os.path import dirname, basename, abspath
from subprocess import Popen,PIPE

def read_annotation(dbfile):
    exp_dict=dict()
    if dbfile=='-':
        lines=sys.stdin.read().splitlines()
    else:
        fp=open(dbfile)
        lines=fp.read().splitlines()
        fp.close()
    for line in lines:
        sacc,GOterm_list=line.split('\t')
        exp_dict[sacc]=GOterm_list.split(',')
    return exp_dict

def parse_blast(infile):
    blast_dict=dict()
    fp=open(infile)
    stdout=fp.read()
    fp.close()
    for line in stdout.splitlines():
        items=line.split('\t')
        if len(items)!=8:
            continue
        qacc,qlen,sacc,slen,evalue,bitscore,length,nident=items
        qlen=float(qlen)
        slen=float(slen)
        bitscore=float(bitscore)
        nident=float(nident)
        length=float(length)
        evalue=float(evalue)
        if nident<=1:
            nident*=length
        if not qacc in blast_dict:
            blast_dict[qacc]=[]
            if len(blast_dict) % 1000 == 0:
                print("parse %d %s"%(len(blast_dict),qacc))
        ID=nident/max((qlen,slen))
        if length<1:
            length=1
        blast_dict[qacc].append([sacc,
            bitscore,
            bitscore * ID,
            nident / qlen,
            nident / slen,
            nident / length,
            ID,
            1.0,
        ])
    return blast_dict

def write_output(blast_dict,exp_dict,outfile,scoring):
    print("writing "+outfile)
    txt=''
    for t,target in enumerate(sorted(blast_dict.keys())):
        if t%1000==0:
            print("predict %d %s"%(t+1,target))
        predict_dict=dict()
        denominator=0
        for items in blast_dict[target]:
            sacc=items[0]
            score=items[scoring % 10]
            if not sacc in exp_dict:
                continue
            denominator+=score
            for GOterm in exp_dict[sacc]:
                if not GOterm in predict_dict:
                    predict_dict[GOterm]=0
                if scoring<10:
                    predict_dict[GOterm]+=score
                else:
                    predict_dict[GOterm]=max((score,predict_dict[GOterm]))
                    
        for cscore,GOterm in sorted([(predict_dict[GOterm],
            GOterm) for GOterm in predict_dict],reverse=True):
            if scoring<10 and cscore>0:
                cscore/=denominator
            cscore="%.3f"%cscore
            if cscore=="0.000":
                break
            txt+="%s\t%s\t%s\n"%(target,GOterm,cscore)
    fp=open(outfile,'w')
    fp.write(txt)
    fp.close()
    return

if __name__=="__main__":
    if len(sys.argv)!=5:
        sys.stderr.write(docstring)
        exit()

    infile =    sys.argv[1]
    dbfile =    sys.argv[2]
    outfile=    sys.argv[3]
    scoring=int(sys.argv[4])
    exp_dict=read_annotation(dbfile)
    blast_dict=parse_blast(infile)
    write_output(blast_dict,exp_dict,outfile,scoring)
