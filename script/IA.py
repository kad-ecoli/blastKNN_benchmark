#!/usr/bin/python2
docstring='''
IA.py Aspect is_a.csv train_terms.3 IA.txt
    calculate information content (a.k.a. information accretion)
for GO terms

Input:
    Aspect        - GO aspect to consider: F, P, C, or FPC
    is_a.csv      - is_a relation
                    [1] GO term
                    [2] Aspect
                    [3] Direct parent
                    [4] indirect parent (ignored)
    train_terms.3 - GO annotations
                    [1] accession
                    [2] comma separate list of GO term

Output:
    IA.txt        - information content
                    [1] GO term, sorted in numerical order
                    [2] information content, calculated by
                        IA(q)=-log2((1+N(q))/(1+N_p(q))
                        Here, N(q) is the number of accession
                        with GO term q; N_p(q) is the number of
                        accession with all parents of GO term q

'''
import sys
from math import log
ln2log2=1./log(2)

def read_isafile(isafile,Aspect):
    isa_dict=dict()
    fp=open(isafile)
    for line in fp.read().splitlines():
        items=line.split('\t')
        if not items[1] in Aspect:
            continue
        GOterm=items[0]
        parents=items[2]
        if not parents:
            isa_dict[GOterm]=[]
        else:
            isa_dict[GOterm]=parents.split(',')
    fp.close()
    return isa_dict

def read_annotation(dbfile,isa_dict):
    GOdict=dict()
    fp=open(dbfile)
    for line in fp.read().splitlines():
        accession,GOterms=line.split('\t')
        for GOterm in GOterms.split(','):
            if not GOterm in isa_dict:
                continue
            if not GOterm in GOdict:
                GOdict[GOterm]=[accession]
            else:
                GOdict[GOterm].append(accession)
    fp.close()
    return GOdict

def calculate_ia(GOdict,isa_dict):
    ia_dict=dict()
    for GOterm in isa_dict:
        if len(isa_dict[GOterm])==0:
            ia_dict[GOterm]=0
            continue
        NumGOterm=0
        if GOterm in GOdict:
            NumGOterm=len(GOdict[GOterm])
        accession_set=set()
        for p,parent in enumerate(isa_dict[GOterm]):
            parent_set=set()
            if parent in GOdict:
                parent_set=set(GOdict[parent])
            if p==0:
                accession_set=parent_set
            else:
                accession_set=accession_set.intersection(parent_set)
        NumParent=len(accession_set)
        ia_dict[GOterm]=-log((1.+NumGOterm)/(1.+NumParent))*ln2log2
    return ia_dict

if __name__=="__main__":
    if len(sys.argv)!=5:
        sys.stderr.write(docstring)
        exit()

    Aspect =sys.argv[1]
    isafile=sys.argv[2]
    dbfile =sys.argv[3]
    outfile=sys.argv[4]

    isa_dict=read_isafile(isafile,Aspect)
    GOdict=read_annotation(dbfile,isa_dict)
    ia_dict=calculate_ia(GOdict,isa_dict)

    txt=''
    for GOterm in sorted(ia_dict.keys()):
        line="%s\t%.15f\t%s\n"%(GOterm,ia_dict[GOterm],Aspect)
        if line.endswith("-0.000000000000000\t%s\n"%Aspect):
            line="%s\t0.000000000000000\t%s\n"%(GOterm,Aspect)
        txt+=line
    fp=open(outfile,'w')
    fp.write(txt)
    fp.close()
