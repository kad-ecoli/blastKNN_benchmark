#!/usr/bin/env python3
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import os
import numpy as np
bindir=os.path.dirname(os.path.abspath(__file__))
rootdir=os.path.dirname(bindir)

method_list=[
    "blastp",
    "mmseqs",
    "diamond",
    "phmmer",
    "jackhmmer",
    "psiblast",
    "hhblits",
    ]

method_dict={
    "blastp"   :"BLASTp",
    "mmseqs"   :"MMseqs2",
    "diamond"  :"DIAMOND",
    "phmmer"   :"phmmer",
    "jackhmmer":"jackhmmer",
    "psiblast" :"PSI-BLAST",
    "hhblits"  :"HHblits",
    }

color_list=[
    "black",
    "white",
    "black",
    "black",
    "white",
    "black",
    "black",
    "black",
    "white",
    "black",
    "white",
    ]

scoring_list=[
    ( 1,"red",      "weight: bitscore"),
    ( 2,"black",    "weight: bitscore * nident / max(qlen, tlen)"),
    ( 3,"pink",     "weight: nident / qlen "),
    ( 4,"green",    "weight: nident / tlen"),
    ( 5,"blue",     "weight: nident / alnlen"),
    ( 6,"orange",   "weight: nident / max(qlen, tlen)"),
    ( 7,"magenta",  "weight: 1"),
    (13,"lightgrey","max: nident / qlen"),
    (14,"purple",   "max: nident / tlen"),
    (15,"yellow",   "max: nident / alnlen"),
    (16,"grey",     "max: nident / max(qlen, tlen)"),
    ]


fontsize=9
width=1./len(scoring_list)
for metric in ["Fmax","wFmax"]:
    if metric=="wFmax":
        plt.figure(figsize=(7.87,1.2*len(method_list)))
    else:
        plt.figure(figsize=(7.87,1.7*len(method_list)))
    for m,method in enumerate(method_list):
        ax=plt.subplot(len(method_list),1,m+1)
        for s,(suffix,color,label) in enumerate(scoring_list):
            filename="%s/result/%s_validate.%d.txt"%(rootdir,method,suffix)
            if not os.path.isfile(filename):
                continue
            value_list=[]
            sem_list=[]
            fp=open(filename,'r')
            for line in fp.read().splitlines():
                if not line.startswith('#'):
                    items=line.split('\t')
                    if metric=="Fmax":
                        value_list.append(float(items[1]))
                        sem_list.append(float(items[3]))
                    elif metric=="wFmax":
                        value_list.append(float(items[7]))
                        sem_list.append(float(items[9]))
            fp.close()
            plt.bar(np.arange(3)+width*(s+0.5), value_list,
                width=width,facecolor=color, label=label)
            for v,value in enumerate(value_list):
                plt.text(v+width*(s+0.5),
                    value+sem_list[v]+0.02,"%.3f"%value,va="bottom",ha="center",
                    fontsize=fontsize,rotation=70)
                plt.plot([v+width*(s+0.5),v+width*(s+0.5)],
                    [value,value+sem_list[v]],'k-')
                plt.text(v+width*(s+0.5),0.01,"S%d"%(s+1),fontsize=fontsize,
                    color=color_list[s],va="bottom",ha="center")
        plt.yticks([0,0.2,0.4,0.6,0.8],
            ["0","0.2","0.4","0.6","0.8"],fontsize=fontsize)
        plt.ylabel("%s for\n%s"%(metric,method_dict[method]),
            labelpad=0.2,fontsize=fontsize)
        if m+1== len(method_list):
            plt.xticks([0.5,1.5,2.5],[
                "Molecular Function (MF)",
                "Biological Process (BP)",
                "Cellular Component (CC)"],fontsize=fontsize)
        else:
            plt.xticks([])
        plt.plot([1,1],[0,1],'k--')
        plt.plot([2,2],[0,1],'k--')
        plt.axis([0,3,0,1])
        ax.tick_params(axis='y', which='both',
            pad=0.5, length=2, direction='in')
        ax.tick_params(axis='x', which='both',
            pad=0.5, width=0, length=1, direction='out')
    plt.tight_layout(pad=0.1,w_pad=-1.3,h_pad=-1.3)
    plt.savefig("%s/result/default_%s.png"%(rootdir,metric),dpi=300)
