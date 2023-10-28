#!/usr/bin/env python3
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import os
import numpy as np
bindir=os.path.dirname(os.path.abspath(__file__))
rootdir=os.path.dirname(bindir)

method_list=[
    ("BLASTp (-evalue 0.1 -max_target_seqs 25)","red","blastp_0.1_25_validate.txt"),
    ("DIAMOND (--evalue 1 --ultra-sensitive)","green","diamond_1_25_validate.txt"),
    ("MMseqs2 (-s 7.5 --max-seqs 50)","black","mmseqs_7.5_100_validate.txt"),
    ("BLASTp (default)","pink","blastp_validate.2.txt"),
    ("DIAMOND (default)","lightgreen","diamond_validate.2.txt"),
    ("MMseqs2 (default)","lightgrey","mmseqs_validate.2.txt"),
    ]

fontsize=10
width=1./len(method_list)
for c,metric in enumerate(["Fmax","wFmax"]):
    plt.figure(figsize=(7.87,1.9))
    ax=plt.subplot(1,1,1)
    for m,items in enumerate(method_list):
        method=items[0]
        color=items[1]
        filename=rootdir+"/result/"+items[2]
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
        label=method.replace('\n',' ')
        plt.bar(np.arange(3)+width*(m+0.5), value_list,
            width=width,facecolor=color,#edgecolor='k',
            label=label)
        for v,value in enumerate(value_list):
            plt.text(v+width*(m+0.5),
                value+sem_list[v]+0.02 #+0.18*(v==1)*(m==3)
                ,"%.3f"%value,va="bottom",ha="center",
                fontsize=fontsize,rotation=10)
            plt.plot([v+width*(m+0.5),v+width*(m+0.5)],
                [value,value+sem_list[v]],'k-')
            #plt.text(v+width*(m+0.5),0.01,method,
                #va="bottom",ha="center",
                #fontsize=fontsize,rotation=90)
    plt.yticks([0,0.2,0.4,0.6,0.8,1],[0,0.2,0.4,0.6,0.8,1],fontsize=fontsize)
    plt.ylabel("%s"%(metric),fontsize=fontsize)
    if c==1:
        plt.xticks([0.5,1.5,2.5],[
                "Molecular Function (MF)",
                "Biological Process (BP)",
                "Cellular Component (CC)"],fontsize=fontsize)
    else:
        plt.xticks([])
    plt.legend(loc="upper left",fontsize=fontsize,ncol=2,
        borderpad=0.1,labelspacing=0.1,handlelength=1,
        handletextpad=0.1,borderaxespad=0.2,columnspacing=0.1)
    plt.plot([1,1],[0,1],'k--')
    plt.plot([2,2],[0,1],'k--')
    plt.axis([0,3,0,1])
    ax.tick_params(axis='y', which='both',
        pad=0.5, length=2, direction='in')
    ax.tick_params(axis='x', which='both',
        pad=0.5, width=0, length=1, direction='out')
    plt.tight_layout(pad=0.1,w_pad=0.1,h_pad=0.1)
    plt.savefig("%s/result/best_%s.png"%(rootdir,metric),dpi=300)
