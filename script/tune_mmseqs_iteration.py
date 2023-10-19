#!/usr/bin/env python3
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import os
import numpy as np
bindir=os.path.dirname(os.path.abspath(__file__))
rootdir=os.path.dirname(bindir)

method_list=["25","50","100"]

scoring_list=[
    ("1","pink"  ),
    ("2","green" ),
    ("3","blue"  ),
    ("4","orange"),
    ]


fontsize=9
width=1./len(scoring_list)
txt=''
for metric in ["Fmax","wFmax"]:
    txt+='####'+metric+'####\n'
    data=np.zeros((len(method_list),len(scoring_list)*3))
    txt+='#Aspect\t'
    for a,Aspect in enumerate(("MF","BP","CC")):
        for s,(label,color) in enumerate(scoring_list):
            txt+=Aspect+'\t'
    txt=txt[:-1]+'\n'
    txt+="#evalue\t"
    for a,Aspect in enumerate(("MF","BP","CC")):
        for s,(label,color) in enumerate(scoring_list):
            txt+=label+'\t'
    txt=txt[:-1]+'\n'
    plt.figure(figsize=(7.87,2*(len(method_list))))
    for m,method in enumerate(method_list):
        ax=plt.subplot(len(method_list),1,m+1)
        for s,(label,color) in enumerate(scoring_list):
            filename="%s/result/mmseqs_7.5_%s_%s_validate.txt"%(rootdir,method,label)
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
                    data[m][(len(value_list)-1)*len(scoring_list)+s]=value_list[-1]
            fp.close()
            plt.bar(np.arange(3)+width*(s+0.5), value_list,
                width=width,facecolor=color, label="--num-iterations "+label)
            for v,value in enumerate(value_list):
                plt.text(v+width*(s+0.5),
                    value+sem_list[v]+0.02,"%.3f"%value,va="bottom",ha="center",
                    fontsize=fontsize,rotation=0)
                plt.plot([v+width*(s+0.5),v+width*(s+0.5)],
                    [value,value+sem_list[v]],'k-')
        plt.yticks(fontsize=fontsize)
        plt.ylabel("%s for --max-seqs %s"%(metric,method),fontsize=fontsize)
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
        ax.tick_params(axis='both', which='both',
            pad=0.5, length=1, direction='in')
        if m==0:
            plt.legend(loc="upper center",fontsize=fontsize,ncol=1,
                borderpad=0.1,labelspacing=0.1,handlelength=1,
                handletextpad=0.1,borderaxespad=0.1,columnspacing=0.5)
    plt.tight_layout(pad=0.9,w_pad=0.1,h_pad=0.1)
    plt.savefig("%s/result/tune_mmseqs_iteration_%s.png"%(rootdir,metric),dpi=300)

    for m,method in enumerate(method_list):
        txt+=method+'\t'
        for s in range(len(data[m])):
            txt+="%.3f\t"%data[m][s]
        txt=txt[:-1]+'\n'
fp=open("%s/result/tune_mmseqs_iteration.tsv"%(rootdir),'w')
fp.write(txt)
fp.close()
