#!/usr/bin/env python3
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import os
import numpy as np
bindir=os.path.dirname(os.path.abspath(__file__))
rootdir=os.path.dirname(bindir)

hit_num_list=["10","25","50","100",#"200",
    "300","500",#"1000"
]
evalue_list=[#"0.0001",
    "0.001","0.01","0.1","1","10",#"100"
    ]
s_list=["1.0","4.0","5.7","7.5","10"]
method_dict=dict(
    blastp="BLASTp",
    diamond="DIAMOND",
    mmseqs="MMseqs2",
)
aspect_dict=dict(
    MF="Molecular Function (MF)",
    BP="Biological Process (BP)",
    CC="Cellular Componenet (CC)"
    )

fontsize=10
plt.figure(figsize=(7.87,5))
for a,Aspect in enumerate(("MF","BP","CC")):
    for m,method in enumerate(("mmseqs","diamond","blastp")):
        ax=plt.subplot(3,3,3*m+a+1)
        if m==0:
            y_list=s_list
        else:
            y_list=evalue_list
        data=np.zeros((len(y_list),len(hit_num_list)))
        for e,evalue in enumerate(y_list):
            for h,hit_num in enumerate(hit_num_list):        
                filename="%s/result/%s_%s_%s_validate.txt"%(
                    rootdir,method,evalue,hit_num)
                if not os.path.isfile(filename):
                    continue
                fp=open(filename,'r')
                lines=fp.read().splitlines()
                fp.close()
                wF=float(lines[1+a].split()[7])
                data[e][h]=wF

        if a==0:
            vmin=0.45
            vmax=0.55
        elif a==1:
            vmin=0.3
            vmax=0.4
        elif a==2:
            vmin=0.6
            vmax=0.7
        plt.pcolormesh(data,vmin=vmin,vmax=vmax)
        if m==0:
            ticks=[vmin,vmin+0.05,vmax]
            label=[("%.2f"%c).rstrip('0') for c in ticks]
            cbar=plt.colorbar(pad=0.05,ticks=ticks,shrink=0.85,
                format="%.2f",location='top')
            cbar.ax.tick_params(labelsize=fontsize,pad=0.1)
        for e,evalue in enumerate(y_list):
            for h,hit_num in enumerate(hit_num_list):        
                plt.text(h+0.5,e+0.5,("%.3f"%data[e][h]).replace("0.","."),
                    fontsize=fontsize,ha="center",va="center")
        ax.tick_params(axis='both', width=0, pad=0.1)
        if m==2:
            plt.xticks(0.5+np.arange(len(hit_num_list)),
                hit_num_list,fontsize=fontsize)
        else:
            plt.xticks([])
        if m==0:
            plt.plot([4,5,5,4,4],[2,2,3,3,2],'k-')
            plt.plot([3,4,4,3,3],[3,3,4,4,3],'r--')
        elif m==1:
            plt.plot([1,2,2,1,1],[3,3,4,4,3],'r--')
            plt.plot([1,2,2,1,1],[0.02,0.02,1,1,0.02],'k-')
        elif m==2:
            plt.plot([1,2,2,1,1],[2,2,3,3,2],'r--')
            plt.plot([5,5.98,5.98,5,5],[4,4,4.98,4.98,4],'k-')
        if a==0:
            plt.yticks(0.5+np.arange(len(y_list)),y_list,fontsize=fontsize)
            if m==0:
                plt.ylabel(method_dict[method]+"\nsensitivity",
                    labelpad=0.1,fontsize=fontsize)
            else:
                plt.ylabel(method_dict[method]+" E-value",
                    labelpad=0.1,fontsize=fontsize)
        else:
            plt.yticks([])
        if m==2:
            #plt.xlabel("Hit number for\n"+aspect_dict[Aspect],fontsize=fontsize)
            plt.xlabel("Hit number for "+Aspect,fontsize=fontsize)

plt.tight_layout(pad=0.1,w_pad=0.1,h_pad=0.1)
plt.savefig("%s/result/tune.png"%(rootdir),dpi=300)
