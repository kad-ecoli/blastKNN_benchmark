#!/usr/bin/python2
docstring='''assess_result.py IA.txt input.tsv test.3.is_a output.txt
    calculate Fmax, coverage, weighted Fmax and Smin by full mode for prediction 

Input:
    IA.txt      - information content in the following format
                  [1] GO term
                  [2] information accretion, i.e., information content
                  [3] GO aspect
                  all subsequent columns are ignored
    input.txt   - GO prediction in the following format
                  [1] accession
                  [2] GO term
                  [3] cscore ranging from 0 (exclusive) to 1 (inclusive)
    test.3.is_a - ground truth GO annotation in the following format
                  [1] accession
                  [2] comma separated list of GO terms
Output:
    output.txt  - assessment result
                  [1] Aspect
                  [2] Fmax
                  [3] cutoff for Fmax
                  [4] SEM for Fmax
                  [5] Smin
                  [6] cutoff for Smin
                  [7] SEM for Smin
                  [8] wFmax
                  [9] cutoff for wFmax
                  [10] SEM for wFmax
                  [11] coverage
'''
import sys
from os.path import dirname, basename, abspath
from math import sqrt
root_terms={ "GO:0005575", "GO:0003674", "GO:0008150"}
Aspect_list="FPC"


def read_information_content(iafile):
    ic_dict=dict()
    Aspect_dict=dict()
    fp=open(iafile)
    for line in fp.read().splitlines():
        if line.startswith('#'):
            continue
        GOterm,ic,Aspect=line.split('\t')[:3]
        if GOterm in root_terms:
            continue
        #if float(ic)==0: # This makes incorrect Fmax but is consistent with
        #    continue     # https://github.com/BioComputingUP/CAFA-evaluator
        ic_dict[GOterm]=float(ic)
        Aspect_dict[GOterm]=Aspect
    fp.close()
    return ic_dict,Aspect_dict

def read_label(Aspect_dict,labelfile):
    label_dict=dict(F=dict(),P=dict(),C=dict())
    fp=open(labelfile)
    for line in fp.read().splitlines():
        target,GOterms=line.split('\t')
        for GOterm in GOterms.split(','):
            if not GOterm in Aspect_dict:
                continue
            Aspect=Aspect_dict[GOterm]
            if not target in label_dict[Aspect]:
                label_dict[Aspect][target]=[]
            label_dict[Aspect][target].append(GOterm)
    fp.close()
    for Aspect in Aspect_list:
        print("%d %s labeled targets"%(len(label_dict[Aspect]),Aspect))
        for target in label_dict[Aspect]:
            label_dict[Aspect][target]=set(label_dict[Aspect][target])
    return label_dict

def read_prediction(infile,label_dict,Aspect_dict):
    predict_dict=dict(F=dict(),P=dict(),C=dict())
    fp=open(infile)
    for line in fp.read().splitlines():
        if line.startswith('#'):
            continue
        items=line.split('\t')
        target=items[0]
        GOterm=items[1]
        cscore=items[-1]
        if not GOterm in Aspect_dict:
            continue
        Aspect=Aspect_dict[GOterm]
        if not target in label_dict[Aspect]:
            continue
        if not target in predict_dict[Aspect]:
            predict_dict[Aspect][target]=[]
        predict_dict[Aspect][target].append((float(cscore),GOterm))
    fp.close()
    for Aspect in predict_dict:
        print("%s %s predicted targets"%(len(predict_dict[Aspect]),Aspect))
        for target in predict_dict[Aspect]:
            predict_dict[Aspect][target]=[(GOterm,cscore) for cscore,GOterm in sorted(
            predict_dict[Aspect][target],reverse=True)]
    return predict_dict

def sum_ic(GOterm_list,ic_dict):
    return sum([ic_dict[GOterm] for GOterm in GOterm_list if GOterm in ic_dict])
    
def dot_product(list1,list2):
    if len(list1)!=len(list2):
        sys.stderr.write("ERROR! unequal list length\n")
        exit(1)
    return sum([list1[i]*list2[i] for i in range(len(list1))])

def assess_result(label_dict,ic_dict,predict_dict,outfile):
    txt="#Aspect\tFmax\tCutoff\tFsem\tSmin\tCutoff\tSsem\twFmax\tCutoff\twFsem\tCoverage\n"
    Fmax_list=[]
    Cutoff_F_list=[]
    Fsem_list=[]
    Smin_list=[]
    Cutoff_S_list=[]
    Ssem_list=[]
    wFmax_list=[]
    Cutoff_wF_list=[]
    wFsem_list=[]
    Coverage_list=[]
    for Aspect in Aspect_list:
        Fmax=0
        Cutoff_F=0
        F_list=[]
        Smin=0
        Cutoff_S=1
        S_list=[]
        wFmax=0
        wF_list=[]
        Cutoff_wF=0
        total_label=len(label_dict[Aspect])
        if total_label==0:
            continue
        print("evaluating %s"%Aspect)
        Coverage=1.*len(predict_dict[Aspect])/total_label
        cscore_list=[]
        label_ic_dict=dict()
        target_ic_list=[]
        for target in label_dict[Aspect]:
            label_ic_dict[target]=sum_ic(label_dict[Aspect][target],ic_dict)
            target_ic_list.append(label_ic_dict[target])
        sum_target_ic=sum(target_ic_list)
        Smin=sum_target_ic/total_label
        for target,GOterms in predict_dict[Aspect].items():
            cscore_list+=[cscore for GOterm,cscore in GOterms]
        cscore_list=sorted(set(cscore_list))
        for cutoff in cscore_list:
            total_precision_list=[]
            total_recall_list=[]
            total_predict=0
            total_wprecision_list=[]
            total_wrecall_list=[]
            total_ru_list=[]
            total_mi_list=[]
            for target in label_dict[Aspect]:
                if not target in predict_dict[Aspect]:
                    total_precision_list.append(0)
                    total_recall_list.append(0)
                    total_wprecision_list.append(0)
                    total_wrecall_list.append(0)
                    total_ru_list.append(label_ic_dict[target])
                    total_mi_list.append(0)
                    continue
                predict_list=[]
                for GOterm,cscore in predict_dict[Aspect][target]:
                    if cscore<cutoff:
                        break
                    predict_list.append(GOterm)
                if len(predict_list)==0:
                    total_precision_list.append(0)
                    total_recall_list.append(0)
                    total_wprecision_list.append(0)
                    total_wrecall_list.append(0)
                    total_ru_list.append(label_ic_dict[target])
                    total_mi_list.append(0)
                    continue
                label_set=label_dict[Aspect][target]
                total_predict+=1
                wpredict=sum_ic(predict_list,ic_dict)
                tp=label_set.intersection(predict_list)
                wtp=sum_ic(tp,ic_dict)
                precision=1.*len(tp)/len(predict_list)
                wprecision=0 if wpredict==0 else wtp/wpredict
                recall=1.*len(tp)/len(label_set)
                wrecall=0 if label_ic_dict[target]==0 else wtp/label_ic_dict[target]
                total_precision_list.append(precision)
                total_recall_list.append(recall)
                total_wprecision_list.append(wprecision)
                total_wrecall_list.append(wrecall)
                total_mi=0
                for GOterm in predict_list:
                    if not GOterm in label_set:
                        total_mi+=ic_dict[GOterm]
                total_mi_list.append(total_mi)
                total_ru=0
                for GOterm in label_set:
                    if not GOterm in predict_list:
                        total_ru+=ic_dict[GOterm]
                total_ru_list.append(total_ru)
            total_label=len(label_dict[Aspect])
            precision =0 if total_predict==0 else sum(total_precision_list)/total_predict
            wprecision=0 if total_predict==0 else sum(total_wprecision_list)/total_predict
            recall   =sum(total_recall_list)/total_label
            wrecall  =sum(total_wrecall_list)/total_label
            divide   =total_label
            #divide   =total_predict # This makses incorrect Smin but is consistent with https://github.com/BioComputingUP/CAFA-evaluator
            mi       =0 if divide==0 else sum(total_mi_list)/divide
            ru       =0 if divide==0 else sum(total_ru_list)/divide

            F  =0 if  precision * recall ==0 else 2/(1/ precision +1/ recall )
            wF =0 if wprecision *wrecall ==0 else 2/(1/wprecision +1/wrecall )
            S  =sqrt(ru *ru +mi *mi )
            if F>=Fmax:
                Fmax=F
                Cutoff_F=cutoff
            if S<=Smin:
                Smin=S
                Cutoff_S=cutoff
            if wF>=wFmax:
                wFmax=wF
                Cutoff_wF=cutoff
        target_list=[]
        for target in label_dict[Aspect]:
            target_list.append(target)
            if not target in predict_dict[Aspect]:
                F_list.append(0)
                wF_list.append(0)
                S_list.append(label_ic_dict[target])
                continue
            else:
                F_predict_list=[]
                wF_predict_list=[]
                S_predict_list=[]
                if target in predict_dict[Aspect]:
                    for GOterm,cscore in predict_dict[Aspect][target]:
                        if cscore>=Cutoff_F:
                            F_predict_list.append(GOterm)
                        if cscore>=Cutoff_wF:
                            wF_predict_list.append(GOterm)
                        if cscore>=Cutoff_S:
                            S_predict_list.append(GOterm)
                label_set=label_dict[Aspect][target]
                wpredict=sum_ic(wF_predict_list,ic_dict)
                tp=label_set.intersection(F_predict_list)
                wtp=sum_ic(label_set.intersection(wF_predict_list),ic_dict)
                precision=1.*len(tp)/len(F_predict_list) if len(F_predict_list) else 0
                wprecision=wtp/wpredict if wpredict>0 else 0
                recall=1.*len(tp)/len(label_set)
                wrecall=wtp/label_ic_dict[target] if label_ic_dict[target]>0 else 0
                F_list.append(2./(1./precision+1./recall) if precision*recall>0 else 0)
                wF_list.append(2./(1./wprecision+1./wrecall) if wprecision*wrecall>0 else 0)
                mi=sum([ic_dict[GOterm] for GOterm in S_predict_list if not GOterm in label_set])
                ru=sum([ic_dict[GOterm] for GOterm in label_set if not GOterm in S_predict_list])
                S_list.append(sqrt(mi*mi+ru*ru))
        targetNum=1.*len(F_list)
        txt_Aspect="#target\tF\tS\twF\n"
        for t in range(len(target_list)):
            txt_Aspect+="%s\t%.3f\t%.3f\t%.3f\n"%(target_list[t],F_list[t],S_list[t],wF_list[t])
        txt_Aspect+="#mean\t%.3f\t%.3f\t%.3f\n"%(sum(F_list)/targetNum,
            sum(S_list)/targetNum,sum(wF_list)/targetNum)
        fp=open(outfile+'.'+Aspect,'w')
        fp.write(txt_Aspect)
        fp.close()
        Fmean =sum(F_list)/targetNum
        wFmean=sum(wF_list)/targetNum
        Smean =sum(S_list)/targetNum
        Fsem  =sqrt(sum([( F- Fmean)*( F- Fmean) for  F in  F_list]))/targetNum
        wFsem =sqrt(sum([(wF-wFmean)*(wF-wFmean) for wF in wF_list]))/targetNum
        Ssem  =sqrt(sum([( S- Smean)*( S- Smean) for  S in  S_list]))/targetNum
        txt+="%s\t%.4f\t%.3f\t%.4f\t%.4f\t%.3f\t%.4f\t%.4f\t%.3f\t%.4f\t%.4f\n"%(
            Aspect,Fmax,Cutoff_F,Fsem,Smin,Cutoff_S,Ssem,wFmax,Cutoff_wF,wFsem,Coverage)
        Fmax_list.append(Fmax)
        Cutoff_F_list.append(Cutoff_F)
        Fsem_list.append(Fsem)
        Smin_list.append(Smin)
        Cutoff_S_list.append(Cutoff_S)
        Ssem_list.append(Ssem)
        wFmax_list.append(wFmax)
        Cutoff_wF_list.append(Cutoff_wF)
        wFsem_list.append(wFsem)
        Coverage_list.append(Coverage)

    divide=len(Fmax_list)
    if divide:
        txt+="#mean\t%.4f\t%.3f\t%.4f\t%.4f\t%.3f\t%.4f\t%.4f\t%.3f\t%.4f\t%.4f\n"%(
            sum(Fmax_list)/divide, sum(Cutoff_F_list)/divide, sum(Fsem_list)/divide,
            sum(Smin_list)/divide, sum(Cutoff_S_list)/divide, sum(Ssem_list)/divide,
            sum(wFmax_list)/divide, sum(Cutoff_wF_list)/divide, sum(wFsem_list)/divide,
            sum(Coverage_list)/divide)
    fp=open(outfile,'w')
    fp.write(txt)
    fp.close()
    return

if __name__=="__main__":
    if len(sys.argv)!=5:
        sys.stderr.write(docstring)
        exit()

    iafile   =sys.argv[1]
    infile   =sys.argv[2]
    labelfile=sys.argv[3]
    outfile  =sys.argv[4]

    ic_dict,Aspect_dict=read_information_content(iafile)
    label_dict=read_label(Aspect_dict,labelfile)
    predict_dict=read_prediction(infile,label_dict,Aspect_dict)
    assess_result(label_dict,ic_dict,predict_dict,outfile)
