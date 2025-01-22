#! /usr/bin/env python3
# -*- coding:utf-8 -*-

import sys
import numpy as np


def get_exon(F):
    gene_d = {}
    with open(F,'r') as f:
        for line in f:
            if line[0] != '#':
                line = line.strip().split('\t')
                if line[2] == 'mRNA':
                    tid = line[8][line[8].find('ID'):].split('=')[1].split(';')[0]
                    gene_d[tid] = []
                else:
                    if line[2] == 'exon':
                        pid = line[8][line[8].find('Parent'):].split('=')[1].split(';')[0]
                        gene_d[pid].append(int(line[4])-int(line[3])+1)


    return gene_d

def N(lis,sum,num):
    total = 0
    ID = 'N'+str(num)
    for i in lis:
        total += i
        if total >= float('0.'+str(num))*sum:
            ID = i
            break
    return ID



def main():
    g_lis = get_exon(sys.argv[1])

    single_exon = []
    for k,v in g_lis.items():
        if len(v) == 1:
            single_exon.append(v[0])

    single_exon = sorted(single_exon,reverse=False)
    exon_lis = np.array(single_exon) 

    sum_l = np.sum(exon_lis)

    N50=N(single_exon,sum_l,5)
    N10=N(single_exon,sum_l,1)

    #print(np.max(exon_lis),np.min(exon_lis),np.percentile(exon_lis,1),N50,sep='\t')

    print('max:',str(np.max(exon_lis)),sep='\t')
    print('min:',str(np.min(exon_lis)),sep='\t')
    print('1%:',str(np.percentile(exon_lis,1)),sep='\t')
    print('5%:',str(np.percentile(exon_lis,5)),sep='\t')
    print('N50:',str(N50),sep='\t')
    print('N10:',str(N10),sep='\t')



if __name__ == '__main__':
    main()



