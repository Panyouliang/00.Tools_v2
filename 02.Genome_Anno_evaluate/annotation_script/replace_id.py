#! /usr/bin/env python3
# -*- coding:utf-8 -*-


import sys



def read_old(F):
    g_t_p = {}
    with open(F,'r') as f:
        for line in f:
            if line.strip():
                if line[0] != '#':
                    line = line.strip().split('\t')
                    if line[2] == 'mRNA':
                        tid = line[8][line[8].find('ID'):].split('=')[1].split(';')[0]
                        gid = line[8][line[8].find('Parent'):].split('=')[1].split(';')[0]
                        g_t_p[gid] = tid

    return g_t_p


def replace(F,g_t_p):
    mRNA_d = {}
    feature = {}
    with open(F,'r') as f:
        for line in f:
            if line.strip():
                if line[0] != '#':
                    line = line.strip().split('\t')
                    gid = line[8].split('=')[1].split(';')[0]
                    if line[2] == 'mRNA':

                        lis = line[0:8].copy()
                        lis[1] = 'ReName'
                        des = ['ID='+g_t_p[gid]+';Parent='+gid+';']
                        lis += des
                        mRNA_d[gid] = lis
                        feature[gid] = []

                    else:
                        flis = line[0:8].copy()
                        pdes = ['Parent='+g_t_p[gid]+';']
                        flis += pdes

                        feature[gid].append(flis)

    feature_sort = {}
    for k,v in feature.items():
        if mRNA_d[k][6] == '+':
            v.sort(key=lambda x:int(x[3]),reverse=False)
        if mRNA_d[k][6] == '-':
            v.sort(key=lambda x:int(x[3]),reverse=True)

    return mRNA_d,feature 


def main():
    gname = read_old(sys.argv[1])
    mRNA,feature = replace(sys.argv[2],gname)
    for k,v in feature.items():
        print('\t'.join(x for x in mRNA[k]))
        for line in v:
            print('\t'.join(x for x in line))


if __name__ == '__main__':
    main()

