#! /usr/bin/env python3
# -*- coding:utf-8 -*-


import sys,gzip


def stat_gL(F):
    if F[-3:] == '.gz':
        f = gzip.open(F,'rt')
    else:
        f = open(F,'r')

    d = {}
    T = 0
    for line in f:
        if line.strip():
            if line[0] == '>':
                ids = line[1:].strip().split()[0]
                l = 0
            else:
                l += len(line.strip())
                T += len(line.strip())
                d[ids] = l

    d['genome'] = T
    return d

def stat_dp(F):
    if F[-3:] == '.gz':
        f = gzip.open(F,'rt')
    else:
        f = open(F,'r')

    d = {}
    for line in f:
        if line.strip():
            line = line.strip().split()
            if int(line[2]) >= 5:
                if line[0] in d.keys():
                    d[line[0]].append(int(line[2]))
                else:
                    d[line[0]] = [int(line[2])]
    return d


genomel = stat_gL(sys.argv[1])
depth = stat_dp(sys.argv[2])


print('chr_num','length','coveRate','depth avg(×)',sep='\t')

cover_bp,depth_read = 0,0
for k,v in depth.items():
    cover_bp += len(v)
    depth_read += sum(v)
    print(k,genomel[k],len(v)/genomel[k],sum(v)/genomel[k],sep='\t')

print('genome',genomel['genome'],cover_bp/genomel['genome'],depth_read/genomel['genome'],sep='\t')


