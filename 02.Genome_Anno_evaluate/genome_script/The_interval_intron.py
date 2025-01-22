#! /usr/bin/env python3 
# -*- coding:utf-8 -*-
#

import sys,gzip,re

F = sys.argv[1]

if F[-3:] == '.gz':
    f = gzip.open(F,'rt')
    t = 'gene'
else:
    f = open(F,'r')
    t = 'mRNA'

dict_g = {}
for line in f:
    if line[0] != '#':
        line = line.strip().split('\t')
        if line[2] == t:
            lstr = re.split(r';',line[8])
            name = re.split(r'=',lstr[0])[1]
            exon = []
        elif line[2] == 'CDS':
            cds = [int(line[3]),int(line[4]),line[6]]
            exon.append(cds)
            dict_g[name] = exon

print('gene_id','Species','intron',sep='\t')
for k,v in dict_g.items():
    if len(v) > 1:
        for i in range(len(v)-1):
            if v[i][2] == '+':
                intron_l = abs(v[i][1]-v[i+1][0])
                print(k,sys.argv[2],intron_l,sep='\t')
            elif v[i][2] == '-':
                intron_2 = abs(v[i][0]-v[i+1][1])
                print(k,sys.argv[2],intron_2,sep='\t')



