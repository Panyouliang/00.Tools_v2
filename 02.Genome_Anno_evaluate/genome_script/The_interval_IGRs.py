#1 /usr/bin/env python3
# -*- coding:utf-8 -*-
#
import sys,gzip,re,math
import numpy as np

F = sys.argv[1]

if F[-3:] == '.gz':
    f = gzip.open(F,'rt')
else:
    f = open(F,'r')


list_s = []
for line in f:
    if line[0] != '#':
        line = line.strip().split('\t')
        if line[2] == 'mRNA':
            list_s.append([line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8]])

list_sort = sorted(list_s,key = lambda x:[x[0],int(x[3])],reverse=False)


dict_s = {}
for i in list_sort:
    key = i[0]+i[3]+i[6]+i[4]
    dict_s[key] = i


dict_g = {}
seq = ''
for k,line in dict_s.items():
    name = line[8].split('=')[1].split(';')[0]
    if line[0] != seq:
        seq = line[0]
        d = [[int(line[3]),int(line[4]),str(name)]]
        dict_g[seq] = d
    elif line[0] == seq:
        seq = line[0]
        d.append([int(line[3]),int(line[4]),str(name)])
        dict_g[seq] = d


#print('Chromosome','geneID','Species','length',sep='\t')

IGRlis = []

for k,v in dict_g.items():
    if len(v) > 1:
        for i in range(len(v)-1):
            length = abs(int(v[i][1])-int(v[i+1][0]))
            IGRlis.append(length)
            #print(k,v[i][2],sys.argv[2],length,sep='\t')


print('mean: '+str(np.mean(IGRlis)),'median: '+str(np.median(IGRlis)),sep='\t')
