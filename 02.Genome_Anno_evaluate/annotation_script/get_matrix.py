import sys
import numpy as np


# feature_count matrix
gcount = {}
for line in open(sys.argv[1],'r'):
    if line.strip():
        line = line.strip().split('\t')
        if line[2] == 'mRNA':
            ids = line[8].split('=')[1].split(';')[0]
            gcount[ids] = 1

with open(sys.argv[2],'r') as f:
    for line in f:
        if line.strip():
            if line[0] != '#':
                line = line.strip().split()
                total = 1
                if line[0] != 'Geneid':
                    for num in line[6:]:
                        total += int(num)
                    gcount[line[0]] = total


# overlap pair genes

lis = []
overlap_p = []
with open(sys.argv[3],'r') as f:
    for line in f:
        if line.strip():
            if line[0] != '#':
                line = line.strip().split()
                g1,g2 = line[0],line[1]
                overlap_p.append([g1,g2])
                value = abs(gcount[g1]-gcount[g2])
                lis.append(value)


lis  = sorted(lis,reverse=True)

total = sum(lis)

num = 0
for i in lis:
    num += int(i)
    if num >= total*0.99:
        n95 = i
        break


for v in overlap_p:
    gd = {}
    g1,g2 = v[0],v[1]
    gd[g1] = gcount[g1]
    gd[g2] = gcount[g2]
    if gcount[g1]/gcount[g2] >3 or gcount[g1]/gcount[g2] < 0.333:
        if abs(gcount[g1]-gcount[g2]) > n95:
            print(min(gd,key=gd.get))
