import sys
import numpy as np

lis = []
with open(sys.argv[1],'r') as f:
    for line in f:
        line = line.strip().split()
        if line[0] != 'gene_id':
            lis.append(int(line[2]))


lis = sorted(lis)

lis = np.array(lis)

percentile1 = np.percentile(lis,90)
percentile2 = np.percentile(lis,95)
percentile3 = np.percentile(lis,98)
percentile4 = np.percentile(lis,99)
maxnum = np.percentile(lis,100)

print('specie','90','95','98','99','max',sep='\t')
print(sys.argv[2],percentile1,percentile2,percentile3,percentile4,maxnum,sep='\t')
