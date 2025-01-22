import sys
import numpy as np

lis = []
with open(sys.argv[1],'r') as f:
    for line in f:
        if line.strip():
            if line[0] != '#':
                line = line.strip().split('\t')
                if line[2] == 'transcript':
                    idx = line[8].find('cov')
                    TPM = line[8][idx:].split('"')[1]
                    lis.append(float(TPM))



lis = sorted(lis)
lis = np.array(lis)

media = np.median(lis)
mean = np.mean(lis)

p25 = np.percentile(lis,25)

values, counts = np.unique(lis, return_counts=True)
mode = values[np.argmax(counts)]


print('media: '+str(format(media,'.2f')),'p25: '+str(format(p25,'.2f')),'mean: '+str(format(mean,'.2f')),'mode: '+str(format(mode,'.2f')),sep='\t')

