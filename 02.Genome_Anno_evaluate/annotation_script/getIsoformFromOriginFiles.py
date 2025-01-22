import sys




def getGeneN(F):
    geneset = []
    with open(F,'r') as f:
        for line in f:
            line = line.strip().split()
            if line[2] == 'mRNA':
                geneN = line[8][line[8].find('Parent'):].split('=')[1].split(';')[0]
                geneset.append(geneN)

    return geneset

def readgff(F):
    genedict = {}
    genetagline = {}
    with open(F,'r') as f:
        for line in f:
            if line.startswith('#'):
                pass
            else:
                line = line.strip().split()
                if line[2] == 'gene':
                    geneN = line[8][line[8].find('ID'):].split('=')[1].split(';')[0]
                    genedict[geneN] = {}
                    genetagline[geneN] = line
                elif line[2] == 'mRNA':
                    rnaN  = line[8][line[8].find('ID'):].split('=')[1].split(';')[0]
                    genedict[geneN][rnaN] = [line]
                else:
                    prnaN = line[8][line[8].find('Parent'):].split('=')[1].split(';')[0]
                    genedict[geneN][prnaN].append(line)
    return genedict, genetagline


if __name__ == '__main__':
    geneset = getGeneN(sys.argv[1])     #input the filter gff
    genedict,tag = readgff(sys.argv[2]) #input the old gff

    for gene in geneset:
        tag[gene][1] = 'GeMoMa'
        print('\t'.join(x for x in tag[gene][0:8]),';'.join(x for x in tag[gene][8].split(';')[0:2])+';',sep='\t')
        for rnaN,transcript in genedict[gene].items():
            for line in transcript:
                line[1] = 'GeMoMa'
                if line[2] == 'mRNA':
                    print('\t'.join(x for x in line[0:8]),'ID='+rnaN+';'+'Parent='+gene+';',sep='\t')
                else:
                    print('\t'.join(x for x in line[0:8]),'Parent='+rnaN+';',sep='\t')


