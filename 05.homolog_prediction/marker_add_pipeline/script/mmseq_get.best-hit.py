import sys


def getblast(F):
    dicts = {}
    with open(F,'r') as f:
        for line in f:
            line = line.strip().split()
            if line[0] in dicts.keys():
                dicts[line[0]][line[1]] = float(line[11])
            else:
                dicts[line[0]] = {line[1]:float(line[11])}
    return dicts

def getop(dicts):
    existgene = []
    for k,v in dicts.items():
        mc = max(v, key=v.get)
        ids = k.split('-')[0]
        if ids[0:4] in mc:
            print(k,mc,sep='\t')
            existgene.append(mc)
    uniqhit = list(set(existgene))
    return uniqhit

def getgene(F):
    genelis = []
    with open(F,'r') as f:
        for line in f:
            if line[0] == '>':
                name = line[1:].strip()
                genelis.append(name)
    return genelis


def main():
    blasthit = getblast(sys.argv[1])
    genelis = getgene(sys.argv[2])

    uniqhit = getop(blasthit)

    nonhit = open('non-existing-predict-gene.lis','w')
    for i in genelis:
        if i[0:4] not in '\t'.join(x for x in uniqhit):
            nonhit.write(i+'\n')


if __name__ == '__main__':
    main()
