import sys


filter_id = []
with open(sys.argv[1],'r') as f:
    for line in f:
        if line[0] != '#':
            line = line.strip()
            filter_id.append(line)



with open(sys.argv[2],'r') as f:
    gene_d = {}
    for line in f:
        if line.strip():
            if line[0] != '#':
                line = line.strip().split('\t')
                if line[2] == 'mRNA':
                    ID = line[8][line[8].find('ID'):].split('=')[1].split(';')[0]
                    #PID = line[8][line[8].find('Parent'):].split('=')[1].split(';')[0]
                    gene_d[ID] = [line]
                else:
                    gene_d[ID].append(line)



for k,v in gene_d.items():
    if k not in filter_id:
        for line in v:
            print('\t'.join(x for x in line))
