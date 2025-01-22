import sys

dicts = {}
with open(sys.argv[1],'r') as f:
    for line in f:
        if line.strip():
            if line[0] != '#':
                line = line.strip().split('\t')
                if line[2] == 'mRNA':
                    idx = line[8].find('Parent')
                    mdx = line[8].find('ID')
                    geneID = line[8][idx:].split('=')[1].split(';')[0]
                    mrnaID = line[8][mdx:].split('=')[1].split(';')[0]
                    dicts[mrnaID] = geneID

overlocus,overcds = [],[]
with open(sys.argv[2],'r') as f:
    for line in f:
        line = line.strip().split('\t')
        if line[0] != '#id1':
            if float(line[5]) > 0.5 or float(line[6]) > 0.5:
                if dicts[line[0]] != dicts[line[1]]:
                    overlocus.append(dicts[line[0]])

            if float(line[10]) > 0.5 or float(line[11]) > 0.5:
                if dicts[line[0]] != dicts[line[1]]:
                    overcds.append(dicts[line[0]])


overlocus = list(set(overlocus))

overcds = list(set(overcds))


print('locus overlap: '+str(len(overlocus)),'CDS overlap: '+str(len(overcds)),sep='\t')


