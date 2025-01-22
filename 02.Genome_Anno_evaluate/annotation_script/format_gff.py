import sys




with open(sys.argv[1],'r') as f:
    for line in f:
        if line.strip():
            if line[0] != '#':
                line = line.strip().split('\t')
                if line[2] == 'mRNA':
                    ID = line[8][line[8].find('ID'):].split('=')[1].split(';')[0]
                    geneID = line[8][line[8].find('Parent'):].split('=')[1].split(';')[0] if 'Parent' in line[8] else line[8][line[8].find('gene_id'):].split('=')[1].split(';')[0]
                    print('\t'.join(x for x in line[0:8]),'ID='+ID+';Parent='+geneID,sep='\t')
                if line[2] == 'CDS':
                    print('\t'.join(x for x in line[0:8]),'Parent='+ID+';',sep='\t')

                if line[2] == 'UTR5':
                    print('\t'.join(x for x in line[0:8]),'Parent='+ID+';',sep='\t')

                if line[2] == 'UTR3':
                    print('\t'.join(x for x in line[0:8]),'Parent='+ID+';',sep='\t')

