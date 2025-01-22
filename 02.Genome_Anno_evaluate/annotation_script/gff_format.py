import sys



with open(sys.argv[1],'r') as f:
    for line in f:
        if line[0] != '#':
            line = line.strip().split('\t')
            if line[2] == 'gene':
                gid = line[8].split('=')[1].split(';')[0]
            if line[2] == 'mRNA':
                tid = line[8].split('=')[1].split(';')[0]
                print('\t'.join(x for x in line[0:8]),'ID='+gid+';',sep='\t')
            else:
                line[1] = 'ReName'
                if line[2] == 'UTR5' or line[2] == 'five_prime_UTR':
                    print('\t'.join(x for x in line[0:2]),'UTR5','\t'.join(x for x in line[3:8]),'Parent='+gid+';',sep='\t')
                if line[2] == 'UTR3' or line[2] == 'three_prime_UTR':
                    print('\t'.join(x for x in line[0:2]),'UTR3','\t'.join(x for x in line[3:8]),'Parent='+gid+';',sep='\t')
                if line[2] == 'CDS':
                    print('\t'.join(x for x in line[0:8]),'Parent='+gid+';',sep='\t')
