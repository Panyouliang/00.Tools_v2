import sys


for line in open(sys.argv[1],'r'):
    if line[0] != '#':
        line = line.strip().split('\t')
        if line[2] == 'mRNA':
            ID = line[8].split('=')[1].split(';')[0]
            align = line[8].split('=')[-1].split(';')[0]
            if float(align) > 0:
                print(ID)
