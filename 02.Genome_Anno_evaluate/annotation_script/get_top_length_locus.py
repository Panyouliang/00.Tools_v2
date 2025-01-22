import sys


def get_pair(F):
    locus_p = []
    with open(F,'r') as f:
        for line in f:
            if line[0] != '#':
                line = line.strip().split()
                locus_p.append([line[0],line[1]])

    return locus_p


def read_gff(F):
    gene_d = {}
    gene_l = {}
    with open(F,'r') as f:
        for line in f:
            if line[0] != '#':
                line = line.strip().split()
                if line[2] == 'mRNA':
                    ID = line[8][line[8].find('ID'):].split('=')[1].split(';')[0]
                    gene_d[ID] = [line]
                    g_len = 0
                else:
                    PID= line[8][line[8].find('Parent'):].split('=')[1].split(';')[0]
                    gene_d[PID].append(line)
                    if line[2] == 'CDS':
                        g_len += int(line[4])-int(line[3]) + 1
                        gene_l[PID] = g_len
    return gene_d,gene_l



def main():
    over_p = get_pair(sys.argv[1])
    gene_d,gene_l = read_gff(sys.argv[2])
    ids = []
    for pair in over_p:
        if pair[0] in gene_l.keys() and pair[1] in gene_l.keys():
            if gene_l[pair[0]] >= gene_l[pair[1]]:
                ids.append(pair[1])
            if gene_l[pair[0]] < gene_l[pair[1]]:
                ids.append(pair[0])

    for k,v in gene_d.items():
        if k not in ids:
            for line in v:
                print('\t'.join(x for x in line))


if __name__ == '__main__':
    main()



