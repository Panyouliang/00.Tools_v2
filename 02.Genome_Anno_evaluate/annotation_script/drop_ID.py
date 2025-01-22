#! /usr/bin/env python3
# -*- coding:utf-8 -*-

import sys


tidlist = []
for line in open(sys.argv[1],'r'):
    if line.strip():
        line = line.strip().split()
        tidlist.append(line[0])


def readgff(F,lis):
    gene_d,all_d = {},{}
    with open(F,'r') as f:
        for line in f:
            if line[0] != '#':
                if line.strip():
                    line = line.strip().split('\t')
                    if line[2] == 'gene':
                        gid = line[8][line[8].find('ID'):].split('=')[1].split(';')[0]
                        gene_d[gid] = line
                    else:
                        if line[2] == 'mRNA':
                            tid = line[8][line[8].find('ID'):].split('=')[1].split(';')[0]
                            if 'Parent' in line[8]:
                                gid = line[8][line[8].find('Parent'):].split('=')[1].split(';')[0]
                            else:
                                gid = line[8][line[8].find('gene_id'):].split('=')[1].split(';')[0]

                            if tid not in lis:
                                all_d[gid] = [line]
                            else:
                                pass
                        else:
                            ptid = line[8][line[8].find('Parent'):].split('=')[1].split(';')[0]
                            if ptid not in lis:
                                all_d[gid].append(line)
    return gene_d,all_d



def main():
    gene_d,all_d = readgff(sys.argv[2],tidlist)
    for k,v in all_d.items():
        print('\t'.join(x for x in gene_d[k]))
        for line in v:
            print('\t'.join(x for x in line))



if __name__ == '__main__':
    main()



