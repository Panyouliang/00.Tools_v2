#! /usr/bin/env python3
# -*- coding:utf-8 -*-


import sys


def readgff(F):
    featured = {}
    gened = {}
    with open(F,'r') as f:
        for line in f:
            if line[0] == '#':
                pass
            else:
                line = line.strip().split()
                if line[2] == 'mRNA':
                    geneN = line[8][line[8].find('Parent'):].split('=')[1].split(';')[0]
                    rnaN  = line[8][line[8].find('ID'):].split('=')[1].split(';')[0]
                    featured[rnaN] = []
                    gened[rnaN] = [line]
                elif line[2] == 'gene':
                    pass
                else:
                    featured[rnaN].append([int(line[3]),int(line[4])])
                    gened[rnaN].append(line)

    feature_adj = {}
    for k,v in featured.items():
        vsort = sorted(v, key=lambda x:int(x[1]))
        feature_adj[k] = merge_adjacent_coordinates(vsort)

    return feature_adj, gened


def merge_adjacent_coordinates(coordinates):
    merged_coordinates = []
    if not coordinates:
        return merged_coordinates
    current_start, current_end = coordinates[0]
    for coord in coordinates[1:]:
        if coord[0] <= current_end + 1:
            current_end = max(current_end, coord[1])
        else:
            merged_coordinates.append([current_start, current_end])
            current_start, current_end = coord
    merged_coordinates.append([current_start, current_end])
    return merged_coordinates



if __name__ == '__main__':

    feature, gened = readgff(sys.argv[1])
    featureold, genedold = readgff(sys.argv[2])

    ## extract single exon gene ID
    single = []
    for k,v in feature.items():
        if len(v) == 1:
            single.append(k)


    ## get suited conditions geneID
    suited = []
    for k,v in genedold.items():
        if k in single:
            for line in v:
                if line[2] == 'mRNA':
                    rnaN  = line[8][line[8].find('ID='):].split('=')[1].split(';')[0]
                    aa    = line[8][line[8].find(';aa='):].split('=')[1].split(';')[0]
                    score = line[8][line[8].find(';score='):].split('=')[1].split(';')[0]
                    start = line[8][line[8].find(';start='):].split('=')[1].split(';')[0]
                    stop  = line[8][line[8].find(';stop='):].split('=')[1].split(';')[0]

                    if float(int(score)/int(aa)) >= float(sys.argv[3]):
                        if start == 'M' and stop == '*':
                        #suited.append(rnaN)
                            print('\t'.join(x for x in line))

    ## get single exon genes gff file
    #for k,v in gened.items():
    #    if k in suited:
    #        for line in v:
    #            print('\t'.join(x for x in line))








