#! /usr/bin/env python3
# -*- coding:utf-8 -*-

import sys
from itertools import zip_longest


def readgff(F):

    geneline,transcript_dict,genefeature_dict,exon_coord = {},{},{},{}

    with open(F,'r') as f:
        for line in f:
            if line.strip():
                if line[0] != '#':
                    line = line.strip().split('\t')
                    if line[2] == 'mRNA':
                        rnaID = line[8].split('=')[1].split(';')[0]
                        idx = line[8].find('Parent') if 'Parent' in line[8] else line[8].find('gene_id')
                        geneID = line[8][idx:].split('=')[1].split(';')[0]+' '+line[6]

                        genetag = line[0:8].copy()
                        genetag[1] = 'ReName'

                        geneline[geneID] = genetag.copy()
                        geneline[geneID][2] = 'gene'
                        transcript_dict[geneID] = {rnaID:genetag.copy()}
                        transcript_dict[geneID][rnaID][2] = 'mRNA'
                        genefeature_dict[geneID] = {rnaID:[]}
                        exon_coord[geneID] = {rnaID:[]}
                    elif line[2] == 'gene':
                        pass

                    else:
                        if line[2] != 'exon':
                            prnaID = line[8].split('=')[1].split(';')[0]
                            genefeature_dict[geneID][prnaID].append(line[0:8])
                            exon_coord[geneID][prnaID].append([int(line[3]),int(line[4])])

    exon_coord_merge = {}
    for gene_id,v in exon_coord.items():
        exon_coord_merge[gene_id] = {}
        for rna_id,coord in v.items():
            if gene_id.endswith('-'):
                exon_coord_merge[gene_id][rna_id] = sorted(merge_adjacent_coordinates(coord),key=lambda x:x[0],reverse=True)
            else:
                exon_coord_merge[gene_id][rna_id] = sorted(merge_adjacent_coordinates(coord),key=lambda x:x[0],reverse=False)

    genefeature_dict_sort = {}
    for geneid,v in genefeature_dict.items():
        genefeature_dict_sort[geneid] = {}
        for rnaid,tag in v.items():
            if geneid.endswith('-'):
                genefeature_dict_sort[geneid][rnaid] = sorted(tag, key=lambda x:int(x[3]),reverse=True)
            else:
                genefeature_dict_sort[geneid][rnaid] = sorted(tag, key=lambda x:int(x[3]),reverse=False)

    return geneline,transcript_dict,genefeature_dict_sort,exon_coord_merge



def merge_adjacent_coordinates(coordinates):
    merged_coordinates = []
    if not coordinates:
        return merged_coordinates

    current_start, current_end = coordinates[0]
    for coord in coordinates[1:]:
        if coord[0] == current_end + 1:
            current_end = max(current_end, coord[1])
        else:
            merged_coordinates.append([current_start, current_end])
            current_start, current_end = coord

    merged_coordinates.append([current_start, current_end])

    return merged_coordinates


def main():
    geneline,transcript_dict,feature,exon_coord = readgff(sys.argv[1])

    for geneid,coord in exon_coord.items():
        gid = geneid.split()[0]
        print('\t'.join(x for x in geneline[geneid]),'ID='+gid+';',sep='\t')

        for rnaid in coord.keys():
            print('\t'.join(x for x in transcript_dict[geneid][rnaid]),'ID='+rnaid+';'+'Parent='+gid+';',sep='\t')
            gffline = transcript_dict[geneid][rnaid].copy()
            gffline[2] = 'exon'

            exon_num,cds_num,utr5_num,utr3_num = 0,0,0,0
            for cds,exon in zip_longest(feature[geneid][rnaid], coord[rnaid], fillvalue=None):
                if exon != None:
                    exon_num += 1
                    exonline = gffline.copy()
                    exonline[3] = exon[0]
                    exonline[4] = exon[1]
                    print('\t'.join(str(x) for x in exonline),'Parent='+rnaid+';',sep='\t')
                if cds != None:
                    if cds[2] == 'CDS':
                        cds_num += 1 
                        print('\t'.join(str(x) for x in cds),'Parent='+rnaid+';',sep='\t')
                    if cds[2] == 'UTR5':
                        utr5_num += 1
                        print('\t'.join(str(x) for x in cds),'Parent='+rnaid+';',sep='\t')
                    if cds[2] == 'UTR3':
                        utr3_num += 1
                        print('\t'.join(str(x) for x in cds),'Parent='+rnaid+';',sep='\t')



if __name__ == '__main__':
    main()

