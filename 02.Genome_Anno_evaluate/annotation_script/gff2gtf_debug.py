#! /usr/bin/env python3
# -*- coding:utf-8 -*-

import sys
from itertools import zip_longest


def readgff(F):

    geneline,transcript_dict,genefeature_dict,exon_coord = {},{},{},{}
    gname = {}

    with open(F,'r') as f:
        for line in f:
            if line.strip():
                if line[0] != '#':
                    line = line.strip().split('\t')
                    if line[2] == 'mRNA':
                        rnaID = line[8].split('=')[1].split(';')[0]
                        geneID = line[8][line[8].find('Parent'):].split('=')[1].split(';')[0]+' '+line[6]
                        genetag = line[0:8].copy()
                        genetag[1] = 'ReName'
                        if 'Name' in line[8]:
                            name = line[8].split('Name=')[1].split(';')[0]
                            if name != 'Unknown':
                                gname[geneID.split()[0]] = name
                            else:
                                gname[geneID.split()[0]] = geneID.split()[0]

                        geneline[geneID] = genetag.copy()
                        geneline[geneID][2] = 'gene'
                        transcript_dict[geneID] = {rnaID:genetag.copy()}
                        transcript_dict[geneID][rnaID][2] = 'transcript'
                        genefeature_dict[geneID] = {rnaID:[]}
                        exon_coord[geneID] = {rnaID:[]}

                    else:
                        if line[2] != 'exon' and line[2] != 'gene':
                            prnaID = line[8][line[8].find('Parent'):].split('=')[1].split(';')[0]
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

    return geneline,transcript_dict,genefeature_dict_sort,exon_coord_merge,gname



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
    geneline,transcript_dict,feature,exon_coord,gname = readgff(sys.argv[1])

    utr_5,utr_3 = ['five_prime_utr', 'five_prime_utr', 'UTR5'],['three_prime_utr','three_prime_UTR','UTR3']

    for geneid,coord in exon_coord.items():
        gid = geneid.split()[0]
        print('\t'.join(x for x in geneline[geneid]),'gene_id "'+gid+'"; gene_name "'+gid+'"; gene_biotype "protein-coding";',sep='\t')

        for rnaid in coord.keys():
            print('\t'.join(x for x in transcript_dict[geneid][rnaid]),'gene_id "'+gid+'"; transcript_id "'+rnaid+'"; gene_name "'+gid+'"; gene_biotype "protein-coding";',sep='\t')
            gffline = transcript_dict[geneid][rnaid].copy()
            gffline[2] = 'exon'

            exon_num,cds_num,utr5_num,utr3_num = 0,0,0,0
            for cds,exon in zip_longest(feature[geneid][rnaid], coord[rnaid], fillvalue=None):
                if exon != None:
                    exon_num += 1
                    exonline = gffline.copy()
                    exonline[3] = exon[0]
                    exonline[4] = exon[1]
                    print('\t'.join(str(x) for x in exonline),'gene_id "'+gid+'"; transcript_id "'+rnaid+'"; gene_name "'+gid+'"; exon_id "'+rnaid+'.exon.'+str(exon_num)+'"; transcript_name "'+rnaid+'"; gene_biotype "protein-coding";',sep='\t')
                if cds != None:
                    if cds[2] == 'CDS':
                        cds_num += 1 
                        print('\t'.join(str(x) for x in cds),'gene_id "'+gid+'"; transcript_id "'+rnaid+'"; gene_name "'+gid+'"; exon_id "'+rnaid+'.cds.'+str(cds_num)+'"; transcript_name "'+rnaid+'"; gene_biotype "protein-coding";',sep='\t')
                    if cds[2] in utr_5:
                        utr5_num += 1
                        print('\t'.join(str(x) for x in cds),'gene_id "'+gid+'"; transcript_id "'+rnaid+'"; gene_name "'+gid+'"; exon_id "'+rnaid+'.utr5.'+str(utr5_num)+'"; transcript_name "'+rnaid+'"; gene_biotype "protein-coding";',sep='\t')
                    if cds[2] in utr_3:
                        utr3_num += 1
                        print('\t'.join(str(x) for x in cds),'gene_id "'+gid+'"; transcript_id "'+rnaid+'"; gene_name "'+gid+'"; exon_id "'+rnaid+'.utr3.'+str(utr3_num)+'"; transcript_name "'+rnaid+'"; gene_biotype "protein-coding";',sep='\t')



if __name__ == '__main__':
    main()

