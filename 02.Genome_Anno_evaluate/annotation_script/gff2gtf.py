#! /usr/bin/env python3
# -*- coding:utf-8 -*-

import sys

def read_gff(F):
    if F[-3:] == '.gz':
        f = gzip.open(F,'rt')
    else:
        f = open(F,'r')
    gene_d,UTR5_d,CDS_d,exon_d,UTR3_d = {},{},{},{},{}
    g_exon = {}
    allfea = {}
    for line in f.readlines():
        if line.strip():
            if line[0] != '#':
                line = line.strip().split('\t')
                if line[2] == 'mRNA':
                    ID = line[8][line[8].find('ID'):].split('=')[1].split(';')[0]+' '+line[6]
                    if 'geneID' in line[8]:
                        pID = line[8][line[8].find('geneID'):].split('=')[1].split(';')[0]
                    if 'Parent' in line[8]:
                        pID = line[8][line[8].find('Parent'):].split('=')[1].split(';')[0]

                    gene_d[pID+' '+ID] = line
                else:
                    fID = line[8][line[8].find('Parent'):].split('=')[1].split(';')[0]+' '+line[6]
                    ppID = pID+' '+fID
                    if ppID in allfea.keys():
                        allfea[ppID].append(line)
                    else:
                        allfea[ppID] = [line]

                    if line[2] == 'exon':
                        g_exon[ppID] = 'exon'

    sort_fea = {}
    for k,v in allfea.items():
        if k.endswith('+'):
            sv = sorted(v, key=lambda x: int(x[3]),reverse=False)
            sort_fea[k] = sv
        else:
            sv = sorted(v, key=lambda x: int(x[3]),reverse=True)
            sort_fea[k] = sv

    for k,v in sort_fea.items():
        pID = k.split()[0]
        for line in v:
            if line[2] == 'UTR5':
                if pID not in UTR5_d.keys():
                    UTR5_d[pID] = [line]
                else:
                    UTR5_d[pID].append(line)
                if line[2] == 'CDS':
                    if pID not in CDS_d.keys():
                        CDS_d[pID] = [line]
                    else:
                        CDS_d[pID].append(line)
                if line[2] == 'UTR3':
                    if pID not in UTR3_d.keys():
                        UTR3_d[pID] = [line]
                    else:
                        UTR3_d[pID].append(line)

    return sort_fea,gene_d,g_exon




def main():

    d,gd,g_exon = read_gff(sys.argv[1])

    for k,v in d.items():
        gid = k.split()[0]
        tid = k.split()[1]
        exon,cds = 0,0
        gd[k][1] = 'ReName'
        print('\t'.join(x for x in gd[k][0:2]),'gene','\t'.join(x for x in gd[k][3:8]),'gene_id "'+gid+'"; gene_name "'+gid+'"; gene_biotype "protein-coding";',sep='\t')
        print('\t'.join(x for x in gd[k][0:2]),'transcript','\t'.join(x for x in gd[k][3:8]),'gene_id "'+gid+'"; transcript_id "'+tid+'"; gene_name "'+gid+'"; gene_biotype "protein-coding";',sep='\t')
        for line in v:
            line[1] = 'ReName'
            #print('\t'.join(x for x in line[0:2]),'exon','\t'.join(x for x in line[3:8]),'gene_id "'+gid+'"; transcript_id "'+gid+'"; gene_name "'+gid+'"; exon_id "'+gid+'.exon.'+str(exon)+'"; transcript_name "'+gid+'"; gene_biotype "protein-coding";',sep='\t') 
            if k in g_exon.keys():
                if g_exon[k] == 'exon':
                    if line[2] == 'exon':
                        exon += 1
                        print('\t'.join(x for x in line[0:8]),'gene_id "'+gid+'"; transcript_id "'+tid+'"; gene_name "'+gid+'"; exon_id "'+tid+'.exon.'+str(exon)+'"; transcript_name "'+tid+'"; gene_biotype "protein-coding";',sep='\t')
                    if line[2] == 'CDS':
                        cds += 1
                        print('\t'.join(x for x in line[0:8]),'gene_id "'+gid+'"; transcript_id "'+tid+'"; gene_name "'+gid+'"; exon_id "'+tid+'.cds.'+str(cds)+'"; transcript_name "'+tid+'"; gene_biotype "protein-coding";',sep='\t')
            else:
                exon += 1
                cds += 1
                print('\t'.join(x for x in line[0:2]),'exon','\t'.join(x for x in line[3:8]),'gene_id "'+gid+'"; transcript_id "'+tid+'"; gene_name "'+gid+'"; exon_id "'+tid+'.exon.'+str(exon)+'"; transcript_name "'+tid+'"; gene_biotype "protein-coding";',sep='\t')
                print('\t'.join(x for x in line[0:2]),'CDS','\t'.join(x for x in line[3:8]),'gene_id "'+gid+'"; transcript_id "'+tid+'"; gene_name "'+gid+'"; exon_id "'+tid+'.exon.'+str(cds)+'"; transcript_name "'+tid+'"; gene_biotype "protein-coding";',sep='\t')



if __name__ == '__main__':
    main()

