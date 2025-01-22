import sys



def is_overlap(interval1, interval2):
    return interval1[0] <= interval2[1] and interval2[0] <= interval1[1]


def overlap_length(interval1, interval2):
    if is_overlap(interval1, interval2):
        return min(interval1[1], interval2[1]) - max(interval1[0], interval2[0])
    return 0


def get_locus_cluster(anno_gene_site,isoform_gene_site):
    gene_cluster = {}
    for chrs, gene_loci in anno_gene_site.items():
        for geneID, site in gene_loci.items():
            gene_cluster[geneID] = [geneID]
            for gene_name, loci in isoform_gene_site[chrs].items():
                t_overlap_l = 0
                t_q_l = site[1] - site[0]
                t_t_l = loci[1] - loci[0]
                if is_overlap(site,loci):
                    t_overlap_l += overlap_length(site,loci)
                else:
                    pass

                q_rate = t_overlap_l/t_q_l
                t_rate = t_overlap_l/t_t_l
                if q_rate > 0.5 and t_rate > 0.5:
                    if geneID not in gene_cluster:
                        gene_cluster[geneID] = [gene_name]
                    else:
                        gene_cluster[geneID].append(gene_name)

    return gene_cluster



# 构造函数，计算一个locus的exon区域是否和其他locus存在overlap
def search_exon_overlap_rate(iso_sites,target_sites):

    query_length = sum(interval[1]-interval[0] for interval in iso_sites)
    target_length= sum(interval[1]-interval[0] for interval in target_sites)
    total_overlap_length = 0

    for query_loc in iso_sites:
        for target_loc in target_sites:
            if is_overlap(query_loc,target_loc):
                total_overlap_length += overlap_length(query_loc,target_loc)

    rate_query = total_overlap_length / query_length
    rate_target= total_overlap_length / target_length

    return rate_query,rate_target

def filter_overlap_isos(gene_cluster,cds_dict1,cds_dict2):
    filter_lis = []
    for GeneID, isogeneNlis in gene_cluster.items():
        for qisoN, qcds_site in cds_dict1[GeneID].items():
            for isogeneN in isogeneNlis:
                if isogeneN in cds_dict2:
                    for isornaN,iso_site in cds_dict2[isogeneN].items():
                        rate_q,rate_t = search_exon_overlap_rate(qcds_site,iso_site)
                        if rate_q > 0.95 and rate_t > 0.95:
                            filter_lis.append(isornaN)

    return filter_lis



def readgff(F,spes):
    gff_dicts = {}
    gene_loci = {}
    isoform_loci = {}
    exon_dict = {}
    with open(F,'r') as f:
        for line in f:
            if line.strip():
                if line[0] != '#':
                    line = line.strip().split('\t')
                    if line[2] == 'mRNA':
                        chrs = line[0]
                        geneID = line[8].split('Parent=')[1].split(';')[0]+'#'+spes+'|'+line[6]
                        name = geneID+'|'+line[8].split(';')[0].split('=')[1]+'|'+spes+'|'+line[6]
                        gff_dicts[name] = [line]

                        if chrs not in gene_loci:
                            gene_loci[chrs] = {geneID:[int(line[3]),int(line[4])]}
                        else:
                            if geneID not in gene_loci[chrs]:
                                gene_loci[chrs][geneID] = [int(line[3]),int(line[4])]
                            else:
                                gene_loci[chrs][geneID].extend([int(line[3]),int(line[4])])

                        if geneID not in isoform_loci:
                            isoform_loci[geneID] = {name:[int(line[3]),int(line[4])]}
                        else:
                            isoform_loci[geneID][name] = [int(line[3]),int(line[4])]

                    else:
                        if (line[2] != 'gene'):
                            gff_dicts[name].append(line)
                            if line[2] != 'exon':
                                if geneID not in exon_dict:
                                    exon_dict[geneID] = {name:[[int(line[3]),int(line[4])]]}
                                else:
                                    if name not in exon_dict[geneID]:
                                        exon_dict[geneID][name] = [[int(line[3]),int(line[4])]]
                                    else:
                                        exon_dict[geneID][name].append([int(line[3]),int(line[4])])

    gene_loci_ = {}
    for chrs,gene_loc in gene_loci.items():
        gene_loci_[chrs] = {}
        for geneID,loc in gene_loc.items():
            gene_loci_[chrs][geneID] = (min(loc),max(loc))

    gene_loci_sort = {}
    for chrs,gene_loc in gene_loci_.items():
        gene_loci_sort[chrs] = dict(sorted(gene_loc.items(), key=lambda item:item[1]))

    return gff_dicts, gene_loci_sort, isoform_loci, exon_dict



def merge_two_level_dicts(dict1, dict2):
    merged_dict = dict1.copy()
    for outer_key, inner_dict in dict2.items():
        if outer_key in merged_dict:
            merged_dict[outer_key].update(inner_dict)
        else:
            merged_dict[outer_key] = inner_dict
    merged_dict_sort = {}
    for chrs,gene_locs in merged_dict.items():
        merged_dict_sort[chrs] = dict(sorted(gene_locs.items(), key=lambda item:item[1]))

    return merged_dict_sort


def main():
    gff_dicts,gene_loci,isoform_loci,exon_dicts = readgff(sys.argv[1],'q')     # GFF
    gff_dicts1,gene_loci1,isoform_loci1,exon_dicts1 = readgff(sys.argv[2],'t') # isoforms.gff

    gff_dicts.update(gff_dicts1)
    isoform_loci.update(isoform_loci1)

    gene_2_gene = get_locus_cluster(gene_loci, gene_loci1)
    filter_lis = filter_overlap_isos(gene_2_gene,exon_dicts,exon_dicts1)
    for geneID,IDlis in gene_2_gene.items():
        isonum = 0
        for ID in IDlis:
            for iso in isoform_loci[ID].keys():
                if iso not in filter_lis:
                    isonum += 1
                    for line in gff_dicts[iso]:
                        geneID = geneID.split('#')[0]
                        line[1] = 'Pogona_vitticeps'
                        if line[2] == 'mRNA':
                            print('\t'.join(line[0:8]),'ID='+geneID+'.'+str(isonum)+';Parent='+geneID+';',sep='\t')
                        else:
                            print('\t'.join(line[0:8]),'Parent='+geneID+'.'+str(isonum)+';',sep='\t')



if __name__ == '__main__':
    main()

