import sys



def is_overlap(interval1, interval2):
    return interval1[0] <= interval2[1] and interval2[0] <= interval1[1]


def overlap_length(interval1, interval2):
    if is_overlap(interval1, interval2):
        return min(interval1[1], interval2[1]) - max(interval1[0], interval2[0])
    return 0

def cluster_obit_site(gene_dicts):
    gene_cluster_dicts = {}
    for chrs,gene_loci in gene_dicts.items():
        for geneID, site in gene_loci.items():
            gene_cluster_dicts[geneID] = [geneID]
            for name, loci in gene_loci.items():
                if (geneID != name) and (geneID.split('|')[-1] == name.split('|')[-1]):
                    total_overlap_length = 0
                    total_small_length = site[1]-site[0]
                    total_large_length = loci[1]-loci[0]
                    if is_overlap(site,loci):
                        total_overlap_length += overlap_length(site,loci)
                    else:
                        pass

                    if total_overlap_length > 0:
                        overlap_rate_small = total_overlap_length / total_small_length
                        overlap_rate_large = total_overlap_length / total_large_length
                        if (overlap_rate_small >0.5) and (overlap_rate_large > 0.5):
                            if geneID not in gene_cluster_dicts:
                                gene_cluster_dicts[geneID] = [name]
                            else:
                                gene_cluster_dicts[geneID].append(name)
            #print(geneID,gene_cluster_dicts[geneID])

    return gene_cluster_dicts



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



def readgff(F,spes):
    gff_dicts = {}
    gene_loci = {}
    isoform_loci = {}
    exon_dict = {}
    with open(F,'r') as f:
        for line in f:
            if line[0] != '#':
                line = line.strip().split('\t')
                if line[2] == 'transcript':
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
                    if line[2] == 'exon':
                        gff_dicts[name].append(line)
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



def readcsv(files):
    gff_all_dicts = {}
    gene_locus = {}
    isoform_locus = {}
    exon_dicts = {}
    num = 0
    for file in files:
        num += 1
        gff_dicts,gene_loci,isoform_loci,exon_dict = readgff(file,str(num))

        gff_all_dicts.update(gff_dicts)
        gene_locus = merge_two_level_dicts(gene_locus, gene_loci)
        isoform_locus.update(isoform_loci)
        exon_dicts.update(exon_dict)

    return gff_all_dicts,gene_locus,isoform_locus,exon_dicts



def stat_isoforms_overlap(gene_clusterd, exon_dicts, cutrate):

    filter_ids = []

    for qgeneID,geneID_lis in gene_clusterd.items():
        for tgeneID in geneID_lis:
            if qgeneID != tgeneID:
                for qisoID,qisoloc in exon_dicts[qgeneID].items():
                    if qisoID not in filter_ids:
                        for tisoID,tisoloc in exon_dicts[tgeneID].items():
                            if tisoID not in filter_ids:
                                q_rate,t_rate = search_exon_overlap_rate(qisoloc,tisoloc)
                                if q_rate > cutrate and t_rate > cutrate:
                                    filter_ids.append(tisoID)

    return filter_ids








def print_contain_isoform(gene_clusterd,exon_dicts,gff_dicts,filter_ids):
    geneID_lis = []
    novel_gene = 0

    gene_isos = {}

    for geneID,IDlis in gene_clusterd.items():

        IDs = geneID.split('#')[0]
        if IDs not in gene_isos:
            gene_isos[IDs] = 0

        if 'novel' not in geneID:
            if geneID not in geneID_lis:
                for gID in IDlis:
                    geneID_lis.append(gID)
                    for iso_name in exon_dicts[gID].keys():
                        if iso_name not in filter_ids:
                            gene_isos[IDs] = gene_isos[IDs] + 1
                            for line in gff_dicts[iso_name]:
                                gname = geneID.split('#')[0]
                                if line[2] == 'transcript':
                                    print('\t'.join(line[0:8]),'ID='+gname+'.'+str(gene_isos[IDs])+';Parent='+gname+';',sep='\t')
                                else:
                                    print('\t'.join(line[0:8]),'Parent='+gname+'.'+str(gene_isos[IDs])+';',sep='\t')
        else:
            novel_gene += 1
            if geneID not in geneID_lis:
                for gID in IDlis:
                    geneID_lis.append(gID)
                    for iso_name in exon_dicts[gID].keys():
                        if iso_name not in filter_ids:
                            gene_isos[IDs] = gene_isos[IDs] + 1
                            for line in gff_dicts[iso_name]:
                                gname = 'novel_gene'+str(novel_gene)
                                if line[2] == 'transcript':
                                    print('\t'.join(line[0:8]),'ID='+gname+'.'+str(gene_isos[IDs])+';Parent='+gname+';',sep='\t')
                                else:
                                    print('\t'.join(line[0:8]),'Parent='+gname+'.'+str(gene_isos[IDs])+';',sep='\t')




def readfiles(F):
    files = []
    with open(F,'r') as f:
        for line in f:
            line = line.strip()
            files.append(line)
    return files

def main():
    files = readfiles(sys.argv[1])
    cutrate = float(sys.argv[2])
    gff_all_dicts,gene_loci,isoform_loci,exon_dicts = readcsv(files)
    gene_clusterd = cluster_obit_site(gene_loci)
    filter_ids = stat_isoforms_overlap(gene_clusterd,exon_dicts,cutrate)
    print_contain_isoform(gene_clusterd, exon_dicts, gff_all_dicts, filter_ids)

if __name__ == '__main__':
    main()

