import sys,gzip
import numpy  as np
import matplotlib.pyplot as plt


#input fasta specieName sampleName


def readfa(F,species,sample):
    total = 0
    reads, length, rlist = 0, 0, []
    G_count, C_count, N_count = 0, 0, 0

    read_len = {}

    line_num = 0

    k10,k20,k30,k40,k50,k80 = 0,0,0,0,0,0

    if F[-3:] == '.gz':
        f = gzip.open(F,'rt')
    else:
        f = open(F,'r')

    for line in f:
        line = line.strip()
        if line[0] == '>':
            name = line[1:]
            readL = 0
        else:
            if len(line) >= 0:
                total += len(line)
                readL += len(line)

                read_len[name] = readL

                reads += 1
                length += len(line)
                rlist.append(len(line))

                G_count += line.count('G')
                C_count += line.count('C')
                N_count += line.count('N')
                if len(line) > 10000:
                    k10 += len(line)
                if len(line) > 20000:
                    k20 += len(line)
                if len(line) > 30000:
                    k30 += len(line)
                if len(line) > 40000:
                    k40 += len(line)
                if len(line) > 50000:
                    k50 += len(line)
                if len(line) > 80000:
                    k80 += len(line)

#=====#statistic value #=====#

    read_len_sort = sorted(read_len.items(),key=lambda x:x[1],reverse=True)

    outf = open(sample+'.readLen.lis','w')
    for cut in read_len_sort:
        outf.write(str(cut[0])+'\t'+str(cut[1])+'\n')


    GC_rate = (G_count + C_count)/(length - N_count)

    max_read = max(rlist)
    min_read = min(rlist)

    k10rate,k20rate,k30rate,k40rate,k50rate,k80rate = format(k10/length,'.2f'), format(k20/length, '.2f'), format(k30/length,'.2f'), format(k40/length), format(k50/length,'.2f'), format(k80/length,'.2f')

    sortlis = rlist.copy()

    sortlis.sort(reverse=True)

    N50_L = total/2
    Nnum = 0
    for rL in sortlis:
        Nnum += rL
        if Nnum >= N50_L:
            readN50 = rL
            break


    nprlis = np.array(rlist)
    rmean = np.mean(nprlis)
    rmedian = np.median(nprlis)

    print(species,sample,reads,max_read,min_read,k10rate,k20rate,k30rate,k40rate,k50rate,k80rate,readN50,format(rmean,'.2f'),format(rmedian,'.2f'),length,format(GC_rate,'.2f'),sep='\t')

    #plt.hist(nprlis, bins=5000, edgecolor='black')
    #plt.xlim((-1,10000))
    #plt.xlabel('read length')
    #plt.ylabel('Frequency')
    #plt.title('Histogram of Cyclone')

    #plt.savefig(sample+'.pdf', bbox_inches='tight')




readfa(sys.argv[1],sys.argv[2],sys.argv[3])
