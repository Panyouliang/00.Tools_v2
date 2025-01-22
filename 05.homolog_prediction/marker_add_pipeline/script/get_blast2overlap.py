import sys




def id2id(F):
    name2name = {}
    with open(F,'r') as f:
        for line in f:
            line = line.strip().split()
            name2name[line[0]] = line[1]
    return name2name


def overlap(F,name2name,alghit):
    existhit = []
    with open(F,'r') as f:
        for line in f:
            if line[0] != '#':
                line = line.strip().split()
                gwid,refid = line[0],line[1]
                overlapr1,overlapr2 = float(line[10]),float(line[11])
                gname = gwid.split('-')[0]

                if gwid in name2name.keys():
                    if overlapr1 > 0.5 or overlapr2 >0.5:
                        print(gwid,name2name[gwid],'\t'.join(x for x in alghit[gwid]),refid,overlapr1,overlapr2,sep='\t')
                        existhit.append(gwid)
    nonhit = []
    for k,v in name2name.items():
        if k in existhit:
            pass
        else:
            nonhit.append(k)

    return nonhit

def getalg(F):
    hitd = {}
    with open(F,'r') as f:
        for line in f:
            if line[0] != '#':
                line = line.strip().split()
                hitd[line[0]] = [line[5],line[8],line[9],line[10]]
    return hitd


def main():
    name2name = id2id(sys.argv[1])
    alghitd = getalg(sys.argv[2])
    nonhit = overlap(sys.argv[3],name2name,alghitd)


    nonid = open('existing-but-non-overlap.genes','w')
    for i in nonhit:
        nonid.write(name2name[i]+'\n')
        print(i,name2name[i],'\t'.join(x for x in alghitd[i]),sep='\t')


if __name__ == '__main__':
    print('#GwisID\tGeneName\tChromosome\tAlignRate\tScore\tIdentity\tSubjectID\tPer1\tPer2')
    main()


