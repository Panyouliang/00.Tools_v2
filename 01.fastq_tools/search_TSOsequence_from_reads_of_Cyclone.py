import sys,gzip


def readfa(F, sample, TSOF, TSOR, OligoF, OligoR, max_errors):
    whole, chimeric, caps, tails, miss, total_reads = 0, 0, 0, 0, 0, 0
    name,seq = '',[]

    if F[-3:] == '.gz':
        lines = gzip.open(F,'rt')
    else:
        lines = open(F,'r')

    out = open(sample+'.complete.reads.fa','w')
    for line in lines:
        line = line.strip()
        if line[0] == '>':
            total_reads += 1
            if name != '':
                reads = ''.join(seq)
                matcheTSOF = find_matches(reads, TSOF, max_errors)
                matcheTSOR = find_matches(reads, TSOR, max_errors)
                matcheOligoF = find_matches(reads, OligoF, max_errors)
                matcheOligoR = find_matches(reads, OligoR, max_errors)
                w, c, t, m, q, pos1, pos2 = judge(name,reads,matcheTSOF,matcheTSOR,matcheOligoF,matcheOligoR)
                whole += w
                chimeric += q
                caps += c
                tails += t
                miss += m
                if w == 1:
                    if (pos2 > pos1 > 0) and ((pos2-pos1+1) >= 150):
                        out.write(f">{name}\n")
                        out.write(f"{reads[pos1:pos2]}\n")
            else:
                pass
            name = line[1:]
            seq = []
        else:
            seq.append(line)
    lines.close()
    reads = ''.join(seq)
    matcheTSOF = find_matches(reads, TSOF, max_errors)
    matcheTSOR = find_matches(reads, TSOR, max_errors)
    matcheOligoF = find_matches(reads, OligoF, max_errors)
    matcheOligoR = find_matches(reads, OligoR, max_errors)
    w, c, t, m, q, pos1, pos2 = judge(name,reads,matcheTSOF,matcheTSOR,matcheOligoF,matcheOligoR)
    whole += w
    chimeric += q
    caps += c
    tails += t
    miss += m
    if w == 1:
        if (pos2 > pos1 > 0) and ((pos2-pos1+1)>= 150):
            out.write(f">{name}\n")
            out.write(f"{reads[pos1:pos2]}\n")

    out.close()

    return whole, chimeric, caps, tails, miss, total_reads


def find_matches(long_read, query_sequence, max_errors):
    matches = []
    for i in range(len(long_read) - len(query_sequence) + 1):
        errors = 0
        for j in range(len(query_sequence)):
            if long_read[i + j] != query_sequence[j]:
                errors += 1
                if errors > max_errors:
                    break
        if errors <= max_errors:
            matches.append((i, i + len(query_sequence) - 1))
    return matches


def judge(name,reads,TSOF,TSOR,OliF,OliR):
    whole, cap, tail, miss, chimeric = 0, 0, 0, 0, 0
    pos1,pos2 = 0,0

    if (TSOF or TSOR) and (OliF or OliR):
        if (TSOF and OliR) or (TSOR and OliF):
            if TSOF and OliR:
                Lvalue = min([max(pair) for pair in TSOF])
                Rvalue = max([min(pair) for pair in OliR])
                if Rvalue > Lvalue:
                    pos1,pos2 = Lvalue,Rvalue
                    whole += 1
                else:
                    chimeric += 1
            if TSOR and OliF:
                Lvalue = max([min(pair) for pair in TSOR])
                Rvalue = min([max(pair) for pair in OliF])
                if Rvalue < Lvalue:
                    pos1,pos2 = Rvalue,Lvalue
                    whole += 1
                else:
                    chimeric += 1
        else:
            chimeric += 1


    else:
        if not (TSOF or TSOR) and (OliF or OliR):
            tail += 1
        if not (OliF or OliR) and (TSOF or TSOR):
            cap += 1
        if not (TSOF or TSOR) and  not (OliF or OliR):
            miss += 1

    return whole, cap, tail, miss, chimeric, pos1, pos2


def main():

    TSOF_ = "AAGCAGTGGTATCAACGCAGAGTACATGGG"
    TSOR_ = "CCCATGTACTCTGCGTTGATACCACTGCTT"
    OligoF_="AAGCAGTGGTATCAACGCAGAGTACTTTTTT"
    OligoR_="AAAAAAGTACTCTGCGTTGATACCACTGCTT"
    max_errors = 2

    fasta, name = sys.argv[1], sys.argv[2]  ## input the fasta file and name!

    whole, chimeric, caps, tails, miss, total_reads = readfa(fasta, name, TSOF_,TSOR_, OligoF_, OligoR_, max_errors)

    sys.stdout = sys.__stdout__
    sys.stderr.write(f"complete:{whole}\t{whole/total_reads}\nchimeric:{chimeric}\t{chimeric/total_reads}\
            \ncaps:{caps}\t{caps/total_reads}\ntails:{tails}\t{tails/total_reads}\nmissing:{miss}\t{miss/total_reads}\
            \ntotal_reads:{total_reads}\n")


if __name__ == "__main__":
    main()
