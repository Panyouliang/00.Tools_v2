import gzip,sys



def readfile(F, bitch_size=4):

    fq1 = gzip.open(F.split('.')[0]+'.fq1.gz','wt')
    fq2 = gzip.open(F.split('.')[0]+'.fq2.gz','wt')


    with gzip.open(F,'rt') as f:
        lines = []
        rd1,rd2 = 0,0
        for line in f:
            lines.append(line)
            if len(lines) == bitch_size:
                if lines[0].strip()[-3:] == './1':
                    rd1 += 1
                    read1 = lines[0].strip().split('.')[0]
                    fq1.write(read1+'.'+str(rd1)+'./1\n')
                    for later in lines[1:]:
                        fq1.write(later)

                if lines[0].strip()[-3:] == './2':
                    rd2 += 1
                    read2 = lines[0].strip().split('.')[0]
                    fq2.write(read2+'.'+str(rd2)+'./2\n')
                    for later in lines[1:]:
                        fq2.write(later)

                lines = []



readfile(sys.argv[1])
