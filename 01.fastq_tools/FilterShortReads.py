import sys
import re
import pysam




if ("--help" in sys.argv) or ("-h" in sys.argv):
    sys.stderr.write("usage: python3 FilterShortReads.py  <bam> <cut length below>\n")

try:
    f1=sys.argv[1]
    p2=int(sys.argv[2])
except:
    exit("usage: python3 FilterShortReads.py <bam>  <cut length below>")

out_file='.'.join(f1.split('.')[:-1])+'.lengthcut.bam'

samfile = pysam.AlignmentFile(f1, "rb")
out_file = pysam.AlignmentFile(out_file, "wb", template=samfile)

flag=0
input_reads = 0
for read in samfile:
    input_reads += 1
    cigar = read.cigar
    #print(cigar)
    RL=0
    for i in cigar:
        if i[0]==0:#meaning match in pysam format:
            RL+=i[1]
    if RL>=p2:
        out_file.write(read)
        flag+=1

print('input reads = ' + str(input_reads))
print('retained >= '+str(p2)+': '+str(flag))



