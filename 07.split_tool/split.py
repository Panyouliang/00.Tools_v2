import sys,os


if ('--help' in sys.argv) or ('-h' in sys.argv):
    sys.stderr.write("usage: python3 script.py <query split file> <split num> <spacies name>\n")
try:
    dirs = sys.argv[1]
except:
    exit("Please enter need to split file")
try:
    split_num = sys.argv[2]
except:
    exit("Please enter split num")
try:
    species = sys.argv[3]
except:
    exit("Please enter new file name")


f2 = open(dirs,'r')


filed = {}
seq,ids = [],''
for line in f2:
    line = line.strip()
    if '>' in line:
        filed[ids] = ''.join(seq)
        ids = line
        seq = []
    else:
        seq.append(line)

filed[ids] = ''.join(seq)



def wrap_sequence(value, line_length=60):
    wrapped_sequence = ""
    for i in range(0, len(value), line_length):
        wrapped_sequence += value[i:i + line_length] + '\n'
    return wrapped_sequence



lines = len(filed.keys())
n = int(split_num)

os.system('mkdir split')

if lines%n == 0:
    qty = lines//n
else:
    qty = lines//n + 1


lines,file_num = 0,0
for k,v in filed.items():
    if k == '':
        pass
    else:
        lines += 1
        if lines % qty == 1:
            file_num += 1
            f = open('split/'+species+str(file_num)+'.fasta','w')
            f.write(k+'\n')
            f.write(wrap_sequence(v))
        else:
            f.write(k+'\n')
            f.write(wrap_sequence(v))
