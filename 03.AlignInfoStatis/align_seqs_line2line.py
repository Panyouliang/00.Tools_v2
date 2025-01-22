import sys



def readf(F):
    seqd = {}
    name,seq = '',[]
    lis = []
    with open(F,'r') as f:
        for line in f:
            line = line.strip()
            if line[0] == '>':
                if name != '':
                    seqd[name] = seq
                name = line[1:]
                lis.append(name)
                seq = []
            else:
                seq.append(line.upper())
    seqd[name] = seq

    return seqd,lis


def print_lines(input_str, line_length=60):
    seq = []
    for i in range(0, len(input_str), line_length):
        seq.append(input_str[i:i+line_length])
    return '\n'.join(x for x in seq)

def links(length):
        return '|' * length



if __name__ == '__main__':
    seqd,lis = readf(sys.argv[1])
    for i in range(len(seqd[lis[0]])):
        print(f"{seqd[lis[0]][i]}")
        print(links(len(seqd[lis[0]][i])))
        print(f"{seqd[lis[1]][i]}")
        marker = []
        for base1,base2 in zip(seqd[lis[0]][i], seqd[lis[1]][i]):
            if base1 == base2:
                marker.append('.')
            else:
                marker.append('*')

        print(''.join(marker))
        print(' ')






