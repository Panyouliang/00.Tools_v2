import sys



def read(F):
    gened,feature = {},{}
    with open(F,'r') as f:
        for line in f:
            line = line.strip().split()
            if line[2] == 'gene':
                pass
            elif line[2] == 'mRNA':
                name = line[8].split('=')[1].split(';')[0]
                feature[name] = []
                gened[name] = [line]
            else:
                pname = line[8][line[8].find('Parent'):].split('=')[1].split(';')[0]
                feature[pname].append([int(line[3]),int(line[4])])
                gened[name].append(line)

    feature_adj = {}
    for k,v in feature.items():
        vsort = sorted(v, key=lambda x:int(x[1]))
        feature_adj[k] = merge_adjacent_coordinates(vsort)

    return feature_adj, gened




def merge_adjacent_coordinates(coordinates):
    merged_coordinates = []
    if not coordinates:
        return merged_coordinates
    current_start, current_end = coordinates[0]
    for coord in coordinates[1:]:
        if coord[0] <= current_end + 1:
            current_end = max(current_end, coord[1])
        else:
            merged_coordinates.append([current_start, current_end])
            current_start, current_end = coord
    merged_coordinates.append([current_start, current_end])
    return merged_coordinates



if __name__ == '__main__':
    cdsd,gened = read(sys.argv[1])

    for k,v in cdsd.items():
        if len(v) > 1:
            for line in gened[k]:
                print('\t'.join(line))









