#! /usr/bin/env python3
# -*- coding:utf-8 -*-
import sys

f = open(sys.argv[1])

g_id = ''
list1 = []
for line in f:
    line = line.strip().split()
    if line[0][0] != '#':
        if line[0] == g_id:
            pass
        else:
            g_id = line[0]
            d = [x for x in line]
            list1.append(d)
    else:
        print('\t'.join(x for x in line))

#print(list1)
for i in list1:
    print('\t'.join(x for x in i))
#list1 = sorted(list1, key=lambda x: (x[0],x[11]), reverse = False)




