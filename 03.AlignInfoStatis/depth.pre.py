#! /usr/bin/env python3
# -*- coding:utf-8 -*-

import sys

F = open(sys.argv[1])

fileiter = (x.strip().split() for x in F)

dicts = {}
test = ''
cutoff = int(sys.argv[2])
for line in fileiter:
    length = int(line[3])
    if line[0] == test:
        if int(line[1]) >= cutoff:
            rate += float(line[4])
            dep += int(line[1])*int(line[2])
            dicts[line[0]+'\t'+'depth>='+str(cutoff)]=[rate,dep]
            test = line[0]
        else:
            if int(line[1]) != 0:
                rat += float(line[4])
                de += int(line[1])*int(line[2])
                dicts[line[0]+'\t'+'depth<'+str(cutoff)]=[rat,de]
    else:
        test = line[0]
        rate=rat=dep=de=0
        if int(line[1]) >= cutoff:
            rate += float(line[4])
            dep += int(line[1])*int(line[2])
            dicts[line[0]+'\t'+'depth>='+str(cutoff)]=[rate,dep]
            test = line[0]
        else:
            if int(line[1]) != 0:
                rat += float(line[4])
                de += int(line[1])*int(line[2])
                dicts[line[0]+'\t'+'depth<'+str(cutoff)]=[rat,de]


for k,v in dicts.items():
    print(k,length,format(v[0],'.4f'),format(float(v[1]/length),'.4f'),sep='\t')

