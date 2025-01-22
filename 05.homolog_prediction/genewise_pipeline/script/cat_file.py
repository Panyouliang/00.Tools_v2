#! /usr/bin/env python3
# -*- coding:utf-8 -*-



import os
import sys
import glob
import argparse

parser = argparse.ArgumentParser(description='cat tools')
parser.add_argument('-p', '--path', type=str, help='you need input path', required=True)
parser.add_argument('-s', '--suffix', type=str, help='input the suffix', required=True)
parser.add_argument('-n', '--outname', type=str, help='input the output name', required=True)
parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose output')
args = parser.parse_args()

fold=args.path
tail=args.suffix
name=args.outname

fout = open(name+tail,'w')
fout.close()


for root, ds, fs in os.walk(fold):
    for f in fs:
        if f.endswith(tail):
            full_name=os.path.join(root, f)
            p = open(full_name)
            line = p.read()
            p.close()
            fout = open(name+tail,'a')
            fout.write(line)
            fout.close()

