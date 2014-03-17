#!/usr/bin/python

import screed
from fastp_iter import *
import sys

infile = sys.argv[1]
class PInfo:
    def __init__(self, pid):
        self.pid = pid
        self.count = 0
        self.g = 0
        self.c = 0
        self.a = 0
        self.t = 0
    def add_nuc(self, nuc):
        nuc.upper()
        if nuc == 'C':
            self.c += 1
        elif nuc == 'G':
            self.g += 1
        elif nuc == 'A':
            self.a += 1
        else:
            self.t += 1
    def calc_gc(self):
        return float(self.c + self.g) / float(self.g + self.c + self.a + self.t)
    
data = {}
total = 0
print 'iterating through {}'.format(infile)
for name, pid, seq in fastp_iter(infile):
    if total % 100000 == 0:
        print 'processed {} transcripts'.format(total)
    if pid in data:
        for c in seq:
            data[pid].add_nuc(c)
    else:
        tinfo = PInfo(pid)
        for c in seq:
            tinfo.add_nuc(c)
        data[pid] = tinfo
    total += 1
print '{} transcripts total, {} partitions'.format(total, len(data))
print 'sorting for output...'
datalist = data.values()
datalist.sort(key=lambda pinfo: pinfo.pid)
outfile = infile + '.pgc'
print 'calculating writing data to {}'.format(outfile)
with open(outfile, 'wb') as outfp:
    for p in datalist:
        outfp.write('{},{}\n'.format(p.pid, p.calc_gc()))
