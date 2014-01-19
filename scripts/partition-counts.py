#!/usr/bin/python

import screed
from fastp_iter import *
import sys

infile = sys.argv[1]
class PInfo:
    def __init__(self, pid):
        self.pid = pid
        self.count = 0

data = {}
total = 0
print 'iterating through {}'.format(infile)
for name, pid, seq in fastp_iter(infile):
    if total % 100000 == 0:
        print 'processed {} transcripts'.format(total)
    if pid in data:
        data[pid].count += 1
    else:
        tinfo = PInfo(pid)
        tinfo.count += 1
        data[pid] = tinfo
    total += 1
print '{} transcripts total, {} partitions'.format(total, len(data))
print 'sorting for output...'
datalist = data.values()
datalist.sort(key=lambda pinfo: pinfo.pid)
outfile = infile + '.pcounts'
print 'writing data to {}'.format(outfile)
with open(outfile, 'wb') as outfp:
    for p in datalist:
        outfp.write('{},{}\n'.format(p.pid, p.count))
