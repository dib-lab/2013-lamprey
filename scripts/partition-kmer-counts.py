#!/usr/bin/python

import screed
from fastp_iter import *
import sys
import khmer

infile = sys.argv[1]
class PInfo:
    def __init__(self, pid):
        self.pid = pid
        self.seqs = []
    def add_seq(self, seq):
        self.seqs.append(seq)

def kmers(seq, K):
    n = (len(seq) - K) + 1
    for i in range(0, n):
        yield seq[i:i+K]

ht = khmer.new_counting_hash(25, 1e9, 4)

data = {}
total = 0
print 'iterating through {}'.format(infile)
for name, pid, seq in fastp_iter(infile):
    if total % 100000 == 0:
        print 'processed {} transcripts'.format(total)
    if pid in data:
        data[pid].add_seq(seq)
    else:
        tinfo = PInfo(pid)
        tinfo.add_seq(seq)
        data[pid] = tinfo
    total += 1
print '{} transcripts total, {} partitions'.format(total, len(data))
print 'sorting for output...'
datalist = data.values()
datalist.sort(key=lambda pinfo: pinfo.pid)
outfile = infile + '.kcounts'
print 'writing data to {}'.format(outfile)
all_kmers = 0
with open(outfile, 'wb') as outfp:
    for p in datalist:
        t_ht = khmer.new_counting_hash(25, 4e5, 4)
        t_count = 0
        for seq in p.seqs:
            for kmer in kmers(seq, 25):
                if not t_ht.get(kmer):
                    t_count += 1
                t_ht.count(kmer)
                if not ht.get(kmer):
                    all_kmers += 1
                ht.count(kmer)
        outfp.write('{},{},{}\n'.format(p.pid, t_count, all_kmers))
