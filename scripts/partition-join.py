#!/usr/bin/python

import screed
from fastp_iter import *
import sys
import khmer

file1 = sys.argv[1]
file2 = sys.argv[2]

class PInfo:
    def __init__(self, pid):
        self.pid = pid
        self.seqs = []
    def add_seq(self, seq):
        self.seqs.append(seq)

class Seq:
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

    def kmers(K):
        n = (len(self.seq) - K) + 1
        for i in range(0, n):
            yield self.seq[i:i+K]

def load_data(filename):
    data = {}
    print 'iterating through {}'.format(filename)
    total = 0
    for name, pid, seq in fastp_iter(filename):
        if total % 100000 == 0:
            print 'processed {} transcripts'.format(total)
        S = Seq(name, seq)
        if pid in data:
            data[pid].add_seq(S)
        else:
            tinfo = PInfo(pid)
            tinfo.add_seq(S)
            data[pid] = tinfo
        total += 1
    print '{} transcripts total, {} partitions'.format(total, len(data))
    return data

data1 = load_data(file1)
data2 = load_data(file2)

keys1 = data1.viewkeys()
keys2 = data2.viewkeys()

print >>sys.stderr, 'joining files on partition id'
join_keys = keys1 & keys2
print >>sys.stderr, 'joined on {} partitions'.format(len(join_keys))
print >>sys.stderr, 'outputting partitioned-joined sequences to files'
with open(file1 + '.joined.fa', 'wb') as outfp:
    for key in data1:
        if key in join_keys:
            for S in data1[key].seqs:
                outfp.write('>{}\t{}\n{}\n'.format(S.name, key, S.seq))

with open(file2 + '.joined.fa', 'wb') as outfp:
    for key in data2:
        if key in join_keys:
            for S in data2[key].seqs:
                outfp.write('>{}\t{}\n{}\n'.format(S.name, key, S.seq))

