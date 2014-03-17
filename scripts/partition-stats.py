#!/usr/bin/python

import screed
from fastp_iter import *
import sys
import khmer
import scipy.stats as sps


infile = sys.argv[1]

class PInfo:
    def __init__(self, pid):
        self.pid = pid
        self.seqs = []
    def add_seq(self, seq):
        self.seqs.append(seq)

    def stats(self):
        lens = [float(len(S.seq)) for S in self.seqs]
        size, minmax, mean, var, _, _ = sps.describe(lens)
        return [self.pid, size, minmax[0], minmax[1], mean, var]

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

data = load_data(infile)

with open(infile + '.stats', 'wb') as outfp:
    for key in data:
        stats = data[key].stats()
        outfp.write('{},{},{},{},{},{}\n'.format(stats[0], stats[1], stats[2], stats[3], stats[4], stats[5]))
