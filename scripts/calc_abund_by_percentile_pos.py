#! /usr/bin/env python

import khmer
import sys
import utils
from fastp_iter import *
import numpy as np
from os.path import basename as base
import argparse

def perc_bin_abunds(seq, ht, func):
    abunds = get_abunds(seq, ht)
    L = len(abunds) 
    tmp = [list() for i in xrange(100)]
    for i, x in enumerate(abunds):
        bucket = int(np.floor(float(i) / L * 100))
        tmp[bucket].append(x)
    return np.array(map(func, tmp))

def get_abunds(seq, ht):
    abunds = np.array([ht.get(kmer) for kmer in utils.kmers(seq, ht.ksize())])
    return abunds

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--function', dest='function', choices=['min', 'max', 'mean', 'var'], default='min')
    parser.add_argument('ht')
    parser.add_argument('infiles', nargs='+')
    args = parser.parse_args()

    if args.function == 'min':
        func = np.min
    elif args.function == 'max':
        func = np.max
    elif args.function == 'mean':
        func = np.mean
    else:
        func = np.var

    print >>sys.stderr, 'loading counting hash from {f}...'.format(f=args.ht)
    ht = khmer.load_counting_hash(args.ht)

    for infile in args.infiles:
        print >>sys.stderr, 'processing {i}...'.format(i=infile)
        n = 0
        with open('{i}.{f}.csv'.format(i=base(infile),f=args.function), 'wb') as outfp:
            for name, pid, sequence in fastp_iter(infile):
                n += 1
                if n % 10000 == 0:
                    print >>sys.stderr, 'processed {n} sequences...'.format(n=n)
                
                vals = perc_bin_abunds(sequence, ht, func)
                outfp.write('{pid},'.format(pid=pid) + ','.join([str(x) for x in vals]) + '\n')

if __name__ == '__main__':
    main()
