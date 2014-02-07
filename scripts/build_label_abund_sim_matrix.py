import pandas as pd
import numpy as np
import itertools
import argparse
import sys

'''
mkunion taken from:
http://stackoverflow.com/questions/3438140/how-to-create-the-union-of-many-sets-using-a-generator-expression
'''
def mkunion(*args):
    return frozenset(itertools.chain.from_iterable(args))

def build_label_intersect_mat(df, binsize, func, out):

    arr_s = int(np.max(df.nlabels))/binsize + 1
    M = np.zeros( (arr_s, arr_s) )

    # assign each sequence a bin based on label abundance
    bins = df.nlabels/binsize
    bins = bins.astype(int)
    df.insert(3, 'bin', bins)
    
    # for each abundance bin, create a set of labels
    bin_labels = []
    for b in xrange(arr_s):
        bin_labels.append(mkunion(*df[df['bin'] == b].labels.tolist()))

    # get size of set op between abundance bins
    for x in xrange(arr_s):
        for y in xrange(arr_s):
            print >>sys.stderr, '({x},{y})'.format(x=x,y=y)
            common = float(len(func(bin_labels[x], bin_labels[y])))/len(bin_labels[x])
            M[x,y] = common

    np.savetxt(out, M, delimiter=',')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', dest='binsize', default=5, type=int)
    parser.add_argument('-o', dest='outfile', 
        type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('--func', dest='func', default='intersect',
        choices=['intersect', 'union', 'diff', 'symdiff'])
    parser.add_argument('infile')
    args = parser.parse_args()
    
    df = pd.read_csv(args.infile, sep=';', header=None, 
        names=['sid','name', 'nlabels', 'labels'],
        dtype={'sid':np.int, 'name':str, 'nlabels':np.float64, 
                'labels':set}, converters={'labels':lambda x: set(x.split(','))})

    if args.func == 'difference':
        func = frozenset.difference
    elif args.func == 'union':
        func=frozenset.union
    elif args.func == 'intersect':
        func=frozenset.intersection
    else:
        func = frozenset.symmetric_difference

    build_label_intersect_mat(df, args.binsize, func, args.outfile)

if __name__ == '__main__':
    main()
