#! /usr/bin/env python

import khmer
import screed
import argparse
import sys

def output_labels(infile, out, lh):
    for n, record in enumerate(screed.open(infile)):
        if n % 25000 == 0:
            print >>sys.stderr, '{n}...'.format(n=n)
        seq = record.sequence
        labels = lh.sweep_label_neighborhood(seq, 0)
        out.write('{sid};{name};{nl};{labels}\n'.format(sid=n, name=record.name.split()[0],
            nl=len(labels), labels=','.join([str(l) for l in labels])))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', dest='ksize', type=int, default=25)
    parser.add_argument('-x', dest='hashsize', type=float, default=1e9)
    parser.add_argument('-N', dest='nhashes', type=int, default=4)
    parser.add_argument('-o', dest='outfile', 
        type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('infile')
    args = parser.parse_args()
    
    lh = khmer.LabelHash(args.ksize, args.hashsize, args.nhashes)
    lh.consume_fasta_and_tag_with_labels(args.infile)
    
    output_labels(args.infile, args.outfile, lh)

if __name__ == '__main__':
    main()
