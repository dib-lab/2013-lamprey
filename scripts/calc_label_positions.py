#! /usr/bin/env python

import khmer
import screed
import argparse
import sys

def calc_label_positions(infile, out, lh):

    for n, record in enumerate(screed.open(infile)):
        if n % 25000 == 0:
            print >>sys.stderr, '{n}...'.format(n=n)
        seq = record.sequence
        L = len(seq)
        tags = lh.sweep_tag_neighborhood(seq, 0)
        X = []
        for tag in tags:
            labels = lh.get_tag_labels(tag)
            kmer = khmer.reverse_hash(tag, lh.ksize())
            pos = seq.find(kmer)
            X.append((pos,len(labels)))
        q1, q1l = 0,0
        q2, q2l = 0,0
        q3, q3l = 0,0
        q4, q4l = 0,0
        for x,l in X:
            p = float(x)/L
            if p < .25:
                q1+=1
                q1l+=l
            elif p < .5:
                q2+=1
                q2l+=l
            elif p < .75:
                q3+=1
                q3l+=l
            else:
                q4+=1
                q4l+=l
        try:
            q1l = float(q1l)/q1
        except ZeroDivisionError:
            q1l = 0.0
        try:
            q2l = float(q2l)/q2
        except ZeroDivisionError:
            q2l = 0.0
        try:
            q3l = float(q3l)/q3
        except ZeroDivisionError:
            q3l = 0.0
        try:
            q4l = float(q4l)/q4
        except ZeroDivisionError:
            q4l = 0.0
        
        out.write('{n},{q1},{q1l},{q2},{q2l},{q3},{q3l},{q4},{q4l}\n'.format(
                    n=n, q1=q1, q1l=q1l, q2=q2, q2l=q2l, 
                    q3=q3, q3l=q3l, q4=q4, q4l=q4l))
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
    
    calc_label_positions(args.infile, args.outfile, lh)

if __name__ == '__main__':
    main()
