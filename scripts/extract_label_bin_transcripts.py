import pandas as pd
import numpy as np
import itertools
import argparse
import sys
import screed

def extract_label_bin_transcripts(df, binsize, ebins, in_fp, out_a, out_b):

    arr_s = int(np.max(df.nlabels))/binsize + 1
    M = np.zeros( (arr_s, arr_s) )

    # assign each sequence a bin based on label abundance
    bins = df.nlabels/binsize
    bins = bins.astype(int)
    df.insert(3, 'bin', bins)
    
    for n, record in enumerate(in_fp):
        if n % 25000 == 0 and n:
            print >>sys.stderr, '{n}...'.format(n=n)
        try:
            b = df.iat[n,3]
            df_n = df.iat[n,1]
        except:
            print >>sys.stderr, 'ERROR in df lookup'
            sys.exit()

        #print df[df.sid == n], b, df_n
        if record.name.split('\t')[0] != df_n:
            raise LookupError('ERROR: mismatch between dataframe and fasta\n \
                               {n1} =/= {n2}'.format(n1=record.name, n2=df_n))
        if b in ebins:
            out_b.write('>{name} {b}\n{seq}\n'.format(name=record.name,
                seq=record.sequence, b=b))
        else:
            out_a.write('>{name} {b}\n{seq}\n'.format(name=record.name,
                seq=record.sequence, b=b))
    

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', dest='binsize', default=5, type=int)
    parser.add_argument('-e', dest='extract', type=int, nargs='+')
    parser.add_argument('-o', dest='outpref')
    parser.add_argument('in_labels')
    parser.add_argument('in_fasta')
    args = parser.parse_args()
    
    try:
        e_ln = '-'.join(map(str,args.extract)) 
        out_a = open('{o}.remaining.{e}.fa'.format(
            o=args.outpref, e=e_ln), 'wb')
        out_b = open('{o}.extracted.{e}.fa'.format(
            o=args.outpref, e=e_ln), 'wb')
    except:
        sys.exit()

    ebins = map(lambda x: x/args.binsize, args.extract)

    df = pd.read_csv(args.in_labels, sep=';', header=None, 
        names=['sid','name', 'nlabels', 'labels'],
        dtype={'sid':np.int, 'name':str, 'nlabels':np.float64, 
                'labels':set}, converters={'labels':lambda x: set(x.split(','))})
    try:
        in_fp = screed.open(args.in_fasta)
    except IOError as e:
        print >>sys.stderr, e
        sys.exit()
    else:
        try:
            extract_label_bin_transcripts(df, args.binsize, 
                                          ebins, in_fp, out_a, out_b)
        except LookupError as e:
            print >>sys.stderr, e
            sys.exit()

if __name__ == '__main__':
    main()
