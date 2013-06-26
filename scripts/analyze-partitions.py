import argparse
import os
import hambone as hb
'''
CS Welcher
MSU GED

Read in partitioned fasta files and corresponding blast files,
then do ~things~ with them.

Usage:
analyze-partitions.py -i <partitioned fasta 1> <partitioned fasta 2> .. -j <fasta 1 vs X> <fasta 2 vs X> ..
'''

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--out_name', dest='out_name')
    parser.add_argument('-i', '--fastp', dest='input_fastp', nargs='+')
    parser.add_argument('-j', '--blast', dest='input_blast_csv', nargs='+')
    args = parser.parse_args()
    
    
    output_basename = args.out_name
    files = zip(args.input_fastp, args.input_blast_csv)
    sets = {}
    
    # in the future, this will do all pairs for any number of files
    # for now, just assuming two
    for fastp, blast in files:
        ps = hb.PartitionSet.associateHits(fastp, blast)
        fn = os.path.basename(fastp)
        sets[fn] = ps
    
    fastp_A = sets.keys()[0]
    fastp_B = sets.keys()[1]
    
    print "Performing intersect..."
    A_intersect_B = sets[fastp_A] & sets[fastp_B]
    print "Performing A-B..."
    A_only = sets[fastp_A] - sets[fastp_B]
    print "Performing B-A..."
    B_only = sets[fastp_B] - sets[fastp_A]
    print "Performing union..."
    A_union_B = sets[fastp_A] + sets[fastp_B]
    
    print '*' * 40
    print "Intersect:\n", A_intersect_B
    print "{} - {}\n".format(fastp_A, fastp_B), A_only
    print "{} - {}\n".format(fastp_B, fastp_A), B_only
    print "Union:\n", A_union_B
    
    for fn in sets:
        print sets[fn]
        with open(fn + '.hist', 'wb') as histfp:
            hist = sets[fn].histogram()
            for row in hist:
                histfp.write(row)

if __name__ == '__main__':
    main()
