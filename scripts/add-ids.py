import screed
import sys
from fastp_iter import *

input_file = sys.argv[1]
tag = sys.argv[2]

with open(input_file + '.id.fa', 'wb') as outfp:
    for n, record in enumerate(fastp_iter(input_file)):
        name, pid, sequence = record
        if n % 100000 == 0:
            print >>sys.stderr, '...processed {} sequences'.format(n)
        outfp.write('>{}_{}\t{}\n{}\n'.format(tag, name, pid, sequence))
