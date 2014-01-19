import khmer
import screed
import sys

input_file = sys.argv[1]

five_mers = {}

def kmers(seq, K):
    n = (len(seq) - K) + 1
    for i in range(0, n):
        yield seq[i:i+K]

for n, record in enumerate(screed.open(input_file)):
    if n % 100000 == 0:
        print >>sys.stderr, '...processed {num} sequences'.format(num=n)
    for kmer in kmers(record.sequence, 5):
        five_mers[kmer] = five_mers.get(kmer, 0) + 1

mer_list = zip(five_mers.keys(), five_mers.values())
mer_list.sort(key=lambda i: i[1])

with open(input_file + '.mers', 'wb') as outfp:
    for row in mer_list:
        outfp.write('{mer}, {c}\n'.format(mer=row[0], c=row[1]))
