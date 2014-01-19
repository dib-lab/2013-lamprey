import screed
import sys

in_file = sys.argv[1]
out_file = in_file + '.fa'

with open(out_file, 'wb') as outfp:
    for record in screed.open(in_file):
        name_line = record.name.split()
        name = name_line[0]
        seq = record.sequence
        outfp.write('>{}\n{}\n'.format(name, seq))
