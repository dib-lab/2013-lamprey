#!/usr/bin/python
import screed

def fastp_iter(filename):
    for record in screed.open(filename, parse_description=False):
        name = record.name
        try:
            name, partition_id = name.rsplit('\t', 1)
        except ValueError:
            print "Derp! Is this file partitioned?"
            sys.exit(1)
        # convert name to blast format if necessary
        nname = name.split('|', 2)
        if len(nname) >= 2:
            name = nname[2]
        name = name.split(' ')[0]
        yield name, int(partition_id), record.sequence
