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

class Sequence:

    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

def partition_iter(filename):
    cur_pid = -1
    part = []
    for name, pid, seq in fastp_iter(filename):
        if pid != cur_pid:
            ret_pid = cur_pid
            ret_part = part
            cur_pid = pid
            part = []
            yield ret_pid, ret_part
        else:
            part.append(Sequence(name, sequence))
    yield cur_pid, part

