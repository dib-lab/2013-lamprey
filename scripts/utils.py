#!/usr/bin/python
import screed
import csv

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

def kmers(seq, K):
    n = (len(seq) - K) + 1
    for i in xrange(n):
        yield seq[i:i+K]

def make_table(grid):
    max_cols = [max(out) for out in map(list, zip(*[[len(item) for item in row] for row in grid]))]
    rst = table_div(max_cols, 1)

    for i, row in enumerate(grid):
        header_flag = False
        if i == 0 or i == len(grid)-1: header_flag = True
        rst += normalize_row(row,max_cols)
        rst += table_div(max_cols, header_flag )
    return rst

def table_div(max_cols, header_flag=1):
    out = ""
    if header_flag == 1:
        style = "="
    else:
        style = "-"

    for max_col in max_cols:
        out += max_col * style + " "

    out += "\n"
    return out


def normalize_row(row, max_cols):
    r = ""
    for i, max_col in enumerate(max_cols):
        r += row[i] + (max_col  - len(row[i]) + 1) * " "

    return r + "\n"

def iter_distinct_hits(blast_file, header="qseqid sacc length qstart qend sstart send bitscore evalue"):
    header = header.split()
    
    with open(blast_file, 'rb') as csv_fp:
        reader = csv.DictReader(csv_fp, skipinitialspace=True, fieldnames=header)
        
        best_hit = reader.next()
        best_query = best_hit['qseqid']
        for n, this_hit in enumerate(reader):
            this_query = this_hit['qseqid']
            if this_query == best_query:
                # Maximize score and minimize expect
                # (score in the alignments section is the sum of all HSP bitscores, so already converted for information content)
                # Little bit of convoluted logic for when expect is 0; favor the 0 expect, or choose by score if they're both 0
                
                this_hit_e = float(this_hit['evalue'])
                this_hit_score = float(this_hit['bitscore'])
                
                best_hit_e = float(best_hit['evalue'])
                best_hit_score = float(best_hit['bitscore'])
                
                if this_hit_e == 0 or best_hit_e == 0:
                    if this_hit_e < best_hit_e:
                        best_hit = this_hit
                    elif this_hit_e == best_hit_e:
                        if this_hit_score > best_hit_score:
                            best_hit = this_hit
                            best_query = this_query
                elif this_hit_score / this_hit_e > best_hit_score / best_hit_e:
                    best_hit = this_hit
                    best_query = this_query
            elif this_query != best_query:
            
                yield best_hit

                best_hit = this_hit
                best_query = this_query
            if n % 10000 == 0:
                print >>sys.stderr, "hit: {}, {}".format(n, this_hit)


