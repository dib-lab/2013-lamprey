'''
It's not pea soup without a hambone!

Collection of classes and functions used in peasoup, including:
-blast parsers
-classes for homologies
-classes for scaffolding
...etc

From now on all utility classes and functions will live here,
or if it gets too big, other similary named modules

CS Welcher
GED Lab
Michigan State University
'''

import screed
import csv
import khmer
import sys

class BlastHit:

    def __init__(self, query, subject, length, nident, qstart, qend, sstart, send, evalue, bitscore):
        self.query = query
        self.subject = subject
        self.length = int(length)
        self.nident = int(nident)
        self.qstart = int(qstart)
        self.qend = int(qend)
        self.sstart = int(sstart)
        self.send = int(send)
        self.bitscore = float(bitscore)
        self.evalue = float(evalue)
    
    def __len__(self):
        return self.length
    
    # maximizes bitscore and minimizes evalue
    # favors a 0 evalue, or the bitscore if both are 0
    # otherwise divides bitscore by evalue
    '''
    def __gt__(self, other):
        if self.evalue == 0.0 or other.evalue == 0.0:
            if self.evalue < other.evalue:
                return True
            elif self.evalue == other.evalue:
                return self.bitscore > other.bitscore
        else:
            return self.bitscore/self.evalue > other.bitscore/other.evalue
    '''
    def __cmp__(self, other):
        if self.evalue == other.evalue:
            if self.bitscore == other.bitscore:
                return 0
            elif self.bitscore > other.bitscore:
                return 1
            else:
                return -1
        elif self.evalue == 0.0:
            return 1
        elif other.evalue == 0.0:
            return -1
        else:
            thisScore = self.bitscore/self.evalue
            otherScore = other.bitscore/other.evalue
            if thisScore > otherScore:
                return 1
            elif thisScore == otherScore:
                return 0
            else:
                return -1
    
    def __str__(self):
        return ','.join([self.query, self.subject, str(self.length), str(self.evalue), str(self.qstart), str(self.qend)])
    
    def __repr__(self):
        return self.query + '--' + self.subject

class ScaffoldTNode:

    def __init__(self, name):
        self.name = name
        self.children = []
    
    def addChild(self, gnode):
        self.children.append(gnode)

class ScaffoldGNode:

    def __init__(self, name):
        self.name = name
        self.hits = []
    
    
    

'''
Collection of BlastHits. Reads a file with the given column names,
filters by length or bitscore, and aggregates by query name
'''
class BlastDB:
    
    def __init__(self, blast_file, cols="qseqid sseqid length nident qstart qend sstart send bitscore evalue qframe sframe"):
        self.hits = []
        
        cols = cols.split()
        with open(blast_file, 'rb') as bfp:
            print "Parsing hits..."
            reader = csv.DictReader(bfp, skipinitialspace=True, fieldnames=cols)
            for n, row in enumerate(reader):
                hit = BlastHit(row['qseqid'], row['sseqid'], row['length'], row['nident'], \
                               row['qstart'], row['qend'], row['sstart'], row['send'], row['evalue'], row['bitscore'])
                if n % 10000 == 0:
                    print "Parsed", n, "rows..."
                    print hit
                self.hits.append(hit)
        print "...done parsing {} hits.".format(len(self.hits))
    
    def filterByLength(self, min_length):
        print "filtering out hits less than {} bases...".format(min_length)
        newHits = []
        for n, hit in enumerate(self.hits):
            if n % 10000 == 0:
                print "...looked through", n, "hits"
                print hit
            if len(hit) >= min_length:
                newHits.append(hit)
        self.hits = newHits
        print "...done filtering, {} hits remaining".format(len(self.hits))
    
    def filterByBitscore(self, min_bitscore):
        newHits = []
        for hit in self.hits:
            if hit.bitscore > min_bitscore:
                newHits.append(hit)
        self.hits = newHits

    '''
    Creates a dict of lists of the hits, where each sublist
    is the hits from a single query
    
    Optionally remove groups under min_size hits
    '''
    def partitionByQuery(self, min_size=1):
        self.hitsByQuery = {}
        print "sorting hits by query..."
        self.hits.sort(key=lambda hit: hit.query)
        print "sorted {} hits".format(len(self.hits))
        
        cur_hit = self.hits[0]
        tmplst = [cur_hit]
        
        print "partitioning hits by query..."
        for hit in self.hits[1:]:
            if hit.query == cur_hit.query:
                tmplst.append(hit)
                cur_hit = hit
            else:
                if len(tmplst) >= min_size:
                    self.hitsByQuery[cur_hit.query] = tmplst
                tmplst = [hit]
                cur_hit = hit
        self.hitsByQuery[cur_hit.query] = tmplst
        print "...done!"
        print len(self.hitsByQuery), "queries over", len(self.hits), "hits"
        return self.hitsByQuery

def parse_blast_distinctly(blast_file, header="query subject length qstart qend sstart send bitscore expect"):

    header = header.split()
    distinct_hits = {}
    
    with open(blast_file, 'rb') as csvfp:
        reader = csv.DictReader(csvfp, skipinitialspace=True, fieldnames=header)
        
        best_hit = reader.next()
        best_query = best_hit['query']
        for n, this_hit in enumerate(reader):
            this_query = this_hit['query']
            if this_query == best_query:
                # Maximize score and minimize expect
                # (score in the alignments section is the sum of all HSP bitscores, so already converted for information content)
                # Little bit of convoluted logic for when expect is 0; favor the 0 expect, or choose by score if they're both 0
                
                this_hit_e = float(this_hit['expect'])
                this_hit_score = float(this_hit['bitscore'])
                
                best_hit_e = float(best_hit['expect'])
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
            
                distinct_hits[best_query] = [best_hit[k] for k in header]

                best_hit = this_hit
                best_query = this_query
            if n % 10000 == 0:
                print "hit: {}, {}".format(n, this_hit)
        distinct_hits[best_query] = [best_hit[k] for k in header]
    return distinct_hits

class Seq(object):

    def __init__(self, name, seq, **kwds):
        self.name = name
        self.seq = seq
        self.color = 0
        super(Seq, self).__init__(**kwds)

    def __len__(self):
        return len(self.seq)
    
    def kmers(self, K):
        n = (len(self.seq) - K) + 1
        for i in range(0, n):
            yield self.seq[i:i+K]
    
    def paint(self, pColor):
        self.color = pColor

class AlignedSeq(Seq):

    def __init__(self, **kwds):
        self.hits = []
        super(AlignedSeq, self).__init__(**kwds)
    
    def addHit(self, hit):
        self.hits.append(hit)
    
    def addHits(self, hits):
        self.hits.extend(hits)
    
    def numHits(self):
        return len(self.hits)
    
    
class Partition:
    
    def __init__(self, pid):
        self.pid = pid
        self.seqs = []
    
    def addSeq(self, seq):
        self.seqs.append(seq)
    
    def addSeqs(self, newSeqs):
        self.seqs.extend(newSeqs)
    
    def numSeqs(self):
        return len(self.seqs)
    
    def numHits(self):
        return sum([s.numHits() for s in self.seqs])
        
    def consume(self, other):
        # we're consuming another Partition with the same pid, or should be
        assert self.pid == other.pid
        for s in other.seqs:
            self.addSeq(s)
    
    def __iter__(self):
        for seq in self.seqs:
            yield seq
    
    def __str__(self):
        info = "Partition {}: {} sequences, {} hits".format(self.numSeqs, self.numHits)
        return info
    
    def __repr__(self):
        return "{}  {}".format(self.pid, self.numSeqs())
    
    # Construct a new Partition from two other partitions
    # There is no integrity check on the pid, so it's up to the
    # user to make sure she doesn't introduce conflicts
    @classmethod
    def merge(cls, pA, pB, newPid):
        newPartition = cls(newPid)
        newPartition.addSeqs(pA.seqs)
        newPartition.addSeqs(pB.seqs)
        
        assert newPartition.numSeqs() == pA.numSeqs()+pB.numSeqs()
        assert newPartition.numHits() == pA.numHits()+pB.numHits()
        
        return newPartition
    
    @classmethod
    def new(cls, other):
        newP = cls(other.pid)
        for s in other:
            newP.addSeq(s)
        return newP

    # reads a sequence file with tab-separated partition ids
    @staticmethod
    def fastaPartitionIter(filename):
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
    
class PartitionSet:

    def __init__(self):
        self.partitions = {}
    
    # add a partition to the set;
    # if it exists, it will be consumed by the *existing* partition
    def addPartition(self, P):
        if P.pid in self.partitions:
            self.partitions[P.pid].consume(P)
        else:
            self.partitions[P.pid] = P

    # add a sequence to the set; if the partition exists,
    # add it, otherwise construct a new partition for it
    def addSeq(self, seq, pid):
        if pid in self.partitions:
            self.partitions[pid].addSeq(seq)
        else:
            newP = Partition(pid)
            newP.addSeq(seq)
            self.partitions[pid] = newP
    
    def numSeqs(self):
        n = 0
        for p in self.partitions.viewvalues():
            n += p.numSeqs()
        return n
    
    def numHits(self):
        n = 0
        for p in self.partitions.viewvalues():
            n += p.numHits()
        return n
    
    def avgSize(self):
        m = float(len(self))
        n = float(self.numSeqs())
        return n/m
    
    def histogram(self):
        print "Generating histogram..."
        h = self.partitions.values()
        h.sort(key=lambda p: p.numSeqs())
        hist = ["{},{}\n".format(p.pid, p.numSeqs()) for p in h]
        print "...done!"
        return hist
    
    # union on two partition sets
    def __add__(self, other):
        newPS = PartitionSet()
        for p in self:
            newPS.addPartition(Partition.new(p))
        for p in other:
            newPS.addPartition(Partition.new(p))
        return newPS
    
    # one-way difference, ie things only in this
    def __sub__(self, other):
        newPS = PartitionSet()
        for p in self:
            if p not in other:
                newPS.addPartition(Partition.new(p))
        return newPS
    
    # intersect
    def __and__(self, other):
        newPS = PartitionSet()
        for p in self:
            if p in other:
                newPS.addPartition(Partition.new(p))
        return newPS
    
    def __xor__(self, other):
        return (self + other) - (self & other)
    
    # iterator over the partition objects
    def __iter__(self):
        for p in self.partitions.viewvalues():
            yield p
    
    def __len__(self):
        return len(self.partitions)
    
    def __repr__(self):
        info = '''# PartitionSet containing:
        # {} partitions
        # comprising {} seqs
        # with {} seqs/partition avg\n'''.format(len(self), self.numSeqs(), round(self.avgSize(),3))
        return info
    
    # does the set contain a partition with this 
    # partition *object's* pid
    def __contains__(self, other):
        if other.pid in self.partitions:
            return True
        return False
    
    # given a partitioned fasta and a corresponding blast file,
    # construct a new PartitionSet with associated hits
    @classmethod
    def associateHits(cls, fastp_file, blast_file):
        ps = cls()
        bdp = BlastDB(blast_file)
        hits = bdp.partitionByQuery()
        print "Reading sequences and associating hits..."
        for name, pid, seq in Partition.fastaPartitionIter(fastp_file):
            s = AlignedSeq(name=name, seq=seq)
            if name in hits:
                s.addHits(hits[name])
            ps.addSeq(s, pid)
        print "...done!"
        return ps


'''
Kitchen-sink class for holding partitioned sequences, doing set ops,
associating homologies, etc
'''
# TODO: This is a fucking travesty. Subclass shit, fix the style conventions,
# maybe just delete it and start again...
class SeqSet:

    db = {}
    dbkeys = db.viewkeys()
    hits = {}
    partCounts = {}
    pids = partCounts.viewkeys()
    partHitCounts = {}
    
    def __init__(self, db={}, hits={}, partCounts={}, partHitCounts={}):
        self.db = db
        self.dbkeys = self.db.viewkeys()
        self.hits = hits
        self.partCounts = partCounts
        self.pids = self.partCounts.viewkeys()
        self.partHitCounts = partHitCounts
        self.checkIntegrity()
    
    def checkIntegrity(self):
        print "checking basic SeqSet integrity..."
        assert len(self.db) == len(self.dbkeys)
        assert len(self.hits) <= len(self.dbkeys)
        assert len(self.partCounts) == len(self.partHitCounts)
        for pid in self.pids:
            assert self.partCounts[pid] >= self.partHitCounts[pid]
        for pid in self.partHitCounts:
            assert pid in self.pids
        print "...done with integrity check."
    
    @classmethod
    def init_fromFiles(cls, seqFile, blastcsv):
        hits = parse_blast_distinctly(blastcsv)
        db = {}
        partCounts = {}
        partHitCounts = {}
        
        print "loading {} and {} into seqSet...".format(seqFile, blastcsv)
        for n, name, pid, seq in cls.partitionedSeqIter(seqFile):
            if n % 10000 == 0:
                print "loaded {} seqs in {}...".format(n, seqFile)
            if pid not in partHitCounts:
                partHitCounts[pid] = 0
            # check if the seq has a homology
            hit = False
            if name in hits:
                hit = True
                partHitCounts[pid] = partHitCounts.get(pid, 0) + 1
            # update the partition counts
            partCounts[pid] = partCounts.get(pid, 0) + 1
            db[name] = {'id':n, 'pid':pid, 'seq':seq, 'hit':hit}

        newhits = {k: hits[k] for k in db.viewkeys()}
        hits = newhits
        newSeqSet = cls(db, hits, partCounts, partHitCounts)
                
        print "...finished loading {} seqs and {} hits from {} and {}".format(len(newSeqSet), newSeqSet.numHits(), seqFile, blastcsv)
        return newSeqSet
    
    '''
    Builds a new SeqSet from the partitions shared between
    this SeqSet and otherSeqSet
    '''
    def init_partIntersect(self, otherSeqSet):
        newSeqSet = SeqSet({},{},{},{})
        intersectPids = self.pids & otherSeqSet.pids
        for key, record in self.iter_seqsByPidSet(intersectPids):
            newSeqSet.copySeq(key, self)
        for key, record in otherSeqSet.iter_seqsByPidSet(intersectPids):
            newSeqSet.copySeq(key, otherSeqSet)
        return newSeqSet
    '''
    Builds a new SeqSet from partitions found
    by the difference between self and otherSeqSet
    '''
    def init_partDifference(self, otherSeqSet):
        newSeqSet = SeqSet({},{},{},{})
        for key, record in self:
            pid = record['pid']
            if pid not in otherSeqSet.pids:
                newSeqSet.copySeq(key, self)
        return newSeqSet
    '''
    Builds the union between self and otherSeqSet
    '''
    def init_union(self, otherSeqSet):
        newSeqSet = SeqSet({},{},{},{})
        for key, record in self:
            newSeqSet.copySeq(key, self)
        for key, record in otherSeqSet:
            newSeqSet.copySeq(key, otherSeqSet)
        newSeqSet.checkIntegrity()
        return newSeqSet
    
    '''
    copies the sequence record, hit info, and partition info for a sequence into self
    '''
    def copySeq(self, fromKey, fromSeqSet):
        record = fromSeqSet.db[fromKey]
        self.db[fromKey] = record
        pid = record['pid']
        self.partCounts[pid] = self.partCounts.get(pid, 0) + 1
        if pid not in self.partHitCounts:
            self.partHitCounts[pid] = 0
        if record['hit'] == True:
            self.hits[fromKey] = fromSeqSet.hits[fromKey]
            self.partHitCounts[pid] = self.partHitCounts.get(pid, 0) + 1
    
    def dumpPartCounts(self, outname):
        with open(outname, 'wb') as out_fp:
            for pid in self.partCounts:
                out_fp.write('{}\n'.format(self.partCounts[pid]))
    
    # Dump seqs in outkeys to a fasta    
    def dumpSubsetFa(self, outname, outkeys):
        with open(outname, 'wb') as out_fp:
            for key, record in self:
                if key in outkeys:
                    out_fp.write('>{}\t{}\n{}\n'.format(key, record['pid'], record['seq']))
    
    # Dump seqs with hits to fasta
    def dumpSeqsHitsFa(self, outname):
        with open(outname, 'wb') as out_fp:
            for key, record in self:
                if key in self.hits:
                    out_fp.write('>{}\t{}\n{}\n'.format(key, record['pid'], record['seq']))
    
    # Dump seqs without hits to fasta
    def dumpSeqsPartHitsFa(self, outname):
        with open(outname, 'wb') as out_fp:
            for key, record in self:
                pid = self.db[key]['pid']
                if self.partHitCounts.get(pid, 0) > 0:
                    out_fp.write('>{}\t{}\n{}\n'.format(key, record['pid'], record['seq']))
    
    def dumpSeqsNoPartHitsFa(self, outname):
        with open(outname, 'wb') as out_fp:
            for key, record in self:
                pid = self.db[key]['pid']
                if self.partHitCounts.get(pid, 0) == 0:
                    out_fp.write('>{}\t{}\n{}\n'.format(key, record['pid'], record['seq']))
    
    # Dump seqs without hits to fasta
    def dumpSeqsNoHitsFa(self, outname):
        with open(outname, 'wb') as out_fp:
            for key, record in self:
                if key not in self.hits:
                    out_fp.write('>{}\t{}\n{}\n'.format(key, record['pid'], record['seq']))
    
    # Dump all sequences to outname
    def dumpSeqsFa(self, outname):
        with open(outname, 'wb') as out_fp:
            for key, record in self:
                out_fp.write('>{}\t{}\n{}\n'.format(key, record['pid'], record['seq']))
    
    # iterator over sequences with blast hits
    def iter_seqsWithHits(self):
        for key in self.dbkeys:
            if self.db[key]['hit'] == True:
                yield key, self.db[key]
    # iterator over sequences with the given partition id
    def iter_seqsByPid(self, pid):
        for key in self.dbkeys:
            if self.db[key]['pid'] == pid:
                yield key, self.db[key]
    
    # iterator over the sequences with any of the partition ids in pids
    def iter_seqsByPidSet(self, pids):
        for key in self.dbkeys:
            if self.db[key]['pid'] in pids:
                yield key, self.db[key]
                
    def iter_hits(self):
        for key in self.hits:
            yield key, self.hits[key]
    
    
    # store a counting hash for reference use
    def attach_countingHash(self, ht):
        self.refht = ht
    
    # build a counting hash from my own sequences
    def build_countingHash(self, ksize, htsize):
        self.ht = khmer.new_counting_hash(ksize, htsize, 4)
        for key, record in self:
            self.ht.consume(record['seq'])
    
    # reads a sequence file with tab-separated partition ids
    @staticmethod
    def partitionedSeqIter(filename):
        for n, record in enumerate(screed.open(filename, parse_description=False)):

            name = record.name
            name, partition_id = name.rsplit('\t', 1)

            # convert name to blast format if necessary
            nname = name.split('|', 2)
            if len(nname) >= 2:
                name = nname[2]
            name = name.split(' ')[0]    
            
            if n % 10000 == 0:
                print "seq: {},{}".format(n, name)
                
            yield n, name, int(partition_id), record.sequence
    
    def getPartCounts(self):
        return list(self.partCounts.viewvalues())
    
    # danger: needs testing...
    def getPartHitCounts(self):
        return zip(self.partCounts.viewvalues(), self.partHitCounts.viewvalues())
    
    def numParts(self):
        return len(self.pids)
    
    def numPartHits(self):
        i = 0
        for pid in self.pids:
            if self.partHitCounts[pid] > 0:
                i += 1
        return i
    
    def numHits(self):
        return len(self.hits)
    
    def printStats(self, name, header=False):
        if header:
            print "setName, numSeqs, numHits, numParts, numPartHits" 
        print "{}, {}, {}, {}, {}".format(name, len(self), self.numHits(), self.numParts(), self.numPartHits())
    
    def about(self, name):
        print '#' * 40
        self.printStats(name, True)
        print "\nSeqs:"
        for key, record in self:
            print key, record
        print "\nHits:"
        for key, hit in self.iter_hits():
            print key, hit
        print '#' * 40
    
    # returns stats as a string for file ops
    def __str__(self):
        return "{}, {}, {}, {}".format(len(self), self.numHits(), self.numParts(), self.numPartHits())
    
    # number of seqs in set
    def __len__(self):
        return len(self.dbkeys)
    
    # iterator over sequences
    def __iter__(self):
        for key in self.dbkeys:
            yield key, self.db[key]
    
    # method allowing for membership tests
    def __contains__(self, item):
        if item in self.dbkeys:
            return True
        return False
