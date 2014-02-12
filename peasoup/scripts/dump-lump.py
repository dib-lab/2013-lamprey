import sys
sys.path.append('/w/peasoup/scripts/')
import hambone as hb

infile = sys.argv[1]
out_basename = sys.argv[2]

lump_ps = hb.PartitionSet.fromFastp(infile)
dist = lump_ps.histogram()

print dist[0]
print lump_ps.numSeqs(), len(lump_ps), lump_ps.avgSize()

print "Dumping largest partition..."
with open(out_basename + '_lump.fp', 'wb') as outfp:
    lump_ps.writeFastp(outfp, dist[0][0])
with open(out_basename + '_nonlump.fp', 'wb') as outfp:
    for p in dist[1:]:
        lump_ps.writeFastp(outfp, p[0])
