include ../Makefile.inc

# Everything to take the assembly from raw to processed and ready for analysis
# For now, filenames are lamprey-specific

all: petMar_mrna.fp
	cp petMar_mrna.fp ../

clean:
	rm -f *.stoptags *.part *.ht *.fp *.stopfilt *.merged *.info petMar_mrna.fa

nodeps:
	touch petMar_mrna.fa
	touch petMar_part.ht
	touch petMar_part.info
	touch petMar_part.pmap.merged
	touch petMar_mrna.fa.part
	touch petMar_mrna_lump.fp
	touch petMar_lump.ht
	touch petMar_lump.stoptags
	touch petMar_lump.info
	touch petMar_mrna_lump.fp.stopfilt
	touch lumpfilt.ht
	touch lumpfilt.info
	touch lumpfilt.pmap.merged
	touch lumpfilt.fp.stopfilt.part
	touch petMar_mrna.fp

READS=*.fq.gz


petMar_reads_part.ht: $(READS)
	$(KHMER)/scripts/load-graph.py -k 25 -N 4 -x 16e9 petMar_reads_part $(READS)
	
petMar_reads_part.info: petMar_reads_part.ht
	$(KHMER)/scripts/partition-graph.py -T 7 -s 1e6 petMar_reads_part

petMar_reads_part.pmap.merged: petMar_reads_part.info
	$(KHMER)/scripts/merge-partitions.py -k 25 petMar_reads_part

petMar_mrna.fa.part: petMar_reads_part.pmap.merged
	$(KHMER)/scripts/annotate-partitions.py -k 25 petMar_reads_part petMar_mrna.fa
	
petMar_mrna_lump.fp: petMar_mrna.fa.part
	python ../scripts/dump-lump.py petMar_mrna.fa.part petMar_mrna
	
petMar_lump.ht: petMar_mrna_lump.fp
	$(KHMER)/scripts/load-graph.py -x 16e9 -k 25 petMar_lump petMar_mrna_lump.fp

petMar_lump.stoptags: petMar_lump.ht
	$(KHMER)/scripts/make-initial-stoptags.py petMar_lump

petMar_lump.info: petMar_lump.stoptags
	$(KHMER)/scripts/partition-graph.py -T 8 --stoptags petMar_lump.stoptags petMar_lump

# Unfortunately, needs two commands in one rule because of filename conventions
petMar_mrna_lump.fp.stopfilt: petMar_lump.info
	$(KHMER)/scripts/find-knots.py -x 2e8 -N 4 petMar_lump
	$(KHMER)/scripts/filter-stoptags.py -k 25 petMar_lump.stoptags petMar_mrna_lump.fp

lumpfilt.ht: petMar_mrna_lump.fp.stopfilt
	$(KHMER)/scripts/load-graph.py -x 8e9 -k 25 lumpfilt petMar_mrna_lump.fp.stopfilt

lumpfilt.info: lumpfilt.ht
	$(KHMER)/scripts/partition-graph.py -T 8 lumpfilt

lumpfilt.pmap.merged: lumpfilt.info
	$(KHMER)/scripts/merge-partitions.py -k 25 lumpfilt

lumpfilt.fp.stopfilt.part: lumpfilt.pmap.merged
	$(KHMER)/scripts/annotate-partitions.py -k 25 lumpfilt petMar_mrna_lump.fp.stopfilt

petMar_mrna.fp: lumpfilt.fp.stopfilt.part
	mv petMar_mrna_lump.fp.stopfilt.part petMar_mrna_delumped.fp
	cat petMar_mrna_nonlump.fp petMar_mrna_delumped.fp > petMar_mrna.fp

petMar_mrna.fp.nhr: petMar_mrna.fp
	makeblastdb -in petMar_mrna.fp -dbtype nucl
