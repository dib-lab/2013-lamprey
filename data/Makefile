include ../Makefile.inc

# Everything to take the assembly from raw to processed and ready for analysis
# For now, filenames are lamprey-specific

all: petMar_mrna.fp
	cp petMar_mrna.fp ../

lamp3.fasta: 
#	get <url to raw assembly, or include it>


lamp3.fasta.dedupe.fa: lamp3.fasta
	$(CDHITEST) -i lamp3.fasta -o lamp3.dedupe.fa -c .98 -n 9 -g 1 -T 4

petMar_mrna.fa: lamp3.fasta.dedupe.fa
	$(SEQCLEAN) lamp3.fasta.dedupe.fa -c 6 -v ../db/univec_core.fa -o petMar_mrna.fa

petMar_mrna.fp: petMar_mrna.fa

petMar_part.ht: petMar_mrna.fa
	$(KHMER)/scripts/load-graph.py -k 25 -N 4 -x 16e9 petMar_part petMar_mrna.fa
	
petMar_part.info: petMar_part.ht
	$(KHMER)/scripts/partition-graph.py -T 7 -s 1e6 petMar_part

petMar_part.pmap.merged: petMar_part.info
	$(KHMER)/scripts/merge-partitions.py -k 25 petMar_part

petMar_mrna.fa.part: petMar_part.pmap.merged
	$(KHMER)/scripts/annotate-partitions.py -k 25 petMar_part petMar_mrna.fa
	
petMar_mrna_lump.fp: petMar_part.fa.part
	python ../scripts/dump-lump.py petMar_mrna.fa.part petMar_mrna
	
petMar_lump.ht: petMar_mrna_lump.fp
	$(KHMER)/scripts/load-graph.py -x 16e9 -k 25 petMar_lump petMar_mrna_lump.fp

petMar_lump.stoptags: petMar_lump.ht
	$(KHMER)/scripts/make-initial-stoptags.py petMar_lump

petMar_lump.info: petMar_lump.stoptags
	$(KHMER)/scripts/partition-graph.py -T 8 --stoptags petMar_lump.stoptags petMar_lump

	$(KHMER)/scripts/find-knots.py -x 2e8 -N 4 petMar_lump
#TODO: check these deps once internet comes back up...
petMar_mrna_lump.fp.stopfilt:
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
