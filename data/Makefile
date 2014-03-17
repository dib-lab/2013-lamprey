include ../Makefile.inc

# Everything to take the assembly from raw to processed and ready for analysis
# For now, filenames are lamprey-specific

all: $(OUTPRF).fp
	cp $(OUTPRF).fp ../

clean:
	rm -f *.tagset *.stoptags *.part *.ht *.fp *.stopfilt *.merged *.info $(OUTPRF).fa

nodeps:
	touch $(OUTPRF).fa
	touch $(OUTPRF)_part.ht
	touch $(OUTPRF)_part.info
	touch $(OUTPRF)_part.pmap.merged
	touch $(OUTPRF).fp
	touch $(OUTPRF)_lump.fp
	touch $(OUTPRF)_lump.ht
	touch $(OUTPRF)_lump.stoptags
	touch $(OUTPRF)_lump.info
	touch $(OUTPRF)_lump.fp.stopfilt
	touch lumpfilt.ht
	touch lumpfilt.info
	touch lumpfilt.pmap.merged
	touch lumpfilt.fp.stopfilt.part
	touch $(OUTPRF).fp

lamp0.fasta.nhr:
	makeblastdb -in lamp0.fasta -dbtype nucl

$(OUTPRF).fa: Trinity.fasta
	ln -s $< $@

$(OUTPRF)_part.ht: $(OUTPRF).fa
	$(KHMER)/scripts/load-graph.py -k 25 -N 4 -x $(HT) $(OUTPRF)_part $<
	
$(OUTPRF)_part.info: $(OUTPRF)_part.ht
	$(KHMER)/scripts/partition-graph.py -T 7 -s 1e6 $(OUTPRF)_part

$(OUTPRF)_part.pmap.merged: $(OUTPRF)_part.info
	$(KHMER)/scripts/merge-partitions.py -k 25 $(OUTPRF)_part

$(OUTPRF).fp: $(OUTPRF)_part.pmap.merged
	$(KHMER)/scripts/annotate-partitions.py -k 25 $(OUTPRF)_part $(OUTPRF).fa
	mv $(OUTPRF).fa.part $@

$(OUTPRF).renamed.fp.gz: $(OUTPRF).fp
	python ../eel-pond/rename-with-partitions.py $(OUTPRF) $<
	mv $<.renamed.fasta.gz $(OUTPRF).renamed.fp.gz

#
# Begin delumping process
#
	
$(OUTPRF)_lump.fp: $(OUTPRF).fp
	python ../scripts/dump-lump.py $< $(OUTPRF)

$(OUTPRF)_lump.ht: $(OUTPRF)_lump.fp
	$(KHMER)/scripts/load-graph.py -x $(HT) -k 25 $(OUTPRF)_lump $(OUTPRF)_lump.fp

$(OUTPRF)_lump.stoptags: $(OUTPRF)_lump.ht
	$(KHMER)/scripts/make-initial-stoptags.py $(OUTPRF)_lump

$(OUTPRF)_lump.info: $(OUTPRF)_lump.stoptags
	$(KHMER)/scripts/partition-graph.py -T 8 --stoptags $(OUTPRF)_lump.stoptags $(OUTPRF)_lump

# Unfortunately, needs two commands in one rule because of filename conventions
$(OUTPRF)_lump.fp.stopfilt: $(OUTPRF)_lump.info
	$(KHMER)/scripts/find-knots.py -x 2e8 -N 4 $(OUTPRF)_lump
	$(KHMER)/scripts/filter-stoptags.py -k 25 $(OUTPRF)_lump.stoptags $(OUTPRF)_lump.fp

lumpfilt.ht: $(OUTPRF)_lump.fp.stopfilt
	$(KHMER)/scripts/load-graph.py -x 8e9 -k 25 lumpfilt $(OUTPRF)_lump.fp.stopfilt

lumpfilt.info: lumpfilt.ht
	$(KHMER)/scripts/partition-graph.py -T 8 lumpfilt

lumpfilt.pmap.merged: lumpfilt.info
	$(KHMER)/scripts/merge-partitions.py -k 25 lumpfilt

lumpfilt.fp.stopfilt.part: lumpfilt.pmap.merged
	$(KHMER)/scripts/annotate-partitions.py -k 25 lumpfilt $(OUTPRF)_lump.fp.stopfilt

$(OUTPRF)_delumped.fp: lumpfilt.fp.stopfilt.part
	mv $(OUTPRF)_lump.fp.stopfilt.part $(OUTPRF)_delumped.fp; cat $(OUTPRF)_nonlump.fp $(OUTPRF)_delumped.fp > $(OUTPRF).fp

$(OUTPRF).fp.nhr: $(OUTPRF)_delumped.fp
	makeblastdb -in $(OUTPRF)_delumped.fp -dbtype nucl

$(OUTPRF).fa.nhr: $(OUTPRF).fa
	makeblastdb -in $< -dbtype nucl
