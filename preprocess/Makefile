
include Makefile.inc


#
# First we run FastQC on all the raw files for baseline quality infos
# Keep the general rules out here, with all the invidual rules in the Makefile.inc
#



fastqc1: $(READS_FASTQC)

fastqc1-pbs: $(patsubst %.zip, %.zip-pbs, $READS_FASTQC)
	echo Submitting PBS jobs for $(patsubst %.zip, %.zip-pbs, $READS_FASTQC)

# Run Trimmomatic on all the raw files



trim: $(READS_TRIM)

trim-pbs: $(patsubst %.trim.fq, %.trim.fq-pbs, $(READS_TRIM))
	echo Submitting PBS jobs for $(patsubst %.trim.fq, %.trim.fq-pbs, $(READS_TRIM))
parallel-dg50-pbs: $(patsubst, %, %-pbs, READS_DG50)

parallel-dg50: READS_DG50

diginorm.c50.ht-pbs:
	echo make diginorm.c50.ht | cat diginorm.inc.sub - | qsub -l $(DGRES) -N lamp4-diginorm.c50


diginorm.c50.ht: $(READS_DG50)
	$(TIME) -o diginorm.c50.time $(KHMER)/scripts/normalize-by-median.py -x 4000000000.0 -N 4 -k 20 -C 50 -s diginorm.c50.ht -R diginorm.c50.report \
$(READS_DG50)
	for F in $(READS_DG50); do 	mv $(F) $(patsubst %.trim.fq.keep.keep, %.trim.keep.keep.fq); \
	done

parallel.c50.ht: $(READS_DG50)
	$(TIME) -o parallel.c50.ht.time $(KHMER)/scripts/load-into-counting.py -x 1e10 -N 4 -k 20 parallel.c50.ht $(READS_DG50)

parallel.c50.ht-pbs: $(DG50)
	echo make parallel.c50.ht | cat diginorm.inc.sub - | qsub -l $(DGRES) -N load-parallel.c50.ht

abundfilt: abundfilt.time

abundfilt-pbs:
	echo make abundfilt.time | cat diginorm.inc.sub - | qsub -l $(DGRES) -N lamp4-abundfilt


diginorm.c20.ht: abundfilt.time
	$(TIME) -o diginorm.c20.time $(KHMER)/scripts/normalize-by-median.py -C 20 -x 10000000000.0 -k 20 -s diginorm.c20.ht -R diginorm.c20.report $(READS_ABF)

diginorm.c20.ht-pbs:
	echo make diginorm.c20.ht | cat diginorm.inc.sub - | qsub -l $(DGRES) -N {prj}-diginorm.c20

clean:
	rm -rf *.time *.clean *.keep *.trim.fq *.report *.log *.trimlog *.ht cleaning_*/ *.cidx *.cln all.c20.fa all.c20.fq *.sort *.abundfilt
