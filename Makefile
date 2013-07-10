include Makefile.inc

# make all for non-cluster runs
# make pbs-blast && make petMar_mrna.fp for cluster rns
all:  petMar_mrna.fp x.genome.csv x.zebrafish.csv x.amph.csv \
	x.mouse.csv x.myx.csv x.est.csv x.lamp0.csv x.cdna.csv

db: force
	cd db; $(MAKE) all

data: force
	

# touch prereq files in order to skip dep building steps
nodep:
	touch lamp3.fasta
	touch lamp3.fasta.dedupe.fa
	touch petMar_mrna.fa
	touch petMar_mrna.fp


# these could be changed to if statements checking the HPCC_MODULES_PATH
# and the PBS_O_WORKDIR to see if we're on the cluster, and then if
# we've been submitted. (or, something more general to clusters? checking
# if the qsub command exists?)
pbs-blast: x.lamp0-pbs x.genome-pbs x.amph-pbs x.mouse-pbs x.zebrafish-pbs x.myx-pbs x.est-pbs
x.lamp0-pbs:
	qsub -l $(BLASTRES) lamp0_blast.pbs
x.genome-pbs:
	qsub -l $(BLASTRES) genome_blast.pbs	
x.amph-pbs:
	qsub -l $(BLASTRES) amph_blast.pbs	
x.mouse-pbs:
	qsub -l $(BLASTRES) mouse_blast.pbs	
x.zebrafish-pbs:
	qsub -l $(BLASTRES) zebrafish_blast.pbs	
x.myx-pbs:
	qsub -l $(BLASTRES) myx_blast.pbs	
x.est-pbs:
	qsub -l $(BLASTRES) est_blast.pbs
x.cdna-pbs:
	qsub -l $(BLASTRES) cdna_blast.pbs
	
petMar_mrna.fp: force
	cd db; $(MAKE) univec_core.fa.nhr
	cd data; $(MAKE) petMar_mrna.fp

x.genome.csv: petMar_mrna.fp
	cd db; $(MAKE) petMar2.fa.masked.phr
	blastn -query petMar_mrna.fp -db db/petMar2.fa.masked -out x.genome_hard.csv -num_threads $(BLAST_CPUS) -evalue 0.0000001 -outfmt "10 $(CSV)"

x.zebrafish.csv: petMar_mrna.fp
	cd db; $(MAKE) zebrafish.protein.fa.phr
	blastx -query petMar_mrna.fp -db db/zebrafish.protein.fa -out x.zebrafish.csv -num_threads $(BLAST_CPUS) -evalue 0.0000001 -outfmt "10 $(CSV)"

x.myx.csv: petMar_mrna.fp
	cd db; $(MAKE) myxinidae.protein.fa.phr
	blastx -query petMar_mrna.fp -db db/myxinidae.protein.fa -out x.myx.csv -num_threads $(BLAST_CPUS) -evalue 0.0000001 -outfmt "10 $(CSV)"

x.amph.csv: petMar_mrna.fp
	cd db; $(MAKE) amphioxus.protein.fa.phr
	blastx -query petMar_mrna.fp -db db/amphioxus.protein.fa -out x.amph.csv -num_threads $(BLAST_CPUS) -evalue 0.0000001 -outfmt "10 $(CSV)"
	
x.mouse.csv: petMar_mrna.fp
	cd db; $(MAKE) mouse.protein.fa.phr
	blastx -query petMar_mrna.fp -db db/mouse.protein.fa -out x.mouse.csv -num_threads $(BLAST_CPUS) -evalue 0.0000001 -outfmt "10 $(CSV)"

x.est.csv: petMar_mrna.fp
	cd db; $(MAKE) est.fa.nhr
	blastn -query petMar_mrna.fp -db db/est.fa -out x.est.csv -num_threads $(BLAST_CPUS) -evalue 0.0000001 -outfmt "10 $(CSV)"

x.cdna.csv: petMar_mrna.fp
	cd db; $(MAKE) cdna.all.fa.nhr
	blastn -query petMar_mrna.fp -db db/cdna.all.fa -out x.cdna.csv -num_threads $(BLAST_CPUS) -evalue 0.0000001 -outfmt "10 $(CSV)"

x.lamp0.csv: petMar_mrna.fp
	cd data; $(MAKE) lamp0.fasta.nhr
	blastn -query petMar_mrna.fp -db data/lamp0.fasta -out x.lamp0.csv -num_threads $(BLAST_CPUS) -evalue 0.0000001 -outfmt "10 $(CSV)"


force:
	true
