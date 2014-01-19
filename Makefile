include Makefile.inc

# make all for non-cluster runs
# make pbs-blast && make $(OUTPRF)_mrna.fp for cluster rns
all:  db data x.genome.csv x.zebrafish.csv x.amph.csv \
	x.mouse.csv x.myx.csv x.est.csv x.lamp0.csv x.cdna.csv

db: force
	cd db; $(MAKE) all

data: force
	cd data; $(MAKE) all

# TODO: add variable for making with cluster
pbs-blast: x.lamp0-pbs x.genome-pbs x.amph-pbs x.mouse-pbs x.zebrafish-pbs x.myx-pbs x.est-pbs
x.lamp0-pbs:
	echo make x.lamp0.csv | cat pbs/blast.sub - | qsub -l $(BLASTRES) -N x.lamp0.csv
x.genome-pbs:
	echo make x.genome.csv | cat pbs/blast.sub - | qsub -l $(BLASTRES) -N x.genome.csv	
x.amph-pbs:
	echo make x.amph.csv | cat pbs/blast.sub - | qsub -l $(BLASTRES) -N x.amph.csv	
x.mouse-pbs:
	echo make x.mouse.csv | cat pbs/blast.sub - | qsub -l $(BLASTRES) -N x.mouse.csv
x.zebrafish-pbs:
	echo make x.zebrafish.csv | cat pbs/blast.pbs - | qsub -l $(BLASTRES) -N x.zebrafish.csv	
x.myx-pbs:
	echo x.myx.csv | cat pbs/blast.sub - | qsub -l $(BLASTRES) -N x.myx.csv
x.est-pbs:
	echo x.est.csv | cat pbs/blast.sub - | qsub -l $(BLASTRES) -N x.est.csv
x.cdna-pbs:
	echo x.cdna.csv | cat pbs/blast.sub - | qsub -l $(BLASTRES) -N x.cdna.csv
	
$(OUTPRF)_mrna.fp:
	cd db; $(MAKE) univec_core.fa.nhr
	cd data; $(MAKE) all

x.genome.csv: $(OUTPRF).fp
	cd db; $(MAKE) petMar2.fa.masked.phr
	blastn -query $(OUTPRF)_mrna.fp -db db/petMar2.fa.masked -out x.genome_hard.csv -num_threads $(BLAST_CPUS) -evalue $(BLAST_EVALUE) -outfmt "10 $(CSV)"

x.zebrafish.csv: $(OUTPRF).fp
	cd db; $(MAKE) zebrafish.protein.fa.phr
	blastx -query $(OUTPRF)_mrna.fp -db db/zebrafish.protein.fa -out x.zebrafish.csv -num_threads $(BLAST_CPUS) -evalue $(BLAST_EVALUE) -outfmt "10 $(CSV)"

x.myx.csv: $(OUTPRF).fp
	cd db; $(MAKE) myxinidae.protein.fa.phr
	blastx -query $(OUTPRF)_mrna.fp -db db/myxinidae.protein.fa -out x.myx.csv -num_threads $(BLAST_CPUS) -evalue $(BLAST_EVALUE) -outfmt "10 $(CSV)"

x.amph.csv: $(OUTPRF).fp
	cd db; $(MAKE) amphioxus.protein.fa.phr
	blastx -query $(OUTPRF)_mrna.fp -db db/amphioxus.protein.fa -out x.amph.csv -num_threads $(BLAST_CPUS) -evalue $(BLAST_EVALUE) -outfmt "10 $(CSV)"
	
x.mouse.csv: $(OUTPRF).fp
	cd db; $(MAKE) mouse.protein.fa.phr
	blastx -query $(OUTPRF)_mrna.fp -db db/mouse.protein.fa -out x.mouse.csv -num_threads $(BLAST_CPUS) -evalue $(BLAST_EVALUE) -outfmt "10 $(CSV)"

x.est.csv: $(OUTPRF).fp
	cd db; $(MAKE) est.fa.nhr
	blastn -query $(OUTPRF)_mrna.fp -db db/est.fa -out x.est.csv -num_threads $(BLAST_CPUS) -evalue $(BLAST_EVALUE) -outfmt "10 $(CSV)"

x.cdna.csv: $(OUTPRF).fp
	cd db; $(MAKE) cdna.all.fa.nhr
	blastn -query $(OUTPRF)_mrna.fp -db db/cdna.all.fa -out x.cdna.csv -num_threads $(BLAST_CPUS) -evalue $(BLAST_EVALUE) -outfmt "10 $(CSV)"

x.lamp0.csv: $(OUTPRF).fp
	cd data; $(MAKE) lamp0.fasta.nhr
	blastn -query $(OUTPRF)_mrna.fp -db data/lamp0.fasta -out x.lamp0.csv $(BLAST_OPTS) 

x.human.protein.csv: $(OUTPRF).fp
	cd db; $(MAKE) human.protein.fa.phr
	blastx -query $< -db db/human.protein.fa -out x.human.protein.csv $(BLAST_OPTS)

x.human.cdna.csv: $(OUTPRF).fp
	cd db; $(MAKE) human.cdna.fa.nhr
	blastn -query $< -db db/human.cdna.fa -out x.human.cdna.csv $(BLAST_OPTS)
 

force:
	true
