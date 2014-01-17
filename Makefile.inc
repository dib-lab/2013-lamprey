KHMER=~/tools/khmer
PEASOUP=/w/peasoup/scripts
SEQCLEAN=seqclean
CDHITEST=cd-hit-est

CSV=qseqid sseqid length nident qstart qend sstart send bitscore evalue
BLAST_EVALUE=0.001
BLAST_CPUS=4
BLASTRES=nodes=1:ppn=9,walltime=48:00:00,mem=8gb
BLAST_OPTS=-num_threads $(BLAST_CPUS) -evalue $(BLAST_EVALUE) -outfmt "10 $(CSV)" 

define newline

endef
