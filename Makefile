KHMER=~/tools/khmer
PEASOUP=/w/peasoup/scripts
CSV="qseqid sseqid length nident qstart qend sstart send bitscore evalue"

genome:
	curl http://hgdownload.soe.ucsc.edu/goldenPath/petMar2/bigZips/petMar2.fa.gz -o db/petMar2.fa.gz
	curl http://hgdownload.soe.ucsc.edu/goldenPath/petMar2/bigZips/petMar2.fa.out.gz -o db/petMar2.fa.out.gz
	curl http://hgdownload.soe.ucsc.edu/goldenPath/petMar2/bigZips/petMar2.fa.masked.gz -o db/petMar2.fa.masked.gz
	curl http://hgdownload.soe.ucsc.edu/goldenPath/petMar2/bigZips/est.fa.gz -o db/est.fa.gz
	curl http://hgdownload.soe.ucsc.edu/goldenPath/petMar2/bigZips/mrna.fa.gz -o db/mrna.fa.gz
	
	gunzip db/petMar2.fa.gz
	gunzip db/petMar2.fa.masked.gz
	gunzip db/est.fa.gz
	gunzip db/mrna.fa.gz
	
proteins:
	curl ftp://ftp.ensembl.org/pub/release-71/fasta/mus_musculus/pep/Mus_musculus.GRCm38.71.pep.all.fa.gz -o db/Mus_musculus.GRCm38.71.pep.all.fa.gz
	curl ftp://ftp.ensembl.org/pub/release-71/fasta/danio_rerio/pep/Danio_rerio.Zv9.71.pep.all.fa.gz -o db/Danio_rerio.Zv9.71.pep.all.fa.gz
	
	perl scripts/get_uniprot_fam.pl Myxinidae
	mv Myxinidae.fasta db/myxinidae.protein.fa
	
	perl scripts/get_uniprot_org.pl 7739
	mv 7739.fasta db/amphioxus.protein.fa
	
	gunzip -c db/Mus_musculus.GRCm38.71.pep.all.fa.gz > db/mouse.protein.fa
	gunzip -c db/Danio_rerio.Zv9.71.pep.all.fa.gz > db/zebrafish.protein.fa

blastdb:
	makeblastdb -in db/myxinidae.protein.fa -dbtype prot
	makeblastdb -in db/amphioxus.protein.fa -dbtype prot
	makeblastdb -in db/zebrafish.protein.fa -dbtype prot
	makeblastdb -in db/mouse.protein.fa -dbtype prot
	
	makeblastdb -in db/petMar2.fa -dbtype nucl
	makeblastdb -in db/petMar2.fa.masked -dbtype nucl
	makeblastdb -in db/est.fa -dbtype nucl
	makeblastdb -in db/mrna.fa -dbtype nucl

lamp3.fasta.dedupe.fa: data/lamp3.fasta
	cd-hit-est -i data/lamp3.fasta -o data/lamp3.dedupe.fa -c .98 -n 9 -g 1 -T 4

petMar_mrna.fa: lamp3.fasta.dedupe.fa
	seqclean data/lamp3.fasta.dedupe.fa -c 6 -v db/UniVec.fa -o data/petMar_mrna.fa

petMar_mrna.fp: petMar_mrna.fa
	$KHMER/scripts/load-graph.py -k 25 -N 4 -x 16e9 petMar_part petMar_mrna.fa
	$KHMER/scripts/partition-graph.py -T 7 -s 1e6 petMar_part
	$KHMER/scripts/merge-partitions.py -k 25 petMar_part
	$KHMER/scripts/annotate-partitions.py -k 25 petMar_part petMar_mrna.fa
	
	python scripts/dump-lump.py petMar_mrna.fa.part petMar_mrna
	
	$KHMER/scripts/load-graph.py -x 16e9 -k 25 petMar_lump petMar_mrna_lump.fp
	$KHMER/scripts/make-initial-stoptags.py petMar_lump
	$KHMER/scripts/partition-graph.py -T 8 --stoptags petMar_lump.stoptags petMar_lump
	$KHMER/scripts/find-knots.py -x 2e8 -N 4 petMar_lump
	$KHMER/scripts/filter-stoptags.py -k 25 petMar_lump.stoptags petMar_mrna_lump.fp
	$KHMER/scripts/load-graph.py -x 8e9 -k 25 lumpfilt petMar_mrna_lump.fp.stopfilt
	$KHMER/scripts/partition-graph.py -T 8 lumpfilt
	$KHMER/scripts/merge-partitions.py -k 25 lumpfilt
	$KHMER/scripts/annotate-partitions.py -k 25 lumpfilt petMar_mrna_lump.fp.stopfilt
	
	mv petMar_mrna_lump.fp.stopfilt.part petMar_mrna_delumped.fp
	cat petMar_mrna_nonlump.fp petMar_mrna_delumped.fp > petMar_mrna.fp

x.genome.csv: petMar_mrna.fp
	blastn -query petMar_mrna.fp -db db/petMar2.fa.masked -out x.genome_hard.csv -num_threads 8 -evalue 0.0000001 -outfmt \"10 $CSV\"

x.zebrafish.csv: petMar_mrna.fp
	blastx -query petMar_mrna.fp -db db/zebrafish.protein.fa -out x.zebrafish.csv -num_threads 8 -evalue 0.0000001 -outfmt \"10 $CSV\"

x.myx.csv: petMar_mrna.fp
	blastx -query petMar_mrna.fp -db db/myxinidae.protein.fa -out x.myx.csv -num_threads 8 -evalue 0.0000001 -outfmt \"10 $CSV\"

x.amph.csv: petMar_mrna.fp
	blastx -query petMar_mrna.fp -db db/amphioxus.protein.fa -out x.amph.csv -num_threads 8 -evalue 0.0000001 -outfmt \"10 $CSV\"

x.est.csv: petMar_mrna.fp
	blastn -query petMar_mrna.fp -db db/est.fa -out x.EST.csv -num_threads 8 -evalue 0.0000001 -outfmt \"10 $CSV\"

x.lamp0.csv: petMar_mrna.fp
	blastn -query petMar_mrna.fp -db data/lamp0.fasta -out x.lamp0.csv -num_threads 8 -evalue 0.0000001 -outfmt \"10 $CSV\"
