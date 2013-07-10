include ../Makefile.inc

OBJS = petMar2.fa.nhr petMar2.fa.masked.nhr est.fa.nhr cdna.all.fa.nhr \
	mrna.all.fa.nhr mouse.protein.fa.phr zebrafish.protein.fa.phr \
	myxinidae.protein.fa.phr amphioxus.protein.fa.phr

all: $(OBJS)

clean:
	echo Running clean on db directory
	rm -f *.phr *.nhr *.fa .fa.masked

nodeps:
	touch petMar2.fa.nhr
	touch petMar2.fa.masked.nhr 
	touch est.fa.nhr 
	touch cdna.all.fa.nhr
	touch mrna.all.fa.nhr
	touch mouse.protein.fa.phr
	touch zebrafish.protein.fa.phr
	touch myxinidae.protein.fa.phr
	touch amphioxus.protein.fa.phr
	touch univec_core.fa.phr

petMar2.fa:
	curl http://hgdownload.soe.ucsc.edu/goldenPath/petMar2/bigZips/petMar2.fa.gz -o petMar2.fa.gz
	gunzip petMar2.fa.gz

petMar2.fa.nhr: petMar2.fa
	makeblastdb -in petMar2.fa -dbtype nucl

petMar2.fa.out:
	curl http://hgdownload.soe.ucsc.edu/goldenPath/petMar2/bigZips/petMar2.fa.out.gz -o petMar2.fa.out.gz
	gunzip petMar2.fa.out.gz

petMar2.fa.masked:
	curl http://hgdownload.soe.ucsc.edu/goldenPath/petMar2/bigZips/petMar2.fa.masked.gz -o petMar2.fa.masked.gz
	gunzip petMar2.fa.masked.gz

petMar2.fa.masked.nhr: petMar2.fa.masked
	makeblastdb -in petMar2.fa.masked -dbtype nucl

est.fa:
	curl http://hgdownload.soe.ucsc.edu/goldenPath/petMar2/bigZips/est.fa.gz -o est.fa.gz
	gunzip est.fa.gz

est.fa.nhr: est.fa
	makeblastdb -in est.fa -dbtype nucl

mrna.fa:
	curl http://hgdownload.soe.ucsc.edu/goldenPath/petMar2/bigZips/mrna.fa.gz -o mrna.fa.gz
	gunzip mrna.fa.gz

mrna.fa.nhr: mrna.fa
	makeblastdb -in mrna.fa -dbtype nucl

cdna.all.fa:
	curl -O ftp://ftp.ensembl.org/pub/release-72/fasta/petromyzon_marinus/cdna/Petromyzon_marinus.Pmarinus_7.0.72.cdna.all.fa.gz
	gunzip Petromyzon_marinus.Pmarinus_7.0.72.cdna.all.fa.gz
	mv Petromyzon_marinus.Pmarinus_7.0.72.cdna.all.fa cdna.all.fa

cdna.all.fa.nhr: cdna.all.fa
	makeblastdb -in cdna.all.fa -dbtype nucl

mouse.protein.fa:
	curl ftp://ftp.ensembl.org/pub/release-71/fasta/mus_musculus/pep/Mus_musculus.GRCm38.71.pep.all.fa.gz -o Mus_musculus.GRCm38.71.pep.all.fa.gz
	gunzip -c Mus_musculus.GRCm38.71.pep.all.fa.gz > mouse.protein.fa

mouse.protein.fa.phr: mouse.protein.fa
	makeblastdb -in mouse.protein.fa -dbtype prot

zebrafish.protein.fa:
	curl ftp://ftp.ensembl.org/pub/release-71/fasta/danio_rerio/pep/Danio_rerio.Zv9.71.pep.all.fa.gz -o Danio_rerio.Zv9.71.pep.all.fa.gz
	gunzip -c Danio_rerio.Zv9.71.pep.all.fa.gz > zebrafish.protein.fa

zebrafish.protein.fa.phr: zebrafish.protein.fa
	makeblastdb -in zebrafish.protein.fa -dbtype prot

myxinidea.protein.fa:	
	perl ../scripts/get_uniprot_fam.pl Myxinidae
	mv Myxinidae.fasta myxinidae.protein.fa

myxinidea.protein.fa.phr: myxinidea.protein.fa
	makeblastdb -in myxinidae.protein.fa -dbtype prot
	
amphioxus.protein.fa:
	perl scripts/get_uniprot_org.pl 7739
	mv 7739.fasta amphioxus.protein.fa
	
amphioxus.protein.fa.phr: amphioxus.protein.fa
	makeblastdb -in amphioxus.protein.fa -dbtype prot

univec_core.fa:
	curl ftp://ftp.ncbi.nih.gov/pub/UniVec/UniVec_Core -o univec_core.fa

univec_core.fa.nhr: univec_core.fa
	makeblastdb -in univec_core.fa -dbtype nucl

force:
	true
