
KHMER=~/bin/khmer
SEQCLEAN=seqclean
TRIM=~/bin/Trimmomatic-0.25/trimmomatic-0.25.jar
TIME=/usr/bin/time -f %E%U%K%W%w%I%O
KMERGENIE=~/bin/kmergenie-1.5397/kmergenie

KMERRES=nodes=1:ppn=4,walltime=96:00:00,mem=32gb
TRIMRES=nodes=1:ppn=4,walltime=12:00:00,mem=4gb
QCRES=nodes=1:ppn=1,walltime=2:00:00,mem=4gb
DGRES=nodes=1:ppn=1,walltime=72:00:00,mem=42gb

ABYSSRES=nodes=1:ppn=24,walltime=3:00:00:00,mem=384gb
MIRARES=nodes=1:ppn=16,walltime=7:00:00:00,mem744gb

2_ACAGTG_L001_R1_001.fastqc1:
	fastqc --noextract 2_ACAGTG_L001_R1_001.fastq.gz

2_ACAGTG_L001_R1_001.fastqc1-pbs:
	qsub -l $(QCRES) 2_ACAGTG_L001_R1_001.fastqc1.sub

2_ACAGTG_L001_R2_001.fastqc1:
	fastqc --noextract 2_ACAGTG_L001_R2_001.fastq.gz

2_ACAGTG_L001_R2_001.fastqc1-pbs:
	qsub -l $(QCRES) 2_ACAGTG_L001_R2_001.fastqc1.sub

2_ACAGTG_L002_R1_001.fastqc1:
	fastqc --noextract 2_ACAGTG_L002_R1_001.fastq.gz

2_ACAGTG_L002_R1_001.fastqc1-pbs:
	qsub -l $(QCRES) 2_ACAGTG_L002_R1_001.fastqc1.sub

2_ACAGTG_L002_R2_001.fastqc1:
	fastqc --noextract 2_ACAGTG_L002_R2_001.fastq.gz

2_ACAGTG_L002_R2_001.fastqc1-pbs:
	qsub -l $(QCRES) 2_ACAGTG_L002_R2_001.fastqc1.sub

4_GCCAAT_L001_R1_001.fastqc1:
	fastqc --noextract 4_GCCAAT_L001_R1_001.fastq.gz

4_GCCAAT_L001_R1_001.fastqc1-pbs:
	qsub -l $(QCRES) 4_GCCAAT_L001_R1_001.fastqc1.sub

4_GCCAAT_L001_R2_001.fastqc1:
	fastqc --noextract 4_GCCAAT_L001_R2_001.fastq.gz

4_GCCAAT_L001_R2_001.fastqc1-pbs:
	qsub -l $(QCRES) 4_GCCAAT_L001_R2_001.fastqc1.sub

4_GCCAAT_L002_R1_001.fastqc1:
	fastqc --noextract 4_GCCAAT_L002_R1_001.fastq.gz

4_GCCAAT_L002_R1_001.fastqc1-pbs:
	qsub -l $(QCRES) 4_GCCAAT_L002_R1_001.fastqc1.sub

4_GCCAAT_L002_R2_001.fastqc1:
	fastqc --noextract 4_GCCAAT_L002_R2_001.fastq.gz

4_GCCAAT_L002_R2_001.fastqc1-pbs:
	qsub -l $(QCRES) 4_GCCAAT_L002_R2_001.fastqc1.sub

5_CTTGTA_L001_R1_001.fastqc1:
	fastqc --noextract 5_CTTGTA_L001_R1_001.fastq.gz

5_CTTGTA_L001_R1_001.fastqc1-pbs:
	qsub -l $(QCRES) 5_CTTGTA_L001_R1_001.fastqc1.sub

5_CTTGTA_L001_R2_001.fastqc1:
	fastqc --noextract 5_CTTGTA_L001_R2_001.fastq.gz

5_CTTGTA_L001_R2_001.fastqc1-pbs:
	qsub -l $(QCRES) 5_CTTGTA_L001_R2_001.fastqc1.sub

5_CTTGTA_L002_R1_001.fastqc1:
	fastqc --noextract 5_CTTGTA_L002_R1_001.fastq.gz

5_CTTGTA_L002_R1_001.fastqc1-pbs:
	qsub -l $(QCRES) 5_CTTGTA_L002_R1_001.fastqc1.sub

5_CTTGTA_L002_R2_001.fastqc1:
	fastqc --noextract 5_CTTGTA_L002_R2_001.fastq.gz

5_CTTGTA_L002_R2_001.fastqc1-pbs:
	qsub -l $(QCRES) 5_CTTGTA_L002_R2_001.fastqc1.sub

6_GTGAAA_L001_R1_001.fastqc1:
	fastqc --noextract 6_GTGAAA_L001_R1_001.fastq.gz

6_GTGAAA_L001_R1_001.fastqc1-pbs:
	qsub -l $(QCRES) 6_GTGAAA_L001_R1_001.fastqc1.sub

6_GTGAAA_L001_R2_001.fastqc1:
	fastqc --noextract 6_GTGAAA_L001_R2_001.fastq.gz

6_GTGAAA_L001_R2_001.fastqc1-pbs:
	qsub -l $(QCRES) 6_GTGAAA_L001_R2_001.fastqc1.sub

6_GTGAAA_L002_R1_001.fastqc1:
	fastqc --noextract 6_GTGAAA_L002_R1_001.fastq.gz

6_GTGAAA_L002_R1_001.fastqc1-pbs:
	qsub -l $(QCRES) 6_GTGAAA_L002_R1_001.fastqc1.sub

6_GTGAAA_L002_R2_001.fastqc1:
	fastqc --noextract 6_GTGAAA_L002_R2_001.fastq.gz

6_GTGAAA_L002_R2_001.fastqc1-pbs:
	qsub -l $(QCRES) 6_GTGAAA_L002_R2_001.fastqc1.sub

fastqc1-pbs: 2_ACAGTG_L001_R1_001.fastqc1-pbs 2_ACAGTG_L001_R2_001.fastqc1-pbs 2_ACAGTG_L002_R1_001.fastqc1-pbs 2_ACAGTG_L002_R2_001.fastqc1-pbs 4_GCCAAT_L001_R1_001.fastqc1-pbs 4_GCCAAT_L001_R2_001.fastqc1-pbs 4_GCCAAT_L002_R1_001.fastqc1-pbs 4_GCCAAT_L002_R2_001.fastqc1-pbs 5_CTTGTA_L001_R1_001.fastqc1-pbs 5_CTTGTA_L001_R2_001.fastqc1-pbs 5_CTTGTA_L002_R1_001.fastqc1-pbs 5_CTTGTA_L002_R2_001.fastqc1-pbs 6_GTGAAA_L001_R1_001.fastqc1-pbs 6_GTGAAA_L001_R2_001.fastqc1-pbs 6_GTGAAA_L002_R1_001.fastqc1-pbs 6_GTGAAA_L002_R2_001.fastqc1-pbs
2_ACAGTG_L001_R1_001.trim.fq: 2_ACAGTG_L001_R1_001.fastq.gz
		$(TIME) -o 2_ACAGTG_L001_R1_001.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 2_ACAGTG_L001_R1_001.trimlog -threads 4 2_ACAGTG_L001_R1_001.fastq.gz 2_ACAGTG_L001_R1_001.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

2_ACAGTG_L001_R1_001.trim.fq-pbs:
	qsub -l $(TRIMRES) 2_ACAGTG_L001_R1_001.trim.sub

2_ACAGTG_L001_R2_001.trim.fq: 2_ACAGTG_L001_R2_001.fastq.gz
		$(TIME) -o 2_ACAGTG_L001_R2_001.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 2_ACAGTG_L001_R2_001.trimlog -threads 4 2_ACAGTG_L001_R2_001.fastq.gz 2_ACAGTG_L001_R2_001.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

2_ACAGTG_L001_R2_001.trim.fq-pbs:
	qsub -l $(TRIMRES) 2_ACAGTG_L001_R2_001.trim.sub

2_ACAGTG_L002_R1_001.trim.fq: 2_ACAGTG_L002_R1_001.fastq.gz
		$(TIME) -o 2_ACAGTG_L002_R1_001.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 2_ACAGTG_L002_R1_001.trimlog -threads 4 2_ACAGTG_L002_R1_001.fastq.gz 2_ACAGTG_L002_R1_001.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

2_ACAGTG_L002_R1_001.trim.fq-pbs:
	qsub -l $(TRIMRES) 2_ACAGTG_L002_R1_001.trim.sub

2_ACAGTG_L002_R2_001.trim.fq: 2_ACAGTG_L002_R2_001.fastq.gz
		$(TIME) -o 2_ACAGTG_L002_R2_001.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 2_ACAGTG_L002_R2_001.trimlog -threads 4 2_ACAGTG_L002_R2_001.fastq.gz 2_ACAGTG_L002_R2_001.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

2_ACAGTG_L002_R2_001.trim.fq-pbs:
	qsub -l $(TRIMRES) 2_ACAGTG_L002_R2_001.trim.sub

4_GCCAAT_L001_R1_001.trim.fq: 4_GCCAAT_L001_R1_001.fastq.gz
		$(TIME) -o 4_GCCAAT_L001_R1_001.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 4_GCCAAT_L001_R1_001.trimlog -threads 4 4_GCCAAT_L001_R1_001.fastq.gz 4_GCCAAT_L001_R1_001.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

4_GCCAAT_L001_R1_001.trim.fq-pbs:
	qsub -l $(TRIMRES) 4_GCCAAT_L001_R1_001.trim.sub

4_GCCAAT_L001_R2_001.trim.fq: 4_GCCAAT_L001_R2_001.fastq.gz
		$(TIME) -o 4_GCCAAT_L001_R2_001.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 4_GCCAAT_L001_R2_001.trimlog -threads 4 4_GCCAAT_L001_R2_001.fastq.gz 4_GCCAAT_L001_R2_001.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

4_GCCAAT_L001_R2_001.trim.fq-pbs:
	qsub -l $(TRIMRES) 4_GCCAAT_L001_R2_001.trim.sub

4_GCCAAT_L002_R1_001.trim.fq: 4_GCCAAT_L002_R1_001.fastq.gz
		$(TIME) -o 4_GCCAAT_L002_R1_001.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 4_GCCAAT_L002_R1_001.trimlog -threads 4 4_GCCAAT_L002_R1_001.fastq.gz 4_GCCAAT_L002_R1_001.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

4_GCCAAT_L002_R1_001.trim.fq-pbs:
	qsub -l $(TRIMRES) 4_GCCAAT_L002_R1_001.trim.sub

4_GCCAAT_L002_R2_001.trim.fq: 4_GCCAAT_L002_R2_001.fastq.gz
		$(TIME) -o 4_GCCAAT_L002_R2_001.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 4_GCCAAT_L002_R2_001.trimlog -threads 4 4_GCCAAT_L002_R2_001.fastq.gz 4_GCCAAT_L002_R2_001.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

4_GCCAAT_L002_R2_001.trim.fq-pbs:
	qsub -l $(TRIMRES) 4_GCCAAT_L002_R2_001.trim.sub

5_CTTGTA_L001_R1_001.trim.fq: 5_CTTGTA_L001_R1_001.fastq.gz
		$(TIME) -o 5_CTTGTA_L001_R1_001.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 5_CTTGTA_L001_R1_001.trimlog -threads 4 5_CTTGTA_L001_R1_001.fastq.gz 5_CTTGTA_L001_R1_001.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

5_CTTGTA_L001_R1_001.trim.fq-pbs:
	qsub -l $(TRIMRES) 5_CTTGTA_L001_R1_001.trim.sub

5_CTTGTA_L001_R2_001.trim.fq: 5_CTTGTA_L001_R2_001.fastq.gz
		$(TIME) -o 5_CTTGTA_L001_R2_001.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 5_CTTGTA_L001_R2_001.trimlog -threads 4 5_CTTGTA_L001_R2_001.fastq.gz 5_CTTGTA_L001_R2_001.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

5_CTTGTA_L001_R2_001.trim.fq-pbs:
	qsub -l $(TRIMRES) 5_CTTGTA_L001_R2_001.trim.sub

5_CTTGTA_L002_R1_001.trim.fq: 5_CTTGTA_L002_R1_001.fastq.gz
		$(TIME) -o 5_CTTGTA_L002_R1_001.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 5_CTTGTA_L002_R1_001.trimlog -threads 4 5_CTTGTA_L002_R1_001.fastq.gz 5_CTTGTA_L002_R1_001.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

5_CTTGTA_L002_R1_001.trim.fq-pbs:
	qsub -l $(TRIMRES) 5_CTTGTA_L002_R1_001.trim.sub

5_CTTGTA_L002_R2_001.trim.fq: 5_CTTGTA_L002_R2_001.fastq.gz
		$(TIME) -o 5_CTTGTA_L002_R2_001.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 5_CTTGTA_L002_R2_001.trimlog -threads 4 5_CTTGTA_L002_R2_001.fastq.gz 5_CTTGTA_L002_R2_001.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

5_CTTGTA_L002_R2_001.trim.fq-pbs:
	qsub -l $(TRIMRES) 5_CTTGTA_L002_R2_001.trim.sub

6_GTGAAA_L001_R1_001.trim.fq: 6_GTGAAA_L001_R1_001.fastq.gz
		$(TIME) -o 6_GTGAAA_L001_R1_001.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 6_GTGAAA_L001_R1_001.trimlog -threads 4 6_GTGAAA_L001_R1_001.fastq.gz 6_GTGAAA_L001_R1_001.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

6_GTGAAA_L001_R1_001.trim.fq-pbs:
	qsub -l $(TRIMRES) 6_GTGAAA_L001_R1_001.trim.sub

6_GTGAAA_L001_R2_001.trim.fq: 6_GTGAAA_L001_R2_001.fastq.gz
		$(TIME) -o 6_GTGAAA_L001_R2_001.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 6_GTGAAA_L001_R2_001.trimlog -threads 4 6_GTGAAA_L001_R2_001.fastq.gz 6_GTGAAA_L001_R2_001.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

6_GTGAAA_L001_R2_001.trim.fq-pbs:
	qsub -l $(TRIMRES) 6_GTGAAA_L001_R2_001.trim.sub

6_GTGAAA_L002_R1_001.trim.fq: 6_GTGAAA_L002_R1_001.fastq.gz
		$(TIME) -o 6_GTGAAA_L002_R1_001.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 6_GTGAAA_L002_R1_001.trimlog -threads 4 6_GTGAAA_L002_R1_001.fastq.gz 6_GTGAAA_L002_R1_001.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

6_GTGAAA_L002_R1_001.trim.fq-pbs:
	qsub -l $(TRIMRES) 6_GTGAAA_L002_R1_001.trim.sub

6_GTGAAA_L002_R2_001.trim.fq: 6_GTGAAA_L002_R2_001.fastq.gz
		$(TIME) -o 6_GTGAAA_L002_R2_001.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 6_GTGAAA_L002_R2_001.trimlog -threads 4 6_GTGAAA_L002_R2_001.fastq.gz 6_GTGAAA_L002_R2_001.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

6_GTGAAA_L002_R2_001.trim.fq-pbs:
	qsub -l $(TRIMRES) 6_GTGAAA_L002_R2_001.trim.sub

trim-pbs: 2_ACAGTG_L001_R1_001.trim.fq-pbs 2_ACAGTG_L001_R2_001.trim.fq-pbs 2_ACAGTG_L002_R1_001.trim.fq-pbs 2_ACAGTG_L002_R2_001.trim.fq-pbs 4_GCCAAT_L001_R1_001.trim.fq-pbs 4_GCCAAT_L001_R2_001.trim.fq-pbs 4_GCCAAT_L002_R1_001.trim.fq-pbs 4_GCCAAT_L002_R2_001.trim.fq-pbs 5_CTTGTA_L001_R1_001.trim.fq-pbs 5_CTTGTA_L001_R2_001.trim.fq-pbs 5_CTTGTA_L002_R1_001.trim.fq-pbs 5_CTTGTA_L002_R2_001.trim.fq-pbs 6_GTGAAA_L001_R1_001.trim.fq-pbs 6_GTGAAA_L001_R2_001.trim.fq-pbs 6_GTGAAA_L002_R1_001.trim.fq-pbs 6_GTGAAA_L002_R2_001.trim.fq-pbs
2_ACAGTG_L001_R1_001.trim.fq.fastqc2: 2_ACAGTG_L001_R1_001.trim.fq
	fastqc --noextract 2_ACAGTG_L001_R1_001.trim.fq

2_ACAGTG_L001_R1_001.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) 2_ACAGTG_L001_R1_001.trim.fq.fastqc2.sub

2_ACAGTG_L001_R2_001.trim.fq.fastqc2: 2_ACAGTG_L001_R2_001.trim.fq
	fastqc --noextract 2_ACAGTG_L001_R2_001.trim.fq

2_ACAGTG_L001_R2_001.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) 2_ACAGTG_L001_R2_001.trim.fq.fastqc2.sub

2_ACAGTG_L002_R1_001.trim.fq.fastqc2: 2_ACAGTG_L002_R1_001.trim.fq
	fastqc --noextract 2_ACAGTG_L002_R1_001.trim.fq

2_ACAGTG_L002_R1_001.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) 2_ACAGTG_L002_R1_001.trim.fq.fastqc2.sub

2_ACAGTG_L002_R2_001.trim.fq.fastqc2: 2_ACAGTG_L002_R2_001.trim.fq
	fastqc --noextract 2_ACAGTG_L002_R2_001.trim.fq

2_ACAGTG_L002_R2_001.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) 2_ACAGTG_L002_R2_001.trim.fq.fastqc2.sub

4_GCCAAT_L001_R1_001.trim.fq.fastqc2: 4_GCCAAT_L001_R1_001.trim.fq
	fastqc --noextract 4_GCCAAT_L001_R1_001.trim.fq

4_GCCAAT_L001_R1_001.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) 4_GCCAAT_L001_R1_001.trim.fq.fastqc2.sub

4_GCCAAT_L001_R2_001.trim.fq.fastqc2: 4_GCCAAT_L001_R2_001.trim.fq
	fastqc --noextract 4_GCCAAT_L001_R2_001.trim.fq

4_GCCAAT_L001_R2_001.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) 4_GCCAAT_L001_R2_001.trim.fq.fastqc2.sub

4_GCCAAT_L002_R1_001.trim.fq.fastqc2: 4_GCCAAT_L002_R1_001.trim.fq
	fastqc --noextract 4_GCCAAT_L002_R1_001.trim.fq

4_GCCAAT_L002_R1_001.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) 4_GCCAAT_L002_R1_001.trim.fq.fastqc2.sub

4_GCCAAT_L002_R2_001.trim.fq.fastqc2: 4_GCCAAT_L002_R2_001.trim.fq
	fastqc --noextract 4_GCCAAT_L002_R2_001.trim.fq

4_GCCAAT_L002_R2_001.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) 4_GCCAAT_L002_R2_001.trim.fq.fastqc2.sub

5_CTTGTA_L001_R1_001.trim.fq.fastqc2: 5_CTTGTA_L001_R1_001.trim.fq
	fastqc --noextract 5_CTTGTA_L001_R1_001.trim.fq

5_CTTGTA_L001_R1_001.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) 5_CTTGTA_L001_R1_001.trim.fq.fastqc2.sub

5_CTTGTA_L001_R2_001.trim.fq.fastqc2: 5_CTTGTA_L001_R2_001.trim.fq
	fastqc --noextract 5_CTTGTA_L001_R2_001.trim.fq

5_CTTGTA_L001_R2_001.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) 5_CTTGTA_L001_R2_001.trim.fq.fastqc2.sub

5_CTTGTA_L002_R1_001.trim.fq.fastqc2: 5_CTTGTA_L002_R1_001.trim.fq
	fastqc --noextract 5_CTTGTA_L002_R1_001.trim.fq

5_CTTGTA_L002_R1_001.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) 5_CTTGTA_L002_R1_001.trim.fq.fastqc2.sub

5_CTTGTA_L002_R2_001.trim.fq.fastqc2: 5_CTTGTA_L002_R2_001.trim.fq
	fastqc --noextract 5_CTTGTA_L002_R2_001.trim.fq

5_CTTGTA_L002_R2_001.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) 5_CTTGTA_L002_R2_001.trim.fq.fastqc2.sub

6_GTGAAA_L001_R1_001.trim.fq.fastqc2: 6_GTGAAA_L001_R1_001.trim.fq
	fastqc --noextract 6_GTGAAA_L001_R1_001.trim.fq

6_GTGAAA_L001_R1_001.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) 6_GTGAAA_L001_R1_001.trim.fq.fastqc2.sub

6_GTGAAA_L001_R2_001.trim.fq.fastqc2: 6_GTGAAA_L001_R2_001.trim.fq
	fastqc --noextract 6_GTGAAA_L001_R2_001.trim.fq

6_GTGAAA_L001_R2_001.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) 6_GTGAAA_L001_R2_001.trim.fq.fastqc2.sub

6_GTGAAA_L002_R1_001.trim.fq.fastqc2: 6_GTGAAA_L002_R1_001.trim.fq
	fastqc --noextract 6_GTGAAA_L002_R1_001.trim.fq

6_GTGAAA_L002_R1_001.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) 6_GTGAAA_L002_R1_001.trim.fq.fastqc2.sub

6_GTGAAA_L002_R2_001.trim.fq.fastqc2: 6_GTGAAA_L002_R2_001.trim.fq
	fastqc --noextract 6_GTGAAA_L002_R2_001.trim.fq

6_GTGAAA_L002_R2_001.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) 6_GTGAAA_L002_R2_001.trim.fq.fastqc2.sub

fastqc2-pbs: 2_ACAGTG_L001_R1_001.trim.fq.fastqc2-pbs 2_ACAGTG_L001_R2_001.trim.fq.fastqc2-pbs 2_ACAGTG_L002_R1_001.trim.fq.fastqc2-pbs 2_ACAGTG_L002_R2_001.trim.fq.fastqc2-pbs 4_GCCAAT_L001_R1_001.trim.fq.fastqc2-pbs 4_GCCAAT_L001_R2_001.trim.fq.fastqc2-pbs 4_GCCAAT_L002_R1_001.trim.fq.fastqc2-pbs 4_GCCAAT_L002_R2_001.trim.fq.fastqc2-pbs 5_CTTGTA_L001_R1_001.trim.fq.fastqc2-pbs 5_CTTGTA_L001_R2_001.trim.fq.fastqc2-pbs 5_CTTGTA_L002_R1_001.trim.fq.fastqc2-pbs 5_CTTGTA_L002_R2_001.trim.fq.fastqc2-pbs 6_GTGAAA_L001_R1_001.trim.fq.fastqc2-pbs 6_GTGAAA_L001_R2_001.trim.fq.fastqc2-pbs 6_GTGAAA_L002_R1_001.trim.fq.fastqc2-pbs 6_GTGAAA_L002_R2_001.trim.fq.fastqc2-pbs

diginorm.c5: 2_ACAGTG_L001_R1_001.trim.fq 2_ACAGTG_L001_R2_001.trim.fq 2_ACAGTG_L002_R1_001.trim.fq 2_ACAGTG_L002_R2_001.trim.fq 4_GCCAAT_L001_R1_001.trim.fq 4_GCCAAT_L001_R2_001.trim.fq 4_GCCAAT_L002_R1_001.trim.fq 4_GCCAAT_L002_R2_001.trim.fq 5_CTTGTA_L001_R1_001.trim.fq 5_CTTGTA_L001_R2_001.trim.fq 5_CTTGTA_L002_R1_001.trim.fq 5_CTTGTA_L002_R2_001.trim.fq 6_GTGAAA_L001_R1_001.trim.fq 6_GTGAAA_L001_R2_001.trim.fq 6_GTGAAA_L002_R1_001.trim.fq 6_GTGAAA_L002_R2_001.trim.fq
	$(TIME) -o diginorm.c5.time $(KHMER)/scripts/normalize-by-median.py -x 10000000000.0 -N 4 -k 20 -C 5 -s diginorm.c5.ht -R diginorm.c5.report 2_ACAGTG_L001_R1_001.trim.fq 2_ACAGTG_L001_R2_001.trim.fq 2_ACAGTG_L002_R1_001.trim.fq 2_ACAGTG_L002_R2_001.trim.fq 4_GCCAAT_L001_R1_001.trim.fq 4_GCCAAT_L001_R2_001.trim.fq 4_GCCAAT_L002_R1_001.trim.fq 4_GCCAAT_L002_R2_001.trim.fq 5_CTTGTA_L001_R1_001.trim.fq 5_CTTGTA_L001_R2_001.trim.fq 5_CTTGTA_L002_R1_001.trim.fq 5_CTTGTA_L002_R2_001.trim.fq 6_GTGAAA_L001_R1_001.trim.fq 6_GTGAAA_L001_R2_001.trim.fq 6_GTGAAA_L002_R1_001.trim.fq 6_GTGAAA_L002_R2_001.trim.fq

diginorm.c5-pbs:
	qsub -l $(DGRES) diginorm.c5.sub

clean:
	rm -rf *.time *.clean *.keep *.trim.fq *.report *.log *.trimlog *.ht cleaning_*/ *.cidx *.cln all.c20.fa all.c20.fq *.sort


kmergenie: 2_ACAGTG_L001_R1_001.trim.fq 2_ACAGTG_L001_R2_001.trim.fq 2_ACAGTG_L002_R1_001.trim.fq 2_ACAGTG_L002_R2_001.trim.fq 4_GCCAAT_L001_R1_001.trim.fq 4_GCCAAT_L001_R2_001.trim.fq 4_GCCAAT_L002_R1_001.trim.fq 4_GCCAAT_L002_R2_001.trim.fq 5_CTTGTA_L001_R1_001.trim.fq 5_CTTGTA_L001_R2_001.trim.fq 5_CTTGTA_L002_R1_001.trim.fq 5_CTTGTA_L002_R2_001.trim.fq 6_GTGAAA_L001_R1_001.trim.fq 6_GTGAAA_L001_R2_001.trim.fq 6_GTGAAA_L002_R1_001.trim.fq 6_GTGAAA_L002_R2_001.trim.fq
	$(TIME) -o kmergenie.time $(KMERGENIE) reads.txt  --diploid -k 100 -l 20 

kmergenie-pbs:
	qsub -l $(KMERRES) kmergenie.sub


dg.kmergenie: 2_ACAGTG_L001_R1_001.trim.fq.keep 2_ACAGTG_L001_R2_001.trim.fq.keep 2_ACAGTG_L002_R1_001.trim.fq.keep 2_ACAGTG_L002_R2_001.trim.fq.keep 4_GCCAAT_L001_R1_001.trim.fq.keep 4_GCCAAT_L001_R2_001.trim.fq.keep 4_GCCAAT_L002_R1_001.trim.fq.keep 4_GCCAAT_L002_R2_001.trim.fq.keep 5_CTTGTA_L001_R1_001.trim.fq.keep 5_CTTGTA_L001_R2_001.trim.fq.keep 5_CTTGTA_L002_R1_001.trim.fq.keep 5_CTTGTA_L002_R2_001.trim.fq.keep 6_GTGAAA_L001_R1_001.trim.fq.keep 6_GTGAAA_L001_R2_001.trim.fq.keep 6_GTGAAA_L002_R1_001.trim.fq.keep 6_GTGAAA_L002_R2_001.trim.fq.keep
	$(TIME) -o kmergenie.time $(KMERGENIE) dg.reads.txt  --diploid -k 100 -l 20 

dg.kmergenie-pbs:
	qsub -l $(KMERRES) dg.kmergenie.sub

abyss-assemble-31/abyss-lamprey-germline-31-scaffolds.fa: 2_ACAGTG_L001_R1_001.trim.fq.keep 2_ACAGTG_L001_R2_001.trim.fq.keep 2_ACAGTG_L002_R1_001.trim.fq.keep 2_ACAGTG_L002_R2_001.trim.fq.keep 4_GCCAAT_L001_R1_001.trim.fq.keep 4_GCCAAT_L001_R2_001.trim.fq.keep 4_GCCAAT_L002_R1_001.trim.fq.keep 4_GCCAAT_L002_R2_001.trim.fq.keep 5_CTTGTA_L001_R1_001.trim.fq.keep 5_CTTGTA_L001_R2_001.trim.fq.keep 5_CTTGTA_L002_R1_001.trim.fq.keep 5_CTTGTA_L002_R2_001.trim.fq.keep 6_GTGAAA_L001_R1_001.trim.fq.keep 6_GTGAAA_L001_R2_001.trim.fq.keep 6_GTGAAA_L002_R1_001.trim.fq.keep 6_GTGAAA_L002_R2_001.trim.fq.keep
	WORK=abyss-assemble-31
	mkdir -p $(WORK)
	ln -s *.trim.fq.keep $(WORK)/
	abyss-pe v=-v k=31 -C $(WORK) name=abyss-lamprey-germline-31 np=$PBS_NUM_PPN lib="2_ACAGTG_L001_ 2_ACAGTG_L002_ 4_GCCAAT_L001_ 4_GCCAAT_L002_ 5_CTTGTA_L001_ 5_CTTGTA_L002_ 6_GTGAAA_L001_ 6_GTGAAA_L002_" 2_ACAGTG_L001_="2_ACAGTG_L001_R1_001.trim.fq.keep 2_ACAGTG_L001_R2_001.trim.fq.keep" 2_ACAGTG_L002_="2_ACAGTG_L002_R1_001.trim.fq.keep 2_ACAGTG_L002_R2_001.trim.fq.keep" 4_GCCAAT_L001_="4_GCCAAT_L001_R1_001.trim.fq.keep 4_GCCAAT_L001_R2_001.trim.fq.keep" 4_GCCAAT_L002_="4_GCCAAT_L002_R1_001.trim.fq.keep 4_GCCAAT_L002_R2_001.trim.fq.keep" 5_CTTGTA_L001_="5_CTTGTA_L001_R1_001.trim.fq.keep 5_CTTGTA_L001_R2_001.trim.fq.keep" 5_CTTGTA_L002_="5_CTTGTA_L002_R1_001.trim.fq.keep 5_CTTGTA_L002_R2_001.trim.fq.keep" 6_GTGAAA_L001_="6_GTGAAA_L001_R1_001.trim.fq.keep 6_GTGAAA_L001_R2_001.trim.fq.keep" 6_GTGAAA_L002_="6_GTGAAA_L002_R1_001.trim.fq.keep 6_GTGAAA_L002_R2_001.trim.fq.keep" | tee abyss-lamprey-germline-31-`date +%F_%T`.out

abyss-pbs:
	qsub -l $(ABYSSRES) abyss-assemble.sub
