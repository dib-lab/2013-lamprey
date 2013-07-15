
KHMER=~/tools/khmer
SEQCLEAN=seqclean
TRIM=~/tools/Trimmomatic-0.22/trimmomatic-0.22.jar
TRIMRES=nodes=1:ppn=4,walltime=4:00:00,mem=2gb
DG1RES=nodes=1:ppn=4,walltime=4:00:00,mem=8gb
DG2RES=nodes=1:ppn=4,walltime=4:00:00,mem=8gb

614GMAAXX_1_1_pf.trim.fq: 614GMAAXX_1_1_pf.qc.fq.gz
	ln -s 614GMAAXX_1_1_pf.qc.fq.gz /mnt/home/welcherc/lamprey/reads/all/614GMAAXX_1_1_pf.qc.fq.gz
	time -o 614GMAAXX_1_1_pf.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 614GMAAXX_1_1_pf.trimlog -threads 4 614GMAAXX_1_1_pf.qc.fq.gz 614GMAAXX_1_1_pf.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

614GMAAXX_2_1_pf.trim.fq: 614GMAAXX_2_1_pf.qc.fq.gz
	ln -s 614GMAAXX_2_1_pf.qc.fq.gz /mnt/home/welcherc/lamprey/reads/all/614GMAAXX_2_1_pf.qc.fq.gz
	time -o 614GMAAXX_2_1_pf.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 614GMAAXX_2_1_pf.trimlog -threads 4 614GMAAXX_2_1_pf.qc.fq.gz 614GMAAXX_2_1_pf.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

614GMAAXX_3_1_pf.trim.fq: 614GMAAXX_3_1_pf.qc.fq.gz
	ln -s 614GMAAXX_3_1_pf.qc.fq.gz /mnt/home/welcherc/lamprey/reads/all/614GMAAXX_3_1_pf.qc.fq.gz
	time -o 614GMAAXX_3_1_pf.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 614GMAAXX_3_1_pf.trimlog -threads 4 614GMAAXX_3_1_pf.qc.fq.gz 614GMAAXX_3_1_pf.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

614GMAAXX_4_1_pf.trim.fq: 614GMAAXX_4_1_pf.qc.fq.gz
	ln -s 614GMAAXX_4_1_pf.qc.fq.gz /mnt/home/welcherc/lamprey/reads/all/614GMAAXX_4_1_pf.qc.fq.gz
	time -o 614GMAAXX_4_1_pf.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 614GMAAXX_4_1_pf.trimlog -threads 4 614GMAAXX_4_1_pf.qc.fq.gz 614GMAAXX_4_1_pf.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

614GMAAXX_6_1_pf.trim.fq: 614GMAAXX_6_1_pf.qc.fq.gz
	ln -s 614GMAAXX_6_1_pf.qc.fq.gz /mnt/home/welcherc/lamprey/reads/all/614GMAAXX_6_1_pf.qc.fq.gz
	time -o 614GMAAXX_6_1_pf.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 614GMAAXX_6_1_pf.trimlog -threads 4 614GMAAXX_6_1_pf.qc.fq.gz 614GMAAXX_6_1_pf.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

614GMAAXX_7_1_pf.trim.fq: 614GMAAXX_7_1_pf.qc.fq.gz
	ln -s 614GMAAXX_7_1_pf.qc.fq.gz /mnt/home/welcherc/lamprey/reads/all/614GMAAXX_7_1_pf.qc.fq.gz
	time -o 614GMAAXX_7_1_pf.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 614GMAAXX_7_1_pf.trimlog -threads 4 614GMAAXX_7_1_pf.qc.fq.gz 614GMAAXX_7_1_pf.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

614GMAAXX_8_1_pf.trim.fq: 614GMAAXX_8_1_pf.qc.fq.gz
	ln -s 614GMAAXX_8_1_pf.qc.fq.gz /mnt/home/welcherc/lamprey/reads/all/614GMAAXX_8_1_pf.qc.fq.gz
	time -o 614GMAAXX_8_1_pf.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 614GMAAXX_8_1_pf.trimlog -threads 4 614GMAAXX_8_1_pf.qc.fq.gz 614GMAAXX_8_1_pf.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

AB057JABXX_s_7_pe.trim.fq: AB057JABXX_s_7_pe.trim.fq.gz
	ln -s AB057JABXX_s_7_pe.trim.fq.gz /mnt/home/welcherc/lamprey/reads/all/AB057JABXX_s_7_pe.trim.fq.gz
	time -o AB057JABXX_s_7_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog AB057JABXX_s_7_pe.trimlog -threads 4 AB057JABXX_s_7_pe.trim.fq.gz AB057JABXX_s_7_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

AI.trim.fq: AI.fq.gz
	ln -s AI.fq.gz /mnt/home/welcherc/lamprey/reads/all/AI.fq.gz
	time -o AI.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog AI.trimlog -threads 4 AI.fq.gz AI.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

AK.trim.fq: AK.fq.gz
	ln -s AK.fq.gz /mnt/home/welcherc/lamprey/reads/all/AK.fq.gz
	time -o AK.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog AK.trimlog -threads 4 AK.fq.gz AK.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

BB0582ABXX_s_3_pe.trim.fq: BB0582ABXX_s_3_pe.trim.fq.gz
	ln -s BB0582ABXX_s_3_pe.trim.fq.gz /mnt/home/welcherc/lamprey/reads/all/BB0582ABXX_s_3_pe.trim.fq.gz
	time -o BB0582ABXX_s_3_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog BB0582ABXX_s_3_pe.trimlog -threads 4 BB0582ABXX_s_3_pe.trim.fq.gz BB0582ABXX_s_3_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

BB0582ABXX_s_4_pe.trim.fq: BB0582ABXX_s_4_pe.trim.fq.gz
	ln -s BB0582ABXX_s_4_pe.trim.fq.gz /mnt/home/welcherc/lamprey/reads/all/BB0582ABXX_s_4_pe.trim.fq.gz
	time -o BB0582ABXX_s_4_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog BB0582ABXX_s_4_pe.trimlog -threads 4 BB0582ABXX_s_4_pe.trim.fq.gz BB0582ABXX_s_4_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

BB0582ABXX_s_5_pe.trim.fq: BB0582ABXX_s_5_pe.trim.fq.gz
	ln -s BB0582ABXX_s_5_pe.trim.fq.gz /mnt/home/welcherc/lamprey/reads/all/BB0582ABXX_s_5_pe.trim.fq.gz
	time -o BB0582ABXX_s_5_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog BB0582ABXX_s_5_pe.trimlog -threads 4 BB0582ABXX_s_5_pe.trim.fq.gz BB0582ABXX_s_5_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

BB0582ABXX_s_6_pe.trim.fq: BB0582ABXX_s_6_pe.trim.fq.gz
	ln -s BB0582ABXX_s_6_pe.trim.fq.gz /mnt/home/welcherc/lamprey/reads/all/BB0582ABXX_s_6_pe.trim.fq.gz
	time -o BB0582ABXX_s_6_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog BB0582ABXX_s_6_pe.trimlog -threads 4 BB0582ABXX_s_6_pe.trim.fq.gz BB0582ABXX_s_6_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

BB0582ABXX_s_7_pe.trim.fq: BB0582ABXX_s_7_pe.trim.fq.gz
	ln -s BB0582ABXX_s_7_pe.trim.fq.gz /mnt/home/welcherc/lamprey/reads/all/BB0582ABXX_s_7_pe.trim.fq.gz
	time -o BB0582ABXX_s_7_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog BB0582ABXX_s_7_pe.trimlog -threads 4 BB0582ABXX_s_7_pe.trim.fq.gz BB0582ABXX_s_7_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

BB0582ABXX_s_8_pe.trim.fq: BB0582ABXX_s_8_pe.trim.fq.gz
	ln -s BB0582ABXX_s_8_pe.trim.fq.gz /mnt/home/welcherc/lamprey/reads/all/BB0582ABXX_s_8_pe.trim.fq.gz
	time -o BB0582ABXX_s_8_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog BB0582ABXX_s_8_pe.trimlog -threads 4 BB0582ABXX_s_8_pe.trim.fq.gz BB0582ABXX_s_8_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

brain_paired.trim.fq: brain_paired.qc.fq.gz
	ln -s brain_paired.qc.fq.gz /mnt/home/welcherc/lamprey/reads/all/brain_paired.qc.fq.gz
	time -o brain_paired.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog brain_paired.trimlog -threads 4 brain_paired.qc.fq.gz brain_paired.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

FW1-Gill.trim.fq: FW1-Gill.fq.gz
	ln -s FW1-Gill.fq.gz /mnt/home/welcherc/lamprey/reads/all/FW1-Gill.fq.gz
	time -o FW1-Gill.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog FW1-Gill.trimlog -threads 4 FW1-Gill.fq.gz FW1-Gill.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

FW1-Intestine.trim.fq: FW1-Intestine.fq.gz
	ln -s FW1-Intestine.fq.gz /mnt/home/welcherc/lamprey/reads/all/FW1-Intestine.fq.gz
	time -o FW1-Intestine.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog FW1-Intestine.trimlog -threads 4 FW1-Intestine.fq.gz FW1-Intestine.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

FW1-Kidney.trim.fq: FW1-Kidney.fq.gz
	ln -s FW1-Kidney.fq.gz /mnt/home/welcherc/lamprey/reads/all/FW1-Kidney.fq.gz
	time -o FW1-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog FW1-Kidney.trimlog -threads 4 FW1-Kidney.fq.gz FW1-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

G81-a.trim.fq: G81-a.fq.gz
	ln -s G81-a.fq.gz /mnt/home/welcherc/lamprey/reads/all/G81-a.fq.gz
	time -o G81-a.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog G81-a.trimlog -threads 4 G81-a.fq.gz G81-a.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

G81-b.trim.fq: G81-b.fq.gz
	ln -s G81-b.fq.gz /mnt/home/welcherc/lamprey/reads/all/G81-b.fq.gz
	time -o G81-b.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog G81-b.trimlog -threads 4 G81-b.fq.gz G81-b.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

G82-a.trim.fq: G82-a.fq.gz
	ln -s G82-a.fq.gz /mnt/home/welcherc/lamprey/reads/all/G82-a.fq.gz
	time -o G82-a.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog G82-a.trimlog -threads 4 G82-a.fq.gz G82-a.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

G82-b.trim.fq: G82-b.fq.gz
	ln -s G82-b.fq.gz /mnt/home/welcherc/lamprey/reads/all/G82-b.fq.gz
	time -o G82-b.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog G82-b.trimlog -threads 4 G82-b.fq.gz G82-b.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

G83-a.trim.fq: G83-a.fq.gz
	ln -s G83-a.fq.gz /mnt/home/welcherc/lamprey/reads/all/G83-a.fq.gz
	time -o G83-a.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog G83-a.trimlog -threads 4 G83-a.fq.gz G83-a.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

G83-b.trim.fq: G83-b.fq.gz
	ln -s G83-b.fq.gz /mnt/home/welcherc/lamprey/reads/all/G83-b.fq.gz
	time -o G83-b.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog G83-b.trimlog -threads 4 G83-b.fq.gz G83-b.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

Gill.trim.fq: Gill.fq.gz
	ln -s Gill.fq.gz /mnt/home/welcherc/lamprey/reads/all/Gill.fq.gz
	time -o Gill.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Gill.trimlog -threads 4 Gill.fq.gz Gill.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

JV1_liver.trim.fq: JV1_liver.fq.gz
	ln -s JV1_liver.fq.gz /mnt/home/welcherc/lamprey/reads/all/JV1_liver.fq.gz
	time -o JV1_liver.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog JV1_liver.trimlog -threads 4 JV1_liver.fq.gz JV1_liver.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

JVa_liver.trim.fq: JVa_liver.fq.gz
	ln -s JVa_liver.fq.gz /mnt/home/welcherc/lamprey/reads/all/JVa_liver.fq.gz
	time -o JVa_liver.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog JVa_liver.trimlog -threads 4 JVa_liver.fq.gz JVa_liver.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

JVb_liver.trim.fq: JVb_liver.fq.gz
	ln -s JVb_liver.fq.gz /mnt/home/welcherc/lamprey/reads/all/JVb_liver.fq.gz
	time -o JVb_liver.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog JVb_liver.trimlog -threads 4 JVb_liver.fq.gz JVb_liver.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

JV_intest.trim.fq: JV_intest.fq.gz
	ln -s JV_intest.fq.gz /mnt/home/welcherc/lamprey/reads/all/JV_intest.fq.gz
	time -o JV_intest.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog JV_intest.trimlog -threads 4 JV_intest.fq.gz JV_intest.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

Kidney.trim.fq: Kidney.fq.gz
	ln -s Kidney.fq.gz /mnt/home/welcherc/lamprey/reads/all/Kidney.fq.gz
	time -o Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Kidney.trimlog -threads 4 Kidney.fq.gz Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

L81-a.trim.fq: L81-a.fq.gz
	ln -s L81-a.fq.gz /mnt/home/welcherc/lamprey/reads/all/L81-a.fq.gz
	time -o L81-a.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog L81-a.trimlog -threads 4 L81-a.fq.gz L81-a.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

L81-b.trim.fq: L81-b.fq.gz
	ln -s L81-b.fq.gz /mnt/home/welcherc/lamprey/reads/all/L81-b.fq.gz
	time -o L81-b.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog L81-b.trimlog -threads 4 L81-b.fq.gz L81-b.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

L82-a.trim.fq: L82-a.fq.gz
	ln -s L82-a.fq.gz /mnt/home/welcherc/lamprey/reads/all/L82-a.fq.gz
	time -o L82-a.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog L82-a.trimlog -threads 4 L82-a.fq.gz L82-a.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

L82-b.trim.fq: L82-b.fq.gz
	ln -s L82-b.fq.gz /mnt/home/welcherc/lamprey/reads/all/L82-b.fq.gz
	time -o L82-b.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog L82-b.trimlog -threads 4 L82-b.fq.gz L82-b.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

L83-1-a.trim.fq: L83-1-a.fq.gz
	ln -s L83-1-a.fq.gz /mnt/home/welcherc/lamprey/reads/all/L83-1-a.fq.gz
	time -o L83-1-a.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog L83-1-a.trimlog -threads 4 L83-1-a.fq.gz L83-1-a.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

L83-1-b.trim.fq: L83-1-b.fq.gz
	ln -s L83-1-b.fq.gz /mnt/home/welcherc/lamprey/reads/all/L83-1-b.fq.gz
	time -o L83-1-b.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog L83-1-b.trimlog -threads 4 L83-1-b.fq.gz L83-1-b.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

L83-a.trim.fq: L83-a.fq.gz
	ln -s L83-a.fq.gz /mnt/home/welcherc/lamprey/reads/all/L83-a.fq.gz
	time -o L83-a.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog L83-a.trimlog -threads 4 L83-a.fq.gz L83-a.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

L83-b.trim.fq: L83-b.fq.gz
	ln -s L83-b.fq.gz /mnt/home/welcherc/lamprey/reads/all/L83-b.fq.gz
	time -o L83-b.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog L83-b.trimlog -threads 4 L83-b.fq.gz L83-b.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

LI.trim.fq: LI.fq.gz
	ln -s LI.fq.gz /mnt/home/welcherc/lamprey/reads/all/LI.fq.gz
	time -o LI.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog LI.trimlog -threads 4 LI.fq.gz LI.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

Lips.trim.fq: Lips.fq.gz
	ln -s Lips.fq.gz /mnt/home/welcherc/lamprey/reads/all/Lips.fq.gz
	time -o Lips.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Lips.trimlog -threads 4 Lips.fq.gz Lips.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

LK.trim.fq: LK.fq.gz
	ln -s LK.fq.gz /mnt/home/welcherc/lamprey/reads/all/LK.fq.gz
	time -o LK.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog LK.trimlog -threads 4 LK.fq.gz LK.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

M0_liver.trim.fq: M0_liver.fq.gz
	ln -s M0_liver.fq.gz /mnt/home/welcherc/lamprey/reads/all/M0_liver.fq.gz
	time -o M0_liver.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog M0_liver.trimlog -threads 4 M0_liver.fq.gz M0_liver.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

M2_liver.trim.fq: M2_liver.fq.gz
	ln -s M2_liver.fq.gz /mnt/home/welcherc/lamprey/reads/all/M2_liver.fq.gz
	time -o M2_liver.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog M2_liver.trimlog -threads 4 M2_liver.fq.gz M2_liver.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

M5_liver.trim.fq: M5_liver.fq.gz
	ln -s M5_liver.fq.gz /mnt/home/welcherc/lamprey/reads/all/M5_liver.fq.gz
	time -o M5_liver.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog M5_liver.trimlog -threads 4 M5_liver.fq.gz M5_liver.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

Monocytes.trim.fq: Monocytes.fq.gz
	ln -s Monocytes.fq.gz /mnt/home/welcherc/lamprey/reads/all/Monocytes.fq.gz
	time -o Monocytes.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Monocytes.trimlog -threads 4 Monocytes.fq.gz Monocytes.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

mssm_s_1_pair_pe.trim.fq: mssm_s_1_pair_pe.gz
	ln -s mssm_s_1_pair_pe.gz /mnt/home/welcherc/lamprey/reads/all/mssm_s_1_pair_pe.gz
	time -o mssm_s_1_pair_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog mssm_s_1_pair_pe.trimlog -threads 4 mssm_s_1_pair_pe.gz mssm_s_1_pair_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

mssm_s_1_pair_se.trim.fq: mssm_s_1_pair_se.gz
	ln -s mssm_s_1_pair_se.gz /mnt/home/welcherc/lamprey/reads/all/mssm_s_1_pair_se.gz
	time -o mssm_s_1_pair_se.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog mssm_s_1_pair_se.trimlog -threads 4 mssm_s_1_pair_se.gz mssm_s_1_pair_se.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

mssm_s_2_pair_pe.trim.fq: mssm_s_2_pair_pe.gz
	ln -s mssm_s_2_pair_pe.gz /mnt/home/welcherc/lamprey/reads/all/mssm_s_2_pair_pe.gz
	time -o mssm_s_2_pair_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog mssm_s_2_pair_pe.trimlog -threads 4 mssm_s_2_pair_pe.gz mssm_s_2_pair_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

mssm_s_2_pair_se.trim.fq: mssm_s_2_pair_se.gz
	ln -s mssm_s_2_pair_se.gz /mnt/home/welcherc/lamprey/reads/all/mssm_s_2_pair_se.gz
	time -o mssm_s_2_pair_se.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog mssm_s_2_pair_se.trimlog -threads 4 mssm_s_2_pair_se.gz mssm_s_2_pair_se.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

mssm_s_3_pair_pe.trim.fq: mssm_s_3_pair_pe.gz
	ln -s mssm_s_3_pair_pe.gz /mnt/home/welcherc/lamprey/reads/all/mssm_s_3_pair_pe.gz
	time -o mssm_s_3_pair_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog mssm_s_3_pair_pe.trimlog -threads 4 mssm_s_3_pair_pe.gz mssm_s_3_pair_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

mssm_s_3_pair_se.trim.fq: mssm_s_3_pair_se.gz
	ln -s mssm_s_3_pair_se.gz /mnt/home/welcherc/lamprey/reads/all/mssm_s_3_pair_se.gz
	time -o mssm_s_3_pair_se.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog mssm_s_3_pair_se.trimlog -threads 4 mssm_s_3_pair_se.gz mssm_s_3_pair_se.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

mssm_s_5_pair_pe.trim.fq: mssm_s_5_pair_pe.gz
	ln -s mssm_s_5_pair_pe.gz /mnt/home/welcherc/lamprey/reads/all/mssm_s_5_pair_pe.gz
	time -o mssm_s_5_pair_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog mssm_s_5_pair_pe.trimlog -threads 4 mssm_s_5_pair_pe.gz mssm_s_5_pair_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

mssm_s_5_pair_se.trim.fq: mssm_s_5_pair_se.gz
	ln -s mssm_s_5_pair_se.gz /mnt/home/welcherc/lamprey/reads/all/mssm_s_5_pair_se.gz
	time -o mssm_s_5_pair_se.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog mssm_s_5_pair_se.trimlog -threads 4 mssm_s_5_pair_se.gz mssm_s_5_pair_se.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

mssm_s_6_pair_pe.trim.fq: mssm_s_6_pair_pe.gz
	ln -s mssm_s_6_pair_pe.gz /mnt/home/welcherc/lamprey/reads/all/mssm_s_6_pair_pe.gz
	time -o mssm_s_6_pair_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog mssm_s_6_pair_pe.trimlog -threads 4 mssm_s_6_pair_pe.gz mssm_s_6_pair_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

mssm_s_6_pair_se.trim.fq: mssm_s_6_pair_se.gz
	ln -s mssm_s_6_pair_se.gz /mnt/home/welcherc/lamprey/reads/all/mssm_s_6_pair_se.gz
	time -o mssm_s_6_pair_se.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog mssm_s_6_pair_se.trimlog -threads 4 mssm_s_6_pair_se.gz mssm_s_6_pair_se.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

mssm_s_7_pair_pe.trim.fq: mssm_s_7_pair_pe.gz
	ln -s mssm_s_7_pair_pe.gz /mnt/home/welcherc/lamprey/reads/all/mssm_s_7_pair_pe.gz
	time -o mssm_s_7_pair_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog mssm_s_7_pair_pe.trimlog -threads 4 mssm_s_7_pair_pe.gz mssm_s_7_pair_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

mssm_s_7_pair_se.trim.fq: mssm_s_7_pair_se.gz
	ln -s mssm_s_7_pair_se.gz /mnt/home/welcherc/lamprey/reads/all/mssm_s_7_pair_se.gz
	time -o mssm_s_7_pair_se.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog mssm_s_7_pair_se.trimlog -threads 4 mssm_s_7_pair_se.gz mssm_s_7_pair_se.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

mssm_single_s1.trim.fq: mssm_single_s1.trim.fq.gz
	ln -s mssm_single_s1.trim.fq.gz /mnt/home/welcherc/lamprey/reads/all/mssm_single_s1.trim.fq.gz
	time -o mssm_single_s1.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog mssm_single_s1.trimlog -threads 4 mssm_single_s1.trim.fq.gz mssm_single_s1.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

mssm_single_s2.trim.fq: mssm_single_s2.trim.fq.gz
	ln -s mssm_single_s2.trim.fq.gz /mnt/home/welcherc/lamprey/reads/all/mssm_single_s2.trim.fq.gz
	time -o mssm_single_s2.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog mssm_single_s2.trimlog -threads 4 mssm_single_s2.trim.fq.gz mssm_single_s2.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

mssm_single_s3.trim.fq: mssm_single_s3.trim.fq.gz
	ln -s mssm_single_s3.trim.fq.gz /mnt/home/welcherc/lamprey/reads/all/mssm_single_s3.trim.fq.gz
	time -o mssm_single_s3.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog mssm_single_s3.trimlog -threads 4 mssm_single_s3.trim.fq.gz mssm_single_s3.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

mssm_single_s5.trim.fq: mssm_single_s5.trim.fq.gz
	ln -s mssm_single_s5.trim.fq.gz /mnt/home/welcherc/lamprey/reads/all/mssm_single_s5.trim.fq.gz
	time -o mssm_single_s5.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog mssm_single_s5.trimlog -threads 4 mssm_single_s5.trim.fq.gz mssm_single_s5.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

mssm_single_s6.trim.fq: mssm_single_s6.trim.fq.gz
	ln -s mssm_single_s6.trim.fq.gz /mnt/home/welcherc/lamprey/reads/all/mssm_single_s6.trim.fq.gz
	time -o mssm_single_s6.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog mssm_single_s6.trimlog -threads 4 mssm_single_s6.trim.fq.gz mssm_single_s6.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

mssm_single_s7.trim.fq: mssm_single_s7.trim.fq.gz
	ln -s mssm_single_s7.trim.fq.gz /mnt/home/welcherc/lamprey/reads/all/mssm_single_s7.trim.fq.gz
	time -o mssm_single_s7.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog mssm_single_s7.trimlog -threads 4 mssm_single_s7.trim.fq.gz mssm_single_s7.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

Muscle.trim.fq: Muscle.fq.gz
	ln -s Muscle.fq.gz /mnt/home/welcherc/lamprey/reads/all/Muscle.fq.gz
	time -o Muscle.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Muscle.trimlog -threads 4 Muscle.fq.gz Muscle.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

Neutrophils.trim.fq: Neutrophils.fq.gz
	ln -s Neutrophils.fq.gz /mnt/home/welcherc/lamprey/reads/all/Neutrophils.fq.gz
	time -o Neutrophils.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Neutrophils.trimlog -threads 4 Neutrophils.fq.gz Neutrophils.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

OFHS.trim.fq: OFHS.fq.gz
	ln -s OFHS.fq.gz /mnt/home/welcherc/lamprey/reads/all/OFHS.fq.gz
	time -o OFHS.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog OFHS.trimlog -threads 4 OFHS.fq.gz OFHS.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

POF-eye.trim.fq: POF-eye.fq.gz
	ln -s POF-eye.fq.gz /mnt/home/welcherc/lamprey/reads/all/POF-eye.fq.gz
	time -o POF-eye.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog POF-eye.trimlog -threads 4 POF-eye.fq.gz POF-eye.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

POF-ts.trim.fq: POF-ts.fq.gz
	ln -s POF-ts.fq.gz /mnt/home/welcherc/lamprey/reads/all/POF-ts.fq.gz
	time -o POF-ts.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog POF-ts.trimlog -threads 4 POF-ts.fq.gz POF-ts.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

PSMGill.trim.fq: PSMGill.fq.gz
	ln -s PSMGill.fq.gz /mnt/home/welcherc/lamprey/reads/all/PSMGill.fq.gz
	time -o PSMGill.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog PSMGill.trimlog -threads 4 PSMGill.fq.gz PSMGill.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

Rope.trim.fq: Rope.fq.gz
	ln -s Rope.fq.gz /mnt/home/welcherc/lamprey/reads/all/Rope.fq.gz
	time -o Rope.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Rope.trimlog -threads 4 Rope.fq.gz Rope.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

SMGill.trim.fq: SMGill.fq.gz
	ln -s SMGill.fq.gz /mnt/home/welcherc/lamprey/reads/all/SMGill.fq.gz
	time -o SMGill.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog SMGill.trimlog -threads 4 SMGill.fq.gz SMGill.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

SMHS.trim.fq: SMHS.fq.gz
	ln -s SMHS.fq.gz /mnt/home/welcherc/lamprey/reads/all/SMHS.fq.gz
	time -o SMHS.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog SMHS.trimlog -threads 4 SMHS.fq.gz SMHS.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

SN.trim.fq: SN.fq.gz
	ln -s SN.fq.gz /mnt/home/welcherc/lamprey/reads/all/SN.fq.gz
	time -o SN.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog SN.trimlog -threads 4 SN.fq.gz SN.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

SPDI.trim.fq: SPDI.fq.gz
	ln -s SPDI.fq.gz /mnt/home/welcherc/lamprey/reads/all/SPDI.fq.gz
	time -o SPDI.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog SPDI.trimlog -threads 4 SPDI.fq.gz SPDI.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

spinal_cord_paired.trim.fq: spinal_cord_paired.qc.fq.gz
	ln -s spinal_cord_paired.qc.fq.gz /mnt/home/welcherc/lamprey/reads/all/spinal_cord_paired.qc.fq.gz
	time -o spinal_cord_paired.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog spinal_cord_paired.trimlog -threads 4 spinal_cord_paired.qc.fq.gz spinal_cord_paired.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

SPK.trim.fq: SPK.fq.gz
	ln -s SPK.fq.gz /mnt/home/welcherc/lamprey/reads/all/SPK.fq.gz
	time -o SPK.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog SPK.trimlog -threads 4 SPK.fq.gz SPK.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

SPPI.trim.fq: SPPI.fq.gz
	ln -s SPPI.fq.gz /mnt/home/welcherc/lamprey/reads/all/SPPI.fq.gz
	time -o SPPI.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog SPPI.trimlog -threads 4 SPPI.fq.gz SPPI.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

Stg-1-Intestine.trim.fq: Stg-1-Intestine.fq.gz
	ln -s Stg-1-Intestine.fq.gz /mnt/home/welcherc/lamprey/reads/all/Stg-1-Intestine.fq.gz
	time -o Stg-1-Intestine.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-1-Intestine.trimlog -threads 4 Stg-1-Intestine.fq.gz Stg-1-Intestine.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

Stg-1-Kidney.trim.fq: Stg-1-Kidney.fq.gz
	ln -s Stg-1-Kidney.fq.gz /mnt/home/welcherc/lamprey/reads/all/Stg-1-Kidney.fq.gz
	time -o Stg-1-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-1-Kidney.trimlog -threads 4 Stg-1-Kidney.fq.gz Stg-1-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

Stg-2-Intestine.trim.fq: Stg-2-Intestine.fq.gz
	ln -s Stg-2-Intestine.fq.gz /mnt/home/welcherc/lamprey/reads/all/Stg-2-Intestine.fq.gz
	time -o Stg-2-Intestine.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-2-Intestine.trimlog -threads 4 Stg-2-Intestine.fq.gz Stg-2-Intestine.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

Stg-2-Kidney.trim.fq: Stg-2-Kidney.fq.gz
	ln -s Stg-2-Kidney.fq.gz /mnt/home/welcherc/lamprey/reads/all/Stg-2-Kidney.fq.gz
	time -o Stg-2-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-2-Kidney.trimlog -threads 4 Stg-2-Kidney.fq.gz Stg-2-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

Stg-3-Intestine.trim.fq: Stg-3-Intestine.fq.gz
	ln -s Stg-3-Intestine.fq.gz /mnt/home/welcherc/lamprey/reads/all/Stg-3-Intestine.fq.gz
	time -o Stg-3-Intestine.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-3-Intestine.trimlog -threads 4 Stg-3-Intestine.fq.gz Stg-3-Intestine.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

Stg-3-Kidney.trim.fq: Stg-3-Kidney.fq.gz
	ln -s Stg-3-Kidney.fq.gz /mnt/home/welcherc/lamprey/reads/all/Stg-3-Kidney.fq.gz
	time -o Stg-3-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-3-Kidney.trimlog -threads 4 Stg-3-Kidney.fq.gz Stg-3-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

Stg-4-Intestine.trim.fq: Stg-4-Intestine.fq.gz
	ln -s Stg-4-Intestine.fq.gz /mnt/home/welcherc/lamprey/reads/all/Stg-4-Intestine.fq.gz
	time -o Stg-4-Intestine.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-4-Intestine.trimlog -threads 4 Stg-4-Intestine.fq.gz Stg-4-Intestine.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

Stg-4-Kidney.trim.fq: Stg-4-Kidney.fq.gz
	ln -s Stg-4-Kidney.fq.gz /mnt/home/welcherc/lamprey/reads/all/Stg-4-Kidney.fq.gz
	time -o Stg-4-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-4-Kidney.trimlog -threads 4 Stg-4-Kidney.fq.gz Stg-4-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

Stg-5-Intestine.trim.fq: Stg-5-Intestine.fq.gz
	ln -s Stg-5-Intestine.fq.gz /mnt/home/welcherc/lamprey/reads/all/Stg-5-Intestine.fq.gz
	time -o Stg-5-Intestine.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-5-Intestine.trimlog -threads 4 Stg-5-Intestine.fq.gz Stg-5-Intestine.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

Stg-5-Kidney.trim.fq: Stg-5-Kidney.fq.gz
	ln -s Stg-5-Kidney.fq.gz /mnt/home/welcherc/lamprey/reads/all/Stg-5-Kidney.fq.gz
	time -o Stg-5-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-5-Kidney.trimlog -threads 4 Stg-5-Kidney.fq.gz Stg-5-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

Stg-6-Intestine.trim.fq: Stg-6-Intestine.fq.gz
	ln -s Stg-6-Intestine.fq.gz /mnt/home/welcherc/lamprey/reads/all/Stg-6-Intestine.fq.gz
	time -o Stg-6-Intestine.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-6-Intestine.trimlog -threads 4 Stg-6-Intestine.fq.gz Stg-6-Intestine.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

Stg-6-Kidney.trim.fq: Stg-6-Kidney.fq.gz
	ln -s Stg-6-Kidney.fq.gz /mnt/home/welcherc/lamprey/reads/all/Stg-6-Kidney.fq.gz
	time -o Stg-6-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-6-Kidney.trimlog -threads 4 Stg-6-Kidney.fq.gz Stg-6-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

Stg-7-Intestine.trim.fq: Stg-7-Intestine.fq.gz
	ln -s Stg-7-Intestine.fq.gz /mnt/home/welcherc/lamprey/reads/all/Stg-7-Intestine.fq.gz
	time -o Stg-7-Intestine.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-7-Intestine.trimlog -threads 4 Stg-7-Intestine.fq.gz Stg-7-Intestine.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

Stg-7-Kidney.trim.fq: Stg-7-Kidney.fq.gz
	ln -s Stg-7-Kidney.fq.gz /mnt/home/welcherc/lamprey/reads/all/Stg-7-Kidney.fq.gz
	time -o Stg-7-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-7-Kidney.trimlog -threads 4 Stg-7-Kidney.fq.gz Stg-7-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

SW2-Gill.trim.fq: SW2-Gill.fq.gz
	ln -s SW2-Gill.fq.gz /mnt/home/welcherc/lamprey/reads/all/SW2-Gill.fq.gz
	time -o SW2-Gill.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog SW2-Gill.trimlog -threads 4 SW2-Gill.fq.gz SW2-Gill.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

SW2-Intestine-a.trim.fq: SW2-Intestine-a.fq.gz
	ln -s SW2-Intestine-a.fq.gz /mnt/home/welcherc/lamprey/reads/all/SW2-Intestine-a.fq.gz
	time -o SW2-Intestine-a.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog SW2-Intestine-a.trimlog -threads 4 SW2-Intestine-a.fq.gz SW2-Intestine-a.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

SW2-Intestine-b.trim.fq: SW2-Intestine-b.fq.gz
	ln -s SW2-Intestine-b.fq.gz /mnt/home/welcherc/lamprey/reads/all/SW2-Intestine-b.fq.gz
	time -o SW2-Intestine-b.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog SW2-Intestine-b.trimlog -threads 4 SW2-Intestine-b.fq.gz SW2-Intestine-b.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

SW2-Kidney.trim.fq: SW2-Kidney.fq.gz
	ln -s SW2-Kidney.fq.gz /mnt/home/welcherc/lamprey/reads/all/SW2-Kidney.fq.gz
	time -o SW2-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog SW2-Kidney.trimlog -threads 4 SW2-Kidney.fq.gz SW2-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 HEADCROP:9 SLIDINGWINDOW:4:15 MINLEN:31

614GMAAXX_1_1_pf.trim.fq-pbs: 614GMAAXX_1_1_pf.trim.fq.sub
	qsub -l $(TRIMRES) 614GMAAXX_1_1_pf.trim.fq.sub
614GMAAXX_2_1_pf.trim.fq-pbs: 614GMAAXX_2_1_pf.trim.fq.sub
	qsub -l $(TRIMRES) 614GMAAXX_2_1_pf.trim.fq.sub
614GMAAXX_3_1_pf.trim.fq-pbs: 614GMAAXX_3_1_pf.trim.fq.sub
	qsub -l $(TRIMRES) 614GMAAXX_3_1_pf.trim.fq.sub
614GMAAXX_4_1_pf.trim.fq-pbs: 614GMAAXX_4_1_pf.trim.fq.sub
	qsub -l $(TRIMRES) 614GMAAXX_4_1_pf.trim.fq.sub
614GMAAXX_6_1_pf.trim.fq-pbs: 614GMAAXX_6_1_pf.trim.fq.sub
	qsub -l $(TRIMRES) 614GMAAXX_6_1_pf.trim.fq.sub
614GMAAXX_7_1_pf.trim.fq-pbs: 614GMAAXX_7_1_pf.trim.fq.sub
	qsub -l $(TRIMRES) 614GMAAXX_7_1_pf.trim.fq.sub
614GMAAXX_8_1_pf.trim.fq-pbs: 614GMAAXX_8_1_pf.trim.fq.sub
	qsub -l $(TRIMRES) 614GMAAXX_8_1_pf.trim.fq.sub
AB057JABXX_s_7_pe.trim.fq-pbs: AB057JABXX_s_7_pe.trim.fq.sub
	qsub -l $(TRIMRES) AB057JABXX_s_7_pe.trim.fq.sub
AI.trim.fq-pbs: AI.trim.fq.sub
	qsub -l $(TRIMRES) AI.trim.fq.sub
AK.trim.fq-pbs: AK.trim.fq.sub
	qsub -l $(TRIMRES) AK.trim.fq.sub
BB0582ABXX_s_3_pe.trim.fq-pbs: BB0582ABXX_s_3_pe.trim.fq.sub
	qsub -l $(TRIMRES) BB0582ABXX_s_3_pe.trim.fq.sub
BB0582ABXX_s_4_pe.trim.fq-pbs: BB0582ABXX_s_4_pe.trim.fq.sub
	qsub -l $(TRIMRES) BB0582ABXX_s_4_pe.trim.fq.sub
BB0582ABXX_s_5_pe.trim.fq-pbs: BB0582ABXX_s_5_pe.trim.fq.sub
	qsub -l $(TRIMRES) BB0582ABXX_s_5_pe.trim.fq.sub
BB0582ABXX_s_6_pe.trim.fq-pbs: BB0582ABXX_s_6_pe.trim.fq.sub
	qsub -l $(TRIMRES) BB0582ABXX_s_6_pe.trim.fq.sub
BB0582ABXX_s_7_pe.trim.fq-pbs: BB0582ABXX_s_7_pe.trim.fq.sub
	qsub -l $(TRIMRES) BB0582ABXX_s_7_pe.trim.fq.sub
BB0582ABXX_s_8_pe.trim.fq-pbs: BB0582ABXX_s_8_pe.trim.fq.sub
	qsub -l $(TRIMRES) BB0582ABXX_s_8_pe.trim.fq.sub
brain_paired.trim.fq-pbs: brain_paired.trim.fq.sub
	qsub -l $(TRIMRES) brain_paired.trim.fq.sub
FW1-Gill.trim.fq-pbs: FW1-Gill.trim.fq.sub
	qsub -l $(TRIMRES) FW1-Gill.trim.fq.sub
FW1-Intestine.trim.fq-pbs: FW1-Intestine.trim.fq.sub
	qsub -l $(TRIMRES) FW1-Intestine.trim.fq.sub
FW1-Kidney.trim.fq-pbs: FW1-Kidney.trim.fq.sub
	qsub -l $(TRIMRES) FW1-Kidney.trim.fq.sub
G81-a.trim.fq-pbs: G81-a.trim.fq.sub
	qsub -l $(TRIMRES) G81-a.trim.fq.sub
G81-b.trim.fq-pbs: G81-b.trim.fq.sub
	qsub -l $(TRIMRES) G81-b.trim.fq.sub
G82-a.trim.fq-pbs: G82-a.trim.fq.sub
	qsub -l $(TRIMRES) G82-a.trim.fq.sub
G82-b.trim.fq-pbs: G82-b.trim.fq.sub
	qsub -l $(TRIMRES) G82-b.trim.fq.sub
G83-a.trim.fq-pbs: G83-a.trim.fq.sub
	qsub -l $(TRIMRES) G83-a.trim.fq.sub
G83-b.trim.fq-pbs: G83-b.trim.fq.sub
	qsub -l $(TRIMRES) G83-b.trim.fq.sub
Gill.trim.fq-pbs: Gill.trim.fq.sub
	qsub -l $(TRIMRES) Gill.trim.fq.sub
JV1_liver.trim.fq-pbs: JV1_liver.trim.fq.sub
	qsub -l $(TRIMRES) JV1_liver.trim.fq.sub
JVa_liver.trim.fq-pbs: JVa_liver.trim.fq.sub
	qsub -l $(TRIMRES) JVa_liver.trim.fq.sub
JVb_liver.trim.fq-pbs: JVb_liver.trim.fq.sub
	qsub -l $(TRIMRES) JVb_liver.trim.fq.sub
JV_intest.trim.fq-pbs: JV_intest.trim.fq.sub
	qsub -l $(TRIMRES) JV_intest.trim.fq.sub
Kidney.trim.fq-pbs: Kidney.trim.fq.sub
	qsub -l $(TRIMRES) Kidney.trim.fq.sub
L81-a.trim.fq-pbs: L81-a.trim.fq.sub
	qsub -l $(TRIMRES) L81-a.trim.fq.sub
L81-b.trim.fq-pbs: L81-b.trim.fq.sub
	qsub -l $(TRIMRES) L81-b.trim.fq.sub
L82-a.trim.fq-pbs: L82-a.trim.fq.sub
	qsub -l $(TRIMRES) L82-a.trim.fq.sub
L82-b.trim.fq-pbs: L82-b.trim.fq.sub
	qsub -l $(TRIMRES) L82-b.trim.fq.sub
L83-1-a.trim.fq-pbs: L83-1-a.trim.fq.sub
	qsub -l $(TRIMRES) L83-1-a.trim.fq.sub
L83-1-b.trim.fq-pbs: L83-1-b.trim.fq.sub
	qsub -l $(TRIMRES) L83-1-b.trim.fq.sub
L83-a.trim.fq-pbs: L83-a.trim.fq.sub
	qsub -l $(TRIMRES) L83-a.trim.fq.sub
L83-b.trim.fq-pbs: L83-b.trim.fq.sub
	qsub -l $(TRIMRES) L83-b.trim.fq.sub
LI.trim.fq-pbs: LI.trim.fq.sub
	qsub -l $(TRIMRES) LI.trim.fq.sub
Lips.trim.fq-pbs: Lips.trim.fq.sub
	qsub -l $(TRIMRES) Lips.trim.fq.sub
LK.trim.fq-pbs: LK.trim.fq.sub
	qsub -l $(TRIMRES) LK.trim.fq.sub
M0_liver.trim.fq-pbs: M0_liver.trim.fq.sub
	qsub -l $(TRIMRES) M0_liver.trim.fq.sub
M2_liver.trim.fq-pbs: M2_liver.trim.fq.sub
	qsub -l $(TRIMRES) M2_liver.trim.fq.sub
M5_liver.trim.fq-pbs: M5_liver.trim.fq.sub
	qsub -l $(TRIMRES) M5_liver.trim.fq.sub
Monocytes.trim.fq-pbs: Monocytes.trim.fq.sub
	qsub -l $(TRIMRES) Monocytes.trim.fq.sub
mssm_s_1_pair_pe.trim.fq-pbs: mssm_s_1_pair_pe.trim.fq.sub
	qsub -l $(TRIMRES) mssm_s_1_pair_pe.trim.fq.sub
mssm_s_1_pair_se.trim.fq-pbs: mssm_s_1_pair_se.trim.fq.sub
	qsub -l $(TRIMRES) mssm_s_1_pair_se.trim.fq.sub
mssm_s_2_pair_pe.trim.fq-pbs: mssm_s_2_pair_pe.trim.fq.sub
	qsub -l $(TRIMRES) mssm_s_2_pair_pe.trim.fq.sub
mssm_s_2_pair_se.trim.fq-pbs: mssm_s_2_pair_se.trim.fq.sub
	qsub -l $(TRIMRES) mssm_s_2_pair_se.trim.fq.sub
mssm_s_3_pair_pe.trim.fq-pbs: mssm_s_3_pair_pe.trim.fq.sub
	qsub -l $(TRIMRES) mssm_s_3_pair_pe.trim.fq.sub
mssm_s_3_pair_se.trim.fq-pbs: mssm_s_3_pair_se.trim.fq.sub
	qsub -l $(TRIMRES) mssm_s_3_pair_se.trim.fq.sub
mssm_s_5_pair_pe.trim.fq-pbs: mssm_s_5_pair_pe.trim.fq.sub
	qsub -l $(TRIMRES) mssm_s_5_pair_pe.trim.fq.sub
mssm_s_5_pair_se.trim.fq-pbs: mssm_s_5_pair_se.trim.fq.sub
	qsub -l $(TRIMRES) mssm_s_5_pair_se.trim.fq.sub
mssm_s_6_pair_pe.trim.fq-pbs: mssm_s_6_pair_pe.trim.fq.sub
	qsub -l $(TRIMRES) mssm_s_6_pair_pe.trim.fq.sub
mssm_s_6_pair_se.trim.fq-pbs: mssm_s_6_pair_se.trim.fq.sub
	qsub -l $(TRIMRES) mssm_s_6_pair_se.trim.fq.sub
mssm_s_7_pair_pe.trim.fq-pbs: mssm_s_7_pair_pe.trim.fq.sub
	qsub -l $(TRIMRES) mssm_s_7_pair_pe.trim.fq.sub
mssm_s_7_pair_se.trim.fq-pbs: mssm_s_7_pair_se.trim.fq.sub
	qsub -l $(TRIMRES) mssm_s_7_pair_se.trim.fq.sub
mssm_single_s1.trim.fq-pbs: mssm_single_s1.trim.fq.sub
	qsub -l $(TRIMRES) mssm_single_s1.trim.fq.sub
mssm_single_s2.trim.fq-pbs: mssm_single_s2.trim.fq.sub
	qsub -l $(TRIMRES) mssm_single_s2.trim.fq.sub
mssm_single_s3.trim.fq-pbs: mssm_single_s3.trim.fq.sub
	qsub -l $(TRIMRES) mssm_single_s3.trim.fq.sub
mssm_single_s5.trim.fq-pbs: mssm_single_s5.trim.fq.sub
	qsub -l $(TRIMRES) mssm_single_s5.trim.fq.sub
mssm_single_s6.trim.fq-pbs: mssm_single_s6.trim.fq.sub
	qsub -l $(TRIMRES) mssm_single_s6.trim.fq.sub
mssm_single_s7.trim.fq-pbs: mssm_single_s7.trim.fq.sub
	qsub -l $(TRIMRES) mssm_single_s7.trim.fq.sub
Muscle.trim.fq-pbs: Muscle.trim.fq.sub
	qsub -l $(TRIMRES) Muscle.trim.fq.sub
Neutrophils.trim.fq-pbs: Neutrophils.trim.fq.sub
	qsub -l $(TRIMRES) Neutrophils.trim.fq.sub
OFHS.trim.fq-pbs: OFHS.trim.fq.sub
	qsub -l $(TRIMRES) OFHS.trim.fq.sub
POF-eye.trim.fq-pbs: POF-eye.trim.fq.sub
	qsub -l $(TRIMRES) POF-eye.trim.fq.sub
POF-ts.trim.fq-pbs: POF-ts.trim.fq.sub
	qsub -l $(TRIMRES) POF-ts.trim.fq.sub
PSMGill.trim.fq-pbs: PSMGill.trim.fq.sub
	qsub -l $(TRIMRES) PSMGill.trim.fq.sub
Rope.trim.fq-pbs: Rope.trim.fq.sub
	qsub -l $(TRIMRES) Rope.trim.fq.sub
SMGill.trim.fq-pbs: SMGill.trim.fq.sub
	qsub -l $(TRIMRES) SMGill.trim.fq.sub
SMHS.trim.fq-pbs: SMHS.trim.fq.sub
	qsub -l $(TRIMRES) SMHS.trim.fq.sub
SN.trim.fq-pbs: SN.trim.fq.sub
	qsub -l $(TRIMRES) SN.trim.fq.sub
SPDI.trim.fq-pbs: SPDI.trim.fq.sub
	qsub -l $(TRIMRES) SPDI.trim.fq.sub
spinal_cord_paired.trim.fq-pbs: spinal_cord_paired.trim.fq.sub
	qsub -l $(TRIMRES) spinal_cord_paired.trim.fq.sub
SPK.trim.fq-pbs: SPK.trim.fq.sub
	qsub -l $(TRIMRES) SPK.trim.fq.sub
SPPI.trim.fq-pbs: SPPI.trim.fq.sub
	qsub -l $(TRIMRES) SPPI.trim.fq.sub
Stg-1-Intestine.trim.fq-pbs: Stg-1-Intestine.trim.fq.sub
	qsub -l $(TRIMRES) Stg-1-Intestine.trim.fq.sub
Stg-1-Kidney.trim.fq-pbs: Stg-1-Kidney.trim.fq.sub
	qsub -l $(TRIMRES) Stg-1-Kidney.trim.fq.sub
Stg-2-Intestine.trim.fq-pbs: Stg-2-Intestine.trim.fq.sub
	qsub -l $(TRIMRES) Stg-2-Intestine.trim.fq.sub
Stg-2-Kidney.trim.fq-pbs: Stg-2-Kidney.trim.fq.sub
	qsub -l $(TRIMRES) Stg-2-Kidney.trim.fq.sub
Stg-3-Intestine.trim.fq-pbs: Stg-3-Intestine.trim.fq.sub
	qsub -l $(TRIMRES) Stg-3-Intestine.trim.fq.sub
Stg-3-Kidney.trim.fq-pbs: Stg-3-Kidney.trim.fq.sub
	qsub -l $(TRIMRES) Stg-3-Kidney.trim.fq.sub
Stg-4-Intestine.trim.fq-pbs: Stg-4-Intestine.trim.fq.sub
	qsub -l $(TRIMRES) Stg-4-Intestine.trim.fq.sub
Stg-4-Kidney.trim.fq-pbs: Stg-4-Kidney.trim.fq.sub
	qsub -l $(TRIMRES) Stg-4-Kidney.trim.fq.sub
Stg-5-Intestine.trim.fq-pbs: Stg-5-Intestine.trim.fq.sub
	qsub -l $(TRIMRES) Stg-5-Intestine.trim.fq.sub
Stg-5-Kidney.trim.fq-pbs: Stg-5-Kidney.trim.fq.sub
	qsub -l $(TRIMRES) Stg-5-Kidney.trim.fq.sub
Stg-6-Intestine.trim.fq-pbs: Stg-6-Intestine.trim.fq.sub
	qsub -l $(TRIMRES) Stg-6-Intestine.trim.fq.sub
Stg-6-Kidney.trim.fq-pbs: Stg-6-Kidney.trim.fq.sub
	qsub -l $(TRIMRES) Stg-6-Kidney.trim.fq.sub
Stg-7-Intestine.trim.fq-pbs: Stg-7-Intestine.trim.fq.sub
	qsub -l $(TRIMRES) Stg-7-Intestine.trim.fq.sub
Stg-7-Kidney.trim.fq-pbs: Stg-7-Kidney.trim.fq.sub
	qsub -l $(TRIMRES) Stg-7-Kidney.trim.fq.sub
SW2-Gill.trim.fq-pbs: SW2-Gill.trim.fq.sub
	qsub -l $(TRIMRES) SW2-Gill.trim.fq.sub
SW2-Intestine-a.trim.fq-pbs: SW2-Intestine-a.trim.fq.sub
	qsub -l $(TRIMRES) SW2-Intestine-a.trim.fq.sub
SW2-Intestine-b.trim.fq-pbs: SW2-Intestine-b.trim.fq.sub
	qsub -l $(TRIMRES) SW2-Intestine-b.trim.fq.sub
SW2-Kidney.trim.fq-pbs: SW2-Kidney.trim.fq.sub
	qsub -l $(TRIMRES) SW2-Kidney.trim.fq.sub
trim-cluster: 614GMAAXX_1_1_pf.trim.fq-pbs 614GMAAXX_2_1_pf.trim.fq-pbs 614GMAAXX_3_1_pf.trim.fq-pbs 614GMAAXX_4_1_pf.trim.fq-pbs 614GMAAXX_6_1_pf.trim.fq-pbs 614GMAAXX_7_1_pf.trim.fq-pbs 614GMAAXX_8_1_pf.trim.fq-pbs AB057JABXX_s_7_pe.trim.fq-pbs AI.trim.fq-pbs AK.trim.fq-pbs BB0582ABXX_s_3_pe.trim.fq-pbs BB0582ABXX_s_4_pe.trim.fq-pbs BB0582ABXX_s_5_pe.trim.fq-pbs BB0582ABXX_s_6_pe.trim.fq-pbs BB0582ABXX_s_7_pe.trim.fq-pbs BB0582ABXX_s_8_pe.trim.fq-pbs brain_paired.trim.fq-pbs FW1-Gill.trim.fq-pbs FW1-Intestine.trim.fq-pbs FW1-Kidney.trim.fq-pbs G81-a.trim.fq-pbs G81-b.trim.fq-pbs G82-a.trim.fq-pbs G82-b.trim.fq-pbs G83-a.trim.fq-pbs G83-b.trim.fq-pbs Gill.trim.fq-pbs JV1_liver.trim.fq-pbs JVa_liver.trim.fq-pbs JVb_liver.trim.fq-pbs JV_intest.trim.fq-pbs Kidney.trim.fq-pbs L81-a.trim.fq-pbs L81-b.trim.fq-pbs L82-a.trim.fq-pbs L82-b.trim.fq-pbs L83-1-a.trim.fq-pbs L83-1-b.trim.fq-pbs L83-a.trim.fq-pbs L83-b.trim.fq-pbs LI.trim.fq-pbs Lips.trim.fq-pbs LK.trim.fq-pbs M0_liver.trim.fq-pbs M2_liver.trim.fq-pbs M5_liver.trim.fq-pbs Monocytes.trim.fq-pbs mssm_s_1_pair_pe.trim.fq-pbs mssm_s_1_pair_se.trim.fq-pbs mssm_s_2_pair_pe.trim.fq-pbs mssm_s_2_pair_se.trim.fq-pbs mssm_s_3_pair_pe.trim.fq-pbs mssm_s_3_pair_se.trim.fq-pbs mssm_s_5_pair_pe.trim.fq-pbs mssm_s_5_pair_se.trim.fq-pbs mssm_s_6_pair_pe.trim.fq-pbs mssm_s_6_pair_se.trim.fq-pbs mssm_s_7_pair_pe.trim.fq-pbs mssm_s_7_pair_se.trim.fq-pbs mssm_single_s1.trim.fq-pbs mssm_single_s2.trim.fq-pbs mssm_single_s3.trim.fq-pbs mssm_single_s5.trim.fq-pbs mssm_single_s6.trim.fq-pbs mssm_single_s7.trim.fq-pbs Muscle.trim.fq-pbs Neutrophils.trim.fq-pbs OFHS.trim.fq-pbs POF-eye.trim.fq-pbs POF-ts.trim.fq-pbs PSMGill.trim.fq-pbs Rope.trim.fq-pbs SMGill.trim.fq-pbs SMHS.trim.fq-pbs SN.trim.fq-pbs SPDI.trim.fq-pbs spinal_cord_paired.trim.fq-pbs SPK.trim.fq-pbs SPPI.trim.fq-pbs Stg-1-Intestine.trim.fq-pbs Stg-1-Kidney.trim.fq-pbs Stg-2-Intestine.trim.fq-pbs Stg-2-Kidney.trim.fq-pbs Stg-3-Intestine.trim.fq-pbs Stg-3-Kidney.trim.fq-pbs Stg-4-Intestine.trim.fq-pbs Stg-4-Kidney.trim.fq-pbs Stg-5-Intestine.trim.fq-pbs Stg-5-Kidney.trim.fq-pbs Stg-6-Intestine.trim.fq-pbs Stg-6-Kidney.trim.fq-pbs Stg-7-Intestine.trim.fq-pbs Stg-7-Kidney.trim.fq-pbs SW2-Gill.trim.fq-pbs SW2-Intestine-a.trim.fq-pbs SW2-Intestine-b.trim.fq-pbs SW2-Kidney.trim.fq-pbs
614GMAAXX_1_1_pf.trim.fq.keep: 614GMAAXX_1_1_pf.trim.fq
	time -o 614GMAAXX_1_1_pf.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s 614GMAAXX_1_1_pf.trim.fq.ht -R 614GMAAXX_1_1_pf.trim.fq.report 614GMAAXX_1_1_pf.trim.fq

614GMAAXX_2_1_pf.trim.fq.keep: 614GMAAXX_2_1_pf.trim.fq
	time -o 614GMAAXX_2_1_pf.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s 614GMAAXX_2_1_pf.trim.fq.ht -R 614GMAAXX_2_1_pf.trim.fq.report 614GMAAXX_2_1_pf.trim.fq

614GMAAXX_3_1_pf.trim.fq.keep: 614GMAAXX_3_1_pf.trim.fq
	time -o 614GMAAXX_3_1_pf.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s 614GMAAXX_3_1_pf.trim.fq.ht -R 614GMAAXX_3_1_pf.trim.fq.report 614GMAAXX_3_1_pf.trim.fq

614GMAAXX_4_1_pf.trim.fq.keep: 614GMAAXX_4_1_pf.trim.fq
	time -o 614GMAAXX_4_1_pf.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s 614GMAAXX_4_1_pf.trim.fq.ht -R 614GMAAXX_4_1_pf.trim.fq.report 614GMAAXX_4_1_pf.trim.fq

614GMAAXX_6_1_pf.trim.fq.keep: 614GMAAXX_6_1_pf.trim.fq
	time -o 614GMAAXX_6_1_pf.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s 614GMAAXX_6_1_pf.trim.fq.ht -R 614GMAAXX_6_1_pf.trim.fq.report 614GMAAXX_6_1_pf.trim.fq

614GMAAXX_7_1_pf.trim.fq.keep: 614GMAAXX_7_1_pf.trim.fq
	time -o 614GMAAXX_7_1_pf.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s 614GMAAXX_7_1_pf.trim.fq.ht -R 614GMAAXX_7_1_pf.trim.fq.report 614GMAAXX_7_1_pf.trim.fq

614GMAAXX_8_1_pf.trim.fq.keep: 614GMAAXX_8_1_pf.trim.fq
	time -o 614GMAAXX_8_1_pf.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s 614GMAAXX_8_1_pf.trim.fq.ht -R 614GMAAXX_8_1_pf.trim.fq.report 614GMAAXX_8_1_pf.trim.fq

AB057JABXX_s_7_pe.trim.fq.keep: AB057JABXX_s_7_pe.trim.fq
	time -o AB057JABXX_s_7_pe.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s AB057JABXX_s_7_pe.trim.fq.ht -R AB057JABXX_s_7_pe.trim.fq.report AB057JABXX_s_7_pe.trim.fq

AI.trim.fq.keep: AI.trim.fq
	time -o AI.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s AI.trim.fq.ht -R AI.trim.fq.report AI.trim.fq

AK.trim.fq.keep: AK.trim.fq
	time -o AK.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s AK.trim.fq.ht -R AK.trim.fq.report AK.trim.fq

BB0582ABXX_s_3_pe.trim.fq.keep: BB0582ABXX_s_3_pe.trim.fq
	time -o BB0582ABXX_s_3_pe.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s BB0582ABXX_s_3_pe.trim.fq.ht -R BB0582ABXX_s_3_pe.trim.fq.report BB0582ABXX_s_3_pe.trim.fq

BB0582ABXX_s_4_pe.trim.fq.keep: BB0582ABXX_s_4_pe.trim.fq
	time -o BB0582ABXX_s_4_pe.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s BB0582ABXX_s_4_pe.trim.fq.ht -R BB0582ABXX_s_4_pe.trim.fq.report BB0582ABXX_s_4_pe.trim.fq

BB0582ABXX_s_5_pe.trim.fq.keep: BB0582ABXX_s_5_pe.trim.fq
	time -o BB0582ABXX_s_5_pe.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s BB0582ABXX_s_5_pe.trim.fq.ht -R BB0582ABXX_s_5_pe.trim.fq.report BB0582ABXX_s_5_pe.trim.fq

BB0582ABXX_s_6_pe.trim.fq.keep: BB0582ABXX_s_6_pe.trim.fq
	time -o BB0582ABXX_s_6_pe.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s BB0582ABXX_s_6_pe.trim.fq.ht -R BB0582ABXX_s_6_pe.trim.fq.report BB0582ABXX_s_6_pe.trim.fq

BB0582ABXX_s_7_pe.trim.fq.keep: BB0582ABXX_s_7_pe.trim.fq
	time -o BB0582ABXX_s_7_pe.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s BB0582ABXX_s_7_pe.trim.fq.ht -R BB0582ABXX_s_7_pe.trim.fq.report BB0582ABXX_s_7_pe.trim.fq

BB0582ABXX_s_8_pe.trim.fq.keep: BB0582ABXX_s_8_pe.trim.fq
	time -o BB0582ABXX_s_8_pe.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s BB0582ABXX_s_8_pe.trim.fq.ht -R BB0582ABXX_s_8_pe.trim.fq.report BB0582ABXX_s_8_pe.trim.fq

brain_paired.trim.fq.keep: brain_paired.trim.fq
	time -o brain_paired.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s brain_paired.trim.fq.ht -R brain_paired.trim.fq.report brain_paired.trim.fq

FW1-Gill.trim.fq.keep: FW1-Gill.trim.fq
	time -o FW1-Gill.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s FW1-Gill.trim.fq.ht -R FW1-Gill.trim.fq.report FW1-Gill.trim.fq

FW1-Intestine.trim.fq.keep: FW1-Intestine.trim.fq
	time -o FW1-Intestine.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s FW1-Intestine.trim.fq.ht -R FW1-Intestine.trim.fq.report FW1-Intestine.trim.fq

FW1-Kidney.trim.fq.keep: FW1-Kidney.trim.fq
	time -o FW1-Kidney.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s FW1-Kidney.trim.fq.ht -R FW1-Kidney.trim.fq.report FW1-Kidney.trim.fq

G81-a.trim.fq.keep: G81-a.trim.fq
	time -o G81-a.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s G81-a.trim.fq.ht -R G81-a.trim.fq.report G81-a.trim.fq

G81-b.trim.fq.keep: G81-b.trim.fq
	time -o G81-b.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s G81-b.trim.fq.ht -R G81-b.trim.fq.report G81-b.trim.fq

G82-a.trim.fq.keep: G82-a.trim.fq
	time -o G82-a.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s G82-a.trim.fq.ht -R G82-a.trim.fq.report G82-a.trim.fq

G82-b.trim.fq.keep: G82-b.trim.fq
	time -o G82-b.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s G82-b.trim.fq.ht -R G82-b.trim.fq.report G82-b.trim.fq

G83-a.trim.fq.keep: G83-a.trim.fq
	time -o G83-a.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s G83-a.trim.fq.ht -R G83-a.trim.fq.report G83-a.trim.fq

G83-b.trim.fq.keep: G83-b.trim.fq
	time -o G83-b.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s G83-b.trim.fq.ht -R G83-b.trim.fq.report G83-b.trim.fq

Gill.trim.fq.keep: Gill.trim.fq
	time -o Gill.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s Gill.trim.fq.ht -R Gill.trim.fq.report Gill.trim.fq

JV1_liver.trim.fq.keep: JV1_liver.trim.fq
	time -o JV1_liver.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s JV1_liver.trim.fq.ht -R JV1_liver.trim.fq.report JV1_liver.trim.fq

JVa_liver.trim.fq.keep: JVa_liver.trim.fq
	time -o JVa_liver.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s JVa_liver.trim.fq.ht -R JVa_liver.trim.fq.report JVa_liver.trim.fq

JVb_liver.trim.fq.keep: JVb_liver.trim.fq
	time -o JVb_liver.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s JVb_liver.trim.fq.ht -R JVb_liver.trim.fq.report JVb_liver.trim.fq

JV_intest.trim.fq.keep: JV_intest.trim.fq
	time -o JV_intest.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s JV_intest.trim.fq.ht -R JV_intest.trim.fq.report JV_intest.trim.fq

Kidney.trim.fq.keep: Kidney.trim.fq
	time -o Kidney.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s Kidney.trim.fq.ht -R Kidney.trim.fq.report Kidney.trim.fq

L81-a.trim.fq.keep: L81-a.trim.fq
	time -o L81-a.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s L81-a.trim.fq.ht -R L81-a.trim.fq.report L81-a.trim.fq

L81-b.trim.fq.keep: L81-b.trim.fq
	time -o L81-b.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s L81-b.trim.fq.ht -R L81-b.trim.fq.report L81-b.trim.fq

L82-a.trim.fq.keep: L82-a.trim.fq
	time -o L82-a.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s L82-a.trim.fq.ht -R L82-a.trim.fq.report L82-a.trim.fq

L82-b.trim.fq.keep: L82-b.trim.fq
	time -o L82-b.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s L82-b.trim.fq.ht -R L82-b.trim.fq.report L82-b.trim.fq

L83-1-a.trim.fq.keep: L83-1-a.trim.fq
	time -o L83-1-a.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s L83-1-a.trim.fq.ht -R L83-1-a.trim.fq.report L83-1-a.trim.fq

L83-1-b.trim.fq.keep: L83-1-b.trim.fq
	time -o L83-1-b.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s L83-1-b.trim.fq.ht -R L83-1-b.trim.fq.report L83-1-b.trim.fq

L83-a.trim.fq.keep: L83-a.trim.fq
	time -o L83-a.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s L83-a.trim.fq.ht -R L83-a.trim.fq.report L83-a.trim.fq

L83-b.trim.fq.keep: L83-b.trim.fq
	time -o L83-b.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s L83-b.trim.fq.ht -R L83-b.trim.fq.report L83-b.trim.fq

LI.trim.fq.keep: LI.trim.fq
	time -o LI.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s LI.trim.fq.ht -R LI.trim.fq.report LI.trim.fq

Lips.trim.fq.keep: Lips.trim.fq
	time -o Lips.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s Lips.trim.fq.ht -R Lips.trim.fq.report Lips.trim.fq

LK.trim.fq.keep: LK.trim.fq
	time -o LK.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s LK.trim.fq.ht -R LK.trim.fq.report LK.trim.fq

M0_liver.trim.fq.keep: M0_liver.trim.fq
	time -o M0_liver.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s M0_liver.trim.fq.ht -R M0_liver.trim.fq.report M0_liver.trim.fq

M2_liver.trim.fq.keep: M2_liver.trim.fq
	time -o M2_liver.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s M2_liver.trim.fq.ht -R M2_liver.trim.fq.report M2_liver.trim.fq

M5_liver.trim.fq.keep: M5_liver.trim.fq
	time -o M5_liver.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s M5_liver.trim.fq.ht -R M5_liver.trim.fq.report M5_liver.trim.fq

Monocytes.trim.fq.keep: Monocytes.trim.fq
	time -o Monocytes.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s Monocytes.trim.fq.ht -R Monocytes.trim.fq.report Monocytes.trim.fq

mssm_s_1_pair_pe.trim.fq.keep: mssm_s_1_pair_pe.trim.fq
	time -o mssm_s_1_pair_pe.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s mssm_s_1_pair_pe.trim.fq.ht -R mssm_s_1_pair_pe.trim.fq.report mssm_s_1_pair_pe.trim.fq

mssm_s_1_pair_se.trim.fq.keep: mssm_s_1_pair_se.trim.fq
	time -o mssm_s_1_pair_se.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s mssm_s_1_pair_se.trim.fq.ht -R mssm_s_1_pair_se.trim.fq.report mssm_s_1_pair_se.trim.fq

mssm_s_2_pair_pe.trim.fq.keep: mssm_s_2_pair_pe.trim.fq
	time -o mssm_s_2_pair_pe.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s mssm_s_2_pair_pe.trim.fq.ht -R mssm_s_2_pair_pe.trim.fq.report mssm_s_2_pair_pe.trim.fq

mssm_s_2_pair_se.trim.fq.keep: mssm_s_2_pair_se.trim.fq
	time -o mssm_s_2_pair_se.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s mssm_s_2_pair_se.trim.fq.ht -R mssm_s_2_pair_se.trim.fq.report mssm_s_2_pair_se.trim.fq

mssm_s_3_pair_pe.trim.fq.keep: mssm_s_3_pair_pe.trim.fq
	time -o mssm_s_3_pair_pe.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s mssm_s_3_pair_pe.trim.fq.ht -R mssm_s_3_pair_pe.trim.fq.report mssm_s_3_pair_pe.trim.fq

mssm_s_3_pair_se.trim.fq.keep: mssm_s_3_pair_se.trim.fq
	time -o mssm_s_3_pair_se.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s mssm_s_3_pair_se.trim.fq.ht -R mssm_s_3_pair_se.trim.fq.report mssm_s_3_pair_se.trim.fq

mssm_s_5_pair_pe.trim.fq.keep: mssm_s_5_pair_pe.trim.fq
	time -o mssm_s_5_pair_pe.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s mssm_s_5_pair_pe.trim.fq.ht -R mssm_s_5_pair_pe.trim.fq.report mssm_s_5_pair_pe.trim.fq

mssm_s_5_pair_se.trim.fq.keep: mssm_s_5_pair_se.trim.fq
	time -o mssm_s_5_pair_se.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s mssm_s_5_pair_se.trim.fq.ht -R mssm_s_5_pair_se.trim.fq.report mssm_s_5_pair_se.trim.fq

mssm_s_6_pair_pe.trim.fq.keep: mssm_s_6_pair_pe.trim.fq
	time -o mssm_s_6_pair_pe.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s mssm_s_6_pair_pe.trim.fq.ht -R mssm_s_6_pair_pe.trim.fq.report mssm_s_6_pair_pe.trim.fq

mssm_s_6_pair_se.trim.fq.keep: mssm_s_6_pair_se.trim.fq
	time -o mssm_s_6_pair_se.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s mssm_s_6_pair_se.trim.fq.ht -R mssm_s_6_pair_se.trim.fq.report mssm_s_6_pair_se.trim.fq

mssm_s_7_pair_pe.trim.fq.keep: mssm_s_7_pair_pe.trim.fq
	time -o mssm_s_7_pair_pe.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s mssm_s_7_pair_pe.trim.fq.ht -R mssm_s_7_pair_pe.trim.fq.report mssm_s_7_pair_pe.trim.fq

mssm_s_7_pair_se.trim.fq.keep: mssm_s_7_pair_se.trim.fq
	time -o mssm_s_7_pair_se.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s mssm_s_7_pair_se.trim.fq.ht -R mssm_s_7_pair_se.trim.fq.report mssm_s_7_pair_se.trim.fq

mssm_single_s1.trim.fq.keep: mssm_single_s1.trim.fq
	time -o mssm_single_s1.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s mssm_single_s1.trim.fq.ht -R mssm_single_s1.trim.fq.report mssm_single_s1.trim.fq

mssm_single_s2.trim.fq.keep: mssm_single_s2.trim.fq
	time -o mssm_single_s2.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s mssm_single_s2.trim.fq.ht -R mssm_single_s2.trim.fq.report mssm_single_s2.trim.fq

mssm_single_s3.trim.fq.keep: mssm_single_s3.trim.fq
	time -o mssm_single_s3.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s mssm_single_s3.trim.fq.ht -R mssm_single_s3.trim.fq.report mssm_single_s3.trim.fq

mssm_single_s5.trim.fq.keep: mssm_single_s5.trim.fq
	time -o mssm_single_s5.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s mssm_single_s5.trim.fq.ht -R mssm_single_s5.trim.fq.report mssm_single_s5.trim.fq

mssm_single_s6.trim.fq.keep: mssm_single_s6.trim.fq
	time -o mssm_single_s6.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s mssm_single_s6.trim.fq.ht -R mssm_single_s6.trim.fq.report mssm_single_s6.trim.fq

mssm_single_s7.trim.fq.keep: mssm_single_s7.trim.fq
	time -o mssm_single_s7.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s mssm_single_s7.trim.fq.ht -R mssm_single_s7.trim.fq.report mssm_single_s7.trim.fq

Muscle.trim.fq.keep: Muscle.trim.fq
	time -o Muscle.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s Muscle.trim.fq.ht -R Muscle.trim.fq.report Muscle.trim.fq

Neutrophils.trim.fq.keep: Neutrophils.trim.fq
	time -o Neutrophils.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s Neutrophils.trim.fq.ht -R Neutrophils.trim.fq.report Neutrophils.trim.fq

OFHS.trim.fq.keep: OFHS.trim.fq
	time -o OFHS.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s OFHS.trim.fq.ht -R OFHS.trim.fq.report OFHS.trim.fq

POF-eye.trim.fq.keep: POF-eye.trim.fq
	time -o POF-eye.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s POF-eye.trim.fq.ht -R POF-eye.trim.fq.report POF-eye.trim.fq

POF-ts.trim.fq.keep: POF-ts.trim.fq
	time -o POF-ts.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s POF-ts.trim.fq.ht -R POF-ts.trim.fq.report POF-ts.trim.fq

PSMGill.trim.fq.keep: PSMGill.trim.fq
	time -o PSMGill.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s PSMGill.trim.fq.ht -R PSMGill.trim.fq.report PSMGill.trim.fq

Rope.trim.fq.keep: Rope.trim.fq
	time -o Rope.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s Rope.trim.fq.ht -R Rope.trim.fq.report Rope.trim.fq

SMGill.trim.fq.keep: SMGill.trim.fq
	time -o SMGill.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s SMGill.trim.fq.ht -R SMGill.trim.fq.report SMGill.trim.fq

SMHS.trim.fq.keep: SMHS.trim.fq
	time -o SMHS.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s SMHS.trim.fq.ht -R SMHS.trim.fq.report SMHS.trim.fq

SN.trim.fq.keep: SN.trim.fq
	time -o SN.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s SN.trim.fq.ht -R SN.trim.fq.report SN.trim.fq

SPDI.trim.fq.keep: SPDI.trim.fq
	time -o SPDI.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s SPDI.trim.fq.ht -R SPDI.trim.fq.report SPDI.trim.fq

spinal_cord_paired.trim.fq.keep: spinal_cord_paired.trim.fq
	time -o spinal_cord_paired.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s spinal_cord_paired.trim.fq.ht -R spinal_cord_paired.trim.fq.report spinal_cord_paired.trim.fq

SPK.trim.fq.keep: SPK.trim.fq
	time -o SPK.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s SPK.trim.fq.ht -R SPK.trim.fq.report SPK.trim.fq

SPPI.trim.fq.keep: SPPI.trim.fq
	time -o SPPI.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s SPPI.trim.fq.ht -R SPPI.trim.fq.report SPPI.trim.fq

Stg-1-Intestine.trim.fq.keep: Stg-1-Intestine.trim.fq
	time -o Stg-1-Intestine.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s Stg-1-Intestine.trim.fq.ht -R Stg-1-Intestine.trim.fq.report Stg-1-Intestine.trim.fq

Stg-1-Kidney.trim.fq.keep: Stg-1-Kidney.trim.fq
	time -o Stg-1-Kidney.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s Stg-1-Kidney.trim.fq.ht -R Stg-1-Kidney.trim.fq.report Stg-1-Kidney.trim.fq

Stg-2-Intestine.trim.fq.keep: Stg-2-Intestine.trim.fq
	time -o Stg-2-Intestine.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s Stg-2-Intestine.trim.fq.ht -R Stg-2-Intestine.trim.fq.report Stg-2-Intestine.trim.fq

Stg-2-Kidney.trim.fq.keep: Stg-2-Kidney.trim.fq
	time -o Stg-2-Kidney.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s Stg-2-Kidney.trim.fq.ht -R Stg-2-Kidney.trim.fq.report Stg-2-Kidney.trim.fq

Stg-3-Intestine.trim.fq.keep: Stg-3-Intestine.trim.fq
	time -o Stg-3-Intestine.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s Stg-3-Intestine.trim.fq.ht -R Stg-3-Intestine.trim.fq.report Stg-3-Intestine.trim.fq

Stg-3-Kidney.trim.fq.keep: Stg-3-Kidney.trim.fq
	time -o Stg-3-Kidney.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s Stg-3-Kidney.trim.fq.ht -R Stg-3-Kidney.trim.fq.report Stg-3-Kidney.trim.fq

Stg-4-Intestine.trim.fq.keep: Stg-4-Intestine.trim.fq
	time -o Stg-4-Intestine.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s Stg-4-Intestine.trim.fq.ht -R Stg-4-Intestine.trim.fq.report Stg-4-Intestine.trim.fq

Stg-4-Kidney.trim.fq.keep: Stg-4-Kidney.trim.fq
	time -o Stg-4-Kidney.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s Stg-4-Kidney.trim.fq.ht -R Stg-4-Kidney.trim.fq.report Stg-4-Kidney.trim.fq

Stg-5-Intestine.trim.fq.keep: Stg-5-Intestine.trim.fq
	time -o Stg-5-Intestine.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s Stg-5-Intestine.trim.fq.ht -R Stg-5-Intestine.trim.fq.report Stg-5-Intestine.trim.fq

Stg-5-Kidney.trim.fq.keep: Stg-5-Kidney.trim.fq
	time -o Stg-5-Kidney.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s Stg-5-Kidney.trim.fq.ht -R Stg-5-Kidney.trim.fq.report Stg-5-Kidney.trim.fq

Stg-6-Intestine.trim.fq.keep: Stg-6-Intestine.trim.fq
	time -o Stg-6-Intestine.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s Stg-6-Intestine.trim.fq.ht -R Stg-6-Intestine.trim.fq.report Stg-6-Intestine.trim.fq

Stg-6-Kidney.trim.fq.keep: Stg-6-Kidney.trim.fq
	time -o Stg-6-Kidney.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s Stg-6-Kidney.trim.fq.ht -R Stg-6-Kidney.trim.fq.report Stg-6-Kidney.trim.fq

Stg-7-Intestine.trim.fq.keep: Stg-7-Intestine.trim.fq
	time -o Stg-7-Intestine.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s Stg-7-Intestine.trim.fq.ht -R Stg-7-Intestine.trim.fq.report Stg-7-Intestine.trim.fq

Stg-7-Kidney.trim.fq.keep: Stg-7-Kidney.trim.fq
	time -o Stg-7-Kidney.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s Stg-7-Kidney.trim.fq.ht -R Stg-7-Kidney.trim.fq.report Stg-7-Kidney.trim.fq

SW2-Gill.trim.fq.keep: SW2-Gill.trim.fq
	time -o SW2-Gill.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s SW2-Gill.trim.fq.ht -R SW2-Gill.trim.fq.report SW2-Gill.trim.fq

SW2-Intestine-a.trim.fq.keep: SW2-Intestine-a.trim.fq
	time -o SW2-Intestine-a.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s SW2-Intestine-a.trim.fq.ht -R SW2-Intestine-a.trim.fq.report SW2-Intestine-a.trim.fq

SW2-Intestine-b.trim.fq.keep: SW2-Intestine-b.trim.fq
	time -o SW2-Intestine-b.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s SW2-Intestine-b.trim.fq.ht -R SW2-Intestine-b.trim.fq.report SW2-Intestine-b.trim.fq

SW2-Kidney.trim.fq.keep: SW2-Kidney.trim.fq
	time -o SW2-Kidney.trim.fq.keep.time $(KHMER)/scripts/normalize-by-median.py -x 1000000.0 -N 4 -k 23 -C 30 -s SW2-Kidney.trim.fq.ht -R SW2-Kidney.trim.fq.report SW2-Kidney.trim.fq

all.c20.fa: 614GMAAXX_1_1_pf.trim.fq.keep 614GMAAXX_2_1_pf.trim.fq.keep 614GMAAXX_3_1_pf.trim.fq.keep 614GMAAXX_4_1_pf.trim.fq.keep 614GMAAXX_6_1_pf.trim.fq.keep 614GMAAXX_7_1_pf.trim.fq.keep 614GMAAXX_8_1_pf.trim.fq.keep AB057JABXX_s_7_pe.trim.fq.keep AI.trim.fq.keep AK.trim.fq.keep BB0582ABXX_s_3_pe.trim.fq.keep BB0582ABXX_s_4_pe.trim.fq.keep BB0582ABXX_s_5_pe.trim.fq.keep BB0582ABXX_s_6_pe.trim.fq.keep BB0582ABXX_s_7_pe.trim.fq.keep BB0582ABXX_s_8_pe.trim.fq.keep brain_paired.trim.fq.keep FW1-Gill.trim.fq.keep FW1-Intestine.trim.fq.keep FW1-Kidney.trim.fq.keep G81-a.trim.fq.keep G81-b.trim.fq.keep G82-a.trim.fq.keep G82-b.trim.fq.keep G83-a.trim.fq.keep G83-b.trim.fq.keep Gill.trim.fq.keep JV1_liver.trim.fq.keep JVa_liver.trim.fq.keep JVb_liver.trim.fq.keep JV_intest.trim.fq.keep Kidney.trim.fq.keep L81-a.trim.fq.keep L81-b.trim.fq.keep L82-a.trim.fq.keep L82-b.trim.fq.keep L83-1-a.trim.fq.keep L83-1-b.trim.fq.keep L83-a.trim.fq.keep L83-b.trim.fq.keep LI.trim.fq.keep Lips.trim.fq.keep LK.trim.fq.keep M0_liver.trim.fq.keep M2_liver.trim.fq.keep M5_liver.trim.fq.keep Monocytes.trim.fq.keep mssm_s_1_pair_pe.trim.fq.keep mssm_s_1_pair_se.trim.fq.keep mssm_s_2_pair_pe.trim.fq.keep mssm_s_2_pair_se.trim.fq.keep mssm_s_3_pair_pe.trim.fq.keep mssm_s_3_pair_se.trim.fq.keep mssm_s_5_pair_pe.trim.fq.keep mssm_s_5_pair_se.trim.fq.keep mssm_s_6_pair_pe.trim.fq.keep mssm_s_6_pair_se.trim.fq.keep mssm_s_7_pair_pe.trim.fq.keep mssm_s_7_pair_se.trim.fq.keep mssm_single_s1.trim.fq.keep mssm_single_s2.trim.fq.keep mssm_single_s3.trim.fq.keep mssm_single_s5.trim.fq.keep mssm_single_s6.trim.fq.keep mssm_single_s7.trim.fq.keep Muscle.trim.fq.keep Neutrophils.trim.fq.keep OFHS.trim.fq.keep POF-eye.trim.fq.keep POF-ts.trim.fq.keep PSMGill.trim.fq.keep Rope.trim.fq.keep SMGill.trim.fq.keep SMHS.trim.fq.keep SN.trim.fq.keep SPDI.trim.fq.keep spinal_cord_paired.trim.fq.keep SPK.trim.fq.keep SPPI.trim.fq.keep Stg-1-Intestine.trim.fq.keep Stg-1-Kidney.trim.fq.keep Stg-2-Intestine.trim.fq.keep Stg-2-Kidney.trim.fq.keep Stg-3-Intestine.trim.fq.keep Stg-3-Kidney.trim.fq.keep Stg-4-Intestine.trim.fq.keep Stg-4-Kidney.trim.fq.keep Stg-5-Intestine.trim.fq.keep Stg-5-Kidney.trim.fq.keep Stg-6-Intestine.trim.fq.keep Stg-6-Kidney.trim.fq.keep Stg-7-Intestine.trim.fq.keep Stg-7-Kidney.trim.fq.keep SW2-Gill.trim.fq.keep SW2-Intestine-a.trim.fq.keep SW2-Intestine-b.trim.fq.keep SW2-Kidney.trim.fq.keep
	time -o all.c20.fq.time $(KHMER)/scripts/normalize-by-median.py -x 10000000.0 -N 4 -k 23 -C 20 -s all.c20.ht -R all.c20.report 614GMAAXX_1_1_pf.trim.fq.keep 614GMAAXX_2_1_pf.trim.fq.keep 614GMAAXX_3_1_pf.trim.fq.keep 614GMAAXX_4_1_pf.trim.fq.keep 614GMAAXX_6_1_pf.trim.fq.keep 614GMAAXX_7_1_pf.trim.fq.keep 614GMAAXX_8_1_pf.trim.fq.keep AB057JABXX_s_7_pe.trim.fq.keep AI.trim.fq.keep AK.trim.fq.keep BB0582ABXX_s_3_pe.trim.fq.keep BB0582ABXX_s_4_pe.trim.fq.keep BB0582ABXX_s_5_pe.trim.fq.keep BB0582ABXX_s_6_pe.trim.fq.keep BB0582ABXX_s_7_pe.trim.fq.keep BB0582ABXX_s_8_pe.trim.fq.keep brain_paired.trim.fq.keep FW1-Gill.trim.fq.keep FW1-Intestine.trim.fq.keep FW1-Kidney.trim.fq.keep G81-a.trim.fq.keep G81-b.trim.fq.keep G82-a.trim.fq.keep G82-b.trim.fq.keep G83-a.trim.fq.keep G83-b.trim.fq.keep Gill.trim.fq.keep JV1_liver.trim.fq.keep JVa_liver.trim.fq.keep JVb_liver.trim.fq.keep JV_intest.trim.fq.keep Kidney.trim.fq.keep L81-a.trim.fq.keep L81-b.trim.fq.keep L82-a.trim.fq.keep L82-b.trim.fq.keep L83-1-a.trim.fq.keep L83-1-b.trim.fq.keep L83-a.trim.fq.keep L83-b.trim.fq.keep LI.trim.fq.keep Lips.trim.fq.keep LK.trim.fq.keep M0_liver.trim.fq.keep M2_liver.trim.fq.keep M5_liver.trim.fq.keep Monocytes.trim.fq.keep mssm_s_1_pair_pe.trim.fq.keep mssm_s_1_pair_se.trim.fq.keep mssm_s_2_pair_pe.trim.fq.keep mssm_s_2_pair_se.trim.fq.keep mssm_s_3_pair_pe.trim.fq.keep mssm_s_3_pair_se.trim.fq.keep mssm_s_5_pair_pe.trim.fq.keep mssm_s_5_pair_se.trim.fq.keep mssm_s_6_pair_pe.trim.fq.keep mssm_s_6_pair_se.trim.fq.keep mssm_s_7_pair_pe.trim.fq.keep mssm_s_7_pair_se.trim.fq.keep mssm_single_s1.trim.fq.keep mssm_single_s2.trim.fq.keep mssm_single_s3.trim.fq.keep mssm_single_s5.trim.fq.keep mssm_single_s6.trim.fq.keep mssm_single_s7.trim.fq.keep Muscle.trim.fq.keep Neutrophils.trim.fq.keep OFHS.trim.fq.keep POF-eye.trim.fq.keep POF-ts.trim.fq.keep PSMGill.trim.fq.keep Rope.trim.fq.keep SMGill.trim.fq.keep SMHS.trim.fq.keep SN.trim.fq.keep SPDI.trim.fq.keep spinal_cord_paired.trim.fq.keep SPK.trim.fq.keep SPPI.trim.fq.keep Stg-1-Intestine.trim.fq.keep Stg-1-Kidney.trim.fq.keep Stg-2-Intestine.trim.fq.keep Stg-2-Kidney.trim.fq.keep Stg-3-Intestine.trim.fq.keep Stg-3-Kidney.trim.fq.keep Stg-4-Intestine.trim.fq.keep Stg-4-Kidney.trim.fq.keep Stg-5-Intestine.trim.fq.keep Stg-5-Kidney.trim.fq.keep Stg-6-Intestine.trim.fq.keep Stg-6-Kidney.trim.fq.keep Stg-7-Intestine.trim.fq.keep Stg-7-Kidney.trim.fq.keep SW2-Gill.trim.fq.keep SW2-Intestine-a.trim.fq.keep SW2-Intestine-b.trim.fq.keep SW2-Kidney.trim.fq.keep
	cat 614GMAAXX_1_1_pf.trim.fq.keep 614GMAAXX_2_1_pf.trim.fq.keep 614GMAAXX_3_1_pf.trim.fq.keep 614GMAAXX_4_1_pf.trim.fq.keep 614GMAAXX_6_1_pf.trim.fq.keep 614GMAAXX_7_1_pf.trim.fq.keep 614GMAAXX_8_1_pf.trim.fq.keep AB057JABXX_s_7_pe.trim.fq.keep AI.trim.fq.keep AK.trim.fq.keep BB0582ABXX_s_3_pe.trim.fq.keep BB0582ABXX_s_4_pe.trim.fq.keep BB0582ABXX_s_5_pe.trim.fq.keep BB0582ABXX_s_6_pe.trim.fq.keep BB0582ABXX_s_7_pe.trim.fq.keep BB0582ABXX_s_8_pe.trim.fq.keep brain_paired.trim.fq.keep FW1-Gill.trim.fq.keep FW1-Intestine.trim.fq.keep FW1-Kidney.trim.fq.keep G81-a.trim.fq.keep G81-b.trim.fq.keep G82-a.trim.fq.keep G82-b.trim.fq.keep G83-a.trim.fq.keep G83-b.trim.fq.keep Gill.trim.fq.keep JV1_liver.trim.fq.keep JVa_liver.trim.fq.keep JVb_liver.trim.fq.keep JV_intest.trim.fq.keep Kidney.trim.fq.keep L81-a.trim.fq.keep L81-b.trim.fq.keep L82-a.trim.fq.keep L82-b.trim.fq.keep L83-1-a.trim.fq.keep L83-1-b.trim.fq.keep L83-a.trim.fq.keep L83-b.trim.fq.keep LI.trim.fq.keep Lips.trim.fq.keep LK.trim.fq.keep M0_liver.trim.fq.keep M2_liver.trim.fq.keep M5_liver.trim.fq.keep Monocytes.trim.fq.keep mssm_s_1_pair_pe.trim.fq.keep mssm_s_1_pair_se.trim.fq.keep mssm_s_2_pair_pe.trim.fq.keep mssm_s_2_pair_se.trim.fq.keep mssm_s_3_pair_pe.trim.fq.keep mssm_s_3_pair_se.trim.fq.keep mssm_s_5_pair_pe.trim.fq.keep mssm_s_5_pair_se.trim.fq.keep mssm_s_6_pair_pe.trim.fq.keep mssm_s_6_pair_se.trim.fq.keep mssm_s_7_pair_pe.trim.fq.keep mssm_s_7_pair_se.trim.fq.keep mssm_single_s1.trim.fq.keep mssm_single_s2.trim.fq.keep mssm_single_s3.trim.fq.keep mssm_single_s5.trim.fq.keep mssm_single_s6.trim.fq.keep mssm_single_s7.trim.fq.keep Muscle.trim.fq.keep Neutrophils.trim.fq.keep OFHS.trim.fq.keep POF-eye.trim.fq.keep POF-ts.trim.fq.keep PSMGill.trim.fq.keep Rope.trim.fq.keep SMGill.trim.fq.keep SMHS.trim.fq.keep SN.trim.fq.keep SPDI.trim.fq.keep spinal_cord_paired.trim.fq.keep SPK.trim.fq.keep SPPI.trim.fq.keep Stg-1-Intestine.trim.fq.keep Stg-1-Kidney.trim.fq.keep Stg-2-Intestine.trim.fq.keep Stg-2-Kidney.trim.fq.keep Stg-3-Intestine.trim.fq.keep Stg-3-Kidney.trim.fq.keep Stg-4-Intestine.trim.fq.keep Stg-4-Kidney.trim.fq.keep Stg-5-Intestine.trim.fq.keep Stg-5-Kidney.trim.fq.keep Stg-6-Intestine.trim.fq.keep Stg-6-Kidney.trim.fq.keep Stg-7-Intestine.trim.fq.keep Stg-7-Kidney.trim.fq.keep SW2-Gill.trim.fq.keep SW2-Intestine-a.trim.fq.keep SW2-Intestine-b.trim.fq.keep SW2-Kidney.trim.fq.keep > all.c20.fq
	fastq_to_fasta -i all.c20.fq -o all.c20.fa

all.c20.fa.clean: all.c20.fa
	time -o all.c20.clean.time $(SEQCLEAN) all.c20.fa -c 8 -v ../db/univec_core.fa -s adapters.fa -l 31
clean:
	rm -rf *.time *.clean *.keep *.trim.fq *.report *.log *.trimlog *.ht cleaning_*/ *.cidx *.cln all.c20.fa all.c20.fq *.sort
