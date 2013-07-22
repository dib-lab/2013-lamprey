
KHMER=~/bin/khmer
SEQCLEAN=seqclean
TRIM=~/bin/Trimmomatic-0.25/trimmomatic-0.25.jar
TIME=/usr/bin/time -f %E%U%K%W%w%I%O
KMERGENIE=~/bin/kmergenie-1.5397/kmergenie

KMERRES=nodes=1:ppn=4,walltime=18:00:00,mem=32gb
TRIMRES=nodes=1:ppn=4,walltime=18:00:00,mem=4gb
QCRES=nodes=1:ppn=1,walltime=2:00:00,mem=4gb
DGRES=nodes=1:ppn=1,walltime=72:00:00,mem=42gb
CLEANRES=nodes=1:ppn=8,walltime=12:00:00,mem=8gb

614GMAAXX_1_1_pf.fastqc1:
	fastqc --noextract 614GMAAXX_1_1_pf.qc.fq.gz

614GMAAXX_1_1_pf.fastqc1-pbs:
	qsub -l $(QCRES) 614GMAAXX_1_1_pf.fastqc1.sub

614GMAAXX_2_1_pf.fastqc1:
	fastqc --noextract 614GMAAXX_2_1_pf.qc.fq.gz

614GMAAXX_2_1_pf.fastqc1-pbs:
	qsub -l $(QCRES) 614GMAAXX_2_1_pf.fastqc1.sub

614GMAAXX_3_1_pf.fastqc1:
	fastqc --noextract 614GMAAXX_3_1_pf.qc.fq.gz

614GMAAXX_3_1_pf.fastqc1-pbs:
	qsub -l $(QCRES) 614GMAAXX_3_1_pf.fastqc1.sub

614GMAAXX_4_1_pf.fastqc1:
	fastqc --noextract 614GMAAXX_4_1_pf.qc.fq.gz

614GMAAXX_4_1_pf.fastqc1-pbs:
	qsub -l $(QCRES) 614GMAAXX_4_1_pf.fastqc1.sub

614GMAAXX_6_1_pf.fastqc1:
	fastqc --noextract 614GMAAXX_6_1_pf.qc.fq.gz

614GMAAXX_6_1_pf.fastqc1-pbs:
	qsub -l $(QCRES) 614GMAAXX_6_1_pf.fastqc1.sub

614GMAAXX_7_1_pf.fastqc1:
	fastqc --noextract 614GMAAXX_7_1_pf.qc.fq.gz

614GMAAXX_7_1_pf.fastqc1-pbs:
	qsub -l $(QCRES) 614GMAAXX_7_1_pf.fastqc1.sub

614GMAAXX_8_1_pf.fastqc1:
	fastqc --noextract 614GMAAXX_8_1_pf.qc.fq.gz

614GMAAXX_8_1_pf.fastqc1-pbs:
	qsub -l $(QCRES) 614GMAAXX_8_1_pf.fastqc1.sub

AB057JABXX_s_7_pe.fastqc1:
	fastqc --noextract AB057JABXX_s_7_pe.trim.fq.gz

AB057JABXX_s_7_pe.fastqc1-pbs:
	qsub -l $(QCRES) AB057JABXX_s_7_pe.fastqc1.sub

AI.fastqc1:
	fastqc --noextract AI.fq.gz

AI.fastqc1-pbs:
	qsub -l $(QCRES) AI.fastqc1.sub

AK.fastqc1:
	fastqc --noextract AK.fq.gz

AK.fastqc1-pbs:
	qsub -l $(QCRES) AK.fastqc1.sub

BB0582ABXX_s_3_pe.fastqc1:
	fastqc --noextract BB0582ABXX_s_3_pe.trim.fq.gz

BB0582ABXX_s_3_pe.fastqc1-pbs:
	qsub -l $(QCRES) BB0582ABXX_s_3_pe.fastqc1.sub

BB0582ABXX_s_4_pe.fastqc1:
	fastqc --noextract BB0582ABXX_s_4_pe.trim.fq.gz

BB0582ABXX_s_4_pe.fastqc1-pbs:
	qsub -l $(QCRES) BB0582ABXX_s_4_pe.fastqc1.sub

BB0582ABXX_s_5_pe.fastqc1:
	fastqc --noextract BB0582ABXX_s_5_pe.trim.fq.gz

BB0582ABXX_s_5_pe.fastqc1-pbs:
	qsub -l $(QCRES) BB0582ABXX_s_5_pe.fastqc1.sub

BB0582ABXX_s_6_pe.fastqc1:
	fastqc --noextract BB0582ABXX_s_6_pe.trim.fq.gz

BB0582ABXX_s_6_pe.fastqc1-pbs:
	qsub -l $(QCRES) BB0582ABXX_s_6_pe.fastqc1.sub

BB0582ABXX_s_7_pe.fastqc1:
	fastqc --noextract BB0582ABXX_s_7_pe.trim.fq.gz

BB0582ABXX_s_7_pe.fastqc1-pbs:
	qsub -l $(QCRES) BB0582ABXX_s_7_pe.fastqc1.sub

BB0582ABXX_s_8_pe.fastqc1:
	fastqc --noextract BB0582ABXX_s_8_pe.trim.fq.gz

BB0582ABXX_s_8_pe.fastqc1-pbs:
	qsub -l $(QCRES) BB0582ABXX_s_8_pe.fastqc1.sub

brain_paired.fastqc1:
	fastqc --noextract brain_paired.qc.fq.gz

brain_paired.fastqc1-pbs:
	qsub -l $(QCRES) brain_paired.fastqc1.sub

FW1-Gill.fastqc1:
	fastqc --noextract FW1-Gill.fq.gz

FW1-Gill.fastqc1-pbs:
	qsub -l $(QCRES) FW1-Gill.fastqc1.sub

FW1-Intestine.fastqc1:
	fastqc --noextract FW1-Intestine.fq.gz

FW1-Intestine.fastqc1-pbs:
	qsub -l $(QCRES) FW1-Intestine.fastqc1.sub

FW1-Kidney.fastqc1:
	fastqc --noextract FW1-Kidney.fq.gz

FW1-Kidney.fastqc1-pbs:
	qsub -l $(QCRES) FW1-Kidney.fastqc1.sub

G81-a.fastqc1:
	fastqc --noextract G81-a.fq.gz

G81-a.fastqc1-pbs:
	qsub -l $(QCRES) G81-a.fastqc1.sub

G81-b.fastqc1:
	fastqc --noextract G81-b.fq.gz

G81-b.fastqc1-pbs:
	qsub -l $(QCRES) G81-b.fastqc1.sub

G82-a.fastqc1:
	fastqc --noextract G82-a.fq.gz

G82-a.fastqc1-pbs:
	qsub -l $(QCRES) G82-a.fastqc1.sub

G82-b.fastqc1:
	fastqc --noextract G82-b.fq.gz

G82-b.fastqc1-pbs:
	qsub -l $(QCRES) G82-b.fastqc1.sub

G83-a.fastqc1:
	fastqc --noextract G83-a.fq.gz

G83-a.fastqc1-pbs:
	qsub -l $(QCRES) G83-a.fastqc1.sub

G83-b.fastqc1:
	fastqc --noextract G83-b.fq.gz

G83-b.fastqc1-pbs:
	qsub -l $(QCRES) G83-b.fastqc1.sub

Gill.fastqc1:
	fastqc --noextract Gill.fq.gz

Gill.fastqc1-pbs:
	qsub -l $(QCRES) Gill.fastqc1.sub

JV1_liver.fastqc1:
	fastqc --noextract JV1_liver.fq.gz

JV1_liver.fastqc1-pbs:
	qsub -l $(QCRES) JV1_liver.fastqc1.sub

JVa_liver.fastqc1:
	fastqc --noextract JVa_liver.fq.gz

JVa_liver.fastqc1-pbs:
	qsub -l $(QCRES) JVa_liver.fastqc1.sub

JVb_liver.fastqc1:
	fastqc --noextract JVb_liver.fq.gz

JVb_liver.fastqc1-pbs:
	qsub -l $(QCRES) JVb_liver.fastqc1.sub

JV_intest.fastqc1:
	fastqc --noextract JV_intest.fq.gz

JV_intest.fastqc1-pbs:
	qsub -l $(QCRES) JV_intest.fastqc1.sub

Kidney.fastqc1:
	fastqc --noextract Kidney.fq.gz

Kidney.fastqc1-pbs:
	qsub -l $(QCRES) Kidney.fastqc1.sub

L81-a.fastqc1:
	fastqc --noextract L81-a.fq.gz

L81-a.fastqc1-pbs:
	qsub -l $(QCRES) L81-a.fastqc1.sub

L81-b.fastqc1:
	fastqc --noextract L81-b.fq.gz

L81-b.fastqc1-pbs:
	qsub -l $(QCRES) L81-b.fastqc1.sub

L82-a.fastqc1:
	fastqc --noextract L82-a.fq.gz

L82-a.fastqc1-pbs:
	qsub -l $(QCRES) L82-a.fastqc1.sub

L82-b.fastqc1:
	fastqc --noextract L82-b.fq.gz

L82-b.fastqc1-pbs:
	qsub -l $(QCRES) L82-b.fastqc1.sub

L83-1-a.fastqc1:
	fastqc --noextract L83-1-a.fq.gz

L83-1-a.fastqc1-pbs:
	qsub -l $(QCRES) L83-1-a.fastqc1.sub

L83-1-b.fastqc1:
	fastqc --noextract L83-1-b.fq.gz

L83-1-b.fastqc1-pbs:
	qsub -l $(QCRES) L83-1-b.fastqc1.sub

L83-a.fastqc1:
	fastqc --noextract L83-a.fq.gz

L83-a.fastqc1-pbs:
	qsub -l $(QCRES) L83-a.fastqc1.sub

L83-b.fastqc1:
	fastqc --noextract L83-b.fq.gz

L83-b.fastqc1-pbs:
	qsub -l $(QCRES) L83-b.fastqc1.sub

LI.fastqc1:
	fastqc --noextract LI.fq.gz

LI.fastqc1-pbs:
	qsub -l $(QCRES) LI.fastqc1.sub

Lips.fastqc1:
	fastqc --noextract Lips.fq.gz

Lips.fastqc1-pbs:
	qsub -l $(QCRES) Lips.fastqc1.sub

LK.fastqc1:
	fastqc --noextract LK.fq.gz

LK.fastqc1-pbs:
	qsub -l $(QCRES) LK.fastqc1.sub

M0_liver.fastqc1:
	fastqc --noextract M0_liver.fq.gz

M0_liver.fastqc1-pbs:
	qsub -l $(QCRES) M0_liver.fastqc1.sub

M2_liver.fastqc1:
	fastqc --noextract M2_liver.fq.gz

M2_liver.fastqc1-pbs:
	qsub -l $(QCRES) M2_liver.fastqc1.sub

M5_liver.fastqc1:
	fastqc --noextract M5_liver.fq.gz

M5_liver.fastqc1-pbs:
	qsub -l $(QCRES) M5_liver.fastqc1.sub

Monocytes.fastqc1:
	fastqc --noextract Monocytes.fq.gz

Monocytes.fastqc1-pbs:
	qsub -l $(QCRES) Monocytes.fastqc1.sub

mssm_single_s1.fastqc1:
	fastqc --noextract mssm_single_s1.trim.fq.gz

mssm_single_s1.fastqc1-pbs:
	qsub -l $(QCRES) mssm_single_s1.fastqc1.sub

mssm_single_s2.fastqc1:
	fastqc --noextract mssm_single_s2.trim.fq.gz

mssm_single_s2.fastqc1-pbs:
	qsub -l $(QCRES) mssm_single_s2.fastqc1.sub

mssm_single_s3.fastqc1:
	fastqc --noextract mssm_single_s3.trim.fq.gz

mssm_single_s3.fastqc1-pbs:
	qsub -l $(QCRES) mssm_single_s3.fastqc1.sub

mssm_single_s5.fastqc1:
	fastqc --noextract mssm_single_s5.trim.fq.gz

mssm_single_s5.fastqc1-pbs:
	qsub -l $(QCRES) mssm_single_s5.fastqc1.sub

mssm_single_s6.fastqc1:
	fastqc --noextract mssm_single_s6.trim.fq.gz

mssm_single_s6.fastqc1-pbs:
	qsub -l $(QCRES) mssm_single_s6.fastqc1.sub

mssm_single_s7.fastqc1:
	fastqc --noextract mssm_single_s7.trim.fq.gz

mssm_single_s7.fastqc1-pbs:
	qsub -l $(QCRES) mssm_single_s7.fastqc1.sub

Muscle.fastqc1:
	fastqc --noextract Muscle.fq.gz

Muscle.fastqc1-pbs:
	qsub -l $(QCRES) Muscle.fastqc1.sub

Neutrophils.fastqc1:
	fastqc --noextract Neutrophils.fq.gz

Neutrophils.fastqc1-pbs:
	qsub -l $(QCRES) Neutrophils.fastqc1.sub

OFHS.fastqc1:
	fastqc --noextract OFHS.fq.gz

OFHS.fastqc1-pbs:
	qsub -l $(QCRES) OFHS.fastqc1.sub

POF-eye.fastqc1:
	fastqc --noextract POF-eye.fq.gz

POF-eye.fastqc1-pbs:
	qsub -l $(QCRES) POF-eye.fastqc1.sub

POF-ts.fastqc1:
	fastqc --noextract POF-ts.fq.gz

POF-ts.fastqc1-pbs:
	qsub -l $(QCRES) POF-ts.fastqc1.sub

PSMGill.fastqc1:
	fastqc --noextract PSMGill.fq.gz

PSMGill.fastqc1-pbs:
	qsub -l $(QCRES) PSMGill.fastqc1.sub

Rope.fastqc1:
	fastqc --noextract Rope.fq.gz

Rope.fastqc1-pbs:
	qsub -l $(QCRES) Rope.fastqc1.sub

SMGill.fastqc1:
	fastqc --noextract SMGill.fq.gz

SMGill.fastqc1-pbs:
	qsub -l $(QCRES) SMGill.fastqc1.sub

SMHS.fastqc1:
	fastqc --noextract SMHS.fq.gz

SMHS.fastqc1-pbs:
	qsub -l $(QCRES) SMHS.fastqc1.sub

SN.fastqc1:
	fastqc --noextract SN.fq.gz

SN.fastqc1-pbs:
	qsub -l $(QCRES) SN.fastqc1.sub

SPDI.fastqc1:
	fastqc --noextract SPDI.fq.gz

SPDI.fastqc1-pbs:
	qsub -l $(QCRES) SPDI.fastqc1.sub

spinal_cord_paired.fastqc1:
	fastqc --noextract spinal_cord_paired.qc.fq.gz

spinal_cord_paired.fastqc1-pbs:
	qsub -l $(QCRES) spinal_cord_paired.fastqc1.sub

SPK.fastqc1:
	fastqc --noextract SPK.fq.gz

SPK.fastqc1-pbs:
	qsub -l $(QCRES) SPK.fastqc1.sub

SPPI.fastqc1:
	fastqc --noextract SPPI.fq.gz

SPPI.fastqc1-pbs:
	qsub -l $(QCRES) SPPI.fastqc1.sub

Stg-1-Intestine.fastqc1:
	fastqc --noextract Stg-1-Intestine.fq.gz

Stg-1-Intestine.fastqc1-pbs:
	qsub -l $(QCRES) Stg-1-Intestine.fastqc1.sub

Stg-1-Kidney.fastqc1:
	fastqc --noextract Stg-1-Kidney.fq.gz

Stg-1-Kidney.fastqc1-pbs:
	qsub -l $(QCRES) Stg-1-Kidney.fastqc1.sub

Stg-2-Intestine.fastqc1:
	fastqc --noextract Stg-2-Intestine.fq.gz

Stg-2-Intestine.fastqc1-pbs:
	qsub -l $(QCRES) Stg-2-Intestine.fastqc1.sub

Stg-2-Kidney.fastqc1:
	fastqc --noextract Stg-2-Kidney.fq.gz

Stg-2-Kidney.fastqc1-pbs:
	qsub -l $(QCRES) Stg-2-Kidney.fastqc1.sub

Stg-3-Intestine.fastqc1:
	fastqc --noextract Stg-3-Intestine.fq.gz

Stg-3-Intestine.fastqc1-pbs:
	qsub -l $(QCRES) Stg-3-Intestine.fastqc1.sub

Stg-3-Kidney.fastqc1:
	fastqc --noextract Stg-3-Kidney.fq.gz

Stg-3-Kidney.fastqc1-pbs:
	qsub -l $(QCRES) Stg-3-Kidney.fastqc1.sub

Stg-4-Intestine.fastqc1:
	fastqc --noextract Stg-4-Intestine.fq.gz

Stg-4-Intestine.fastqc1-pbs:
	qsub -l $(QCRES) Stg-4-Intestine.fastqc1.sub

Stg-4-Kidney.fastqc1:
	fastqc --noextract Stg-4-Kidney.fq.gz

Stg-4-Kidney.fastqc1-pbs:
	qsub -l $(QCRES) Stg-4-Kidney.fastqc1.sub

Stg-5-Intestine.fastqc1:
	fastqc --noextract Stg-5-Intestine.fq.gz

Stg-5-Intestine.fastqc1-pbs:
	qsub -l $(QCRES) Stg-5-Intestine.fastqc1.sub

Stg-5-Kidney.fastqc1:
	fastqc --noextract Stg-5-Kidney.fq.gz

Stg-5-Kidney.fastqc1-pbs:
	qsub -l $(QCRES) Stg-5-Kidney.fastqc1.sub

Stg-6-Intestine.fastqc1:
	fastqc --noextract Stg-6-Intestine.fq.gz

Stg-6-Intestine.fastqc1-pbs:
	qsub -l $(QCRES) Stg-6-Intestine.fastqc1.sub

Stg-6-Kidney.fastqc1:
	fastqc --noextract Stg-6-Kidney.fq.gz

Stg-6-Kidney.fastqc1-pbs:
	qsub -l $(QCRES) Stg-6-Kidney.fastqc1.sub

Stg-7-Intestine.fastqc1:
	fastqc --noextract Stg-7-Intestine.fq.gz

Stg-7-Intestine.fastqc1-pbs:
	qsub -l $(QCRES) Stg-7-Intestine.fastqc1.sub

Stg-7-Kidney.fastqc1:
	fastqc --noextract Stg-7-Kidney.fq.gz

Stg-7-Kidney.fastqc1-pbs:
	qsub -l $(QCRES) Stg-7-Kidney.fastqc1.sub

SW2-Gill.fastqc1:
	fastqc --noextract SW2-Gill.fq.gz

SW2-Gill.fastqc1-pbs:
	qsub -l $(QCRES) SW2-Gill.fastqc1.sub

SW2-Intestine-a.fastqc1:
	fastqc --noextract SW2-Intestine-a.fq.gz

SW2-Intestine-a.fastqc1-pbs:
	qsub -l $(QCRES) SW2-Intestine-a.fastqc1.sub

SW2-Intestine-b.fastqc1:
	fastqc --noextract SW2-Intestine-b.fq.gz

SW2-Intestine-b.fastqc1-pbs:
	qsub -l $(QCRES) SW2-Intestine-b.fastqc1.sub

SW2-Kidney.fastqc1:
	fastqc --noextract SW2-Kidney.fq.gz

SW2-Kidney.fastqc1-pbs:
	qsub -l $(QCRES) SW2-Kidney.fastqc1.sub

fastqc1-pbs: 614GMAAXX_1_1_pf.fastqc1-pbs 614GMAAXX_2_1_pf.fastqc1-pbs 614GMAAXX_3_1_pf.fastqc1-pbs 614GMAAXX_4_1_pf.fastqc1-pbs 614GMAAXX_6_1_pf.fastqc1-pbs 614GMAAXX_7_1_pf.fastqc1-pbs 614GMAAXX_8_1_pf.fastqc1-pbs AB057JABXX_s_7_pe.fastqc1-pbs AI.fastqc1-pbs AK.fastqc1-pbs BB0582ABXX_s_3_pe.fastqc1-pbs BB0582ABXX_s_4_pe.fastqc1-pbs BB0582ABXX_s_5_pe.fastqc1-pbs BB0582ABXX_s_6_pe.fastqc1-pbs BB0582ABXX_s_7_pe.fastqc1-pbs BB0582ABXX_s_8_pe.fastqc1-pbs brain_paired.fastqc1-pbs FW1-Gill.fastqc1-pbs FW1-Intestine.fastqc1-pbs FW1-Kidney.fastqc1-pbs G81-a.fastqc1-pbs G81-b.fastqc1-pbs G82-a.fastqc1-pbs G82-b.fastqc1-pbs G83-a.fastqc1-pbs G83-b.fastqc1-pbs Gill.fastqc1-pbs JV1_liver.fastqc1-pbs JVa_liver.fastqc1-pbs JVb_liver.fastqc1-pbs JV_intest.fastqc1-pbs Kidney.fastqc1-pbs L81-a.fastqc1-pbs L81-b.fastqc1-pbs L82-a.fastqc1-pbs L82-b.fastqc1-pbs L83-1-a.fastqc1-pbs L83-1-b.fastqc1-pbs L83-a.fastqc1-pbs L83-b.fastqc1-pbs LI.fastqc1-pbs Lips.fastqc1-pbs LK.fastqc1-pbs M0_liver.fastqc1-pbs M2_liver.fastqc1-pbs M5_liver.fastqc1-pbs Monocytes.fastqc1-pbs mssm_single_s1.fastqc1-pbs mssm_single_s2.fastqc1-pbs mssm_single_s3.fastqc1-pbs mssm_single_s5.fastqc1-pbs mssm_single_s6.fastqc1-pbs mssm_single_s7.fastqc1-pbs Muscle.fastqc1-pbs Neutrophils.fastqc1-pbs OFHS.fastqc1-pbs POF-eye.fastqc1-pbs POF-ts.fastqc1-pbs PSMGill.fastqc1-pbs Rope.fastqc1-pbs SMGill.fastqc1-pbs SMHS.fastqc1-pbs SN.fastqc1-pbs SPDI.fastqc1-pbs spinal_cord_paired.fastqc1-pbs SPK.fastqc1-pbs SPPI.fastqc1-pbs Stg-1-Intestine.fastqc1-pbs Stg-1-Kidney.fastqc1-pbs Stg-2-Intestine.fastqc1-pbs Stg-2-Kidney.fastqc1-pbs Stg-3-Intestine.fastqc1-pbs Stg-3-Kidney.fastqc1-pbs Stg-4-Intestine.fastqc1-pbs Stg-4-Kidney.fastqc1-pbs Stg-5-Intestine.fastqc1-pbs Stg-5-Kidney.fastqc1-pbs Stg-6-Intestine.fastqc1-pbs Stg-6-Kidney.fastqc1-pbs Stg-7-Intestine.fastqc1-pbs Stg-7-Kidney.fastqc1-pbs SW2-Gill.fastqc1-pbs SW2-Intestine-a.fastqc1-pbs SW2-Intestine-b.fastqc1-pbs SW2-Kidney.fastqc1-pbs
614GMAAXX_1_1_pf.trim.fq: 614GMAAXX_1_1_pf.qc.fq.gz
		$(TIME) -o 614GMAAXX_1_1_pf.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 614GMAAXX_1_1_pf.trimlog -threads 4 614GMAAXX_1_1_pf.qc.fq.gz 614GMAAXX_1_1_pf.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

614GMAAXX_1_1_pf.trim.fq-pbs:
	qsub -l $(TRIMRES) 614GMAAXX_1_1_pf.trim.sub

614GMAAXX_2_1_pf.trim.fq: 614GMAAXX_2_1_pf.qc.fq.gz
		$(TIME) -o 614GMAAXX_2_1_pf.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 614GMAAXX_2_1_pf.trimlog -threads 4 614GMAAXX_2_1_pf.qc.fq.gz 614GMAAXX_2_1_pf.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

614GMAAXX_2_1_pf.trim.fq-pbs:
	qsub -l $(TRIMRES) 614GMAAXX_2_1_pf.trim.sub

614GMAAXX_3_1_pf.trim.fq: 614GMAAXX_3_1_pf.qc.fq.gz
		$(TIME) -o 614GMAAXX_3_1_pf.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 614GMAAXX_3_1_pf.trimlog -threads 4 614GMAAXX_3_1_pf.qc.fq.gz 614GMAAXX_3_1_pf.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

614GMAAXX_3_1_pf.trim.fq-pbs:
	qsub -l $(TRIMRES) 614GMAAXX_3_1_pf.trim.sub

614GMAAXX_4_1_pf.trim.fq: 614GMAAXX_4_1_pf.qc.fq.gz
		$(TIME) -o 614GMAAXX_4_1_pf.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 614GMAAXX_4_1_pf.trimlog -threads 4 614GMAAXX_4_1_pf.qc.fq.gz 614GMAAXX_4_1_pf.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

614GMAAXX_4_1_pf.trim.fq-pbs:
	qsub -l $(TRIMRES) 614GMAAXX_4_1_pf.trim.sub

614GMAAXX_6_1_pf.trim.fq: 614GMAAXX_6_1_pf.qc.fq.gz
		$(TIME) -o 614GMAAXX_6_1_pf.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 614GMAAXX_6_1_pf.trimlog -threads 4 614GMAAXX_6_1_pf.qc.fq.gz 614GMAAXX_6_1_pf.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

614GMAAXX_6_1_pf.trim.fq-pbs:
	qsub -l $(TRIMRES) 614GMAAXX_6_1_pf.trim.sub

614GMAAXX_7_1_pf.trim.fq: 614GMAAXX_7_1_pf.qc.fq.gz
		$(TIME) -o 614GMAAXX_7_1_pf.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 614GMAAXX_7_1_pf.trimlog -threads 4 614GMAAXX_7_1_pf.qc.fq.gz 614GMAAXX_7_1_pf.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

614GMAAXX_7_1_pf.trim.fq-pbs:
	qsub -l $(TRIMRES) 614GMAAXX_7_1_pf.trim.sub

614GMAAXX_8_1_pf.trim.fq: 614GMAAXX_8_1_pf.qc.fq.gz
		$(TIME) -o 614GMAAXX_8_1_pf.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog 614GMAAXX_8_1_pf.trimlog -threads 4 614GMAAXX_8_1_pf.qc.fq.gz 614GMAAXX_8_1_pf.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

614GMAAXX_8_1_pf.trim.fq-pbs:
	qsub -l $(TRIMRES) 614GMAAXX_8_1_pf.trim.sub

AB057JABXX_s_7_pe.trim.fq: AB057JABXX_s_7_pe.trim.fq.gz
		$(TIME) -o AB057JABXX_s_7_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog AB057JABXX_s_7_pe.trimlog -threads 4 AB057JABXX_s_7_pe.trim.fq.gz AB057JABXX_s_7_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

AB057JABXX_s_7_pe.trim.fq-pbs:
	qsub -l $(TRIMRES) AB057JABXX_s_7_pe.trim.sub

AI.trim.fq: AI.fq.gz
		$(TIME) -o AI.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog AI.trimlog -threads 4 AI.fq.gz AI.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

AI.trim.fq-pbs:
	qsub -l $(TRIMRES) AI.trim.sub

AK.trim.fq: AK.fq.gz
		$(TIME) -o AK.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog AK.trimlog -threads 4 AK.fq.gz AK.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

AK.trim.fq-pbs:
	qsub -l $(TRIMRES) AK.trim.sub

BB0582ABXX_s_3_pe.trim.fq: BB0582ABXX_s_3_pe.trim.fq.gz
		$(TIME) -o BB0582ABXX_s_3_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog BB0582ABXX_s_3_pe.trimlog -threads 4 BB0582ABXX_s_3_pe.trim.fq.gz BB0582ABXX_s_3_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

BB0582ABXX_s_3_pe.trim.fq-pbs:
	qsub -l $(TRIMRES) BB0582ABXX_s_3_pe.trim.sub

BB0582ABXX_s_4_pe.trim.fq: BB0582ABXX_s_4_pe.trim.fq.gz
		$(TIME) -o BB0582ABXX_s_4_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog BB0582ABXX_s_4_pe.trimlog -threads 4 BB0582ABXX_s_4_pe.trim.fq.gz BB0582ABXX_s_4_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

BB0582ABXX_s_4_pe.trim.fq-pbs:
	qsub -l $(TRIMRES) BB0582ABXX_s_4_pe.trim.sub

BB0582ABXX_s_5_pe.trim.fq: BB0582ABXX_s_5_pe.trim.fq.gz
		$(TIME) -o BB0582ABXX_s_5_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog BB0582ABXX_s_5_pe.trimlog -threads 4 BB0582ABXX_s_5_pe.trim.fq.gz BB0582ABXX_s_5_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

BB0582ABXX_s_5_pe.trim.fq-pbs:
	qsub -l $(TRIMRES) BB0582ABXX_s_5_pe.trim.sub

BB0582ABXX_s_6_pe.trim.fq: BB0582ABXX_s_6_pe.trim.fq.gz
		$(TIME) -o BB0582ABXX_s_6_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog BB0582ABXX_s_6_pe.trimlog -threads 4 BB0582ABXX_s_6_pe.trim.fq.gz BB0582ABXX_s_6_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

BB0582ABXX_s_6_pe.trim.fq-pbs:
	qsub -l $(TRIMRES) BB0582ABXX_s_6_pe.trim.sub

BB0582ABXX_s_7_pe.trim.fq: BB0582ABXX_s_7_pe.trim.fq.gz
		$(TIME) -o BB0582ABXX_s_7_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog BB0582ABXX_s_7_pe.trimlog -threads 4 BB0582ABXX_s_7_pe.trim.fq.gz BB0582ABXX_s_7_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

BB0582ABXX_s_7_pe.trim.fq-pbs:
	qsub -l $(TRIMRES) BB0582ABXX_s_7_pe.trim.sub

BB0582ABXX_s_8_pe.trim.fq: BB0582ABXX_s_8_pe.trim.fq.gz
		$(TIME) -o BB0582ABXX_s_8_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog BB0582ABXX_s_8_pe.trimlog -threads 4 BB0582ABXX_s_8_pe.trim.fq.gz BB0582ABXX_s_8_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

BB0582ABXX_s_8_pe.trim.fq-pbs:
	qsub -l $(TRIMRES) BB0582ABXX_s_8_pe.trim.sub

brain_paired.trim.fq: brain_paired.qc.fq.gz
		$(TIME) -o brain_paired.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog brain_paired.trimlog -threads 4 brain_paired.qc.fq.gz brain_paired.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

brain_paired.trim.fq-pbs:
	qsub -l $(TRIMRES) brain_paired.trim.sub

FW1-Gill.trim.fq: FW1-Gill.fq.gz
		$(TIME) -o FW1-Gill.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog FW1-Gill.trimlog -threads 4 FW1-Gill.fq.gz FW1-Gill.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

FW1-Gill.trim.fq-pbs:
	qsub -l $(TRIMRES) FW1-Gill.trim.sub

FW1-Intestine.trim.fq: FW1-Intestine.fq.gz
		$(TIME) -o FW1-Intestine.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog FW1-Intestine.trimlog -threads 4 FW1-Intestine.fq.gz FW1-Intestine.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

FW1-Intestine.trim.fq-pbs:
	qsub -l $(TRIMRES) FW1-Intestine.trim.sub

FW1-Kidney.trim.fq: FW1-Kidney.fq.gz
		$(TIME) -o FW1-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog FW1-Kidney.trimlog -threads 4 FW1-Kidney.fq.gz FW1-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

FW1-Kidney.trim.fq-pbs:
	qsub -l $(TRIMRES) FW1-Kidney.trim.sub

G81-a.trim.fq: G81-a.fq.gz
		$(TIME) -o G81-a.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog G81-a.trimlog -threads 4 G81-a.fq.gz G81-a.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

G81-a.trim.fq-pbs:
	qsub -l $(TRIMRES) G81-a.trim.sub

G81-b.trim.fq: G81-b.fq.gz
		$(TIME) -o G81-b.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog G81-b.trimlog -threads 4 G81-b.fq.gz G81-b.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

G81-b.trim.fq-pbs:
	qsub -l $(TRIMRES) G81-b.trim.sub

G82-a.trim.fq: G82-a.fq.gz
		$(TIME) -o G82-a.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog G82-a.trimlog -threads 4 G82-a.fq.gz G82-a.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

G82-a.trim.fq-pbs:
	qsub -l $(TRIMRES) G82-a.trim.sub

G82-b.trim.fq: G82-b.fq.gz
		$(TIME) -o G82-b.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog G82-b.trimlog -threads 4 G82-b.fq.gz G82-b.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

G82-b.trim.fq-pbs:
	qsub -l $(TRIMRES) G82-b.trim.sub

G83-a.trim.fq: G83-a.fq.gz
		$(TIME) -o G83-a.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog G83-a.trimlog -threads 4 G83-a.fq.gz G83-a.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

G83-a.trim.fq-pbs:
	qsub -l $(TRIMRES) G83-a.trim.sub

G83-b.trim.fq: G83-b.fq.gz
		$(TIME) -o G83-b.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog G83-b.trimlog -threads 4 G83-b.fq.gz G83-b.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

G83-b.trim.fq-pbs:
	qsub -l $(TRIMRES) G83-b.trim.sub

Gill.trim.fq: Gill.fq.gz
		$(TIME) -o Gill.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Gill.trimlog -threads 4 Gill.fq.gz Gill.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Gill.trim.fq-pbs:
	qsub -l $(TRIMRES) Gill.trim.sub

JV1_liver.trim.fq: JV1_liver.fq.gz
		$(TIME) -o JV1_liver.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog JV1_liver.trimlog -threads 4 JV1_liver.fq.gz JV1_liver.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

JV1_liver.trim.fq-pbs:
	qsub -l $(TRIMRES) JV1_liver.trim.sub

JVa_liver.trim.fq: JVa_liver.fq.gz
		$(TIME) -o JVa_liver.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog JVa_liver.trimlog -threads 4 JVa_liver.fq.gz JVa_liver.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

JVa_liver.trim.fq-pbs:
	qsub -l $(TRIMRES) JVa_liver.trim.sub

JVb_liver.trim.fq: JVb_liver.fq.gz
		$(TIME) -o JVb_liver.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog JVb_liver.trimlog -threads 4 JVb_liver.fq.gz JVb_liver.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

JVb_liver.trim.fq-pbs:
	qsub -l $(TRIMRES) JVb_liver.trim.sub

JV_intest.trim.fq: JV_intest.fq.gz
		$(TIME) -o JV_intest.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog JV_intest.trimlog -threads 4 JV_intest.fq.gz JV_intest.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

JV_intest.trim.fq-pbs:
	qsub -l $(TRIMRES) JV_intest.trim.sub

Kidney.trim.fq: Kidney.fq.gz
		$(TIME) -o Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Kidney.trimlog -threads 4 Kidney.fq.gz Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Kidney.trim.fq-pbs:
	qsub -l $(TRIMRES) Kidney.trim.sub

L81-a.trim.fq: L81-a.fq.gz
		$(TIME) -o L81-a.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog L81-a.trimlog -threads 4 L81-a.fq.gz L81-a.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

L81-a.trim.fq-pbs:
	qsub -l $(TRIMRES) L81-a.trim.sub

L81-b.trim.fq: L81-b.fq.gz
		$(TIME) -o L81-b.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog L81-b.trimlog -threads 4 L81-b.fq.gz L81-b.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

L81-b.trim.fq-pbs:
	qsub -l $(TRIMRES) L81-b.trim.sub

L82-a.trim.fq: L82-a.fq.gz
		$(TIME) -o L82-a.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog L82-a.trimlog -threads 4 L82-a.fq.gz L82-a.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

L82-a.trim.fq-pbs:
	qsub -l $(TRIMRES) L82-a.trim.sub

L82-b.trim.fq: L82-b.fq.gz
		$(TIME) -o L82-b.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog L82-b.trimlog -threads 4 L82-b.fq.gz L82-b.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

L82-b.trim.fq-pbs:
	qsub -l $(TRIMRES) L82-b.trim.sub

L83-1-a.trim.fq: L83-1-a.fq.gz
		$(TIME) -o L83-1-a.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog L83-1-a.trimlog -threads 4 L83-1-a.fq.gz L83-1-a.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

L83-1-a.trim.fq-pbs:
	qsub -l $(TRIMRES) L83-1-a.trim.sub

L83-1-b.trim.fq: L83-1-b.fq.gz
		$(TIME) -o L83-1-b.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog L83-1-b.trimlog -threads 4 L83-1-b.fq.gz L83-1-b.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

L83-1-b.trim.fq-pbs:
	qsub -l $(TRIMRES) L83-1-b.trim.sub

L83-a.trim.fq: L83-a.fq.gz
		$(TIME) -o L83-a.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog L83-a.trimlog -threads 4 L83-a.fq.gz L83-a.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

L83-a.trim.fq-pbs:
	qsub -l $(TRIMRES) L83-a.trim.sub

L83-b.trim.fq: L83-b.fq.gz
		$(TIME) -o L83-b.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog L83-b.trimlog -threads 4 L83-b.fq.gz L83-b.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

L83-b.trim.fq-pbs:
	qsub -l $(TRIMRES) L83-b.trim.sub

LI.trim.fq: LI.fq.gz
		$(TIME) -o LI.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog LI.trimlog -threads 4 LI.fq.gz LI.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

LI.trim.fq-pbs:
	qsub -l $(TRIMRES) LI.trim.sub

Lips.trim.fq: Lips.fq.gz
		$(TIME) -o Lips.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Lips.trimlog -threads 4 Lips.fq.gz Lips.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Lips.trim.fq-pbs:
	qsub -l $(TRIMRES) Lips.trim.sub

LK.trim.fq: LK.fq.gz
		$(TIME) -o LK.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog LK.trimlog -threads 4 LK.fq.gz LK.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

LK.trim.fq-pbs:
	qsub -l $(TRIMRES) LK.trim.sub

M0_liver.trim.fq: M0_liver.fq.gz
		$(TIME) -o M0_liver.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog M0_liver.trimlog -threads 4 M0_liver.fq.gz M0_liver.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

M0_liver.trim.fq-pbs:
	qsub -l $(TRIMRES) M0_liver.trim.sub

M2_liver.trim.fq: M2_liver.fq.gz
		$(TIME) -o M2_liver.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog M2_liver.trimlog -threads 4 M2_liver.fq.gz M2_liver.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

M2_liver.trim.fq-pbs:
	qsub -l $(TRIMRES) M2_liver.trim.sub

M5_liver.trim.fq: M5_liver.fq.gz
		$(TIME) -o M5_liver.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog M5_liver.trimlog -threads 4 M5_liver.fq.gz M5_liver.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

M5_liver.trim.fq-pbs:
	qsub -l $(TRIMRES) M5_liver.trim.sub

Monocytes.trim.fq: Monocytes.fq.gz
		$(TIME) -o Monocytes.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Monocytes.trimlog -threads 4 Monocytes.fq.gz Monocytes.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Monocytes.trim.fq-pbs:
	qsub -l $(TRIMRES) Monocytes.trim.sub

mssm_single_s1.trim.fq: mssm_single_s1.trim.fq.gz
		$(TIME) -o mssm_single_s1.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog mssm_single_s1.trimlog -threads 4 mssm_single_s1.trim.fq.gz mssm_single_s1.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

mssm_single_s1.trim.fq-pbs:
	qsub -l $(TRIMRES) mssm_single_s1.trim.sub

mssm_single_s2.trim.fq: mssm_single_s2.trim.fq.gz
		$(TIME) -o mssm_single_s2.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog mssm_single_s2.trimlog -threads 4 mssm_single_s2.trim.fq.gz mssm_single_s2.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

mssm_single_s2.trim.fq-pbs:
	qsub -l $(TRIMRES) mssm_single_s2.trim.sub

mssm_single_s3.trim.fq: mssm_single_s3.trim.fq.gz
		$(TIME) -o mssm_single_s3.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog mssm_single_s3.trimlog -threads 4 mssm_single_s3.trim.fq.gz mssm_single_s3.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

mssm_single_s3.trim.fq-pbs:
	qsub -l $(TRIMRES) mssm_single_s3.trim.sub

mssm_single_s5.trim.fq: mssm_single_s5.trim.fq.gz
		$(TIME) -o mssm_single_s5.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog mssm_single_s5.trimlog -threads 4 mssm_single_s5.trim.fq.gz mssm_single_s5.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

mssm_single_s5.trim.fq-pbs:
	qsub -l $(TRIMRES) mssm_single_s5.trim.sub

mssm_single_s6.trim.fq: mssm_single_s6.trim.fq.gz
		$(TIME) -o mssm_single_s6.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog mssm_single_s6.trimlog -threads 4 mssm_single_s6.trim.fq.gz mssm_single_s6.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

mssm_single_s6.trim.fq-pbs:
	qsub -l $(TRIMRES) mssm_single_s6.trim.sub

mssm_single_s7.trim.fq: mssm_single_s7.trim.fq.gz
		$(TIME) -o mssm_single_s7.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog mssm_single_s7.trimlog -threads 4 mssm_single_s7.trim.fq.gz mssm_single_s7.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

mssm_single_s7.trim.fq-pbs:
	qsub -l $(TRIMRES) mssm_single_s7.trim.sub

Muscle.trim.fq: Muscle.fq.gz
		$(TIME) -o Muscle.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Muscle.trimlog -threads 4 Muscle.fq.gz Muscle.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Muscle.trim.fq-pbs:
	qsub -l $(TRIMRES) Muscle.trim.sub

Neutrophils.trim.fq: Neutrophils.fq.gz
		$(TIME) -o Neutrophils.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Neutrophils.trimlog -threads 4 Neutrophils.fq.gz Neutrophils.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Neutrophils.trim.fq-pbs:
	qsub -l $(TRIMRES) Neutrophils.trim.sub

OFHS.trim.fq: OFHS.fq.gz
		$(TIME) -o OFHS.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog OFHS.trimlog -threads 4 OFHS.fq.gz OFHS.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

OFHS.trim.fq-pbs:
	qsub -l $(TRIMRES) OFHS.trim.sub

POF-eye.trim.fq: POF-eye.fq.gz
		$(TIME) -o POF-eye.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog POF-eye.trimlog -threads 4 POF-eye.fq.gz POF-eye.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

POF-eye.trim.fq-pbs:
	qsub -l $(TRIMRES) POF-eye.trim.sub

POF-ts.trim.fq: POF-ts.fq.gz
		$(TIME) -o POF-ts.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog POF-ts.trimlog -threads 4 POF-ts.fq.gz POF-ts.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

POF-ts.trim.fq-pbs:
	qsub -l $(TRIMRES) POF-ts.trim.sub

PSMGill.trim.fq: PSMGill.fq.gz
		$(TIME) -o PSMGill.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog PSMGill.trimlog -threads 4 PSMGill.fq.gz PSMGill.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

PSMGill.trim.fq-pbs:
	qsub -l $(TRIMRES) PSMGill.trim.sub

Rope.trim.fq: Rope.fq.gz
		$(TIME) -o Rope.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Rope.trimlog -threads 4 Rope.fq.gz Rope.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Rope.trim.fq-pbs:
	qsub -l $(TRIMRES) Rope.trim.sub

SMGill.trim.fq: SMGill.fq.gz
		$(TIME) -o SMGill.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog SMGill.trimlog -threads 4 SMGill.fq.gz SMGill.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

SMGill.trim.fq-pbs:
	qsub -l $(TRIMRES) SMGill.trim.sub

SMHS.trim.fq: SMHS.fq.gz
		$(TIME) -o SMHS.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog SMHS.trimlog -threads 4 SMHS.fq.gz SMHS.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

SMHS.trim.fq-pbs:
	qsub -l $(TRIMRES) SMHS.trim.sub

SN.trim.fq: SN.fq.gz
		$(TIME) -o SN.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog SN.trimlog -threads 4 SN.fq.gz SN.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

SN.trim.fq-pbs:
	qsub -l $(TRIMRES) SN.trim.sub

SPDI.trim.fq: SPDI.fq.gz
		$(TIME) -o SPDI.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog SPDI.trimlog -threads 4 SPDI.fq.gz SPDI.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

SPDI.trim.fq-pbs:
	qsub -l $(TRIMRES) SPDI.trim.sub

spinal_cord_paired.trim.fq: spinal_cord_paired.qc.fq.gz
		$(TIME) -o spinal_cord_paired.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog spinal_cord_paired.trimlog -threads 4 spinal_cord_paired.qc.fq.gz spinal_cord_paired.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

spinal_cord_paired.trim.fq-pbs:
	qsub -l $(TRIMRES) spinal_cord_paired.trim.sub

SPK.trim.fq: SPK.fq.gz
		$(TIME) -o SPK.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog SPK.trimlog -threads 4 SPK.fq.gz SPK.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

SPK.trim.fq-pbs:
	qsub -l $(TRIMRES) SPK.trim.sub

SPPI.trim.fq: SPPI.fq.gz
		$(TIME) -o SPPI.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog SPPI.trimlog -threads 4 SPPI.fq.gz SPPI.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

SPPI.trim.fq-pbs:
	qsub -l $(TRIMRES) SPPI.trim.sub

Stg-1-Intestine.trim.fq: Stg-1-Intestine.fq.gz
		$(TIME) -o Stg-1-Intestine.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-1-Intestine.trimlog -threads 4 Stg-1-Intestine.fq.gz Stg-1-Intestine.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-1-Intestine.trim.fq-pbs:
	qsub -l $(TRIMRES) Stg-1-Intestine.trim.sub

Stg-1-Kidney.trim.fq: Stg-1-Kidney.fq.gz
		$(TIME) -o Stg-1-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-1-Kidney.trimlog -threads 4 Stg-1-Kidney.fq.gz Stg-1-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-1-Kidney.trim.fq-pbs:
	qsub -l $(TRIMRES) Stg-1-Kidney.trim.sub

Stg-2-Intestine.trim.fq: Stg-2-Intestine.fq.gz
		$(TIME) -o Stg-2-Intestine.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-2-Intestine.trimlog -threads 4 Stg-2-Intestine.fq.gz Stg-2-Intestine.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-2-Intestine.trim.fq-pbs:
	qsub -l $(TRIMRES) Stg-2-Intestine.trim.sub

Stg-2-Kidney.trim.fq: Stg-2-Kidney.fq.gz
		$(TIME) -o Stg-2-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-2-Kidney.trimlog -threads 4 Stg-2-Kidney.fq.gz Stg-2-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-2-Kidney.trim.fq-pbs:
	qsub -l $(TRIMRES) Stg-2-Kidney.trim.sub

Stg-3-Intestine.trim.fq: Stg-3-Intestine.fq.gz
		$(TIME) -o Stg-3-Intestine.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-3-Intestine.trimlog -threads 4 Stg-3-Intestine.fq.gz Stg-3-Intestine.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-3-Intestine.trim.fq-pbs:
	qsub -l $(TRIMRES) Stg-3-Intestine.trim.sub

Stg-3-Kidney.trim.fq: Stg-3-Kidney.fq.gz
		$(TIME) -o Stg-3-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-3-Kidney.trimlog -threads 4 Stg-3-Kidney.fq.gz Stg-3-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-3-Kidney.trim.fq-pbs:
	qsub -l $(TRIMRES) Stg-3-Kidney.trim.sub

Stg-4-Intestine.trim.fq: Stg-4-Intestine.fq.gz
		$(TIME) -o Stg-4-Intestine.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-4-Intestine.trimlog -threads 4 Stg-4-Intestine.fq.gz Stg-4-Intestine.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-4-Intestine.trim.fq-pbs:
	qsub -l $(TRIMRES) Stg-4-Intestine.trim.sub

Stg-4-Kidney.trim.fq: Stg-4-Kidney.fq.gz
		$(TIME) -o Stg-4-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-4-Kidney.trimlog -threads 4 Stg-4-Kidney.fq.gz Stg-4-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-4-Kidney.trim.fq-pbs:
	qsub -l $(TRIMRES) Stg-4-Kidney.trim.sub

Stg-5-Intestine.trim.fq: Stg-5-Intestine.fq.gz
		$(TIME) -o Stg-5-Intestine.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-5-Intestine.trimlog -threads 4 Stg-5-Intestine.fq.gz Stg-5-Intestine.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-5-Intestine.trim.fq-pbs:
	qsub -l $(TRIMRES) Stg-5-Intestine.trim.sub

Stg-5-Kidney.trim.fq: Stg-5-Kidney.fq.gz
		$(TIME) -o Stg-5-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-5-Kidney.trimlog -threads 4 Stg-5-Kidney.fq.gz Stg-5-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-5-Kidney.trim.fq-pbs:
	qsub -l $(TRIMRES) Stg-5-Kidney.trim.sub

Stg-6-Intestine.trim.fq: Stg-6-Intestine.fq.gz
		$(TIME) -o Stg-6-Intestine.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-6-Intestine.trimlog -threads 4 Stg-6-Intestine.fq.gz Stg-6-Intestine.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-6-Intestine.trim.fq-pbs:
	qsub -l $(TRIMRES) Stg-6-Intestine.trim.sub

Stg-6-Kidney.trim.fq: Stg-6-Kidney.fq.gz
		$(TIME) -o Stg-6-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-6-Kidney.trimlog -threads 4 Stg-6-Kidney.fq.gz Stg-6-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-6-Kidney.trim.fq-pbs:
	qsub -l $(TRIMRES) Stg-6-Kidney.trim.sub

Stg-7-Intestine.trim.fq: Stg-7-Intestine.fq.gz
		$(TIME) -o Stg-7-Intestine.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-7-Intestine.trimlog -threads 4 Stg-7-Intestine.fq.gz Stg-7-Intestine.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-7-Intestine.trim.fq-pbs:
	qsub -l $(TRIMRES) Stg-7-Intestine.trim.sub

Stg-7-Kidney.trim.fq: Stg-7-Kidney.fq.gz
		$(TIME) -o Stg-7-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog Stg-7-Kidney.trimlog -threads 4 Stg-7-Kidney.fq.gz Stg-7-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-7-Kidney.trim.fq-pbs:
	qsub -l $(TRIMRES) Stg-7-Kidney.trim.sub

SW2-Gill.trim.fq: SW2-Gill.fq.gz
		$(TIME) -o SW2-Gill.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog SW2-Gill.trimlog -threads 4 SW2-Gill.fq.gz SW2-Gill.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

SW2-Gill.trim.fq-pbs:
	qsub -l $(TRIMRES) SW2-Gill.trim.sub

SW2-Intestine-a.trim.fq: SW2-Intestine-a.fq.gz
		$(TIME) -o SW2-Intestine-a.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog SW2-Intestine-a.trimlog -threads 4 SW2-Intestine-a.fq.gz SW2-Intestine-a.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

SW2-Intestine-a.trim.fq-pbs:
	qsub -l $(TRIMRES) SW2-Intestine-a.trim.sub

SW2-Intestine-b.trim.fq: SW2-Intestine-b.fq.gz
		$(TIME) -o SW2-Intestine-b.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog SW2-Intestine-b.trimlog -threads 4 SW2-Intestine-b.fq.gz SW2-Intestine-b.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

SW2-Intestine-b.trim.fq-pbs:
	qsub -l $(TRIMRES) SW2-Intestine-b.trim.sub

SW2-Kidney.trim.fq: SW2-Kidney.fq.gz
		$(TIME) -o SW2-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -trimlog SW2-Kidney.trimlog -threads 4 SW2-Kidney.fq.gz SW2-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

SW2-Kidney.trim.fq-pbs:
	qsub -l $(TRIMRES) SW2-Kidney.trim.sub

trim-pbs: 614GMAAXX_1_1_pf.trim.fq-pbs 614GMAAXX_2_1_pf.trim.fq-pbs 614GMAAXX_3_1_pf.trim.fq-pbs 614GMAAXX_4_1_pf.trim.fq-pbs 614GMAAXX_6_1_pf.trim.fq-pbs 614GMAAXX_7_1_pf.trim.fq-pbs 614GMAAXX_8_1_pf.trim.fq-pbs AB057JABXX_s_7_pe.trim.fq-pbs AI.trim.fq-pbs AK.trim.fq-pbs BB0582ABXX_s_3_pe.trim.fq-pbs BB0582ABXX_s_4_pe.trim.fq-pbs BB0582ABXX_s_5_pe.trim.fq-pbs BB0582ABXX_s_6_pe.trim.fq-pbs BB0582ABXX_s_7_pe.trim.fq-pbs BB0582ABXX_s_8_pe.trim.fq-pbs brain_paired.trim.fq-pbs FW1-Gill.trim.fq-pbs FW1-Intestine.trim.fq-pbs FW1-Kidney.trim.fq-pbs G81-a.trim.fq-pbs G81-b.trim.fq-pbs G82-a.trim.fq-pbs G82-b.trim.fq-pbs G83-a.trim.fq-pbs G83-b.trim.fq-pbs Gill.trim.fq-pbs JV1_liver.trim.fq-pbs JVa_liver.trim.fq-pbs JVb_liver.trim.fq-pbs JV_intest.trim.fq-pbs Kidney.trim.fq-pbs L81-a.trim.fq-pbs L81-b.trim.fq-pbs L82-a.trim.fq-pbs L82-b.trim.fq-pbs L83-1-a.trim.fq-pbs L83-1-b.trim.fq-pbs L83-a.trim.fq-pbs L83-b.trim.fq-pbs LI.trim.fq-pbs Lips.trim.fq-pbs LK.trim.fq-pbs M0_liver.trim.fq-pbs M2_liver.trim.fq-pbs M5_liver.trim.fq-pbs Monocytes.trim.fq-pbs mssm_single_s1.trim.fq-pbs mssm_single_s2.trim.fq-pbs mssm_single_s3.trim.fq-pbs mssm_single_s5.trim.fq-pbs mssm_single_s6.trim.fq-pbs mssm_single_s7.trim.fq-pbs Muscle.trim.fq-pbs Neutrophils.trim.fq-pbs OFHS.trim.fq-pbs POF-eye.trim.fq-pbs POF-ts.trim.fq-pbs PSMGill.trim.fq-pbs Rope.trim.fq-pbs SMGill.trim.fq-pbs SMHS.trim.fq-pbs SN.trim.fq-pbs SPDI.trim.fq-pbs spinal_cord_paired.trim.fq-pbs SPK.trim.fq-pbs SPPI.trim.fq-pbs Stg-1-Intestine.trim.fq-pbs Stg-1-Kidney.trim.fq-pbs Stg-2-Intestine.trim.fq-pbs Stg-2-Kidney.trim.fq-pbs Stg-3-Intestine.trim.fq-pbs Stg-3-Kidney.trim.fq-pbs Stg-4-Intestine.trim.fq-pbs Stg-4-Kidney.trim.fq-pbs Stg-5-Intestine.trim.fq-pbs Stg-5-Kidney.trim.fq-pbs Stg-6-Intestine.trim.fq-pbs Stg-6-Kidney.trim.fq-pbs Stg-7-Intestine.trim.fq-pbs Stg-7-Kidney.trim.fq-pbs SW2-Gill.trim.fq-pbs SW2-Intestine-a.trim.fq-pbs SW2-Intestine-b.trim.fq-pbs SW2-Kidney.trim.fq-pbs
614GMAAXX_1_1_pf.trim.fq.fastqc2: 614GMAAXX_1_1_pf.trim.fq
	fastqc --noextract 614GMAAXX_1_1_pf.trim.fq

614GMAAXX_1_1_pf.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) 614GMAAXX_1_1_pf.trim.fq.fastqc2.sub

614GMAAXX_2_1_pf.trim.fq.fastqc2: 614GMAAXX_2_1_pf.trim.fq
	fastqc --noextract 614GMAAXX_2_1_pf.trim.fq

614GMAAXX_2_1_pf.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) 614GMAAXX_2_1_pf.trim.fq.fastqc2.sub

614GMAAXX_3_1_pf.trim.fq.fastqc2: 614GMAAXX_3_1_pf.trim.fq
	fastqc --noextract 614GMAAXX_3_1_pf.trim.fq

614GMAAXX_3_1_pf.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) 614GMAAXX_3_1_pf.trim.fq.fastqc2.sub

614GMAAXX_4_1_pf.trim.fq.fastqc2: 614GMAAXX_4_1_pf.trim.fq
	fastqc --noextract 614GMAAXX_4_1_pf.trim.fq

614GMAAXX_4_1_pf.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) 614GMAAXX_4_1_pf.trim.fq.fastqc2.sub

614GMAAXX_6_1_pf.trim.fq.fastqc2: 614GMAAXX_6_1_pf.trim.fq
	fastqc --noextract 614GMAAXX_6_1_pf.trim.fq

614GMAAXX_6_1_pf.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) 614GMAAXX_6_1_pf.trim.fq.fastqc2.sub

614GMAAXX_7_1_pf.trim.fq.fastqc2: 614GMAAXX_7_1_pf.trim.fq
	fastqc --noextract 614GMAAXX_7_1_pf.trim.fq

614GMAAXX_7_1_pf.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) 614GMAAXX_7_1_pf.trim.fq.fastqc2.sub

614GMAAXX_8_1_pf.trim.fq.fastqc2: 614GMAAXX_8_1_pf.trim.fq
	fastqc --noextract 614GMAAXX_8_1_pf.trim.fq

614GMAAXX_8_1_pf.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) 614GMAAXX_8_1_pf.trim.fq.fastqc2.sub

AB057JABXX_s_7_pe.trim.fq.fastqc2: AB057JABXX_s_7_pe.trim.fq
	fastqc --noextract AB057JABXX_s_7_pe.trim.fq

AB057JABXX_s_7_pe.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) AB057JABXX_s_7_pe.trim.fq.fastqc2.sub

AI.trim.fq.fastqc2: AI.trim.fq
	fastqc --noextract AI.trim.fq

AI.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) AI.trim.fq.fastqc2.sub

AK.trim.fq.fastqc2: AK.trim.fq
	fastqc --noextract AK.trim.fq

AK.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) AK.trim.fq.fastqc2.sub

BB0582ABXX_s_3_pe.trim.fq.fastqc2: BB0582ABXX_s_3_pe.trim.fq
	fastqc --noextract BB0582ABXX_s_3_pe.trim.fq

BB0582ABXX_s_3_pe.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) BB0582ABXX_s_3_pe.trim.fq.fastqc2.sub

BB0582ABXX_s_4_pe.trim.fq.fastqc2: BB0582ABXX_s_4_pe.trim.fq
	fastqc --noextract BB0582ABXX_s_4_pe.trim.fq

BB0582ABXX_s_4_pe.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) BB0582ABXX_s_4_pe.trim.fq.fastqc2.sub

BB0582ABXX_s_5_pe.trim.fq.fastqc2: BB0582ABXX_s_5_pe.trim.fq
	fastqc --noextract BB0582ABXX_s_5_pe.trim.fq

BB0582ABXX_s_5_pe.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) BB0582ABXX_s_5_pe.trim.fq.fastqc2.sub

BB0582ABXX_s_6_pe.trim.fq.fastqc2: BB0582ABXX_s_6_pe.trim.fq
	fastqc --noextract BB0582ABXX_s_6_pe.trim.fq

BB0582ABXX_s_6_pe.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) BB0582ABXX_s_6_pe.trim.fq.fastqc2.sub

BB0582ABXX_s_7_pe.trim.fq.fastqc2: BB0582ABXX_s_7_pe.trim.fq
	fastqc --noextract BB0582ABXX_s_7_pe.trim.fq

BB0582ABXX_s_7_pe.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) BB0582ABXX_s_7_pe.trim.fq.fastqc2.sub

BB0582ABXX_s_8_pe.trim.fq.fastqc2: BB0582ABXX_s_8_pe.trim.fq
	fastqc --noextract BB0582ABXX_s_8_pe.trim.fq

BB0582ABXX_s_8_pe.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) BB0582ABXX_s_8_pe.trim.fq.fastqc2.sub

brain_paired.trim.fq.fastqc2: brain_paired.trim.fq
	fastqc --noextract brain_paired.trim.fq

brain_paired.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) brain_paired.trim.fq.fastqc2.sub

FW1.fastqc2: FW1
	fastqc --noextract FW1

FW1.fastqc2-pbs:
	qsub -l $(QCRES) FW1.fastqc2.sub

FW1.fastqc2: FW1
	fastqc --noextract FW1

FW1.fastqc2-pbs:
	qsub -l $(QCRES) FW1.fastqc2.sub

FW1.fastqc2: FW1
	fastqc --noextract FW1

FW1.fastqc2-pbs:
	qsub -l $(QCRES) FW1.fastqc2.sub

G81.fastqc2: G81
	fastqc --noextract G81

G81.fastqc2-pbs:
	qsub -l $(QCRES) G81.fastqc2.sub

G81.fastqc2: G81
	fastqc --noextract G81

G81.fastqc2-pbs:
	qsub -l $(QCRES) G81.fastqc2.sub

G82.fastqc2: G82
	fastqc --noextract G82

G82.fastqc2-pbs:
	qsub -l $(QCRES) G82.fastqc2.sub

G82.fastqc2: G82
	fastqc --noextract G82

G82.fastqc2-pbs:
	qsub -l $(QCRES) G82.fastqc2.sub

G83.fastqc2: G83
	fastqc --noextract G83

G83.fastqc2-pbs:
	qsub -l $(QCRES) G83.fastqc2.sub

G83.fastqc2: G83
	fastqc --noextract G83

G83.fastqc2-pbs:
	qsub -l $(QCRES) G83.fastqc2.sub

Gill.trim.fq.fastqc2: Gill.trim.fq
	fastqc --noextract Gill.trim.fq

Gill.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) Gill.trim.fq.fastqc2.sub

JV1_liver.trim.fq.fastqc2: JV1_liver.trim.fq
	fastqc --noextract JV1_liver.trim.fq

JV1_liver.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) JV1_liver.trim.fq.fastqc2.sub

JVa_liver.trim.fq.fastqc2: JVa_liver.trim.fq
	fastqc --noextract JVa_liver.trim.fq

JVa_liver.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) JVa_liver.trim.fq.fastqc2.sub

JVb_liver.trim.fq.fastqc2: JVb_liver.trim.fq
	fastqc --noextract JVb_liver.trim.fq

JVb_liver.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) JVb_liver.trim.fq.fastqc2.sub

JV_intest.trim.fq.fastqc2: JV_intest.trim.fq
	fastqc --noextract JV_intest.trim.fq

JV_intest.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) JV_intest.trim.fq.fastqc2.sub

Kidney.trim.fq.fastqc2: Kidney.trim.fq
	fastqc --noextract Kidney.trim.fq

Kidney.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) Kidney.trim.fq.fastqc2.sub

L81.fastqc2: L81
	fastqc --noextract L81

L81.fastqc2-pbs:
	qsub -l $(QCRES) L81.fastqc2.sub

L81.fastqc2: L81
	fastqc --noextract L81

L81.fastqc2-pbs:
	qsub -l $(QCRES) L81.fastqc2.sub

L82.fastqc2: L82
	fastqc --noextract L82

L82.fastqc2-pbs:
	qsub -l $(QCRES) L82.fastqc2.sub

L82.fastqc2: L82
	fastqc --noextract L82

L82.fastqc2-pbs:
	qsub -l $(QCRES) L82.fastqc2.sub

L83.fastqc2: L83
	fastqc --noextract L83

L83.fastqc2-pbs:
	qsub -l $(QCRES) L83.fastqc2.sub

L83.fastqc2: L83
	fastqc --noextract L83

L83.fastqc2-pbs:
	qsub -l $(QCRES) L83.fastqc2.sub

L83.fastqc2: L83
	fastqc --noextract L83

L83.fastqc2-pbs:
	qsub -l $(QCRES) L83.fastqc2.sub

L83.fastqc2: L83
	fastqc --noextract L83

L83.fastqc2-pbs:
	qsub -l $(QCRES) L83.fastqc2.sub

LI.trim.fq.fastqc2: LI.trim.fq
	fastqc --noextract LI.trim.fq

LI.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) LI.trim.fq.fastqc2.sub

Lips.trim.fq.fastqc2: Lips.trim.fq
	fastqc --noextract Lips.trim.fq

Lips.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) Lips.trim.fq.fastqc2.sub

LK.trim.fq.fastqc2: LK.trim.fq
	fastqc --noextract LK.trim.fq

LK.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) LK.trim.fq.fastqc2.sub

M0_liver.trim.fq.fastqc2: M0_liver.trim.fq
	fastqc --noextract M0_liver.trim.fq

M0_liver.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) M0_liver.trim.fq.fastqc2.sub

M2_liver.trim.fq.fastqc2: M2_liver.trim.fq
	fastqc --noextract M2_liver.trim.fq

M2_liver.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) M2_liver.trim.fq.fastqc2.sub

M5_liver.trim.fq.fastqc2: M5_liver.trim.fq
	fastqc --noextract M5_liver.trim.fq

M5_liver.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) M5_liver.trim.fq.fastqc2.sub

Monocytes.trim.fq.fastqc2: Monocytes.trim.fq
	fastqc --noextract Monocytes.trim.fq

Monocytes.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) Monocytes.trim.fq.fastqc2.sub

mssm_single_s1.trim.fq.fastqc2: mssm_single_s1.trim.fq
	fastqc --noextract mssm_single_s1.trim.fq

mssm_single_s1.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) mssm_single_s1.trim.fq.fastqc2.sub

mssm_single_s2.trim.fq.fastqc2: mssm_single_s2.trim.fq
	fastqc --noextract mssm_single_s2.trim.fq

mssm_single_s2.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) mssm_single_s2.trim.fq.fastqc2.sub

mssm_single_s3.trim.fq.fastqc2: mssm_single_s3.trim.fq
	fastqc --noextract mssm_single_s3.trim.fq

mssm_single_s3.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) mssm_single_s3.trim.fq.fastqc2.sub

mssm_single_s5.trim.fq.fastqc2: mssm_single_s5.trim.fq
	fastqc --noextract mssm_single_s5.trim.fq

mssm_single_s5.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) mssm_single_s5.trim.fq.fastqc2.sub

mssm_single_s6.trim.fq.fastqc2: mssm_single_s6.trim.fq
	fastqc --noextract mssm_single_s6.trim.fq

mssm_single_s6.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) mssm_single_s6.trim.fq.fastqc2.sub

mssm_single_s7.trim.fq.fastqc2: mssm_single_s7.trim.fq
	fastqc --noextract mssm_single_s7.trim.fq

mssm_single_s7.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) mssm_single_s7.trim.fq.fastqc2.sub

Muscle.trim.fq.fastqc2: Muscle.trim.fq
	fastqc --noextract Muscle.trim.fq

Muscle.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) Muscle.trim.fq.fastqc2.sub

Neutrophils.trim.fq.fastqc2: Neutrophils.trim.fq
	fastqc --noextract Neutrophils.trim.fq

Neutrophils.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) Neutrophils.trim.fq.fastqc2.sub

OFHS.trim.fq.fastqc2: OFHS.trim.fq
	fastqc --noextract OFHS.trim.fq

OFHS.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) OFHS.trim.fq.fastqc2.sub

POF.fastqc2: POF
	fastqc --noextract POF

POF.fastqc2-pbs:
	qsub -l $(QCRES) POF.fastqc2.sub

POF.fastqc2: POF
	fastqc --noextract POF

POF.fastqc2-pbs:
	qsub -l $(QCRES) POF.fastqc2.sub

PSMGill.trim.fq.fastqc2: PSMGill.trim.fq
	fastqc --noextract PSMGill.trim.fq

PSMGill.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) PSMGill.trim.fq.fastqc2.sub

Rope.trim.fq.fastqc2: Rope.trim.fq
	fastqc --noextract Rope.trim.fq

Rope.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) Rope.trim.fq.fastqc2.sub

SMGill.trim.fq.fastqc2: SMGill.trim.fq
	fastqc --noextract SMGill.trim.fq

SMGill.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) SMGill.trim.fq.fastqc2.sub

SMHS.trim.fq.fastqc2: SMHS.trim.fq
	fastqc --noextract SMHS.trim.fq

SMHS.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) SMHS.trim.fq.fastqc2.sub

SN.trim.fq.fastqc2: SN.trim.fq
	fastqc --noextract SN.trim.fq

SN.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) SN.trim.fq.fastqc2.sub

SPDI.trim.fq.fastqc2: SPDI.trim.fq
	fastqc --noextract SPDI.trim.fq

SPDI.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) SPDI.trim.fq.fastqc2.sub

spinal_cord_paired.trim.fq.fastqc2: spinal_cord_paired.trim.fq
	fastqc --noextract spinal_cord_paired.trim.fq

spinal_cord_paired.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) spinal_cord_paired.trim.fq.fastqc2.sub

SPK.trim.fq.fastqc2: SPK.trim.fq
	fastqc --noextract SPK.trim.fq

SPK.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) SPK.trim.fq.fastqc2.sub

SPPI.trim.fq.fastqc2: SPPI.trim.fq
	fastqc --noextract SPPI.trim.fq

SPPI.trim.fq.fastqc2-pbs:
	qsub -l $(QCRES) SPPI.trim.fq.fastqc2.sub

Stg.fastqc2: Stg
	fastqc --noextract Stg

Stg.fastqc2-pbs:
	qsub -l $(QCRES) Stg.fastqc2.sub

Stg.fastqc2: Stg
	fastqc --noextract Stg

Stg.fastqc2-pbs:
	qsub -l $(QCRES) Stg.fastqc2.sub

Stg.fastqc2: Stg
	fastqc --noextract Stg

Stg.fastqc2-pbs:
	qsub -l $(QCRES) Stg.fastqc2.sub

Stg.fastqc2: Stg
	fastqc --noextract Stg

Stg.fastqc2-pbs:
	qsub -l $(QCRES) Stg.fastqc2.sub

Stg.fastqc2: Stg
	fastqc --noextract Stg

Stg.fastqc2-pbs:
	qsub -l $(QCRES) Stg.fastqc2.sub

Stg.fastqc2: Stg
	fastqc --noextract Stg

Stg.fastqc2-pbs:
	qsub -l $(QCRES) Stg.fastqc2.sub

Stg.fastqc2: Stg
	fastqc --noextract Stg

Stg.fastqc2-pbs:
	qsub -l $(QCRES) Stg.fastqc2.sub

Stg.fastqc2: Stg
	fastqc --noextract Stg

Stg.fastqc2-pbs:
	qsub -l $(QCRES) Stg.fastqc2.sub

Stg.fastqc2: Stg
	fastqc --noextract Stg

Stg.fastqc2-pbs:
	qsub -l $(QCRES) Stg.fastqc2.sub

Stg.fastqc2: Stg
	fastqc --noextract Stg

Stg.fastqc2-pbs:
	qsub -l $(QCRES) Stg.fastqc2.sub

Stg.fastqc2: Stg
	fastqc --noextract Stg

Stg.fastqc2-pbs:
	qsub -l $(QCRES) Stg.fastqc2.sub

Stg.fastqc2: Stg
	fastqc --noextract Stg

Stg.fastqc2-pbs:
	qsub -l $(QCRES) Stg.fastqc2.sub

Stg.fastqc2: Stg
	fastqc --noextract Stg

Stg.fastqc2-pbs:
	qsub -l $(QCRES) Stg.fastqc2.sub

Stg.fastqc2: Stg
	fastqc --noextract Stg

Stg.fastqc2-pbs:
	qsub -l $(QCRES) Stg.fastqc2.sub

SW2.fastqc2: SW2
	fastqc --noextract SW2

SW2.fastqc2-pbs:
	qsub -l $(QCRES) SW2.fastqc2.sub

SW2.fastqc2: SW2
	fastqc --noextract SW2

SW2.fastqc2-pbs:
	qsub -l $(QCRES) SW2.fastqc2.sub

SW2.fastqc2: SW2
	fastqc --noextract SW2

SW2.fastqc2-pbs:
	qsub -l $(QCRES) SW2.fastqc2.sub

SW2.fastqc2: SW2
	fastqc --noextract SW2

SW2.fastqc2-pbs:
	qsub -l $(QCRES) SW2.fastqc2.sub

fastqc2-pbs: 614GMAAXX_1_1_pf.trim.fq.fastqc2-pbs 614GMAAXX_2_1_pf.trim.fq.fastqc2-pbs 614GMAAXX_3_1_pf.trim.fq.fastqc2-pbs 614GMAAXX_4_1_pf.trim.fq.fastqc2-pbs 614GMAAXX_6_1_pf.trim.fq.fastqc2-pbs 614GMAAXX_7_1_pf.trim.fq.fastqc2-pbs 614GMAAXX_8_1_pf.trim.fq.fastqc2-pbs AB057JABXX_s_7_pe.trim.fq.fastqc2-pbs AI.trim.fq.fastqc2-pbs AK.trim.fq.fastqc2-pbs BB0582ABXX_s_3_pe.trim.fq.fastqc2-pbs BB0582ABXX_s_4_pe.trim.fq.fastqc2-pbs BB0582ABXX_s_5_pe.trim.fq.fastqc2-pbs BB0582ABXX_s_6_pe.trim.fq.fastqc2-pbs BB0582ABXX_s_7_pe.trim.fq.fastqc2-pbs BB0582ABXX_s_8_pe.trim.fq.fastqc2-pbs brain_paired.trim.fq.fastqc2-pbs FW1.fastqc2-pbs FW1.fastqc2-pbs FW1.fastqc2-pbs G81.fastqc2-pbs G81.fastqc2-pbs G82.fastqc2-pbs G82.fastqc2-pbs G83.fastqc2-pbs G83.fastqc2-pbs Gill.trim.fq.fastqc2-pbs JV1_liver.trim.fq.fastqc2-pbs JVa_liver.trim.fq.fastqc2-pbs JVb_liver.trim.fq.fastqc2-pbs JV_intest.trim.fq.fastqc2-pbs Kidney.trim.fq.fastqc2-pbs L81.fastqc2-pbs L81.fastqc2-pbs L82.fastqc2-pbs L82.fastqc2-pbs L83.fastqc2-pbs L83.fastqc2-pbs L83.fastqc2-pbs L83.fastqc2-pbs LI.trim.fq.fastqc2-pbs Lips.trim.fq.fastqc2-pbs LK.trim.fq.fastqc2-pbs M0_liver.trim.fq.fastqc2-pbs M2_liver.trim.fq.fastqc2-pbs M5_liver.trim.fq.fastqc2-pbs Monocytes.trim.fq.fastqc2-pbs mssm_single_s1.trim.fq.fastqc2-pbs mssm_single_s2.trim.fq.fastqc2-pbs mssm_single_s3.trim.fq.fastqc2-pbs mssm_single_s5.trim.fq.fastqc2-pbs mssm_single_s6.trim.fq.fastqc2-pbs mssm_single_s7.trim.fq.fastqc2-pbs Muscle.trim.fq.fastqc2-pbs Neutrophils.trim.fq.fastqc2-pbs OFHS.trim.fq.fastqc2-pbs POF.fastqc2-pbs POF.fastqc2-pbs PSMGill.trim.fq.fastqc2-pbs Rope.trim.fq.fastqc2-pbs SMGill.trim.fq.fastqc2-pbs SMHS.trim.fq.fastqc2-pbs SN.trim.fq.fastqc2-pbs SPDI.trim.fq.fastqc2-pbs spinal_cord_paired.trim.fq.fastqc2-pbs SPK.trim.fq.fastqc2-pbs SPPI.trim.fq.fastqc2-pbs Stg.fastqc2-pbs Stg.fastqc2-pbs Stg.fastqc2-pbs Stg.fastqc2-pbs Stg.fastqc2-pbs Stg.fastqc2-pbs Stg.fastqc2-pbs Stg.fastqc2-pbs Stg.fastqc2-pbs Stg.fastqc2-pbs Stg.fastqc2-pbs Stg.fastqc2-pbs Stg.fastqc2-pbs Stg.fastqc2-pbs SW2.fastqc2-pbs SW2.fastqc2-pbs SW2.fastqc2-pbs SW2.fastqc2-pbs

diginorm.c50: 614GMAAXX_1_1_pf.trim.fq 614GMAAXX_2_1_pf.trim.fq 614GMAAXX_3_1_pf.trim.fq 614GMAAXX_4_1_pf.trim.fq 614GMAAXX_6_1_pf.trim.fq 614GMAAXX_7_1_pf.trim.fq 614GMAAXX_8_1_pf.trim.fq AB057JABXX_s_7_pe.trim.fq AI.trim.fq AK.trim.fq BB0582ABXX_s_3_pe.trim.fq BB0582ABXX_s_4_pe.trim.fq BB0582ABXX_s_5_pe.trim.fq BB0582ABXX_s_6_pe.trim.fq BB0582ABXX_s_7_pe.trim.fq BB0582ABXX_s_8_pe.trim.fq brain_paired.trim.fq FW1 FW1 FW1 G81 G81 G82 G82 G83 G83 Gill.trim.fq JV1_liver.trim.fq JVa_liver.trim.fq JVb_liver.trim.fq JV_intest.trim.fq Kidney.trim.fq L81 L81 L82 L82 L83 L83 L83 L83 LI.trim.fq Lips.trim.fq LK.trim.fq M0_liver.trim.fq M2_liver.trim.fq M5_liver.trim.fq Monocytes.trim.fq mssm_single_s1.trim.fq mssm_single_s2.trim.fq mssm_single_s3.trim.fq mssm_single_s5.trim.fq mssm_single_s6.trim.fq mssm_single_s7.trim.fq Muscle.trim.fq Neutrophils.trim.fq OFHS.trim.fq POF POF PSMGill.trim.fq Rope.trim.fq SMGill.trim.fq SMHS.trim.fq SN.trim.fq SPDI.trim.fq spinal_cord_paired.trim.fq SPK.trim.fq SPPI.trim.fq Stg Stg Stg Stg Stg Stg Stg Stg Stg Stg Stg Stg Stg Stg SW2 SW2 SW2 SW2
	$(TIME) -o diginorm.c50.time $(KHMER)/scripts/normalize-by-median.py -x 10000000000.0 -N 4 -k 20 -C 50 -s diginorm.c50.ht -R diginorm.c50.report 614GMAAXX_1_1_pf.trim.fq 614GMAAXX_2_1_pf.trim.fq 614GMAAXX_3_1_pf.trim.fq 614GMAAXX_4_1_pf.trim.fq 614GMAAXX_6_1_pf.trim.fq 614GMAAXX_7_1_pf.trim.fq 614GMAAXX_8_1_pf.trim.fq AB057JABXX_s_7_pe.trim.fq AI.trim.fq AK.trim.fq BB0582ABXX_s_3_pe.trim.fq BB0582ABXX_s_4_pe.trim.fq BB0582ABXX_s_5_pe.trim.fq BB0582ABXX_s_6_pe.trim.fq BB0582ABXX_s_7_pe.trim.fq BB0582ABXX_s_8_pe.trim.fq brain_paired.trim.fq FW1 FW1 FW1 G81 G81 G82 G82 G83 G83 Gill.trim.fq JV1_liver.trim.fq JVa_liver.trim.fq JVb_liver.trim.fq JV_intest.trim.fq Kidney.trim.fq L81 L81 L82 L82 L83 L83 L83 L83 LI.trim.fq Lips.trim.fq LK.trim.fq M0_liver.trim.fq M2_liver.trim.fq M5_liver.trim.fq Monocytes.trim.fq mssm_single_s1.trim.fq mssm_single_s2.trim.fq mssm_single_s3.trim.fq mssm_single_s5.trim.fq mssm_single_s6.trim.fq mssm_single_s7.trim.fq Muscle.trim.fq Neutrophils.trim.fq OFHS.trim.fq POF POF PSMGill.trim.fq Rope.trim.fq SMGill.trim.fq SMHS.trim.fq SN.trim.fq SPDI.trim.fq spinal_cord_paired.trim.fq SPK.trim.fq SPPI.trim.fq Stg Stg Stg Stg Stg Stg Stg Stg Stg Stg Stg Stg Stg Stg SW2 SW2 SW2 SW2

diginorm.c50-pbs:
	qsub -l $(DGRES) diginorm.c50.sub


all.c20.fa: 614GMAAXX_1_1_pf.trim.fq.keep 614GMAAXX_2_1_pf.trim.fq.keep 614GMAAXX_3_1_pf.trim.fq.keep 614GMAAXX_4_1_pf.trim.fq.keep 614GMAAXX_6_1_pf.trim.fq.keep 614GMAAXX_7_1_pf.trim.fq.keep 614GMAAXX_8_1_pf.trim.fq.keep AB057JABXX_s_7_pe.trim.fq.keep AI.trim.fq.keep AK.trim.fq.keep BB0582ABXX_s_3_pe.trim.fq.keep BB0582ABXX_s_4_pe.trim.fq.keep BB0582ABXX_s_5_pe.trim.fq.keep BB0582ABXX_s_6_pe.trim.fq.keep BB0582ABXX_s_7_pe.trim.fq.keep BB0582ABXX_s_8_pe.trim.fq.keep brain_paired.trim.fq.keep FW1.keep FW1.keep FW1.keep G81.keep G81.keep G82.keep G82.keep G83.keep G83.keep Gill.trim.fq.keep JV1_liver.trim.fq.keep JVa_liver.trim.fq.keep JVb_liver.trim.fq.keep JV_intest.trim.fq.keep Kidney.trim.fq.keep L81.keep L81.keep L82.keep L82.keep L83.keep L83.keep L83.keep L83.keep LI.trim.fq.keep Lips.trim.fq.keep LK.trim.fq.keep M0_liver.trim.fq.keep M2_liver.trim.fq.keep M5_liver.trim.fq.keep Monocytes.trim.fq.keep mssm_single_s1.trim.fq.keep mssm_single_s2.trim.fq.keep mssm_single_s3.trim.fq.keep mssm_single_s5.trim.fq.keep mssm_single_s6.trim.fq.keep mssm_single_s7.trim.fq.keep Muscle.trim.fq.keep Neutrophils.trim.fq.keep OFHS.trim.fq.keep POF.keep POF.keep PSMGill.trim.fq.keep Rope.trim.fq.keep SMGill.trim.fq.keep SMHS.trim.fq.keep SN.trim.fq.keep SPDI.trim.fq.keep spinal_cord_paired.trim.fq.keep SPK.trim.fq.keep SPPI.trim.fq.keep Stg.keep Stg.keep Stg.keep Stg.keep Stg.keep Stg.keep Stg.keep Stg.keep Stg.keep Stg.keep Stg.keep Stg.keep Stg.keep Stg.keep SW2.keep SW2.keep SW2.keep SW2.keep
	$(TIME) -o diginorm.c20.time $(KHMER)/scripts/normalize-by-median.py -x 10000000000.0 -N 4 -k 20 -C 20 -s diginorm.c20.ht -R diginorm.c20.report 614GMAAXX_1_1_pf.trim.fq.keep 614GMAAXX_2_1_pf.trim.fq.keep 614GMAAXX_3_1_pf.trim.fq.keep 614GMAAXX_4_1_pf.trim.fq.keep 614GMAAXX_6_1_pf.trim.fq.keep 614GMAAXX_7_1_pf.trim.fq.keep 614GMAAXX_8_1_pf.trim.fq.keep AB057JABXX_s_7_pe.trim.fq.keep AI.trim.fq.keep AK.trim.fq.keep BB0582ABXX_s_3_pe.trim.fq.keep BB0582ABXX_s_4_pe.trim.fq.keep BB0582ABXX_s_5_pe.trim.fq.keep BB0582ABXX_s_6_pe.trim.fq.keep BB0582ABXX_s_7_pe.trim.fq.keep BB0582ABXX_s_8_pe.trim.fq.keep brain_paired.trim.fq.keep FW1.keep FW1.keep FW1.keep G81.keep G81.keep G82.keep G82.keep G83.keep G83.keep Gill.trim.fq.keep JV1_liver.trim.fq.keep JVa_liver.trim.fq.keep JVb_liver.trim.fq.keep JV_intest.trim.fq.keep Kidney.trim.fq.keep L81.keep L81.keep L82.keep L82.keep L83.keep L83.keep L83.keep L83.keep LI.trim.fq.keep Lips.trim.fq.keep LK.trim.fq.keep M0_liver.trim.fq.keep M2_liver.trim.fq.keep M5_liver.trim.fq.keep Monocytes.trim.fq.keep mssm_single_s1.trim.fq.keep mssm_single_s2.trim.fq.keep mssm_single_s3.trim.fq.keep mssm_single_s5.trim.fq.keep mssm_single_s6.trim.fq.keep mssm_single_s7.trim.fq.keep Muscle.trim.fq.keep Neutrophils.trim.fq.keep OFHS.trim.fq.keep POF.keep POF.keep PSMGill.trim.fq.keep Rope.trim.fq.keep SMGill.trim.fq.keep SMHS.trim.fq.keep SN.trim.fq.keep SPDI.trim.fq.keep spinal_cord_paired.trim.fq.keep SPK.trim.fq.keep SPPI.trim.fq.keep Stg.keep Stg.keep Stg.keep Stg.keep Stg.keep Stg.keep Stg.keep Stg.keep Stg.keep Stg.keep Stg.keep Stg.keep Stg.keep Stg.keep SW2.keep SW2.keep SW2.keep SW2.keep
	cat *.keep.keep > all.c20.fq
	python $(KHMER)/sandbox/fastq-to-fasta.py all.c20.fq

diginorm.c20-pbs:
	qsub -l $(DGRES) diginorm.c20.sub

all.c20.fa.clean: all.c20.fa
	$(TIME) -o all.c20.clean.time $(SEQCLEAN) all.c20.fa -c 8 -v ../db/univec_core.fa -s adapters.fa -l 31
all.c20.fa.clean-pbs: all.c20.fa
	qsub -l $(CLEANRES) seqclean.sub
clean:
	rm -rf *.time *.clean *.keep *.trim.fq *.report *.log *.trimlog *.ht cleaning_*/ *.cidx *.cln all.c20.fa all.c20.fq *.sort
