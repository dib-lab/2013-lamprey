
KHMER=/w/khmer
SEQCLEAN=seqclean
TRIM=~/bin/Trimmomatic-0.25/trimmomatic-0.25.jar
TIME=/usr/bin/time -f "\n%E elapsed\n%U user cpu seconds\n%K average total memory\n%W times swapped out of main mem\n%w total waits\n%I filesystem inputs\n%O filesystem outputs\n"
KMERGENIE=~/bin/kmergenie-1.5397/kmergenie

KMERRES=nodes=1:ppn=4,walltime=72:00:00,mem=32gb
TRIMRES=nodes=1:ppn=4,walltime=18:00:00,mem=4gb
QCRES=nodes=1:ppn=1,walltime=2:00:00,mem=4gb
DGRES=nodes=1:ppn=1,walltime=72:00:00,mem=42gb
CLEANRES=nodes=1:ppn=8,walltime=12:00:00,mem=8gb

614GMAAXX_1_1_pf.fq_fastqc.zip: 614GMAAXX_1_1_pf.qc.fq.gz
	fastqc --noextract 614GMAAXX_1_1_pf.qc.fq.gz

614GMAAXX_1_1_pf.fq_fastqc.zip-pbs:
	echo make 614GMAAXX_1_1_pf.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-614GMAAXX_1_1_pf

614GMAAXX_2_1_pf.fq_fastqc.zip: 614GMAAXX_2_1_pf.qc.fq.gz
	fastqc --noextract 614GMAAXX_2_1_pf.qc.fq.gz

614GMAAXX_2_1_pf.fq_fastqc.zip-pbs:
	echo make 614GMAAXX_2_1_pf.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-614GMAAXX_2_1_pf

614GMAAXX_3_1_pf.fq_fastqc.zip: 614GMAAXX_3_1_pf.qc.fq.gz
	fastqc --noextract 614GMAAXX_3_1_pf.qc.fq.gz

614GMAAXX_3_1_pf.fq_fastqc.zip-pbs:
	echo make 614GMAAXX_3_1_pf.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-614GMAAXX_3_1_pf

614GMAAXX_4_1_pf.fq_fastqc.zip: 614GMAAXX_4_1_pf.qc.fq.gz
	fastqc --noextract 614GMAAXX_4_1_pf.qc.fq.gz

614GMAAXX_4_1_pf.fq_fastqc.zip-pbs:
	echo make 614GMAAXX_4_1_pf.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-614GMAAXX_4_1_pf

614GMAAXX_6_1_pf.fq_fastqc.zip: 614GMAAXX_6_1_pf.qc.fq.gz
	fastqc --noextract 614GMAAXX_6_1_pf.qc.fq.gz

614GMAAXX_6_1_pf.fq_fastqc.zip-pbs:
	echo make 614GMAAXX_6_1_pf.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-614GMAAXX_6_1_pf

614GMAAXX_7_1_pf.fq_fastqc.zip: 614GMAAXX_7_1_pf.qc.fq.gz
	fastqc --noextract 614GMAAXX_7_1_pf.qc.fq.gz

614GMAAXX_7_1_pf.fq_fastqc.zip-pbs:
	echo make 614GMAAXX_7_1_pf.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-614GMAAXX_7_1_pf

614GMAAXX_8_1_pf.fq_fastqc.zip: 614GMAAXX_8_1_pf.qc.fq.gz
	fastqc --noextract 614GMAAXX_8_1_pf.qc.fq.gz

614GMAAXX_8_1_pf.fq_fastqc.zip-pbs:
	echo make 614GMAAXX_8_1_pf.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-614GMAAXX_8_1_pf

AB057JABXX_s_7_pe.fq_fastqc.zip: AB057JABXX_s_7_pe.trim.fq.gz
	fastqc --noextract AB057JABXX_s_7_pe.trim.fq.gz

AB057JABXX_s_7_pe.fq_fastqc.zip-pbs:
	echo make AB057JABXX_s_7_pe.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-AB057JABXX_s_7_pe

AI.fq_fastqc.zip: AI.fq.gz
	fastqc --noextract AI.fq.gz

AI.fq_fastqc.zip-pbs:
	echo make AI.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-AI

AK.fq_fastqc.zip: AK.fq.gz
	fastqc --noextract AK.fq.gz

AK.fq_fastqc.zip-pbs:
	echo make AK.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-AK

BB0582ABXX_s_3_pe.fq_fastqc.zip: BB0582ABXX_s_3_pe.trim.fq.gz
	fastqc --noextract BB0582ABXX_s_3_pe.trim.fq.gz

BB0582ABXX_s_3_pe.fq_fastqc.zip-pbs:
	echo make BB0582ABXX_s_3_pe.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-BB0582ABXX_s_3_pe

BB0582ABXX_s_4_pe.fq_fastqc.zip: BB0582ABXX_s_4_pe.trim.fq.gz
	fastqc --noextract BB0582ABXX_s_4_pe.trim.fq.gz

BB0582ABXX_s_4_pe.fq_fastqc.zip-pbs:
	echo make BB0582ABXX_s_4_pe.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-BB0582ABXX_s_4_pe

BB0582ABXX_s_5_pe.fq_fastqc.zip: BB0582ABXX_s_5_pe.trim.fq.gz
	fastqc --noextract BB0582ABXX_s_5_pe.trim.fq.gz

BB0582ABXX_s_5_pe.fq_fastqc.zip-pbs:
	echo make BB0582ABXX_s_5_pe.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-BB0582ABXX_s_5_pe

BB0582ABXX_s_6_pe.fq_fastqc.zip: BB0582ABXX_s_6_pe.trim.fq.gz
	fastqc --noextract BB0582ABXX_s_6_pe.trim.fq.gz

BB0582ABXX_s_6_pe.fq_fastqc.zip-pbs:
	echo make BB0582ABXX_s_6_pe.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-BB0582ABXX_s_6_pe

BB0582ABXX_s_7_pe.fq_fastqc.zip: BB0582ABXX_s_7_pe.trim.fq.gz
	fastqc --noextract BB0582ABXX_s_7_pe.trim.fq.gz

BB0582ABXX_s_7_pe.fq_fastqc.zip-pbs:
	echo make BB0582ABXX_s_7_pe.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-BB0582ABXX_s_7_pe

BB0582ABXX_s_8_pe.fq_fastqc.zip: BB0582ABXX_s_8_pe.trim.fq.gz
	fastqc --noextract BB0582ABXX_s_8_pe.trim.fq.gz

BB0582ABXX_s_8_pe.fq_fastqc.zip-pbs:
	echo make BB0582ABXX_s_8_pe.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-BB0582ABXX_s_8_pe

brain_paired.fq_fastqc.zip: brain_paired.qc.fq.gz
	fastqc --noextract brain_paired.qc.fq.gz

brain_paired.fq_fastqc.zip-pbs:
	echo make brain_paired.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-brain_paired

FW1-Gill.fq_fastqc.zip: FW1-Gill.fq.gz
	fastqc --noextract FW1-Gill.fq.gz

FW1-Gill.fq_fastqc.zip-pbs:
	echo make FW1-Gill.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-FW1-Gill

FW1-Intestine.fq_fastqc.zip: FW1-Intestine.fq.gz
	fastqc --noextract FW1-Intestine.fq.gz

FW1-Intestine.fq_fastqc.zip-pbs:
	echo make FW1-Intestine.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-FW1-Intestine

FW1-Kidney.fq_fastqc.zip: FW1-Kidney.fq.gz
	fastqc --noextract FW1-Kidney.fq.gz

FW1-Kidney.fq_fastqc.zip-pbs:
	echo make FW1-Kidney.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-FW1-Kidney

G81-a.fq_fastqc.zip: G81-a.fq.gz
	fastqc --noextract G81-a.fq.gz

G81-a.fq_fastqc.zip-pbs:
	echo make G81-a.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-G81-a

G81-b.fq_fastqc.zip: G81-b.fq.gz
	fastqc --noextract G81-b.fq.gz

G81-b.fq_fastqc.zip-pbs:
	echo make G81-b.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-G81-b

G82-a.fq_fastqc.zip: G82-a.fq.gz
	fastqc --noextract G82-a.fq.gz

G82-a.fq_fastqc.zip-pbs:
	echo make G82-a.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-G82-a

G82-b.fq_fastqc.zip: G82-b.fq.gz
	fastqc --noextract G82-b.fq.gz

G82-b.fq_fastqc.zip-pbs:
	echo make G82-b.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-G82-b

G83-a.fq_fastqc.zip: G83-a.fq.gz
	fastqc --noextract G83-a.fq.gz

G83-a.fq_fastqc.zip-pbs:
	echo make G83-a.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-G83-a

G83-b.fq_fastqc.zip: G83-b.fq.gz
	fastqc --noextract G83-b.fq.gz

G83-b.fq_fastqc.zip-pbs:
	echo make G83-b.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-G83-b

Gill.fq_fastqc.zip: Gill.fq.gz
	fastqc --noextract Gill.fq.gz

Gill.fq_fastqc.zip-pbs:
	echo make Gill.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-Gill

JV1_liver.fq_fastqc.zip: JV1_liver.fq.gz
	fastqc --noextract JV1_liver.fq.gz

JV1_liver.fq_fastqc.zip-pbs:
	echo make JV1_liver.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-JV1_liver

JVa_liver.fq_fastqc.zip: JVa_liver.fq.gz
	fastqc --noextract JVa_liver.fq.gz

JVa_liver.fq_fastqc.zip-pbs:
	echo make JVa_liver.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-JVa_liver

JVb_liver.fq_fastqc.zip: JVb_liver.fq.gz
	fastqc --noextract JVb_liver.fq.gz

JVb_liver.fq_fastqc.zip-pbs:
	echo make JVb_liver.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-JVb_liver

JV_intest.fq_fastqc.zip: JV_intest.fq.gz
	fastqc --noextract JV_intest.fq.gz

JV_intest.fq_fastqc.zip-pbs:
	echo make JV_intest.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-JV_intest

Kidney.fq_fastqc.zip: Kidney.fq.gz
	fastqc --noextract Kidney.fq.gz

Kidney.fq_fastqc.zip-pbs:
	echo make Kidney.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-Kidney

L81-a.fq_fastqc.zip: L81-a.fq.gz
	fastqc --noextract L81-a.fq.gz

L81-a.fq_fastqc.zip-pbs:
	echo make L81-a.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-L81-a

L81-b.fq_fastqc.zip: L81-b.fq.gz
	fastqc --noextract L81-b.fq.gz

L81-b.fq_fastqc.zip-pbs:
	echo make L81-b.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-L81-b

L82-a.fq_fastqc.zip: L82-a.fq.gz
	fastqc --noextract L82-a.fq.gz

L82-a.fq_fastqc.zip-pbs:
	echo make L82-a.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-L82-a

L82-b.fq_fastqc.zip: L82-b.fq.gz
	fastqc --noextract L82-b.fq.gz

L82-b.fq_fastqc.zip-pbs:
	echo make L82-b.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-L82-b

L83-1-a.fq_fastqc.zip: L83-1-a.fq.gz
	fastqc --noextract L83-1-a.fq.gz

L83-1-a.fq_fastqc.zip-pbs:
	echo make L83-1-a.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-L83-1-a

L83-1-b.fq_fastqc.zip: L83-1-b.fq.gz
	fastqc --noextract L83-1-b.fq.gz

L83-1-b.fq_fastqc.zip-pbs:
	echo make L83-1-b.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-L83-1-b

L83-a.fq_fastqc.zip: L83-a.fq.gz
	fastqc --noextract L83-a.fq.gz

L83-a.fq_fastqc.zip-pbs:
	echo make L83-a.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-L83-a

L83-b.fq_fastqc.zip: L83-b.fq.gz
	fastqc --noextract L83-b.fq.gz

L83-b.fq_fastqc.zip-pbs:
	echo make L83-b.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-L83-b

LI.fq_fastqc.zip: LI.fq.gz
	fastqc --noextract LI.fq.gz

LI.fq_fastqc.zip-pbs:
	echo make LI.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-LI

Lips.fq_fastqc.zip: Lips.fq.gz
	fastqc --noextract Lips.fq.gz

Lips.fq_fastqc.zip-pbs:
	echo make Lips.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-Lips

LK.fq_fastqc.zip: LK.fq.gz
	fastqc --noextract LK.fq.gz

LK.fq_fastqc.zip-pbs:
	echo make LK.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-LK

M0_liver.fq_fastqc.zip: M0_liver.fq.gz
	fastqc --noextract M0_liver.fq.gz

M0_liver.fq_fastqc.zip-pbs:
	echo make M0_liver.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-M0_liver

M2_liver.fq_fastqc.zip: M2_liver.fq.gz
	fastqc --noextract M2_liver.fq.gz

M2_liver.fq_fastqc.zip-pbs:
	echo make M2_liver.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-M2_liver

M5_liver.fq_fastqc.zip: M5_liver.fq.gz
	fastqc --noextract M5_liver.fq.gz

M5_liver.fq_fastqc.zip-pbs:
	echo make M5_liver.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-M5_liver

Monocytes.fq_fastqc.zip: Monocytes.fq.gz
	fastqc --noextract Monocytes.fq.gz

Monocytes.fq_fastqc.zip-pbs:
	echo make Monocytes.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-Monocytes

mssm_single_s1.fq_fastqc.zip: mssm_single_s1.trim.fq.gz
	fastqc --noextract mssm_single_s1.trim.fq.gz

mssm_single_s1.fq_fastqc.zip-pbs:
	echo make mssm_single_s1.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-mssm_single_s1

mssm_single_s2.fq_fastqc.zip: mssm_single_s2.trim.fq.gz
	fastqc --noextract mssm_single_s2.trim.fq.gz

mssm_single_s2.fq_fastqc.zip-pbs:
	echo make mssm_single_s2.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-mssm_single_s2

mssm_single_s3.fq_fastqc.zip: mssm_single_s3.trim.fq.gz
	fastqc --noextract mssm_single_s3.trim.fq.gz

mssm_single_s3.fq_fastqc.zip-pbs:
	echo make mssm_single_s3.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-mssm_single_s3

mssm_single_s5.fq_fastqc.zip: mssm_single_s5.trim.fq.gz
	fastqc --noextract mssm_single_s5.trim.fq.gz

mssm_single_s5.fq_fastqc.zip-pbs:
	echo make mssm_single_s5.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-mssm_single_s5

mssm_single_s6.fq_fastqc.zip: mssm_single_s6.trim.fq.gz
	fastqc --noextract mssm_single_s6.trim.fq.gz

mssm_single_s6.fq_fastqc.zip-pbs:
	echo make mssm_single_s6.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-mssm_single_s6

mssm_single_s7.fq_fastqc.zip: mssm_single_s7.trim.fq.gz
	fastqc --noextract mssm_single_s7.trim.fq.gz

mssm_single_s7.fq_fastqc.zip-pbs:
	echo make mssm_single_s7.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-mssm_single_s7

Muscle.fq_fastqc.zip: Muscle.fq.gz
	fastqc --noextract Muscle.fq.gz

Muscle.fq_fastqc.zip-pbs:
	echo make Muscle.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-Muscle

Neutrophils.fq_fastqc.zip: Neutrophils.fq.gz
	fastqc --noextract Neutrophils.fq.gz

Neutrophils.fq_fastqc.zip-pbs:
	echo make Neutrophils.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-Neutrophils

OFHS.fq_fastqc.zip: OFHS.fq.gz
	fastqc --noextract OFHS.fq.gz

OFHS.fq_fastqc.zip-pbs:
	echo make OFHS.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-OFHS

POF-eye.fq_fastqc.zip: POF-eye.fq.gz
	fastqc --noextract POF-eye.fq.gz

POF-eye.fq_fastqc.zip-pbs:
	echo make POF-eye.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-POF-eye

POF-ts.fq_fastqc.zip: POF-ts.fq.gz
	fastqc --noextract POF-ts.fq.gz

POF-ts.fq_fastqc.zip-pbs:
	echo make POF-ts.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-POF-ts

PSMGill.fq_fastqc.zip: PSMGill.fq.gz
	fastqc --noextract PSMGill.fq.gz

PSMGill.fq_fastqc.zip-pbs:
	echo make PSMGill.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-PSMGill

Rope.fq_fastqc.zip: Rope.fq.gz
	fastqc --noextract Rope.fq.gz

Rope.fq_fastqc.zip-pbs:
	echo make Rope.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-Rope

SMGill.fq_fastqc.zip: SMGill.fq.gz
	fastqc --noextract SMGill.fq.gz

SMGill.fq_fastqc.zip-pbs:
	echo make SMGill.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-SMGill

SMHS.fq_fastqc.zip: SMHS.fq.gz
	fastqc --noextract SMHS.fq.gz

SMHS.fq_fastqc.zip-pbs:
	echo make SMHS.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-SMHS

SN.fq_fastqc.zip: SN.fq.gz
	fastqc --noextract SN.fq.gz

SN.fq_fastqc.zip-pbs:
	echo make SN.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-SN

SPDI.fq_fastqc.zip: SPDI.fq.gz
	fastqc --noextract SPDI.fq.gz

SPDI.fq_fastqc.zip-pbs:
	echo make SPDI.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-SPDI

spinal_cord_paired.fq_fastqc.zip: spinal_cord_paired.qc.fq.gz
	fastqc --noextract spinal_cord_paired.qc.fq.gz

spinal_cord_paired.fq_fastqc.zip-pbs:
	echo make spinal_cord_paired.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-spinal_cord_paired

SPK.fq_fastqc.zip: SPK.fq.gz
	fastqc --noextract SPK.fq.gz

SPK.fq_fastqc.zip-pbs:
	echo make SPK.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-SPK

SPPI.fq_fastqc.zip: SPPI.fq.gz
	fastqc --noextract SPPI.fq.gz

SPPI.fq_fastqc.zip-pbs:
	echo make SPPI.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-SPPI

Stg-1-Intestine.fq_fastqc.zip: Stg-1-Intestine.fq.gz
	fastqc --noextract Stg-1-Intestine.fq.gz

Stg-1-Intestine.fq_fastqc.zip-pbs:
	echo make Stg-1-Intestine.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-Stg-1-Intestine

Stg-1-Kidney.fq_fastqc.zip: Stg-1-Kidney.fq.gz
	fastqc --noextract Stg-1-Kidney.fq.gz

Stg-1-Kidney.fq_fastqc.zip-pbs:
	echo make Stg-1-Kidney.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-Stg-1-Kidney

Stg-2-Intestine.fq_fastqc.zip: Stg-2-Intestine.fq.gz
	fastqc --noextract Stg-2-Intestine.fq.gz

Stg-2-Intestine.fq_fastqc.zip-pbs:
	echo make Stg-2-Intestine.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-Stg-2-Intestine

Stg-2-Kidney.fq_fastqc.zip: Stg-2-Kidney.fq.gz
	fastqc --noextract Stg-2-Kidney.fq.gz

Stg-2-Kidney.fq_fastqc.zip-pbs:
	echo make Stg-2-Kidney.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-Stg-2-Kidney

Stg-3-Intestine.fq_fastqc.zip: Stg-3-Intestine.fq.gz
	fastqc --noextract Stg-3-Intestine.fq.gz

Stg-3-Intestine.fq_fastqc.zip-pbs:
	echo make Stg-3-Intestine.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-Stg-3-Intestine

Stg-3-Kidney.fq_fastqc.zip: Stg-3-Kidney.fq.gz
	fastqc --noextract Stg-3-Kidney.fq.gz

Stg-3-Kidney.fq_fastqc.zip-pbs:
	echo make Stg-3-Kidney.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-Stg-3-Kidney

Stg-4-Intestine.fq_fastqc.zip: Stg-4-Intestine.fq.gz
	fastqc --noextract Stg-4-Intestine.fq.gz

Stg-4-Intestine.fq_fastqc.zip-pbs:
	echo make Stg-4-Intestine.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-Stg-4-Intestine

Stg-4-Kidney.fq_fastqc.zip: Stg-4-Kidney.fq.gz
	fastqc --noextract Stg-4-Kidney.fq.gz

Stg-4-Kidney.fq_fastqc.zip-pbs:
	echo make Stg-4-Kidney.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-Stg-4-Kidney

Stg-5-Intestine.fq_fastqc.zip: Stg-5-Intestine.fq.gz
	fastqc --noextract Stg-5-Intestine.fq.gz

Stg-5-Intestine.fq_fastqc.zip-pbs:
	echo make Stg-5-Intestine.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-Stg-5-Intestine

Stg-5-Kidney.fq_fastqc.zip: Stg-5-Kidney.fq.gz
	fastqc --noextract Stg-5-Kidney.fq.gz

Stg-5-Kidney.fq_fastqc.zip-pbs:
	echo make Stg-5-Kidney.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-Stg-5-Kidney

Stg-6-Intestine.fq_fastqc.zip: Stg-6-Intestine.fq.gz
	fastqc --noextract Stg-6-Intestine.fq.gz

Stg-6-Intestine.fq_fastqc.zip-pbs:
	echo make Stg-6-Intestine.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-Stg-6-Intestine

Stg-6-Kidney.fq_fastqc.zip: Stg-6-Kidney.fq.gz
	fastqc --noextract Stg-6-Kidney.fq.gz

Stg-6-Kidney.fq_fastqc.zip-pbs:
	echo make Stg-6-Kidney.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-Stg-6-Kidney

Stg-7-Intestine.fq_fastqc.zip: Stg-7-Intestine.fq.gz
	fastqc --noextract Stg-7-Intestine.fq.gz

Stg-7-Intestine.fq_fastqc.zip-pbs:
	echo make Stg-7-Intestine.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-Stg-7-Intestine

Stg-7-Kidney.fq_fastqc.zip: Stg-7-Kidney.fq.gz
	fastqc --noextract Stg-7-Kidney.fq.gz

Stg-7-Kidney.fq_fastqc.zip-pbs:
	echo make Stg-7-Kidney.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-Stg-7-Kidney

SW2-Gill.fq_fastqc.zip: SW2-Gill.fq.gz
	fastqc --noextract SW2-Gill.fq.gz

SW2-Gill.fq_fastqc.zip-pbs:
	echo make SW2-Gill.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-SW2-Gill

SW2-Intestine-a.fq_fastqc.zip: SW2-Intestine-a.fq.gz
	fastqc --noextract SW2-Intestine-a.fq.gz

SW2-Intestine-a.fq_fastqc.zip-pbs:
	echo make SW2-Intestine-a.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-SW2-Intestine-a

SW2-Intestine-b.fq_fastqc.zip: SW2-Intestine-b.fq.gz
	fastqc --noextract SW2-Intestine-b.fq.gz

SW2-Intestine-b.fq_fastqc.zip-pbs:
	echo make SW2-Intestine-b.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-SW2-Intestine-b

SW2-Kidney.fq_fastqc.zip: SW2-Kidney.fq.gz
	fastqc --noextract SW2-Kidney.fq.gz

SW2-Kidney.fq_fastqc.zip-pbs:
	echo make SW2-Kidney.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc-SW2-Kidney

614GMAAXX_1_1_pf.trim.fq: 614GMAAXX_1_1_pf.qc.fq.gz
		$(TIME) -o 614GMAAXX_1_1_pf.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog 614GMAAXX_1_1_pf.trimlog -threads 4 614GMAAXX_1_1_pf.qc.fq.gz 614GMAAXX_1_1_pf.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

614GMAAXX_1_1_pf.trim.fq-pbs:
	echo make 614GMAAXX_1_1_pf.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-614GMAAXX_1_1_pf

614GMAAXX_2_1_pf.trim.fq: 614GMAAXX_2_1_pf.qc.fq.gz
		$(TIME) -o 614GMAAXX_2_1_pf.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog 614GMAAXX_2_1_pf.trimlog -threads 4 614GMAAXX_2_1_pf.qc.fq.gz 614GMAAXX_2_1_pf.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

614GMAAXX_2_1_pf.trim.fq-pbs:
	echo make 614GMAAXX_2_1_pf.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-614GMAAXX_2_1_pf

614GMAAXX_3_1_pf.trim.fq: 614GMAAXX_3_1_pf.qc.fq.gz
		$(TIME) -o 614GMAAXX_3_1_pf.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog 614GMAAXX_3_1_pf.trimlog -threads 4 614GMAAXX_3_1_pf.qc.fq.gz 614GMAAXX_3_1_pf.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

614GMAAXX_3_1_pf.trim.fq-pbs:
	echo make 614GMAAXX_3_1_pf.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-614GMAAXX_3_1_pf

614GMAAXX_4_1_pf.trim.fq: 614GMAAXX_4_1_pf.qc.fq.gz
		$(TIME) -o 614GMAAXX_4_1_pf.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog 614GMAAXX_4_1_pf.trimlog -threads 4 614GMAAXX_4_1_pf.qc.fq.gz 614GMAAXX_4_1_pf.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

614GMAAXX_4_1_pf.trim.fq-pbs:
	echo make 614GMAAXX_4_1_pf.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-614GMAAXX_4_1_pf

614GMAAXX_6_1_pf.trim.fq: 614GMAAXX_6_1_pf.qc.fq.gz
		$(TIME) -o 614GMAAXX_6_1_pf.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog 614GMAAXX_6_1_pf.trimlog -threads 4 614GMAAXX_6_1_pf.qc.fq.gz 614GMAAXX_6_1_pf.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

614GMAAXX_6_1_pf.trim.fq-pbs:
	echo make 614GMAAXX_6_1_pf.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-614GMAAXX_6_1_pf

614GMAAXX_7_1_pf.trim.fq: 614GMAAXX_7_1_pf.qc.fq.gz
		$(TIME) -o 614GMAAXX_7_1_pf.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog 614GMAAXX_7_1_pf.trimlog -threads 4 614GMAAXX_7_1_pf.qc.fq.gz 614GMAAXX_7_1_pf.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

614GMAAXX_7_1_pf.trim.fq-pbs:
	echo make 614GMAAXX_7_1_pf.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-614GMAAXX_7_1_pf

614GMAAXX_8_1_pf.trim.fq: 614GMAAXX_8_1_pf.qc.fq.gz
		$(TIME) -o 614GMAAXX_8_1_pf.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog 614GMAAXX_8_1_pf.trimlog -threads 4 614GMAAXX_8_1_pf.qc.fq.gz 614GMAAXX_8_1_pf.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

614GMAAXX_8_1_pf.trim.fq-pbs:
	echo make 614GMAAXX_8_1_pf.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-614GMAAXX_8_1_pf

AB057JABXX_s_7_pe.trim.fq: AB057JABXX_s_7_pe.trim.fq.gz
		$(TIME) -o AB057JABXX_s_7_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog AB057JABXX_s_7_pe.trimlog -threads 4 AB057JABXX_s_7_pe.trim.fq.gz AB057JABXX_s_7_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

AB057JABXX_s_7_pe.trim.fq-pbs:
	echo make AB057JABXX_s_7_pe.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-AB057JABXX_s_7_pe

AI.trim.fq: AI.fq.gz
		$(TIME) -o AI.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog AI.trimlog -threads 4 AI.fq.gz AI.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

AI.trim.fq-pbs:
	echo make AI.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-AI

AK.trim.fq: AK.fq.gz
		$(TIME) -o AK.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog AK.trimlog -threads 4 AK.fq.gz AK.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

AK.trim.fq-pbs:
	echo make AK.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-AK

BB0582ABXX_s_3_pe.trim.fq: BB0582ABXX_s_3_pe.trim.fq.gz
		$(TIME) -o BB0582ABXX_s_3_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog BB0582ABXX_s_3_pe.trimlog -threads 4 BB0582ABXX_s_3_pe.trim.fq.gz BB0582ABXX_s_3_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

BB0582ABXX_s_3_pe.trim.fq-pbs:
	echo make BB0582ABXX_s_3_pe.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-BB0582ABXX_s_3_pe

BB0582ABXX_s_4_pe.trim.fq: BB0582ABXX_s_4_pe.trim.fq.gz
		$(TIME) -o BB0582ABXX_s_4_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog BB0582ABXX_s_4_pe.trimlog -threads 4 BB0582ABXX_s_4_pe.trim.fq.gz BB0582ABXX_s_4_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

BB0582ABXX_s_4_pe.trim.fq-pbs:
	echo make BB0582ABXX_s_4_pe.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-BB0582ABXX_s_4_pe

BB0582ABXX_s_5_pe.trim.fq: BB0582ABXX_s_5_pe.trim.fq.gz
		$(TIME) -o BB0582ABXX_s_5_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog BB0582ABXX_s_5_pe.trimlog -threads 4 BB0582ABXX_s_5_pe.trim.fq.gz BB0582ABXX_s_5_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

BB0582ABXX_s_5_pe.trim.fq-pbs:
	echo make BB0582ABXX_s_5_pe.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-BB0582ABXX_s_5_pe

BB0582ABXX_s_6_pe.trim.fq: BB0582ABXX_s_6_pe.trim.fq.gz
		$(TIME) -o BB0582ABXX_s_6_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog BB0582ABXX_s_6_pe.trimlog -threads 4 BB0582ABXX_s_6_pe.trim.fq.gz BB0582ABXX_s_6_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

BB0582ABXX_s_6_pe.trim.fq-pbs:
	echo make BB0582ABXX_s_6_pe.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-BB0582ABXX_s_6_pe

BB0582ABXX_s_7_pe.trim.fq: BB0582ABXX_s_7_pe.trim.fq.gz
		$(TIME) -o BB0582ABXX_s_7_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog BB0582ABXX_s_7_pe.trimlog -threads 4 BB0582ABXX_s_7_pe.trim.fq.gz BB0582ABXX_s_7_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

BB0582ABXX_s_7_pe.trim.fq-pbs:
	echo make BB0582ABXX_s_7_pe.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-BB0582ABXX_s_7_pe

BB0582ABXX_s_8_pe.trim.fq: BB0582ABXX_s_8_pe.trim.fq.gz
		$(TIME) -o BB0582ABXX_s_8_pe.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog BB0582ABXX_s_8_pe.trimlog -threads 4 BB0582ABXX_s_8_pe.trim.fq.gz BB0582ABXX_s_8_pe.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

BB0582ABXX_s_8_pe.trim.fq-pbs:
	echo make BB0582ABXX_s_8_pe.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-BB0582ABXX_s_8_pe

brain_paired.trim.fq: brain_paired.qc.fq.gz
		$(TIME) -o brain_paired.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog brain_paired.trimlog -threads 4 brain_paired.qc.fq.gz brain_paired.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

brain_paired.trim.fq-pbs:
	echo make brain_paired.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-brain_paired

FW1-Gill.trim.fq: FW1-Gill.fq.gz
		$(TIME) -o FW1-Gill.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog FW1-Gill.trimlog -threads 4 FW1-Gill.fq.gz FW1-Gill.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

FW1-Gill.trim.fq-pbs:
	echo make FW1-Gill.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-FW1-Gill

FW1-Intestine.trim.fq: FW1-Intestine.fq.gz
		$(TIME) -o FW1-Intestine.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog FW1-Intestine.trimlog -threads 4 FW1-Intestine.fq.gz FW1-Intestine.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

FW1-Intestine.trim.fq-pbs:
	echo make FW1-Intestine.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-FW1-Intestine

FW1-Kidney.trim.fq: FW1-Kidney.fq.gz
		$(TIME) -o FW1-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog FW1-Kidney.trimlog -threads 4 FW1-Kidney.fq.gz FW1-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

FW1-Kidney.trim.fq-pbs:
	echo make FW1-Kidney.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-FW1-Kidney

G81-a.trim.fq: G81-a.fq.gz
		$(TIME) -o G81-a.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog G81-a.trimlog -threads 4 G81-a.fq.gz G81-a.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

G81-a.trim.fq-pbs:
	echo make G81-a.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-G81-a

G81-b.trim.fq: G81-b.fq.gz
		$(TIME) -o G81-b.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog G81-b.trimlog -threads 4 G81-b.fq.gz G81-b.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

G81-b.trim.fq-pbs:
	echo make G81-b.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-G81-b

G82-a.trim.fq: G82-a.fq.gz
		$(TIME) -o G82-a.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog G82-a.trimlog -threads 4 G82-a.fq.gz G82-a.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

G82-a.trim.fq-pbs:
	echo make G82-a.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-G82-a

G82-b.trim.fq: G82-b.fq.gz
		$(TIME) -o G82-b.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog G82-b.trimlog -threads 4 G82-b.fq.gz G82-b.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

G82-b.trim.fq-pbs:
	echo make G82-b.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-G82-b

G83-a.trim.fq: G83-a.fq.gz
		$(TIME) -o G83-a.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog G83-a.trimlog -threads 4 G83-a.fq.gz G83-a.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

G83-a.trim.fq-pbs:
	echo make G83-a.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-G83-a

G83-b.trim.fq: G83-b.fq.gz
		$(TIME) -o G83-b.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog G83-b.trimlog -threads 4 G83-b.fq.gz G83-b.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

G83-b.trim.fq-pbs:
	echo make G83-b.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-G83-b

Gill.trim.fq: Gill.fq.gz
		$(TIME) -o Gill.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog Gill.trimlog -threads 4 Gill.fq.gz Gill.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Gill.trim.fq-pbs:
	echo make Gill.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-Gill

JV1_liver.trim.fq: JV1_liver.fq.gz
		$(TIME) -o JV1_liver.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog JV1_liver.trimlog -threads 4 JV1_liver.fq.gz JV1_liver.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

JV1_liver.trim.fq-pbs:
	echo make JV1_liver.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-JV1_liver

JVa_liver.trim.fq: JVa_liver.fq.gz
		$(TIME) -o JVa_liver.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog JVa_liver.trimlog -threads 4 JVa_liver.fq.gz JVa_liver.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

JVa_liver.trim.fq-pbs:
	echo make JVa_liver.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-JVa_liver

JVb_liver.trim.fq: JVb_liver.fq.gz
		$(TIME) -o JVb_liver.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog JVb_liver.trimlog -threads 4 JVb_liver.fq.gz JVb_liver.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

JVb_liver.trim.fq-pbs:
	echo make JVb_liver.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-JVb_liver

JV_intest.trim.fq: JV_intest.fq.gz
		$(TIME) -o JV_intest.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog JV_intest.trimlog -threads 4 JV_intest.fq.gz JV_intest.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

JV_intest.trim.fq-pbs:
	echo make JV_intest.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-JV_intest

Kidney.trim.fq: Kidney.fq.gz
		$(TIME) -o Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog Kidney.trimlog -threads 4 Kidney.fq.gz Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Kidney.trim.fq-pbs:
	echo make Kidney.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-Kidney

L81-a.trim.fq: L81-a.fq.gz
		$(TIME) -o L81-a.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog L81-a.trimlog -threads 4 L81-a.fq.gz L81-a.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

L81-a.trim.fq-pbs:
	echo make L81-a.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-L81-a

L81-b.trim.fq: L81-b.fq.gz
		$(TIME) -o L81-b.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog L81-b.trimlog -threads 4 L81-b.fq.gz L81-b.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

L81-b.trim.fq-pbs:
	echo make L81-b.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-L81-b

L82-a.trim.fq: L82-a.fq.gz
		$(TIME) -o L82-a.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog L82-a.trimlog -threads 4 L82-a.fq.gz L82-a.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

L82-a.trim.fq-pbs:
	echo make L82-a.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-L82-a

L82-b.trim.fq: L82-b.fq.gz
		$(TIME) -o L82-b.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog L82-b.trimlog -threads 4 L82-b.fq.gz L82-b.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

L82-b.trim.fq-pbs:
	echo make L82-b.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-L82-b

L83-1-a.trim.fq: L83-1-a.fq.gz
		$(TIME) -o L83-1-a.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog L83-1-a.trimlog -threads 4 L83-1-a.fq.gz L83-1-a.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

L83-1-a.trim.fq-pbs:
	echo make L83-1-a.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-L83-1-a

L83-1-b.trim.fq: L83-1-b.fq.gz
		$(TIME) -o L83-1-b.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog L83-1-b.trimlog -threads 4 L83-1-b.fq.gz L83-1-b.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

L83-1-b.trim.fq-pbs:
	echo make L83-1-b.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-L83-1-b

L83-a.trim.fq: L83-a.fq.gz
		$(TIME) -o L83-a.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog L83-a.trimlog -threads 4 L83-a.fq.gz L83-a.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

L83-a.trim.fq-pbs:
	echo make L83-a.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-L83-a

L83-b.trim.fq: L83-b.fq.gz
		$(TIME) -o L83-b.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog L83-b.trimlog -threads 4 L83-b.fq.gz L83-b.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

L83-b.trim.fq-pbs:
	echo make L83-b.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-L83-b

LI.trim.fq: LI.fq.gz
		$(TIME) -o LI.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog LI.trimlog -threads 4 LI.fq.gz LI.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

LI.trim.fq-pbs:
	echo make LI.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-LI

Lips.trim.fq: Lips.fq.gz
		$(TIME) -o Lips.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog Lips.trimlog -threads 4 Lips.fq.gz Lips.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Lips.trim.fq-pbs:
	echo make Lips.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-Lips

LK.trim.fq: LK.fq.gz
		$(TIME) -o LK.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog LK.trimlog -threads 4 LK.fq.gz LK.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

LK.trim.fq-pbs:
	echo make LK.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-LK

M0_liver.trim.fq: M0_liver.fq.gz
		$(TIME) -o M0_liver.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog M0_liver.trimlog -threads 4 M0_liver.fq.gz M0_liver.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

M0_liver.trim.fq-pbs:
	echo make M0_liver.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-M0_liver

M2_liver.trim.fq: M2_liver.fq.gz
		$(TIME) -o M2_liver.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog M2_liver.trimlog -threads 4 M2_liver.fq.gz M2_liver.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

M2_liver.trim.fq-pbs:
	echo make M2_liver.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-M2_liver

M5_liver.trim.fq: M5_liver.fq.gz
		$(TIME) -o M5_liver.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog M5_liver.trimlog -threads 4 M5_liver.fq.gz M5_liver.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

M5_liver.trim.fq-pbs:
	echo make M5_liver.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-M5_liver

Monocytes.trim.fq: Monocytes.fq.gz
		$(TIME) -o Monocytes.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog Monocytes.trimlog -threads 4 Monocytes.fq.gz Monocytes.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Monocytes.trim.fq-pbs:
	echo make Monocytes.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-Monocytes

mssm_single_s1.trim.fq: mssm_single_s1.trim.fq.gz
		$(TIME) -o mssm_single_s1.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog mssm_single_s1.trimlog -threads 4 mssm_single_s1.trim.fq.gz mssm_single_s1.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

mssm_single_s1.trim.fq-pbs:
	echo make mssm_single_s1.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-mssm_single_s1

mssm_single_s2.trim.fq: mssm_single_s2.trim.fq.gz
		$(TIME) -o mssm_single_s2.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog mssm_single_s2.trimlog -threads 4 mssm_single_s2.trim.fq.gz mssm_single_s2.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

mssm_single_s2.trim.fq-pbs:
	echo make mssm_single_s2.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-mssm_single_s2

mssm_single_s3.trim.fq: mssm_single_s3.trim.fq.gz
		$(TIME) -o mssm_single_s3.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog mssm_single_s3.trimlog -threads 4 mssm_single_s3.trim.fq.gz mssm_single_s3.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

mssm_single_s3.trim.fq-pbs:
	echo make mssm_single_s3.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-mssm_single_s3

mssm_single_s5.trim.fq: mssm_single_s5.trim.fq.gz
		$(TIME) -o mssm_single_s5.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog mssm_single_s5.trimlog -threads 4 mssm_single_s5.trim.fq.gz mssm_single_s5.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

mssm_single_s5.trim.fq-pbs:
	echo make mssm_single_s5.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-mssm_single_s5

mssm_single_s6.trim.fq: mssm_single_s6.trim.fq.gz
		$(TIME) -o mssm_single_s6.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog mssm_single_s6.trimlog -threads 4 mssm_single_s6.trim.fq.gz mssm_single_s6.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

mssm_single_s6.trim.fq-pbs:
	echo make mssm_single_s6.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-mssm_single_s6

mssm_single_s7.trim.fq: mssm_single_s7.trim.fq.gz
		$(TIME) -o mssm_single_s7.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog mssm_single_s7.trimlog -threads 4 mssm_single_s7.trim.fq.gz mssm_single_s7.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

mssm_single_s7.trim.fq-pbs:
	echo make mssm_single_s7.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-mssm_single_s7

Muscle.trim.fq: Muscle.fq.gz
		$(TIME) -o Muscle.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog Muscle.trimlog -threads 4 Muscle.fq.gz Muscle.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Muscle.trim.fq-pbs:
	echo make Muscle.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-Muscle

Neutrophils.trim.fq: Neutrophils.fq.gz
		$(TIME) -o Neutrophils.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog Neutrophils.trimlog -threads 4 Neutrophils.fq.gz Neutrophils.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Neutrophils.trim.fq-pbs:
	echo make Neutrophils.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-Neutrophils

OFHS.trim.fq: OFHS.fq.gz
		$(TIME) -o OFHS.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog OFHS.trimlog -threads 4 OFHS.fq.gz OFHS.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

OFHS.trim.fq-pbs:
	echo make OFHS.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-OFHS

POF-eye.trim.fq: POF-eye.fq.gz
		$(TIME) -o POF-eye.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog POF-eye.trimlog -threads 4 POF-eye.fq.gz POF-eye.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

POF-eye.trim.fq-pbs:
	echo make POF-eye.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-POF-eye

POF-ts.trim.fq: POF-ts.fq.gz
		$(TIME) -o POF-ts.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog POF-ts.trimlog -threads 4 POF-ts.fq.gz POF-ts.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

POF-ts.trim.fq-pbs:
	echo make POF-ts.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-POF-ts

PSMGill.trim.fq: PSMGill.fq.gz
		$(TIME) -o PSMGill.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog PSMGill.trimlog -threads 4 PSMGill.fq.gz PSMGill.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

PSMGill.trim.fq-pbs:
	echo make PSMGill.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-PSMGill

Rope.trim.fq: Rope.fq.gz
		$(TIME) -o Rope.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog Rope.trimlog -threads 4 Rope.fq.gz Rope.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Rope.trim.fq-pbs:
	echo make Rope.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-Rope

SMGill.trim.fq: SMGill.fq.gz
		$(TIME) -o SMGill.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog SMGill.trimlog -threads 4 SMGill.fq.gz SMGill.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

SMGill.trim.fq-pbs:
	echo make SMGill.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-SMGill

SMHS.trim.fq: SMHS.fq.gz
		$(TIME) -o SMHS.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog SMHS.trimlog -threads 4 SMHS.fq.gz SMHS.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

SMHS.trim.fq-pbs:
	echo make SMHS.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-SMHS

SN.trim.fq: SN.fq.gz
		$(TIME) -o SN.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog SN.trimlog -threads 4 SN.fq.gz SN.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

SN.trim.fq-pbs:
	echo make SN.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-SN

SPDI.trim.fq: SPDI.fq.gz
		$(TIME) -o SPDI.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog SPDI.trimlog -threads 4 SPDI.fq.gz SPDI.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

SPDI.trim.fq-pbs:
	echo make SPDI.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-SPDI

spinal_cord_paired.trim.fq: spinal_cord_paired.qc.fq.gz
		$(TIME) -o spinal_cord_paired.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog spinal_cord_paired.trimlog -threads 4 spinal_cord_paired.qc.fq.gz spinal_cord_paired.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

spinal_cord_paired.trim.fq-pbs:
	echo make spinal_cord_paired.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-spinal_cord_paired

SPK.trim.fq: SPK.fq.gz
		$(TIME) -o SPK.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog SPK.trimlog -threads 4 SPK.fq.gz SPK.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

SPK.trim.fq-pbs:
	echo make SPK.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-SPK

SPPI.trim.fq: SPPI.fq.gz
		$(TIME) -o SPPI.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog SPPI.trimlog -threads 4 SPPI.fq.gz SPPI.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

SPPI.trim.fq-pbs:
	echo make SPPI.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-SPPI

Stg-1-Intestine.trim.fq: Stg-1-Intestine.fq.gz
		$(TIME) -o Stg-1-Intestine.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog Stg-1-Intestine.trimlog -threads 4 Stg-1-Intestine.fq.gz Stg-1-Intestine.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-1-Intestine.trim.fq-pbs:
	echo make Stg-1-Intestine.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-Stg-1-Intestine

Stg-1-Kidney.trim.fq: Stg-1-Kidney.fq.gz
		$(TIME) -o Stg-1-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog Stg-1-Kidney.trimlog -threads 4 Stg-1-Kidney.fq.gz Stg-1-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-1-Kidney.trim.fq-pbs:
	echo make Stg-1-Kidney.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-Stg-1-Kidney

Stg-2-Intestine.trim.fq: Stg-2-Intestine.fq.gz
		$(TIME) -o Stg-2-Intestine.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog Stg-2-Intestine.trimlog -threads 4 Stg-2-Intestine.fq.gz Stg-2-Intestine.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-2-Intestine.trim.fq-pbs:
	echo make Stg-2-Intestine.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-Stg-2-Intestine

Stg-2-Kidney.trim.fq: Stg-2-Kidney.fq.gz
		$(TIME) -o Stg-2-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog Stg-2-Kidney.trimlog -threads 4 Stg-2-Kidney.fq.gz Stg-2-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-2-Kidney.trim.fq-pbs:
	echo make Stg-2-Kidney.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-Stg-2-Kidney

Stg-3-Intestine.trim.fq: Stg-3-Intestine.fq.gz
		$(TIME) -o Stg-3-Intestine.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog Stg-3-Intestine.trimlog -threads 4 Stg-3-Intestine.fq.gz Stg-3-Intestine.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-3-Intestine.trim.fq-pbs:
	echo make Stg-3-Intestine.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-Stg-3-Intestine

Stg-3-Kidney.trim.fq: Stg-3-Kidney.fq.gz
		$(TIME) -o Stg-3-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog Stg-3-Kidney.trimlog -threads 4 Stg-3-Kidney.fq.gz Stg-3-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-3-Kidney.trim.fq-pbs:
	echo make Stg-3-Kidney.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-Stg-3-Kidney

Stg-4-Intestine.trim.fq: Stg-4-Intestine.fq.gz
		$(TIME) -o Stg-4-Intestine.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog Stg-4-Intestine.trimlog -threads 4 Stg-4-Intestine.fq.gz Stg-4-Intestine.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-4-Intestine.trim.fq-pbs:
	echo make Stg-4-Intestine.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-Stg-4-Intestine

Stg-4-Kidney.trim.fq: Stg-4-Kidney.fq.gz
		$(TIME) -o Stg-4-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog Stg-4-Kidney.trimlog -threads 4 Stg-4-Kidney.fq.gz Stg-4-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-4-Kidney.trim.fq-pbs:
	echo make Stg-4-Kidney.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-Stg-4-Kidney

Stg-5-Intestine.trim.fq: Stg-5-Intestine.fq.gz
		$(TIME) -o Stg-5-Intestine.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog Stg-5-Intestine.trimlog -threads 4 Stg-5-Intestine.fq.gz Stg-5-Intestine.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-5-Intestine.trim.fq-pbs:
	echo make Stg-5-Intestine.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-Stg-5-Intestine

Stg-5-Kidney.trim.fq: Stg-5-Kidney.fq.gz
		$(TIME) -o Stg-5-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog Stg-5-Kidney.trimlog -threads 4 Stg-5-Kidney.fq.gz Stg-5-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-5-Kidney.trim.fq-pbs:
	echo make Stg-5-Kidney.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-Stg-5-Kidney

Stg-6-Intestine.trim.fq: Stg-6-Intestine.fq.gz
		$(TIME) -o Stg-6-Intestine.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog Stg-6-Intestine.trimlog -threads 4 Stg-6-Intestine.fq.gz Stg-6-Intestine.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-6-Intestine.trim.fq-pbs:
	echo make Stg-6-Intestine.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-Stg-6-Intestine

Stg-6-Kidney.trim.fq: Stg-6-Kidney.fq.gz
		$(TIME) -o Stg-6-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog Stg-6-Kidney.trimlog -threads 4 Stg-6-Kidney.fq.gz Stg-6-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-6-Kidney.trim.fq-pbs:
	echo make Stg-6-Kidney.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-Stg-6-Kidney

Stg-7-Intestine.trim.fq: Stg-7-Intestine.fq.gz
		$(TIME) -o Stg-7-Intestine.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog Stg-7-Intestine.trimlog -threads 4 Stg-7-Intestine.fq.gz Stg-7-Intestine.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-7-Intestine.trim.fq-pbs:
	echo make Stg-7-Intestine.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-Stg-7-Intestine

Stg-7-Kidney.trim.fq: Stg-7-Kidney.fq.gz
		$(TIME) -o Stg-7-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog Stg-7-Kidney.trimlog -threads 4 Stg-7-Kidney.fq.gz Stg-7-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

Stg-7-Kidney.trim.fq-pbs:
	echo make Stg-7-Kidney.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-Stg-7-Kidney

SW2-Gill.trim.fq: SW2-Gill.fq.gz
		$(TIME) -o SW2-Gill.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog SW2-Gill.trimlog -threads 4 SW2-Gill.fq.gz SW2-Gill.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

SW2-Gill.trim.fq-pbs:
	echo make SW2-Gill.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-SW2-Gill

SW2-Intestine-a.trim.fq: SW2-Intestine-a.fq.gz
		$(TIME) -o SW2-Intestine-a.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog SW2-Intestine-a.trimlog -threads 4 SW2-Intestine-a.fq.gz SW2-Intestine-a.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

SW2-Intestine-a.trim.fq-pbs:
	echo make SW2-Intestine-a.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-SW2-Intestine-a

SW2-Intestine-b.trim.fq: SW2-Intestine-b.fq.gz
		$(TIME) -o SW2-Intestine-b.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog SW2-Intestine-b.trimlog -threads 4 SW2-Intestine-b.fq.gz SW2-Intestine-b.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

SW2-Intestine-b.trim.fq-pbs:
	echo make SW2-Intestine-b.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-SW2-Intestine-b

SW2-Kidney.trim.fq: SW2-Kidney.fq.gz
		$(TIME) -o SW2-Kidney.trim.time java -classpath $(TRIM) org.usadellab.trimmomatic.TrimmomaticSE -phred33 -trimlog SW2-Kidney.trimlog -threads 4 SW2-Kidney.fq.gz SW2-Kidney.trim.fq ILLUMINACLIP:adapters.fa:2:40:15 SLIDINGWINDOW:4:15 MINLEN:31

SW2-Kidney.trim.fq-pbs:
	echo make SW2-Kidney.trim.fq | cat trim.inc.sub - | qsub -l $(TRIMRES) -N trim-SW2-Kidney

614GMAAXX_1_1_pf.trim.fq.fq_fastqc.zip: 614GMAAXX_1_1_pf.trim.fq
	fastqc --noextract 614GMAAXX_1_1_pf.trim.fq

614GMAAXX_1_1_pf.trim.fq.fq_fastqc-pbs:
	echo make 614GMAAXX_1_1_pf.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-614GMAAXX_1_1_pf

614GMAAXX_2_1_pf.trim.fq.fq_fastqc.zip: 614GMAAXX_2_1_pf.trim.fq
	fastqc --noextract 614GMAAXX_2_1_pf.trim.fq

614GMAAXX_2_1_pf.trim.fq.fq_fastqc-pbs:
	echo make 614GMAAXX_2_1_pf.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-614GMAAXX_2_1_pf

614GMAAXX_3_1_pf.trim.fq.fq_fastqc.zip: 614GMAAXX_3_1_pf.trim.fq
	fastqc --noextract 614GMAAXX_3_1_pf.trim.fq

614GMAAXX_3_1_pf.trim.fq.fq_fastqc-pbs:
	echo make 614GMAAXX_3_1_pf.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-614GMAAXX_3_1_pf

614GMAAXX_4_1_pf.trim.fq.fq_fastqc.zip: 614GMAAXX_4_1_pf.trim.fq
	fastqc --noextract 614GMAAXX_4_1_pf.trim.fq

614GMAAXX_4_1_pf.trim.fq.fq_fastqc-pbs:
	echo make 614GMAAXX_4_1_pf.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-614GMAAXX_4_1_pf

614GMAAXX_6_1_pf.trim.fq.fq_fastqc.zip: 614GMAAXX_6_1_pf.trim.fq
	fastqc --noextract 614GMAAXX_6_1_pf.trim.fq

614GMAAXX_6_1_pf.trim.fq.fq_fastqc-pbs:
	echo make 614GMAAXX_6_1_pf.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-614GMAAXX_6_1_pf

614GMAAXX_7_1_pf.trim.fq.fq_fastqc.zip: 614GMAAXX_7_1_pf.trim.fq
	fastqc --noextract 614GMAAXX_7_1_pf.trim.fq

614GMAAXX_7_1_pf.trim.fq.fq_fastqc-pbs:
	echo make 614GMAAXX_7_1_pf.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-614GMAAXX_7_1_pf

614GMAAXX_8_1_pf.trim.fq.fq_fastqc.zip: 614GMAAXX_8_1_pf.trim.fq
	fastqc --noextract 614GMAAXX_8_1_pf.trim.fq

614GMAAXX_8_1_pf.trim.fq.fq_fastqc-pbs:
	echo make 614GMAAXX_8_1_pf.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-614GMAAXX_8_1_pf

AB057JABXX_s_7_pe.trim.fq.fq_fastqc.zip: AB057JABXX_s_7_pe.trim.fq
	fastqc --noextract AB057JABXX_s_7_pe.trim.fq

AB057JABXX_s_7_pe.trim.fq.fq_fastqc-pbs:
	echo make AB057JABXX_s_7_pe.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-AB057JABXX_s_7_pe

AI.trim.fq.fq_fastqc.zip: AI.trim.fq
	fastqc --noextract AI.trim.fq

AI.trim.fq.fq_fastqc-pbs:
	echo make AI.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-AI

AK.trim.fq.fq_fastqc.zip: AK.trim.fq
	fastqc --noextract AK.trim.fq

AK.trim.fq.fq_fastqc-pbs:
	echo make AK.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-AK

BB0582ABXX_s_3_pe.trim.fq.fq_fastqc.zip: BB0582ABXX_s_3_pe.trim.fq
	fastqc --noextract BB0582ABXX_s_3_pe.trim.fq

BB0582ABXX_s_3_pe.trim.fq.fq_fastqc-pbs:
	echo make BB0582ABXX_s_3_pe.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-BB0582ABXX_s_3_pe

BB0582ABXX_s_4_pe.trim.fq.fq_fastqc.zip: BB0582ABXX_s_4_pe.trim.fq
	fastqc --noextract BB0582ABXX_s_4_pe.trim.fq

BB0582ABXX_s_4_pe.trim.fq.fq_fastqc-pbs:
	echo make BB0582ABXX_s_4_pe.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-BB0582ABXX_s_4_pe

BB0582ABXX_s_5_pe.trim.fq.fq_fastqc.zip: BB0582ABXX_s_5_pe.trim.fq
	fastqc --noextract BB0582ABXX_s_5_pe.trim.fq

BB0582ABXX_s_5_pe.trim.fq.fq_fastqc-pbs:
	echo make BB0582ABXX_s_5_pe.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-BB0582ABXX_s_5_pe

BB0582ABXX_s_6_pe.trim.fq.fq_fastqc.zip: BB0582ABXX_s_6_pe.trim.fq
	fastqc --noextract BB0582ABXX_s_6_pe.trim.fq

BB0582ABXX_s_6_pe.trim.fq.fq_fastqc-pbs:
	echo make BB0582ABXX_s_6_pe.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-BB0582ABXX_s_6_pe

BB0582ABXX_s_7_pe.trim.fq.fq_fastqc.zip: BB0582ABXX_s_7_pe.trim.fq
	fastqc --noextract BB0582ABXX_s_7_pe.trim.fq

BB0582ABXX_s_7_pe.trim.fq.fq_fastqc-pbs:
	echo make BB0582ABXX_s_7_pe.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-BB0582ABXX_s_7_pe

BB0582ABXX_s_8_pe.trim.fq.fq_fastqc.zip: BB0582ABXX_s_8_pe.trim.fq
	fastqc --noextract BB0582ABXX_s_8_pe.trim.fq

BB0582ABXX_s_8_pe.trim.fq.fq_fastqc-pbs:
	echo make BB0582ABXX_s_8_pe.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-BB0582ABXX_s_8_pe

brain_paired.trim.fq.fq_fastqc.zip: brain_paired.trim.fq
	fastqc --noextract brain_paired.trim.fq

brain_paired.trim.fq.fq_fastqc-pbs:
	echo make brain_paired.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-brain_paired

FW1-Gill.trim.fq.fq_fastqc.zip: FW1-Gill.trim.fq
	fastqc --noextract FW1-Gill.trim.fq

FW1-Gill.trim.fq.fq_fastqc-pbs:
	echo make FW1-Gill.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-FW1-Gill

FW1-Intestine.trim.fq.fq_fastqc.zip: FW1-Intestine.trim.fq
	fastqc --noextract FW1-Intestine.trim.fq

FW1-Intestine.trim.fq.fq_fastqc-pbs:
	echo make FW1-Intestine.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-FW1-Intestine

FW1-Kidney.trim.fq.fq_fastqc.zip: FW1-Kidney.trim.fq
	fastqc --noextract FW1-Kidney.trim.fq

FW1-Kidney.trim.fq.fq_fastqc-pbs:
	echo make FW1-Kidney.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-FW1-Kidney

G81-a.trim.fq.fq_fastqc.zip: G81-a.trim.fq
	fastqc --noextract G81-a.trim.fq

G81-a.trim.fq.fq_fastqc-pbs:
	echo make G81-a.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-G81-a

G81-b.trim.fq.fq_fastqc.zip: G81-b.trim.fq
	fastqc --noextract G81-b.trim.fq

G81-b.trim.fq.fq_fastqc-pbs:
	echo make G81-b.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-G81-b

G82-a.trim.fq.fq_fastqc.zip: G82-a.trim.fq
	fastqc --noextract G82-a.trim.fq

G82-a.trim.fq.fq_fastqc-pbs:
	echo make G82-a.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-G82-a

G82-b.trim.fq.fq_fastqc.zip: G82-b.trim.fq
	fastqc --noextract G82-b.trim.fq

G82-b.trim.fq.fq_fastqc-pbs:
	echo make G82-b.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-G82-b

G83-a.trim.fq.fq_fastqc.zip: G83-a.trim.fq
	fastqc --noextract G83-a.trim.fq

G83-a.trim.fq.fq_fastqc-pbs:
	echo make G83-a.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-G83-a

G83-b.trim.fq.fq_fastqc.zip: G83-b.trim.fq
	fastqc --noextract G83-b.trim.fq

G83-b.trim.fq.fq_fastqc-pbs:
	echo make G83-b.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-G83-b

Gill.trim.fq.fq_fastqc.zip: Gill.trim.fq
	fastqc --noextract Gill.trim.fq

Gill.trim.fq.fq_fastqc-pbs:
	echo make Gill.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-Gill

JV1_liver.trim.fq.fq_fastqc.zip: JV1_liver.trim.fq
	fastqc --noextract JV1_liver.trim.fq

JV1_liver.trim.fq.fq_fastqc-pbs:
	echo make JV1_liver.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-JV1_liver

JVa_liver.trim.fq.fq_fastqc.zip: JVa_liver.trim.fq
	fastqc --noextract JVa_liver.trim.fq

JVa_liver.trim.fq.fq_fastqc-pbs:
	echo make JVa_liver.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-JVa_liver

JVb_liver.trim.fq.fq_fastqc.zip: JVb_liver.trim.fq
	fastqc --noextract JVb_liver.trim.fq

JVb_liver.trim.fq.fq_fastqc-pbs:
	echo make JVb_liver.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-JVb_liver

JV_intest.trim.fq.fq_fastqc.zip: JV_intest.trim.fq
	fastqc --noextract JV_intest.trim.fq

JV_intest.trim.fq.fq_fastqc-pbs:
	echo make JV_intest.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-JV_intest

Kidney.trim.fq.fq_fastqc.zip: Kidney.trim.fq
	fastqc --noextract Kidney.trim.fq

Kidney.trim.fq.fq_fastqc-pbs:
	echo make Kidney.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-Kidney

L81-a.trim.fq.fq_fastqc.zip: L81-a.trim.fq
	fastqc --noextract L81-a.trim.fq

L81-a.trim.fq.fq_fastqc-pbs:
	echo make L81-a.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-L81-a

L81-b.trim.fq.fq_fastqc.zip: L81-b.trim.fq
	fastqc --noextract L81-b.trim.fq

L81-b.trim.fq.fq_fastqc-pbs:
	echo make L81-b.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-L81-b

L82-a.trim.fq.fq_fastqc.zip: L82-a.trim.fq
	fastqc --noextract L82-a.trim.fq

L82-a.trim.fq.fq_fastqc-pbs:
	echo make L82-a.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-L82-a

L82-b.trim.fq.fq_fastqc.zip: L82-b.trim.fq
	fastqc --noextract L82-b.trim.fq

L82-b.trim.fq.fq_fastqc-pbs:
	echo make L82-b.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-L82-b

L83-1-a.trim.fq.fq_fastqc.zip: L83-1-a.trim.fq
	fastqc --noextract L83-1-a.trim.fq

L83-1-a.trim.fq.fq_fastqc-pbs:
	echo make L83-1-a.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-L83-1-a

L83-1-b.trim.fq.fq_fastqc.zip: L83-1-b.trim.fq
	fastqc --noextract L83-1-b.trim.fq

L83-1-b.trim.fq.fq_fastqc-pbs:
	echo make L83-1-b.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-L83-1-b

L83-a.trim.fq.fq_fastqc.zip: L83-a.trim.fq
	fastqc --noextract L83-a.trim.fq

L83-a.trim.fq.fq_fastqc-pbs:
	echo make L83-a.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-L83-a

L83-b.trim.fq.fq_fastqc.zip: L83-b.trim.fq
	fastqc --noextract L83-b.trim.fq

L83-b.trim.fq.fq_fastqc-pbs:
	echo make L83-b.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-L83-b

LI.trim.fq.fq_fastqc.zip: LI.trim.fq
	fastqc --noextract LI.trim.fq

LI.trim.fq.fq_fastqc-pbs:
	echo make LI.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-LI

Lips.trim.fq.fq_fastqc.zip: Lips.trim.fq
	fastqc --noextract Lips.trim.fq

Lips.trim.fq.fq_fastqc-pbs:
	echo make Lips.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-Lips

LK.trim.fq.fq_fastqc.zip: LK.trim.fq
	fastqc --noextract LK.trim.fq

LK.trim.fq.fq_fastqc-pbs:
	echo make LK.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-LK

M0_liver.trim.fq.fq_fastqc.zip: M0_liver.trim.fq
	fastqc --noextract M0_liver.trim.fq

M0_liver.trim.fq.fq_fastqc-pbs:
	echo make M0_liver.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-M0_liver

M2_liver.trim.fq.fq_fastqc.zip: M2_liver.trim.fq
	fastqc --noextract M2_liver.trim.fq

M2_liver.trim.fq.fq_fastqc-pbs:
	echo make M2_liver.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-M2_liver

M5_liver.trim.fq.fq_fastqc.zip: M5_liver.trim.fq
	fastqc --noextract M5_liver.trim.fq

M5_liver.trim.fq.fq_fastqc-pbs:
	echo make M5_liver.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-M5_liver

Monocytes.trim.fq.fq_fastqc.zip: Monocytes.trim.fq
	fastqc --noextract Monocytes.trim.fq

Monocytes.trim.fq.fq_fastqc-pbs:
	echo make Monocytes.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-Monocytes

mssm_single_s1.trim.fq.fq_fastqc.zip: mssm_single_s1.trim.fq
	fastqc --noextract mssm_single_s1.trim.fq

mssm_single_s1.trim.fq.fq_fastqc-pbs:
	echo make mssm_single_s1.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-mssm_single_s1

mssm_single_s2.trim.fq.fq_fastqc.zip: mssm_single_s2.trim.fq
	fastqc --noextract mssm_single_s2.trim.fq

mssm_single_s2.trim.fq.fq_fastqc-pbs:
	echo make mssm_single_s2.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-mssm_single_s2

mssm_single_s3.trim.fq.fq_fastqc.zip: mssm_single_s3.trim.fq
	fastqc --noextract mssm_single_s3.trim.fq

mssm_single_s3.trim.fq.fq_fastqc-pbs:
	echo make mssm_single_s3.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-mssm_single_s3

mssm_single_s5.trim.fq.fq_fastqc.zip: mssm_single_s5.trim.fq
	fastqc --noextract mssm_single_s5.trim.fq

mssm_single_s5.trim.fq.fq_fastqc-pbs:
	echo make mssm_single_s5.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-mssm_single_s5

mssm_single_s6.trim.fq.fq_fastqc.zip: mssm_single_s6.trim.fq
	fastqc --noextract mssm_single_s6.trim.fq

mssm_single_s6.trim.fq.fq_fastqc-pbs:
	echo make mssm_single_s6.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-mssm_single_s6

mssm_single_s7.trim.fq.fq_fastqc.zip: mssm_single_s7.trim.fq
	fastqc --noextract mssm_single_s7.trim.fq

mssm_single_s7.trim.fq.fq_fastqc-pbs:
	echo make mssm_single_s7.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-mssm_single_s7

Muscle.trim.fq.fq_fastqc.zip: Muscle.trim.fq
	fastqc --noextract Muscle.trim.fq

Muscle.trim.fq.fq_fastqc-pbs:
	echo make Muscle.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-Muscle

Neutrophils.trim.fq.fq_fastqc.zip: Neutrophils.trim.fq
	fastqc --noextract Neutrophils.trim.fq

Neutrophils.trim.fq.fq_fastqc-pbs:
	echo make Neutrophils.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-Neutrophils

OFHS.trim.fq.fq_fastqc.zip: OFHS.trim.fq
	fastqc --noextract OFHS.trim.fq

OFHS.trim.fq.fq_fastqc-pbs:
	echo make OFHS.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-OFHS

POF-eye.trim.fq.fq_fastqc.zip: POF-eye.trim.fq
	fastqc --noextract POF-eye.trim.fq

POF-eye.trim.fq.fq_fastqc-pbs:
	echo make POF-eye.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-POF-eye

POF-ts.trim.fq.fq_fastqc.zip: POF-ts.trim.fq
	fastqc --noextract POF-ts.trim.fq

POF-ts.trim.fq.fq_fastqc-pbs:
	echo make POF-ts.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-POF-ts

PSMGill.trim.fq.fq_fastqc.zip: PSMGill.trim.fq
	fastqc --noextract PSMGill.trim.fq

PSMGill.trim.fq.fq_fastqc-pbs:
	echo make PSMGill.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-PSMGill

Rope.trim.fq.fq_fastqc.zip: Rope.trim.fq
	fastqc --noextract Rope.trim.fq

Rope.trim.fq.fq_fastqc-pbs:
	echo make Rope.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-Rope

SMGill.trim.fq.fq_fastqc.zip: SMGill.trim.fq
	fastqc --noextract SMGill.trim.fq

SMGill.trim.fq.fq_fastqc-pbs:
	echo make SMGill.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-SMGill

SMHS.trim.fq.fq_fastqc.zip: SMHS.trim.fq
	fastqc --noextract SMHS.trim.fq

SMHS.trim.fq.fq_fastqc-pbs:
	echo make SMHS.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-SMHS

SN.trim.fq.fq_fastqc.zip: SN.trim.fq
	fastqc --noextract SN.trim.fq

SN.trim.fq.fq_fastqc-pbs:
	echo make SN.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-SN

SPDI.trim.fq.fq_fastqc.zip: SPDI.trim.fq
	fastqc --noextract SPDI.trim.fq

SPDI.trim.fq.fq_fastqc-pbs:
	echo make SPDI.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-SPDI

spinal_cord_paired.trim.fq.fq_fastqc.zip: spinal_cord_paired.trim.fq
	fastqc --noextract spinal_cord_paired.trim.fq

spinal_cord_paired.trim.fq.fq_fastqc-pbs:
	echo make spinal_cord_paired.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-spinal_cord_paired

SPK.trim.fq.fq_fastqc.zip: SPK.trim.fq
	fastqc --noextract SPK.trim.fq

SPK.trim.fq.fq_fastqc-pbs:
	echo make SPK.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-SPK

SPPI.trim.fq.fq_fastqc.zip: SPPI.trim.fq
	fastqc --noextract SPPI.trim.fq

SPPI.trim.fq.fq_fastqc-pbs:
	echo make SPPI.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-SPPI

Stg-1-Intestine.trim.fq.fq_fastqc.zip: Stg-1-Intestine.trim.fq
	fastqc --noextract Stg-1-Intestine.trim.fq

Stg-1-Intestine.trim.fq.fq_fastqc-pbs:
	echo make Stg-1-Intestine.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-Stg-1-Intestine

Stg-1-Kidney.trim.fq.fq_fastqc.zip: Stg-1-Kidney.trim.fq
	fastqc --noextract Stg-1-Kidney.trim.fq

Stg-1-Kidney.trim.fq.fq_fastqc-pbs:
	echo make Stg-1-Kidney.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-Stg-1-Kidney

Stg-2-Intestine.trim.fq.fq_fastqc.zip: Stg-2-Intestine.trim.fq
	fastqc --noextract Stg-2-Intestine.trim.fq

Stg-2-Intestine.trim.fq.fq_fastqc-pbs:
	echo make Stg-2-Intestine.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-Stg-2-Intestine

Stg-2-Kidney.trim.fq.fq_fastqc.zip: Stg-2-Kidney.trim.fq
	fastqc --noextract Stg-2-Kidney.trim.fq

Stg-2-Kidney.trim.fq.fq_fastqc-pbs:
	echo make Stg-2-Kidney.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-Stg-2-Kidney

Stg-3-Intestine.trim.fq.fq_fastqc.zip: Stg-3-Intestine.trim.fq
	fastqc --noextract Stg-3-Intestine.trim.fq

Stg-3-Intestine.trim.fq.fq_fastqc-pbs:
	echo make Stg-3-Intestine.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-Stg-3-Intestine

Stg-3-Kidney.trim.fq.fq_fastqc.zip: Stg-3-Kidney.trim.fq
	fastqc --noextract Stg-3-Kidney.trim.fq

Stg-3-Kidney.trim.fq.fq_fastqc-pbs:
	echo make Stg-3-Kidney.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-Stg-3-Kidney

Stg-4-Intestine.trim.fq.fq_fastqc.zip: Stg-4-Intestine.trim.fq
	fastqc --noextract Stg-4-Intestine.trim.fq

Stg-4-Intestine.trim.fq.fq_fastqc-pbs:
	echo make Stg-4-Intestine.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-Stg-4-Intestine

Stg-4-Kidney.trim.fq.fq_fastqc.zip: Stg-4-Kidney.trim.fq
	fastqc --noextract Stg-4-Kidney.trim.fq

Stg-4-Kidney.trim.fq.fq_fastqc-pbs:
	echo make Stg-4-Kidney.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-Stg-4-Kidney

Stg-5-Intestine.trim.fq.fq_fastqc.zip: Stg-5-Intestine.trim.fq
	fastqc --noextract Stg-5-Intestine.trim.fq

Stg-5-Intestine.trim.fq.fq_fastqc-pbs:
	echo make Stg-5-Intestine.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-Stg-5-Intestine

Stg-5-Kidney.trim.fq.fq_fastqc.zip: Stg-5-Kidney.trim.fq
	fastqc --noextract Stg-5-Kidney.trim.fq

Stg-5-Kidney.trim.fq.fq_fastqc-pbs:
	echo make Stg-5-Kidney.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-Stg-5-Kidney

Stg-6-Intestine.trim.fq.fq_fastqc.zip: Stg-6-Intestine.trim.fq
	fastqc --noextract Stg-6-Intestine.trim.fq

Stg-6-Intestine.trim.fq.fq_fastqc-pbs:
	echo make Stg-6-Intestine.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-Stg-6-Intestine

Stg-6-Kidney.trim.fq.fq_fastqc.zip: Stg-6-Kidney.trim.fq
	fastqc --noextract Stg-6-Kidney.trim.fq

Stg-6-Kidney.trim.fq.fq_fastqc-pbs:
	echo make Stg-6-Kidney.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-Stg-6-Kidney

Stg-7-Intestine.trim.fq.fq_fastqc.zip: Stg-7-Intestine.trim.fq
	fastqc --noextract Stg-7-Intestine.trim.fq

Stg-7-Intestine.trim.fq.fq_fastqc-pbs:
	echo make Stg-7-Intestine.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-Stg-7-Intestine

Stg-7-Kidney.trim.fq.fq_fastqc.zip: Stg-7-Kidney.trim.fq
	fastqc --noextract Stg-7-Kidney.trim.fq

Stg-7-Kidney.trim.fq.fq_fastqc-pbs:
	echo make Stg-7-Kidney.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-Stg-7-Kidney

SW2-Gill.trim.fq.fq_fastqc.zip: SW2-Gill.trim.fq
	fastqc --noextract SW2-Gill.trim.fq

SW2-Gill.trim.fq.fq_fastqc-pbs:
	echo make SW2-Gill.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-SW2-Gill

SW2-Intestine-a.trim.fq.fq_fastqc.zip: SW2-Intestine-a.trim.fq
	fastqc --noextract SW2-Intestine-a.trim.fq

SW2-Intestine-a.trim.fq.fq_fastqc-pbs:
	echo make SW2-Intestine-a.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-SW2-Intestine-a

SW2-Intestine-b.trim.fq.fq_fastqc.zip: SW2-Intestine-b.trim.fq
	fastqc --noextract SW2-Intestine-b.trim.fq

SW2-Intestine-b.trim.fq.fq_fastqc-pbs:
	echo make SW2-Intestine-b.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-SW2-Intestine-b

SW2-Kidney.trim.fq.fq_fastqc.zip: SW2-Kidney.trim.fq
	fastqc --noextract SW2-Kidney.trim.fq

SW2-Kidney.trim.fq.fq_fastqc-pbs:
	echo make SW2-Kidney.trim.fq.fq_fastqc.zip | cat fastqc.inc.sub - | qsub -l $(QCRES) -N fqc2-SW2-Kidney

TRIMMED=614GMAAXX_1_1_pf.trim.fq \
			614GMAAXX_2_1_pf.trim.fq \
			614GMAAXX_3_1_pf.trim.fq \
			614GMAAXX_4_1_pf.trim.fq \
			614GMAAXX_6_1_pf.trim.fq \
			614GMAAXX_7_1_pf.trim.fq \
			614GMAAXX_8_1_pf.trim.fq \
			AB057JABXX_s_7_pe.trim.fq \
			AI.trim.fq \
			AK.trim.fq \
			BB0582ABXX_s_3_pe.trim.fq \
			BB0582ABXX_s_4_pe.trim.fq \
			BB0582ABXX_s_5_pe.trim.fq \
			BB0582ABXX_s_6_pe.trim.fq \
			BB0582ABXX_s_7_pe.trim.fq \
			BB0582ABXX_s_8_pe.trim.fq \
			brain_paired.trim.fq \
			FW1-Gill.trim.fq \
			FW1-Intestine.trim.fq \
			FW1-Kidney.trim.fq \
			G81-a.trim.fq \
			G81-b.trim.fq \
			G82-a.trim.fq \
			G82-b.trim.fq \
			G83-a.trim.fq \
			G83-b.trim.fq \
			Gill.trim.fq \
			JV1_liver.trim.fq \
			JVa_liver.trim.fq \
			JVb_liver.trim.fq \
			JV_intest.trim.fq \
			Kidney.trim.fq \
			L81-a.trim.fq \
			L81-b.trim.fq \
			L82-a.trim.fq \
			L82-b.trim.fq \
			L83-1-a.trim.fq \
			L83-1-b.trim.fq \
			L83-a.trim.fq \
			L83-b.trim.fq \
			LI.trim.fq \
			Lips.trim.fq \
			LK.trim.fq \
			M0_liver.trim.fq \
			M2_liver.trim.fq \
			M5_liver.trim.fq \
			Monocytes.trim.fq \
			mssm_single_s1.trim.fq \
			mssm_single_s2.trim.fq \
			mssm_single_s3.trim.fq \
			mssm_single_s5.trim.fq \
			mssm_single_s6.trim.fq \
			mssm_single_s7.trim.fq \
			Muscle.trim.fq \
			Neutrophils.trim.fq \
			OFHS.trim.fq \
			POF-eye.trim.fq \
			POF-ts.trim.fq \
			PSMGill.trim.fq \
			Rope.trim.fq \
			SMGill.trim.fq \
			SMHS.trim.fq \
			SN.trim.fq \
			SPDI.trim.fq \
			spinal_cord_paired.trim.fq \
			SPK.trim.fq \
			SPPI.trim.fq \
			Stg-1-Intestine.trim.fq \
			Stg-1-Kidney.trim.fq \
			Stg-2-Intestine.trim.fq \
			Stg-2-Kidney.trim.fq \
			Stg-3-Intestine.trim.fq \
			Stg-3-Kidney.trim.fq \
			Stg-4-Intestine.trim.fq \
			Stg-4-Kidney.trim.fq \
			Stg-5-Intestine.trim.fq \
			Stg-5-Kidney.trim.fq \
			Stg-6-Intestine.trim.fq \
			Stg-6-Kidney.trim.fq \
			Stg-7-Intestine.trim.fq \
			Stg-7-Kidney.trim.fq \
			SW2-Gill.trim.fq \
			SW2-Intestine-a.trim.fq \
			SW2-Intestine-b.trim.fq \
			SW2-Kidney.trim.fq


abundfilt.time: diginorm.c50.ht
	$(TIME) -o abundfilt.time python $(KHMER)/scripts/filter-abund.py -V -Z 20 diginorm.c50.ht \
614GMAAXX_1_1_pf.trim.fq.keep \
			614GMAAXX_2_1_pf.trim.fq.keep \
			614GMAAXX_3_1_pf.trim.fq.keep \
			614GMAAXX_4_1_pf.trim.fq.keep \
			614GMAAXX_6_1_pf.trim.fq.keep \
			614GMAAXX_7_1_pf.trim.fq.keep \
			614GMAAXX_8_1_pf.trim.fq.keep \
			AB057JABXX_s_7_pe.trim.fq.keep \
			AI.trim.fq.keep \
			AK.trim.fq.keep \
			BB0582ABXX_s_3_pe.trim.fq.keep \
			BB0582ABXX_s_4_pe.trim.fq.keep \
			BB0582ABXX_s_5_pe.trim.fq.keep \
			BB0582ABXX_s_6_pe.trim.fq.keep \
			BB0582ABXX_s_7_pe.trim.fq.keep \
			BB0582ABXX_s_8_pe.trim.fq.keep \
			brain_paired.trim.fq.keep \
			FW1-Gill.trim.fq.keep \
			FW1-Intestine.trim.fq.keep \
			FW1-Kidney.trim.fq.keep \
			G81-a.trim.fq.keep \
			G81-b.trim.fq.keep \
			G82-a.trim.fq.keep \
			G82-b.trim.fq.keep \
			G83-a.trim.fq.keep \
			G83-b.trim.fq.keep \
			Gill.trim.fq.keep \
			JV1_liver.trim.fq.keep \
			JVa_liver.trim.fq.keep \
			JVb_liver.trim.fq.keep \
			JV_intest.trim.fq.keep \
			Kidney.trim.fq.keep \
			L81-a.trim.fq.keep \
			L81-b.trim.fq.keep \
			L82-a.trim.fq.keep \
			L82-b.trim.fq.keep \
			L83-1-a.trim.fq.keep \
			L83-1-b.trim.fq.keep \
			L83-a.trim.fq.keep \
			L83-b.trim.fq.keep \
			LI.trim.fq.keep \
			Lips.trim.fq.keep \
			LK.trim.fq.keep \
			M0_liver.trim.fq.keep \
			M2_liver.trim.fq.keep \
			M5_liver.trim.fq.keep \
			Monocytes.trim.fq.keep \
			mssm_single_s1.trim.fq.keep \
			mssm_single_s2.trim.fq.keep \
			mssm_single_s3.trim.fq.keep \
			mssm_single_s5.trim.fq.keep \
			mssm_single_s6.trim.fq.keep \
			mssm_single_s7.trim.fq.keep \
			Muscle.trim.fq.keep \
			Neutrophils.trim.fq.keep \
			OFHS.trim.fq.keep \
			POF-eye.trim.fq.keep \
			POF-ts.trim.fq.keep \
			PSMGill.trim.fq.keep \
			Rope.trim.fq.keep \
			SMGill.trim.fq.keep \
			SMHS.trim.fq.keep \
			SN.trim.fq.keep \
			SPDI.trim.fq.keep \
			spinal_cord_paired.trim.fq.keep \
			SPK.trim.fq.keep \
			SPPI.trim.fq.keep \
			Stg-1-Intestine.trim.fq.keep \
			Stg-1-Kidney.trim.fq.keep \
			Stg-2-Intestine.trim.fq.keep \
			Stg-2-Kidney.trim.fq.keep \
			Stg-3-Intestine.trim.fq.keep \
			Stg-3-Kidney.trim.fq.keep \
			Stg-4-Intestine.trim.fq.keep \
			Stg-4-Kidney.trim.fq.keep \
			Stg-5-Intestine.trim.fq.keep \
			Stg-5-Kidney.trim.fq.keep \
			Stg-6-Intestine.trim.fq.keep \
			Stg-6-Kidney.trim.fq.keep \
			Stg-7-Intestine.trim.fq.keep \
			Stg-7-Kidney.trim.fq.keep \
			SW2-Gill.trim.fq.keep \
			SW2-Intestine-a.trim.fq.keep \
			SW2-Intestine-b.trim.fq.keep \
			SW2-Kidney.trim.fq.keep
