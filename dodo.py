from peasoup.tasks import BlastTask, BlastFormatTask, CurlTask, GunzipTask, \
                            UniProtQueryTask, TruncateFastaNameTask

blast_params = '-best_hit_score_edge 0.05 -best_hit_overhang 0.25 -max_target_seqs 1'
blast_threads = 8

assem_url = ('http://athyra.ged.msu.edu/~cswelcher/lamprey/lamp10.fasta', 'lamp03.fasta.gz')
assem = assem_url[1].rstrip('.gz')

db_urls = [
('ftp://ftp.ensembl.org/pub/release-75/fasta/petromyzon_marinus/dna/Petromyzon_marinus.Pmarinus_7.0.75.dna_sm.toplevel.fa.gz','petMar2.fa.gz'),
('http://hgdownload.soe.ucsc.edu/goldenPath/petMar2/bigZips/est.fa.gz', 'petMar2.est.fa.gz'),
('ftp://ftp.ensembl.org/pub/release-75/fasta/petromyzon_marinus/pep/Petromyzon_marinus.Pmarinus_7.0.75.pep.all.fa.gz', 'petMar2.pep.fa.gz'),
('ftp://ftp.ensembl.org/pub/release-75/fasta/petromyzon_marinus/cds/Petromyzon_marinus.Pmarinus_7.0.75.cds.all.fa.gz', 'petMar2.cds.fa.gz'),
('ftp://ftp.ensembl.org/pub/release-75/fasta/petromyzon_marinus/ncrna/Petromyzon_marinus.Pmarinus_7.0.75.ncrna.fa.gz', 'petMar2.ncrna.fa.gz'),
('ftp://ftp.ensembl.org/pub/release-75/fasta/petromyzon_marinus/cdna/Petromyzon_marinus.Pmarinus_7.0.75.cdna.all.fa.gz', 'petMar2.cdna.fa.gz'),
('ftp://ftp.ensembl.org/pub/release-75/fasta/mus_musculus/pep/Mus_musculus.GRCm38.75.pep.all.fa.gz','musMus.pep.fa.gz'),
('ftp://ftp.ensembl.org/pub/release-75/fasta/danio_rerio/pep/Danio_rerio.Zv9.75.pep.all.fa.gz' ,'danRer.pep.fa.gz'),
#('ftp://ftp.ncbi.nih.gov/pub/UniVec/UniVec_Core','uniVec.fa'),
#('ftp://ftp.ensembl.org/pub/release-74/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.74.pep.all.fa.gz','homSap.pep.fa.gz'),
#('ftp://ftp.ensembl.org/pub/release-74/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.74.cdna.all.fa.gz','homSap.cdna.fa.gz'),
]

uniprot_terms = [('organism:7739+AND+keyword:1185', 'braFlo.pep.all.fa.gz'),
                 ('taxonomy:7762', 'Myx.pep.all.fa.gz')]

#
# Get the reference databases
#

get_dbs = CurlTask(db_urls + [assem_url])
uniprot_tasks = [UniProtQueryTask(org_id, fn) for org_id, fn in uniprot_terms]

#
# Gunzip the downloaded databases
#

uniprot_dbs = [gz.outputs().next() for gz in uniprot_tasks]
gunzip_dbs = GunzipTask([(src,src+'.tmp') for src in get_dbs.outputs() \
                    if src.endswith('.gz')] + \
                [(src,src+'.tmp') for src in uniprot_dbs])

print gunzip_dbs.links
print list(gunzip_dbs.outputs())

#
# Fix the names for BLAST+
#

fix_names = TruncateFastaNameTask([(fn,fn.rstrip('.gz.tmp')) for fn in gunzip_dbs.outputs()])

#
# run makeblastdb on downloaded databases
#

mkdb_tasks = []
for src in fix_names.outputs():
    db_type='prot' if 'pep' in src else 'nucl'
    mkdb_tasks.append(BlastFormatTask(src, '{}.db'.format(src), db_type))

def task_prep_databases():
    global databases


    yield get_dbs.tasks()
    for task in uniprot_tasks:
        yield task.tasks()

    yield gunzip_dbs.tasks()

    yield fix_names.tasks()

    for task in mkdb_tasks:
        yield task.tasks()

def task_blast():
    for task in mkdb_tasks:
        db_name, db_fn = task.outputs().next()
        db_type = 'prot' if 'pep' in db_name else 'nucl'
        if db_name != '{}.db'.format(assem):
            if db_type == 'prot':
                yield BlastTask('blastx', assem, db_name, 
                                '{0}.x.{1}.tsv'.format(assem, db_name),
                                num_threads=blast_threads,
                                params=blast_params).tasks()
                yield BlastTask('tblastn', db_name.rstrip('.db'), '{}.db'.format(assem),
                                '{0}.x.{1}.tsv'.format(db_name, assem),
                                num_threads=blast_threads,
                                params=blast_params).tasks()
            else:
                yield BlastTask('blastn', assem, db_name, 
                                '{0}.x.{1}.tsv'.format(assem, db_name),
                                num_threads=blast_threads,
                                params=blast_params).tasks()
                yield BlastTask('blastn', db_name.rstrip('.db'), '{}.db'.format(assem),
                                '{0}.x.{1}.tsv'.format(db_name, assem),
                                num_threads=blast_threads,
                                params=blast_params).tasks()

