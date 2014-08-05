from peasoup.tasks import BlastTask, BlastFormatTask, CurlTask, GunzipTask

blast_params = '-best_hit_score_edge 0.05 -best_hit_overhang 0.25 -max_target_seqs 1'
blast_threads = 8

assem_url = ('athyra.ged.msu.edu/~cswelcher/lamp10.fasta.gz', 'lamp10.fasta.gz')
assem = assem_url[1].rstrip('.gz')

db_urls = [
('http://hgdownload.soe.ucsc.edu/goldenPath/petMar2/bigZips/petMar2.fa.gz ','petMar2.fa.gz'),
#('http://hgdownload.soe.ucsc.edu/goldenPath/petMar2/bigZips/petMar2.fa.masked.gz', 'petMar2.fa.masked.gz'),
('http://hgdownload.soe.ucsc.edu/goldenPath/petMar2/bigZips/est.fa.gz', 'petMar2.est.fa.gz'),
('http://hgdownload.soe.ucsc.edu/goldenPath/petMar2/bigZips/mrna.fa.gz', 'petMar2.mrna.fa.gz'),
('ftp://ftp.ensembl.org/pub/release-72/fasta/petromyzon_marinus/cdna/Petromyzon_marinus.Pmarinus_7.0.72.cdna.all.fa.gz', 'petMar2.cdna.fa.gz'),
('ftp://ftp.ensembl.org/pub/release-71/fasta/mus_musculus/pep/Mus_musculus.GRCm38.71.pep.all.fa.gz','musMus.pep.fa.gz'),
('ftp://ftp.ensembl.org/pub/release-71/fasta/danio_rerio/pep/Danio_rerio.Zv9.71.pep.all.fa.gz','danRer.pep.fa.gz'),
('ftp://ftp.ncbi.nih.gov/pub/UniVec/UniVec_Core','uniVec.fa'),
('f#tp://ftp.ensembl.org/pub/release-74/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.74.pep.all.fa.gz','homSap.pep.fa.gz'),
('ftp://ftp.ensembl.org/pub/release-74/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.74.cdna.all.fa.gz','homSap.cdna.fa.gz'),
]

#
# Get the reference databases
#

get_dbs = CurlTask(db_urls.append(assem_url))

#
# Gunzip the downloaded databases
#

gunzip_dbs = GunzipTask([(src,src.rstrip('.gz')) for src in get_dbs.outputs() \
                    if src.endswith('.gz')])

#
# run makeblastdb on downloaded databases
#

mkdb_tasks = []
for src in gunzip_dbs.outputs():
    db_type='prot' if 'pep' in src else 'nucl'
    mkdb_tasks.append(BlastFormatTask(src, '{}.db'.format(src), db_type))
mkdb_tasks.append(BlastFormatTask(assem, '{}.db'.format(assem), 'nucl'))

def task_prep_databases():
    global databases


    yield get_dbs.tasks()
    yield gunzip_dbs.tasks()

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

