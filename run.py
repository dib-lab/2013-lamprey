from peasoup.tasks import BlastTask, BlastFormatTask, CurlTask, GunzipTask

db_urls = [
#('http://hgdownload.soe.ucsc.edu/goldenPath/petMar2/bigZips/petMar2.fa.gz ','petMar2.fa.gz'),
#('http://hgdownload.soe.ucsc.edu/goldenPath/petMar2/bigZips/petMar2.fa.masked.gz', 'petMar2.fa.masked.gz'),
('http://hgdownload.soe.ucsc.edu/goldenPath/petMar2/bigZips/est.fa.gz', 'petMar2.est.fa.gz'),
('http://hgdownload.soe.ucsc.edu/goldenPath/petMar2/bigZips/mrna.fa.gz ','petMar2.mrna.fa.gz'),
('ftp://ftp.ensembl.org/pub/release-72/fasta/petromyzon_marinus/cdna/Petromyzon_marinus.Pmarinus_7.0.72.cdna.all.fa.gz','petMar2.cdna.fa.gz'),
#('ftp://ftp.ensembl.org/pub/release-71/fasta/mus_musculus/pep/Mus_musculus.GRCm38.71.pep.all.fa.gz','musMus.pep.fa.gz'),
#('ftp://ftp.ensembl.org/pub/release-71/fasta/danio_rerio/pep/Danio_rerio.Zv9.71.pep.all.fa.gz','danRer.pep.fa.gz'),
#('ftp://ftp.ncbi.nih.gov/pub/UniVec/UniVec_Core','uniVec.fa'),
#('f#tp://ftp.ensembl.org/pub/release-74/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.74.pep.all.fa.gz','homSap.pep.fa.gz'),
#('ftp://ftp.ensembl.org/pub/release-74/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.74.cdna.all.fa.gz','homSap.cdna.fa.gz'),
]

#
# Get the reference databases
#
def task_all():
    get_dbs = CurlTask(db_urls)
    gunzip_dbs = GunzipTask([(src,src.rstrip('.gz')) for src in get_dbs.outputs() \
                        if src.endswith('.gz')])
    yield get_dbs.tasks()
    yield gunzip_dbs.tasks()
