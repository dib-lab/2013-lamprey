#!/usr/bin/env python

import os
import sys

from doit import get_var

from peasoup.tasks import BlastTask, BlastFormatTask, CurlTask, GunzipTask, \
                            UniProtQueryTask, TruncateFastaNameTask
from peasoup.configtools import get_cfg

cfg_fn = get_var('metadata')
if not cfg_fn:
    print >>sys.stderr, 'No argument to metadata=<file> specified, using metadata.ini'
    metadata = get_cfg('metadata.ini', 'metadata.spec.ini')
else:
    metadata = get_cfg(cfg_fn, 'metadata.spec.ini')

name = metadata['urls']['assembly']['name']
assem = metadata['urls']['assembly']['dest'].rstrip('.gz')
print >>sys.stderr, 'Running', name, 'pipeline\n', '*' * 40


#
# Get the reference databases
#

get_dbs = CurlTask([(e['url'], e['dest']) for e in metadata['urls'].values()])
uniprot_tasks = [UniProtQueryTask(e['terms'], e['dest']) for e in metadata['queries'].values() \
                    if e['q_type'] == 'uniprot']

#
# Gunzip the downloaded databases
#

uniprot_dbs = [gz.outputs().next() for gz in uniprot_tasks]
gunzip_dbs = GunzipTask([(src,src+'.tmp') for src in get_dbs.outputs() \
                    if src.endswith('.gz')] + \
                [(src,src+'.tmp') for src in uniprot_dbs])

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
    blast_threads = metadata['blast']['threads']
    blast_params = metadata['blast']['params']
    for task in mkdb_tasks:
        db_name, db_fn = task.outputs().next()
        db_type = 'prot' if 'pep' in db_name else 'nucl'
        if not db_name.startswith(assem):
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

