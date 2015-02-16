#!/usr/bin/env python

from peasoup import configtools
from peasoup import task_funcs
from peasoup.task_funcs import run_tasks, apply_run_task
import pandas as pd
import os


# Get metadata from config file
metadata = configtools.get_cfg('metadata.ini', 'metadata.spec.ini')

# Convert into more useful Pandas format
db_metadata = pd.concat([pd.DataFrame(metadata['urls'].dict()).transpose(), pd.DataFrame(metadata['queries'].dict()).transpose()])
db_metadata['fn'] = db_metadata.dest.apply(lambda s: s.rstrip('.gz'))

sample_metadata = pd.DataFrame(metadata['samples'].dict()).transpose()

class Chdir:         
    def __init__( self, newPath ):  
        self.savedPath = os.getcwd()
        os.chdir(newPath)
    def __del__( self ):
        os.chdir( self.savedPath )

#Chdir('work')
os.chdir('work')

print os.getcwd()

tasks = []

@task_funcs.create_task_object
def bowtie2_build_task(row):
    return task_funcs.bowtie2_build_task(row.fn, row.fn.strip('.fasta'))

build_task = bowtie2_build_task(db_metadata.ix['assembly'])
index_basename_fn = build_task.targets[-1]
tasks.extend([build_task])

@task_funcs.create_task_object
def split_pairs_task(row):
    return task_funcs.split_pairs_task(row.filename)

split_tasks = sample_metadata[sample_metadata.paired == True].apply(split_pairs_task, axis=1, reduce=False)
tasks.extend(list(split_tasks))

@task_funcs.create_task_object
def bowtie2_align_task(row, idx=index_basename_fn):
    if row.paired == False:
        return task_funcs.bowtie2_align_task(idx, row.filename + '.x.lamp03', singleton_fn=row.filename)
    else:
        return task_funcs.bowtie2_align_task(idx, row.filename + '.x.lamp03', left_fn=row.filename+'.1',
                                                right_fn=row.filename + '.2')

align_tasks = sample_metadata.apply(bowtie2_align_task, axis=1, reduce=False)
print list(align_tasks)
tasks.extend(list(align_tasks))

run_tasks(tasks, n_cpus=16)

