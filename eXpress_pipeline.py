#!/usr/bin/env python

from peasoup import configtools
from peasoup import task_funcs
from peasoup.task_funcs import run_tasks, apply_run_task
import pandas as pd
import os
import argparse

@task_funcs.create_task_object
def bowtie2_build_task(row):
    return task_funcs.bowtie2_build_task(row.fn, row.fn.split('.fasta')[0])

@task_funcs.create_task_object
def split_pairs_task(row):
    return task_funcs.split_pairs_task(row.filename)

@task_funcs.create_task_object
def bowtie2_align_task(row, idx_fn):
    target = '{sample}.x.{idx}'.format(sample=row.filename, idx=idx_fn)
    if row.paired == False:
        return task_funcs.bowtie2_align_task(idx_fn, target, singleton_fn=row.filename)
    else:
        return task_funcs.bowtie2_align_task(idx_fn, target, left_fn=row.filename+'.1',
                                                right_fn=row.filename + '.2', encoding='--phred64')

@task_funcs.create_task_object
def express_task(hits_fn, transcripts_fn):
    folder = hits_fn.split('.bam')[0]
    return task_funcs.eXpress_task(transcripts_fn, hits_fn, folder)

@task_funcs.create_task_object
def samtools_sort_task(hits_fn):
    return task_funcs.samtools_sort_task(hits_fn)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--metadata', default='metadata.ini')
    parser.add_argument('--metadata-spec', dest='metadata_spec',  default='metadata.spec.ini')
    parser.add_argument('-T', '--threads', default=4, type=int)
    parser.add_argument('--wdir', default='work')
    parser.add_argument('--print-tasks', dest='print_tasks', action='store_true', default=False)
    parser.add_argument('-D', '--dry-run', dest='dry_run', action='store_true', default=False)
    parser.add_argument('--align-only', dest='align_only', action='store_true', default=False)
    parser.add_argument('-C', '--clean', action='store_true', default=False)
    args = parser.parse_args()

    metadata = configtools.get_cfg(args.metadata, args.metadata_spec)
    db_df = configtools.cfg_to_dataframe(metadata['urls'])
    db_df['fn'] = db_df.dest.apply(lambda s: s.rstrip('.gz'))
    sample_df = configtools.cfg_to_dataframe(metadata['samples'])

    old_dir = os.getcwd()
    try:
        if not os.path.exists(args.wdir):
            os.makedirs(args.wdir)
        os.chdir(args.wdir)

        tasks = []
    
        build_task = bowtie2_build_task(db_df.ix['assembly'])
        index_basename_fn = build_task.targets[-1]
        tasks.extend([build_task])

        split_tasks = sample_df[sample_df.paired == True].apply(split_pairs_task, axis=1, reduce=False)
        tasks.extend(list(split_tasks))

        align_tasks = sample_df.apply(bowtie2_align_task, args=(index_basename_fn,), axis=1, reduce=False)
        tasks.extend(list(align_tasks))

        align_files_df = task_funcs.get_task_files_df(align_tasks)
        bam_files_to_sort = align_files_df[align_files_df.apply(lambda row: row.dep.endswith('.1'), axis=1)].target
        bam_files_single = align_files_df[align_files_df.apply(lambda row: row.dep.endswith('.fq.gz'), axis=1)].target
        
        sort_tasks = bam_files_to_sort.apply(samtools_sort_task)
        tasks.extend(sort_tasks)

        if not args.align_only:
            hits_files = bam_files_single.append(sort_tasks.apply(lambda t: t.targets[0]))
            express_tasks = hits_files.apply(express_task, args=(db_df.ix['assembly'].fn,))
            tasks.extend(express_tasks)

        if args.print_tasks:
            for task in tasks:
                print '-----\n', task, task.actions, task.targets, task.file_dep, task.dep_changed

        if not args.dry_run:
            run_tasks(tasks, n_cpus=args.threads, clean=args.clean)
        else:
            print 'Dry run; exiting...'
    finally:
        os.chdir(old_dir)


if __name__ == '__main__':
    main()
