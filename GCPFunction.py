# GCPFunction.py
#

import time
import pandas as pd
from google.cloud import storage
import dalmatian as dm
import numpy as np
import os
import signal
import re
from JKBio import Helper as h


def list_blobs_with_prefix(bucket_name, prefix, delimiter=None):
    """Lists all the blobs in the bucket that begin with the prefix.

    This can be used to list all blobs in a "folder", e.g. "public/".

    The delimiter argument can be used to restrict the results to only the
    "files" in the given "folder". Without the delimiter, the entire tree under
    the prefix is returned. For example, given these blobs:

            /a/1.txt
            /a/b/2.txt

    If you just specify prefix = '/a', you'll get back:

            /a/1.txt
            /a/b/2.txt

    However, if you specify prefix='/a' and delimiter='/', you'll get back:

            /a/1.txt

    """
    storage_client = storage.Client()
    bucket = storage_client.get_bucket(bucket_name)
    ret = []
    blobs = bucket.list_blobs(prefix=prefix, delimiter=delimiter)
    for blob in blobs:
        ret.append(blob.name)
    return(ret)


def mvFiles(files, location):
    """
    move a set of files in parallel (when the set is huge)

    Args:
    ----
            files: gs paths
            location: to move the files to
    """
    by = len(files) if len(files) < 50 else 50
    for sfiles in h.grouped(files, by):
        a = ''
        for val in sfiles:
            a += val + ' '
        code = os.system("gsutil -m mv " + a + location)
        if code != 0:
            print('pressed ctrl+c or command failed')
            break


def lsFiles(files, add='', group=50):
    """
    move a set of files in parallel (when the set is huge)

    Args:
    ----
            files: gs paths
            add: additional params to add
    """
    by = len(files) if len(files) < group else group
    res = []
    for sfiles in h.grouped(files, by):
        a = ''
        for val in sfiles:
            a += val + ' '
        data = os.popen("gsutil -m ls " + add + " " + a)
        if data == signal.SIGINT:
            print('Awakened')
            break
        else:
            res += data.read().split('\n')[:-1]
            if "TOTAL:" in res[-1]:
                res = res[:-1]
    return res


def cpFiles(files, location):
    """
    copy a set of files in parallel (when the set is huge)

    Args:
    ----
            files: gs paths
            location to copy
    """
    by = len(files) if len(files) < 50 else 50
    for sfiles in h.grouped(files, by):
        a = ''
        for val in sfiles:
            a += val + ' '
        code = os.system("gsutil -m cp " + a + location)
        if code != 0:
            print('pressed ctrl+c or command failed')
            break


def rmFiles(files):
    """
    remove a set of files in parallel (when the set is huge)

    Args:
    ----
            files: gs paths
    """
    by = len(files) if len(files) < 50 else 50
    for sfiles in h.grouped(files, by):
        a = ''
        for val in sfiles:
            a += val + ' '
        code = os.system("gsutil -m rm " + a)
        if code != 0:
            print('pressed ctrl+c or command failed')
            break


def recoverFiles(files):
    """
    move a set of files in parallel (when the set is huge)

    Args:
    ----
            files: gs paths
            location: to move the files to
    """
    cmd = ['gsutil mv ' + f + ' ' + f.split('#')[0] for f in files]
    h.parrun(cmd, cores=cores)


def patternRN(rename_dict, location, wildcards, types=[], test=True, cores=1):
    """

    Args:
        wildcards: list[str] can be one of  ['**', '.*', '*.','-.*']
    """
    r = 0
    for k, v in rename_dict.items():
        loc = location
        if '**' in wildcards:
            loc += '**/'
        if '*.' in wildcards or '-.*' in wildcards:
            loc += '*'
        loc += k
        if '.*' in wildcards or '-.*' in wildcards:
            loc += '*'
        res = os.popen('gsutil -m ls ' + loc).read().split('\n')[:-1]
        print('found ' + str(len(res)) + ' files to rename')
        if '-.*' in wildcards:
            cmd = ["gsutil mv " + val + " " + '/'.join(val.split('/')[:-1]) + '/' + v + '.' + '.'.join(val.split('/')[-1].split('.')[1:]) for val in res]
        else:
            cmd = ["gsutil mv " + val + " " + val.replace(k, v) for val in res]
        if test:
            print(cmd)
        else:
            h.parrun(cmd, cores=cores)


def get_all_sizes(folder, suffix='*'):
    """

    will sort and list all the files by their sizes. 

    If some files have the same size, will list them together

    Args:
    ----
            folder: gs folder path
            suffix: of a specific file type

    Returns:
    -------
            dict(sizes:[paths])
    """
    samples = os.popen('gsutil -m ls -al ' + folder + '**.' + suffix).read().split('\n')
    # compute size filepath
    sizes = {'gs://' + val.split('gs://')[1].split('#')[0]: int(re.split("\d{4}-\d{2}-\d{2}", val)[0]) for val in samples[:-2]}
    names = {}
    for k, val in sizes.items():
        if val in names:
            names[val].append(k)
        else:
            names[val] = [k]
    if names == {}:
        # we didn't find any valid file paths
        print("We didn't find any valid file paths in folder: " + str(folder))
    return names


def exists(val):
    """
    tells if a gcp path exists
    """
    return os.popen('gsutil ls ' + val).read().split('\n')[0] == val


def extractSize(val):
    """
    extract the size from the string returned by an ls -l|a command
    """
    return int(re.split("\d{4}-\d{2}-\d{2}", val)[0])


def extractPath(val):
    """
    extract the path from the string returned by an ls -l|a command
    """
    return 'gs://' + val.split('gs://')[1].split('#')[0]


def extractBucket(gspath):
    """
    extract the bucket from the string google storage filepath
    """
    return gspath.split('gs://')[1].split('/')[0]
