# Jeremie Kalfon
# for BroadInsitute
# in 2019

from __future__ import print_function
from matplotlib import pyplot as plt
import json
import os
import sys
import string
import subprocess

import pdb
import ipdb
import pandas as pd
import numpy as np
import itertools
import random

from taigapy import TaigaClient
tc = TaigaClient()


rename_mut = {'contig': 'chr', 'position': 'pos', 'Reference_Allele': 'ref', 'ref_allele': 'ref', 'alt_allele': 'alt',
              'Chromosome': 'chr', 'End_postition': 'end', 'Start_position': 'pos', 'Tumor_Seq_Allele1': "alt"}


def fileToList(filename):
    """
    loads an input file with a\\n b\\n.. into a list [a,b,..]
    """
    with open(filename) as f:
        return [val[:-1] for val in f.readlines()]


def listToFile(l, filename):
    """
    loads a list with [a,b,..] into an input file a\\n b\\n..
    """
    with open(filename, 'w') as f:
        for item in l:
            f.write("%s\n" % item)


def dictToFile(d, filename):
    """
    turn a dict into a json file
    """
    with open(filename, 'w') as json_file:
        json.dump(d, json_file)


def fileToDict(filename):
    """
    loads a json file into a python dict
    """
    with open(filename) as f:
        data = json.load(f)
    return data


def batchMove(l, pattern=['*.', '.*'], folder='', add=''):
    """
    moves a set of files l into a folder:

    Args:
    -----
      l: file list
      pattern: if files are a set of patterns to match
      folder: folder to move file into
      add: some additional mv parameters
    """
    for val in l:
        cmd = 'mv '
        if add:
            cmd += add + ' '
        if '*.' in pattern:
            cmd += '*'
        cmd += val
        if '.*' in pattern:
            cmd += '*'
        cmd += " " + folder
        res = os.system(cmd)
        if res != 0:
            raise Exception("Leave command pressed or command failed")


def batchRename(dt, folder='', sudo=False, doAll=False, add='', dryrun=False):
    """
    Given a dict renames corresponding files in a folder

    Args:
    ----
      dt: dict(currentName:newName) renaming dictionnary
      folder: folder to look into
      add: some additional mv parameters
    """
    cmd = 'ls -R ' + folder if doAll else 'ls ' + folder
    files = os.popen(cmd).read().split('\n')
    if doAll:
        prep=''
        f = []
        for val in files:
            if len(val)==0:
                prep=''
                continue
            if val[0]=='.' and len(val)>3:
                prep=val[:-1]
                continue
            if "." in val:
                f.append(prep+"/"+val)
        files = f
    for k, val in dt.items():
        for f in files:
            if k in f:
                cmd = 'sudo mv ' if sudo else 'mv '
                if add:
                    cmd += add + ' '
                if not doAll:
                    cmd += folder 
                cmd += f
                cmd += ' '
                if not doAll:
                    cmd += folder
                cmd += f.replace(k, val)
                if dryrun:
                    print(cmd)
                else:
                    res = os.system(cmd)
                    if res != 0:
                        raise Exception("Leave command pressed or command failed")


def grouped(iterable, n):
    """
    iterate over element of list 2 at a time python
    s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ...
    """
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, n))
        if not chunk:
            return
        yield chunk


def overlap(interval1, interval2):
    """
    Given [0, 4] and [1, 10] returns [1, 4]
    Given [0, 4] and [8, 10] returns False
    """
    if interval2[0] <= interval1[0] <= interval2[1]:
        start = interval1[0]
    elif interval1[0] <= interval2[0] <= interval1[1]:
        start = interval2[0]
    else:
        return False

    if interval2[0] <= interval1[1] <= interval2[1]:
        end = interval1[1]
    elif interval1[0] <= interval2[1] <= interval1[1]:
        end = interval2[1]
    else:
        return False

    return (start, end)


def union(interval1, interval2):
    """
    Given [0, 4] and [1, 10] returns [0, 10]
    Given [0, 4] and [8, 10] returns False
    """
    if interval1[0] <= interval2[0] <= interval1[1]:
        start = interval1[0]
        end = interval1[1] if interval2[1] <= interval1[1] else interval2[1]
    elif interval1[0] <= interval2[1] <= interval1[1]:
        start = interval2[0] if interval2[0] <= interval1[0] else interval1[0]
        end = interval1[1]
    else:
        return False
    return (start, end)


def nans(df): return df[df.isnull().any(axis=1)]


def createFoldersFor(filepath):
    """
    will recursively create folders if needed until having all the folders required to save the file in this filepath
    """
    prevval = ''
    for val in filepath.split('/')[:-1]:
        prevval += val + '/'
        if not os.path.exists(prevval):
            os.mkdir(prevval)


def randomString(stringLength=6, stype='all', withdigits=True):
    """
    Generate a random string of letters and digits

    Args:
    -----
      stringLength: the amount of char
      stype: one of lowercase, uppercase, all
      withdigits: digits allowed in the string?

    Returns:
    -------
      the string
    """
    if stype == 'lowercase':
        lettersAndDigits = ascii_lowercase
    elif stype == 'uppercase':
        lettersAndDigits = ascii_uppercase
    else:
        lettersAndDigits = string.ascii_letters
    if withdigits:
        lettersAndDigits += string.digits
    return ''.join(random.choice(lettersAndDigits) for i in range(stringLength))


def pdDo(df, op="mean", of="value1", over="value2"):
    """
    apply a function to a panda dataframe WIP
    """
    df = df.sort_values(by=over)
    index = []
    data = df.iloc[0, of]
    ret = []
    prev = df.iloc[0, over]
    j = 0
    for k, val in df.iloc[1:].iterrows():
        if val[over] == prev:
            data.append(val[of])
        else:
            if of == "mean":
                ret[j] = np.mean(data)
            elif of == "sum":
                ret[j] = np.sum(data)
            elif of == "max":
                ret[j] = np.max(data)
            elif of == "min":
                ret[j] = np.min(data)
            index.append(k)
            j += 1
            data = [val[of]]
    return index, ret


def parrun(cmds, cores, add=[]):
    """
    runs a set of commands in parallel using the "&" command

    Args:
    -----
      cmds: the list of commands
      cores: number of parallel execution
      add: an additional list(len(cmds)) of command to run in parallel at the end of each parallel run
    """
    count = 0
    exe = ''
    if len(add) != 0 and len(add) != len(cmds):
        raise ValueError("we would want them to be the same size")
    else:
        addexe = ''
    for i, cmd in enumerate(cmds):
        count += 1
        exe += cmd
        if len(add) != 0:
            addexe += add[i]
        if count < cores and i < len(cmds) - 1:
            exe += ' & '
            if len(add) != 0:
                addexe += ' & '
        else:
            count = 0
            res = subprocess.run(exe, capture_output=True, shell=True)
            if res.returncode != 0:
                raise ValueError('issue with the command: ' + str(res.stderr))
            exe = ''
            if len(add) != 0:
                res = subprocess.run(addexe, capture_output=True, shell=True)
                if res.returncode != 0:
                    raise ValueError(
                        'issue with the command: ' + str(res.stderr))
                addexe = ''


def askif(quest):
    """
    asks a y/n question to the user about something and returns true or false given his answer
    """
    print(quest)
    inp = input()
    if inp in ['yes', 'y', 'Y', 'YES', 'oui', 'si']:
        return 1
    elif inp in ['n', 'no', 'nope', 'non', 'N']:
        return 0
    else:
        return askif('you need to answer by yes or no')


def inttodate(i, lim=1965, unknown='U', sep='-', order="asc", startsatyear=0):
    """
    transforms an int representing days into a date

    Args:
    ----
      i: the int
      lim: the limited year below which we have a mistake
      unknown: what to return when unknown (date is bellow the limited year)
      sep: the sep between your date (e.g. /, -, ...)
      order: if 'asc', do d,m,y else do y,m,d
      startsatyear: when is the year to start counting for this int

    Returns:
      the date or unknown
    """
    a = int(i // 365)
    if a > lim:
        a = str(a + startsatyear)
        r = i % 365
        m = str(int(r // 32)) if int(r // 32) > 0 else str(1)
        r = r % 32
        d = str(int(r)) if int(r) > 0 else str(1)
    else:
        return unknown
    return d + sep + m + sep + a if order == "asc" else a + sep + m + sep + d


def datetoint(dt, split='-', unknown='U', order="des"):
    """
    same as inttodate but in the opposite way;

    starts at 0y,0m,0d

    dt: the date string
    split: the splitter in the string (e.g. /,-,...)
    unknown: maybe the some dates are 'U' or 0 and the program will output 0 for unknown instead of crashing
    order: if 'asc', do d,m,y else do y,m,d

    Returns:
      the date
    """
    arr = np.array(dt[0].split(split) if dt[0] !=
                   unknown else [0, 0, 0]).astype(int)
    if len(dt) > 1:
        for val in dt[1:]:
            arr = np.vstack(
                (arr, np.array(val.split(split) if val != unknown and val.count(split) == 2 else [0, 0, 0]).astype(int)))
        arr = arr.T
    res = arr[2] * 365 + arr[1] * 31 + \
        arr[0] if order == "asc" else arr[0] * 365 + arr[1] * 31 + arr[2]
    return [res] if type(res) is np.int64 else res


def showcount(i, size):
    """
    pretty print of i/size%, to put in a for loop
    """
    print(str(1 + int(100 * (i / size))) + '%', end='\r')


def combin(n, k):
    """
    Nombre de combinaisons de n objets pris k a k
    Number of comabination of n object taken k at a time
    """
    if k > n // 2:
        k = n - k
    x = 1
    y = 1
    i = n - k + 1
    while i <= n:
        x = (x * i) // y
        y += 1
        i += 1
    return x


def dups(lst):
    seen = set()
    # adds all elements it doesn't know yet to seen and all other to seen_twice
    seen_twice = set(x for x in lst if x in seen or seen.add(x))
    # turn the set into a list (as requested)
    return list(seen_twice)

def makeCombinations(size, proba):
    """
    produces probability of X event happening at the same time given binomial proba of event occuring
    
    pretty usefull for cobinding analysis
    """
    sums = {i:0 for i in range(1,size)}
    for i in range(size-1,0,-1):
        print(i)
        if sums[i]> 0:
            continue
        print(combin(size+3,i))
        v=0
        for j in itertools.combinations(proba, i):
            v+=np.prod(j)
        sums[i] = v
    for i in range(size-1,0,-1):
        for j in range(i+1,size):
            icomb = combin(j,i)
            sums[i] -= icomb*sums[j]
    sums[0] = 1-sum(list(sums.values()))
    return sums

def closest(lst, K):
    return lst[min(range(len(lst)), key = lambda i: abs(lst[i]-K))] 


def compareDfs(df1, df2):
    """
    compares df1 to df2
    """
    nmissmatchCols = set(df1.columns)-set(df2.columns)
    omissmatchCols = set(df2.columns)-set(df1.columns)
    nmissmatchInds = set(df1.index)-set(df2.index)
    omissmatchInds = set(df2.index)-set(df1.index)
    newNAs = df1.isna().sum().sum() - df2.isna().sum().sum()
    new0s = (df1 == 0).sum().sum() - (df2 == 0).sum().sum()
    print('FOUND missmatch Columns IN df1: ' + str(nmissmatchCols))
    print('FOUND missmatch Columns NOT IN df1: ' + str(omissmatchCols))
    print('FOUND missmatch Index IN df1: ' + str(nmissmatchInds))
    print('FOUND missmatch Index NOT IN df1: ' + str(omissmatchInds))
    print('FOUND new NAs in df1: ' + str(newNAs))
    print('FOUND new 0s in df1: ' + str(new0s))
    return nmissmatchCols, omissmatchCols, nmissmatchInds, omissmatchInds, newNAs, new0s
