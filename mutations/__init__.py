# Jeremie Kalfon
# for BroadInsitute
# in 2019

from __future__ import print_function

import pdb
import ipdb
import pandas as pd
import numpy as np

from taigapy import TaigaClient
tc = TaigaClient()

def vcf_to_df(path, hasfilter=False, samples=['sample'], additional_unique=[]):
    uniqueargs = ['DB', 'SOMATIC', 'GERMLINE', "OVERLAP", "IN_PON",
                  "STR", "ReverseComplementedAlleles"] + additional_unique

    def read_comments(f):
        fields = {}
        description = {}
        for l in f:
            l = l.decode("utf-8") if type(l) is not str else l
            if l.startswith('##'):
                if 'FORMAT' in l[:20]:
                    res = l.split('ID=')[1].split(',')[0]
                    desc = l.split('Description=')[1][:-2]
                    description.update({res: desc})
                if 'INFO' in l[:20]:
                    res = l.split('ID=')[1].split(',')[0]
                    desc = l.split('Description=')[1][:-2]
                    description.update({res: desc})
                    fields.update({res: []})
            else:
                break
        return fields, description
    if path.endswith('.gz'):
        with gzip.open(path, 'r') as f:
            fields, description = read_comments(f)
    else:
        with open(path, 'r') as f:
            fields, description = read_comments(f)
    names = ['chr', 'pos', 'id', 'ref', 'alt', 'qual']
    names += ['filter'] if hasfilter else ['strand']
    names += ['data', 'format'] + samples
    a = pd.read_csv(path, sep='\t', comment="#", header=None,
                    names=names, index_col=False)
    print(description)
    try:
        for j, val in enumerate(a.data.str.split(';').values.tolist()):
            res = dict([(v, True) if v in uniqueargs else tuple(
                v.split('=')) for v in val])
            for k in fields.keys():
                fields[k].append(res.get(k, None))
    except ValueError:
        print(val)
        raise ValueError('unknown field')
    a = pd.concat([a.drop(columns='data'), pd.DataFrame(
        data=fields, index=a.index)], axis=1)
    for sample in samples:
        sorting = a.format[0].split(':')
        res = a[sample].str.split(':').values.tolist()
        maxcols = max([len(v) for v in res])
        if maxcols - len(sorting) > 0:
            for i in range(maxcols - len(sorting)):
                sorting.append(sorting[-1] + '_' + str(i + 1))
        if len(samples) > 1:
            sorting = [sample + '_' + v for v in sorting]
        a = pd.concat([a.drop(columns=sample), pd.DataFrame(
            data=res, columns=sorting, index=a.index)], axis=1)
    return a.drop(columns='format'), description
