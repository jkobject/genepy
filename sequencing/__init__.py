# Jeremie Kalfon
# for BroadInsitute
# in 2019

from __future__ import print_function
import os
import signal
import re

import pdb
import ipdb
import pandas as pd
import numpy as np

from taigapy import TaigaClient
tc = TaigaClient()

def fromGTF2BED(gtfname, bedname, gtftype='geneAnnot'):
    """
    transforms a  gtf file into a bed file

    Args:
    ----
      gtfname: filepath to gtf file
      bedname: filepath to beddfile
      gtftype: only geneAnnot for now

    Returns:
    --------
      newbed: the bedfile as a pandas.df

    """
    if gtftype == 'geneAnnot':
        gtf = pd.read_csv(gtfname, sep='\t', header=0, names=[
                          "chr", "val", "type", "start", 'stop', 'dot', 'strand', 'loc', 'name'])
        gtf['name'] = [i.split('gene_id "')[-1].split('"; trans')[0]
                       for i in gtf['name']]
        prevname = ''
        newbed = {'chr': [], 'start': [], 'end': [], 'gene': []}
        for i, val in gtf.iterrows():
            showcount(i, len(gtf))
            if val['name'] == prevname:
                newbed['end'][-1] = val['stop']
            else:
                newbed['chr'].append(val['chr'])
                newbed['start'].append(val['start'])
                newbed['end'].append(val['stop'])
                newbed['gene'].append(val['name'])
            prevname = val['name']
        newbed = pd.DataFrame(newbed)
        newbed = newbed[~newbed.chr.str.contains('_fix')]
        newbed.to_csv(bedname + ".bed", sep='\t', index=None)
        newbed.to_csv(bedname + "_genes.bed", sep='\t', index=None)
        return newbed


def getBamDate(bams, split='-', order="des", unknown='U'):
    """
    from bam files (could be in a google bucket) returns their likely sequencing date if available in the header

    Args:
    -----
      bams: the bams file|bucket paths 
      split: the splitter in the output date
      unknown: maybe the some dates can't be found the program will output unknown for them
      order: if 'asc', do d,m,y else do y,m,d

    Returns:
    -------
      a list of likely dates or [unknown]s
    """
    DTs = []
    for i, bam in enumerate(bams):
        print(i / len(bams), end='\r')
        data = os.popen('export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`\
       && samtools view -H ' + bam + ' | grep "^@RG"')
        if data == signal.SIGINT:
            print('Awakened')
            break
        else:
            res = data.read()
            dt = re.findall("(?<=\tDT:).+?\t", res)
        if len(dt) > 1:
            arr = np.array(dt[0].split('T')[0].split(split)).astype(int)
            for val in dt[1:]:
                arr = np.vstack(
                    (arr, np.array(val.split('T')[0].split(split)).astype(int)))
            arr = arr.T
            i = arr[0] * 365 + arr[1] * 31 + \
                arr[2] if order == "asc" else arr[2] * \
                365 + arr[1] * 31 + arr[0]
            DTs.append(dt[np.argsort(i)[0]].split('T')[0])
        elif len(dt) == 1:
            DTs.append(dt[0].split('T')[0])
        else:
            DTs.append(unknown)
    return DTs


def indexBams(bucketpath, cores=4):
    """
    given a bucket path, will index all .bam files without an associated index and return their paths
    """
    files = gcp.lsFiles([bucketpath])
    bams = [val for val in files if '.bam' in val[-4:]]
    unindexed = [val for val in bams if val[:-4]+'.bai' not in files and val[:4] +'.bam.bai' not in files]
    print("found "+str(len(unindexed))+" files to reindex")
    h.parrun(["export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token` && samtools index "+val for val in unindexed], cores)
    return {val: val[:-4]+".bam.bai" for val in unindexed}