# Jeremie Kalfon
# for BroadInsitute
# in 2019

from __future__ import print_function
import os
import signal
import re

import pandas as pd
import numpy as np

from taigapy import TaigaClient
tc = TaigaClient()

size = {"GRCh37": 2864785220,
		"GRCh38": 2913022398}

cmaps = ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
		 'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
		 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']

chroms = {'chr1','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
'chr2','chr20','chr21','chr22','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chrX',
 'chrY','1','10','11','12','13','14','15','16','17','18','19','2','20','21','22','3','4','5','6',
 '7','8','9','X','Y'}

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

		Returns:
		-------
			a dict[str:str] of newly indexed bam files with their assocciated bam index file
    """
    files = gcp.lsFiles([bucketpath])
    bams = [val for val in files if '.bam' in val[-4:]]
    unindexed = [val for val in bams if val[:-4]+'.bai' not in files and val[:4] +'.bam.bai' not in files]
    print("found "+str(len(unindexed))+" files to reindex")
    h.parrun(["export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token` && samtools index "+val for val in unindexed], cores)
    return {val: val[:-4]+".bam.bai" for val in unindexed}



def dropWeirdChromosomes(bedfile, keep=[], skip=0):
	"""
	given a bedfile path, removes chromosomes that are not one of chroms

	Args:
	----
		bedfile: str the filepath to the bedfile
		keep: list[str] of additional chromosomes to keep
	"""
	if skip>=20:
		raise ValueError('too many header lines!')
	try:
		bed = pd.read_csv(bedfile, sep='\t',header=None, skiprows=skip)
	except ParserError:
		dropWeirdChromosomes(bedfile, keep, skip+1)
		return
	except EmptyDataError:
		print("empty bed")
		return
	initlen= len(bed)
	if initlen ==0:
		print("empty bed")
		return
	bed = bed[bed[0].isin(chroms|set(keep))]
	if len(bed) < skip and skip > 5:
		raise ValueError('too many header lines!')
	print("found "+str(skip)+" header line... removing")
	if len(bed) != initlen:
		print('removed '+str(initlen-len(bed))+" lines")
	bed.to_csv(bedfile, sep='\t',header=None,index=None)


def extractPairedSingleEndFrom(folder, sep='-', namepos=2):
	"""
	given a folder, find fastq files and sorts paired and single end based on the R1/R2 patterns
		
	Args:
	-----
		folder: the folder where the fastqs are
		sep: the separator in filename
		namepos: the location of the name in this separated list of name from filepath

	Returns:
	-------
		list of filepath to single end files
		df with R1 and R2 filepath
	"""
	single = []
	paired = {}
	for val in os.listdir(folder):
		if ".fastq" in val or ".fq" in val:
			if 'R1' in val:
				name = val.split(sep)[namepos]
				paired[name] = {'R1': val}
			elif 'R2' in val:
				name = val.split(sep)[namepos]
				paired[name].update({'R2': val})
			else:
				single.append(val)
	return single, pd.DataFrame(paired)


def findReplicatesBams(folder, sep='', filetypes=['.bam'], namings='-r([0-9])', namepos=0):
	"""
	creates a dict of name and replicate files given a regexp naming scheme

	Args:
	-----
		folder: str the folderpath containing the bam files
		sep: str how to separate each filename (to extract the name)
		namepos: int which separated location to use (to extract the name)
		filetypes: list[str] if filename appendix to accept
		namings: str regexp string to extract the replicate number
	"""
	rep = {}
	for val in os.listdir(folder):
		if val[-4:] in filetypes:
			match = re.search(namings, val)
			if match:
				number = match.groups()[0]
				name = val.split(sep)[namepos]
				if name in rep:
					rep[name].append(val)
				else:
					rep[name] = [val]

	return rep


def mergeBams(rep, cores=3):
	"""
	uses samtools to merge a set of replicates considered into one file

	Args:
	-----
		rep: dict[str: list[str]] associating merged bam namepath to list of replicate bam filepath to be merged
		cores: int number of parallel threads to use
	"""
	runs = []
	for i, val in rep.items():
		out1 = i + '.merged.bam'
		for bam in val:
			in1 += ' ' + bam
		runs.append("samtools merge " + out1 + in1)
	h.parrun(runs, cores)
