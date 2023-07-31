import os
import pandas as pd
import numpy as np
from genepy.utils import helper as h
import re
import pyBigWig
from scipy.stats import zscore, fisher_exact
import subprocess
from statsmodels.stats.multitest import multipletests

size = {"GRCh37": 2864785220, "GRCh38": 2913022398}

cmaps = ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
				'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
				'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']

chroms = {'chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
'chr2', 'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chrX',
 'chrY', '1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '2', '20', '21', '22', '3', '4', '5', '6',
 '7', '8', '9', 'X', 'Y'}


def bigWigFrom(bams, folder="", numthreads=8, genome='GRCh37', scaling=None, verbose=1):
	"""
	run the bigwig command line for a set of bam files in a folder

	Can apply some scaling and process files in parallel

	Args:
	----
		bams: list[str] of filepath to bam files
		folder: str folder where the bam files would be stored (will save the bigwigs there as well, else to current path ./)
		numthreadds: int number of threads to process in parallel
		genome: str genome version to use (for genome size) only GRCh37 GRCh38 (see size global variable in this file)
		scaling: list[float] a list of scaling values to apply to the bigwigs. if provided, won't scale to read counts.
		verbose: int verbose level

	"""
	if "bigwig" not in os.listdir(folder if folder else '.'):
		os.mkdir(folder + "bigwig")
	for i, bam in enumerate(bams):
		in1 = folder + bam
		out1 = folder + "bigwig/" + bam.split('/')[-1].split('.')[0] + '.bw'
		cmd = "bamCoverage --effectiveGenomeSize " + str(size[genome]) + " -p " + str(numthreads) +\
			" -b " + in1 + " -o " + out1
		if scaling is not None:
			cmd += ' --scaleFactor ' + str(scaling[i])
		if verbose == 0:
			cmd += ' 2> ' + bam + '.error.log'
		res = subprocess.run(cmd, capture_output=True, shell=True)
		if res.returncode != 0:
			raise ValueError('issue with the command: ' + str(res.stderr))


def ReadRoseSuperEnhancers(roseFolder, containsINPUT=True, name="MV411"):
	"""
	reads ROSE2's output and returns its superenhancer bedfile as a pd dataframe. 

	Can be multiple superenhancers from a set of HRK27ac bam files

	Args:
	-----
		roseFolder: str folderpath to ROSE's output
		containsINPUT: bool whether the bedfile contains INPUT signal as well
		name: str sample name from which we computed superenhancers

	Returns:
	--------
		a dataframe of the bed representation of the superenhancers
	"""
	beds = os.listdir(roseFolder)
	superenhan = pd.DataFrame()
	for i in beds:
		if i.endswith('_SuperEnhancers.table.txt'):
			superenhan = superenhan.append(pd.read_csv(roseFolder+i, sep='\t', skiprows=5)).drop(columns=['enhancerRank','isSuper'])
	data = [superenhan.columns[6]] + superenhan.columns[8:].tolist()
	if containsINPUT:
		inputd = superenhan.columns[7]
	superenhan['name'] = [i.split(name)[1].split('_')[2] for i in superenhan['REGION_ID']]
	superenhan = superenhan.rename(columns={'CONSTITUENT_SIZE':'size','REGION_ID':"id",'CHROM':'chrom','START':'start','STOP':'end',"NUM_LOCI":'num'}).replace(np.nan,0)
	superenhan['foldchange'] = superenhan[data].sum(1)/superenhan[inputd]
	superenhan = superenhan.drop(columns=data+[inputd])
	return superenhan.sort_values(by=['chrom','start','end'])


def loadPeaks(peakFile=None, peakfolder=None, isNF=True, CTFlist=[], skiprows=0):
	"""
	loads 1 to many peak bedfile into one pandas dataframe.

	all og the peaks will be concatenated into one dataframe. this function can 
	work with jkobject|nfcore/chipseq nextflow pipelines output (isNF)

	Args:
	-----
		peakFile: str filepath to a peak bedfile
		peakfolder: str folderpath to a folder congaining beddfile
		isNF: bool if the peaks come from an nfcore/chipseq nextflow pipeline (with VA_names)
		CTFlist: list[str] only load peaks corresponding from this list of names
		skiprows: int if there is rows to skip in the file (can be 1, often)

	Returns:
	-------
		a bed like dataframe containing the concatenated peaks of the bed files
	"""
	if peakfolder:
		bindings = pd.DataFrame()
		for folder in os.listdir(peakfolder):
			if isNF:
				if any(tf in folder for tf in CTFlist) or not CTFlist:
					binding = pd.read_csv(peakfolder + folder + "/NA_peaks.narrowPeak", sep='\t', header=None) if\
					os.exists(peakfolder + folder + "/NA_peaks.narrowPeak") else peakfolder + folder + "/NA_peaks.broadPeak"
					binding['name'] = folder.replace('.narrowPeak', '').replace('.broadPeak','')
					bindings = bindings.append(binding)
			else:
				file = folder
				if file[-10:] in ["narrowPeak",".broadPeak"] and (any(tf in file for tf in CTFlist) or not CTFlist):
					print('reading: '+file)
					binding = pd.read_csv(peakfolder + file, sep='\t', header=None, skiprows=skiprows)
					binding['name'] = file.replace('.narrowPeak', '').replace('.broadPeak','')
					bindings = bindings.append(binding)
	elif peakFile:
		bindings = pd.read_csv(peakFile, sep='\t', header=None, skiprows=skiprows)
		bindings['name'] = peakFile.split('/')[-1].split('.')[0]
	else:
		raise ValueError("need to provide one of peakFile or peakfolder")
	if not len(bindings.columns)==6:
		bindings = bindings.drop(5, 1).drop(4, 1).rename(columns={6:4, 7:5, 8:6, 9:7,})
	bindings = bindings.rename(columns={
		0: "chrom",
		1: 'start',
		2: 'end',
		3: 'peak_number',
		4: "foldchange",
		5: "-log10pvalue",
		6: "-log10qvalue",
		7: 'relative_summit_pos'})
	bindings = bindings.sort_values(by=["chrom", "start", "end"], axis=0)
	bindings.start = bindings.start.astype('int')
	bindings.end = bindings.end.astype('int')
	if not len(bindings.columns)==6:
		bindings['relative_summit_pos'] = bindings.relative_summit_pos.astype(
			float) if 'relative_summit_pos' in bindings.columns else bindings.end - bindings.start
		bindings["-log10pvalue"] = bindings["-log10pvalue"].astype('float')
		bindings['-log10qvalue'] = bindings['-log10qvalue'].astype('float')
		loc = bindings['relative_summit_pos'].isna()
		bindings.loc[bindings[loc].index, 'relative_summit_pos'] = bindings[loc].end - bindings[loc].start
		bindings.relative_summit_pos = bindings.relative_summit_pos.astype(int)
	bindings.foldchange = bindings.foldchange.astype('float')
	bindings = bindings.reset_index(drop=True)
	return bindings


def simpleMergePeaks(peaks, window=0, totpeaknumber=0, maxp=True, mergedFold="mean"):
	"""
	simply merges bedfiles from peak callers. providing a concaneted dataframe of bed-like tables

	will recompute pvalues and foldchange from that.

	Args:
	-----
		peaks: pd.df of concatenated bed-like dataframe with a name column differentiating each peak group
			to merge
		window: int the max size between each peaks under which to still merge the peaks together 
		totpeaknumber: int, set to more than 0 for peaks to be counted from this value 
			(e.g. 10, with 100 peaks will give peak_number of 10-110)
		maxp: bool, set to true to merge pvalue uusing the max pvalue, else will take the product of all
		mergedFold: flag "mean"|"max"|"sum", on how to merge foldchanges

	Returns:
	-------
		df bed-like of the merged peaks
	"""
	peaks = peaks.sort_values(by=['chrom', 'start','end'])
	tfs = list(set(peaks['name']))
	mergedpeaksdict = {}
	remove = []
	peaknumber = 0
	merged_bed = {
		"chrom": [peaks.iloc[0]['chrom']],
		"start": [peaks.iloc[0]['start']],
		"end": [],
		"peak_number": [peaknumber + totpeaknumber],
		"foldchange": [],
		"-log10pvalue": [],
		"-log10qvalue": [],
		"relative_summit_pos": []
	}
	foldchange = [peaks.iloc[0].get('foldchange', 1)]
	log10pvalue = [peaks.iloc[0].get('-log10pvalue', 0)]
	log10qvalue = [peaks.iloc[0].get('-log10qvalue', 0)]
	relative_summit_pos = peaks.iloc[1].get('relative_summit_pos', peaks.iloc[0]['start'])
	# computes overlap by extending a bit the window (100bp?) should be ~readsize
	prev_end = peaks.iloc[0]['end']
	prev_chrom = peaks.iloc[0]['chrom']
	tfmerged = {a: [0] for a in tfs}
	tfmerged[peaks.iloc[0]['name']][-1] = peaks.iloc[0].get('foldchange', 1)
	for i, (pos, peak) in enumerate(peaks.iloc[1:].iterrows()):
		print(str(i / len(peaks)), end="\r")
		if prev_end + window > peak['start'] and prev_chrom == peak['chrom']:
			# can be merged
			if peak.get('foldchange', 1) > max(foldchange):
				relative_summit_pos = peak.get('relative_summit_pos', peaks['start'])
			foldchange.append(peak.get('foldchange', 1))
			log10pvalue.append(peak.get('-log10pvalue', 0))
			log10qvalue.append(peak.get('-log10qvalue', 0))

		else:
			# newpeak
			for k, val in tfmerged.items():
				val.append(0)
			peaknumber += 1
			merged_bed['chrom'].append(peak['chrom'])
			merged_bed['start'].append(peak['start'])
			merged_bed['end'].append(prev_end)
			merged_bed['peak_number'].append(peaknumber + totpeaknumber)
			if mergedFold=="mean":
				merged_bed['foldchange'].append(np.mean(foldchange))
			elif mergedFold=="max":
				merged_bed['foldchange'].append(max(foldchange))
			elif mergedFold=="sum":
				merged_bed['foldchange'].append(sum(foldchange))
			else:
				raise ValueError("mergedFold needs to be one of:")
			merged_bed['-log10pvalue'].append(max(log10pvalue) if maxp else np.prod(log10pvalue))
			merged_bed['-log10qvalue'].append(max(log10qvalue) if maxp else np.prod(log10qvalue))
			merged_bed['relative_summit_pos'].append(relative_summit_pos)
			foldchange = [peak.get('foldchange', 1)]
			log10pvalue = [peak.get('-log10pvalue', 0)]
			log10qvalue = [peak.get('-log10qvalue', 0)]
			relative_summit_pos = peak.get('relative_summit_pos', peak['start'])
		prev_end = peak['end']
		prev_chrom = peak['chrom']
		tfmerged[peak['name']][-1] = peak.get('foldchange', 1)
	merged_bed['end'].append(prev_end)
	if mergedFold=="mean":
		merged_bed['foldchange'].append(np.mean(foldchange))
	elif mergedFold=="max":
		merged_bed['foldchange'].append(max(foldchange))
	elif mergedFold=="sum":
		merged_bed['foldchange'].append(sum(foldchange))
	else:
		raise ValueError("mergedFold needs to be one of:")
	merged_bed['-log10pvalue'].append(max(log10pvalue) if maxp else np.prod(log10pvalue))
	merged_bed['-log10qvalue'].append(max(log10qvalue) if maxp else np.prod(log10qvalue))
	merged_bed['relative_summit_pos'].append(relative_summit_pos)

	merged_bed = pd.DataFrame(merged_bed)
	tfmerged = pd.DataFrame(tfmerged)
	return pd.concat([merged_bed, tfmerged], axis=1, sort=False)

#from numba import jit, float32, int8

#@jit(float32[:](float32[:, 3], float32[:, 4], int8, str))
def putInBed(conscensus, value, window=10, mergetype='mean', columns=['foldchange']):
	""" 
	given a conscensus bed-like dataframe and another one, will merge the second one into the first

	Args:
	-----
		conscensus df[start,end,chrom] the conscensus (first one)
		value df[start, end, chrom]+columns the value one (second one)
		mergetype: flag: mean,first,last, on how to merge ttwo peaks that would fall on the same one
		window: int on max distance of second df peaks' to the first ones, to still merge them
			re.g. 0 is hard overap, 10000000 is nearest peak)

	Returns:
	-------
		np.array of values of the second dataframe to add to the first one 
		(e.g. can do df1[newcol] = returned_array)
	"""
	conscensus = conscensus.sort_values(by=['chrom','start','end']).reset_index(drop=True)
	value = value.sort_values(by=['chrom','start','end']).reset_index(drop=True)
	locinvalue=0
	loc=0
	tot=0
	num = []
	res = np.zeros((len(conscensus), len(columns)))
	not_end=True
	def add(res,num,not_end,loc):
		if len(num)>0:
			if mergetype=='mean':
				res[loc,:] = np.mean(num, 0)
			elif mergetype=='sum':
				res[loc,:]=np.sum(num, 0)
			elif mergetype=='first':
				res[loc,:]=num[0]
			elif mergetype=='last':
				res[loc,:]=num[-1]
			else:
				raise ValueError('must be one of')
			num= []
		loc+=1
		if loc == len(conscensus):
			not_end=False
		return res, num, not_end,loc
	while not_end:
		print(loc/len(conscensus),end="\r")
		try:
			a = conscensus.iloc[loc]
			b = value.iloc[locinvalue]
			# we need to switch chromosome for b
			if b.chrom < a.chrom:
				locinvalue+=1
				if locinvalue == len(value):
					not_end=False
			# we need to switch chromosome for b
			elif b.chrom > a.chrom:
				loc+=1
				if loc == len(conscensus):
					not_end=False
			elif b.start<a.start:
				if b.end+window>a.start:
					tot+=1
					num.append(b[columns])
					if b.end>a.end+window:
						res,num,not_end,loc = add(res,num,not_end,loc)
						continue
				locinvalue+=1
				if locinvalue == len(value):
					not_end=False
			elif b.start<a.end+window:
				tot+=1
				num.append(b[columns])
				if b.end>a.end+window:
					res,num,not_end,loc = add(res,num,not_end,loc)
					continue
				locinvalue+=1
				if locinvalue == len(value):
					not_end=False
			else:
				res,num,not_end,loc = add(res,num,not_end,loc)
		except:
			import pdb; pdb.set_trace()
	print(str(tot)+' were merged into conscensus')
	return res


async def pairwiseOverlap(bedfile, norm=True, bedcol=8, correct=True, docorrelation=True):
	"""
	compute pairwise overlap and correlation on this overlap for a set of peaks mappe to a conscensus 

	with each columns after the 7th one representing the signal of a given ChIP experiment
	over this conscensus. will present overlap of row values in col values

	Args:
	----
		bedfile: df bed-like representing a conscensus set of peaks, and a set of values/foldchanges
			over it
		norm: bool wether to zscore normalise the signal or not
		bedcol: int col where the bed information ends and signal information columns start
		correct: bool whether to correct for fully similar lines/ columns by removing 1 on their last value
		docorrelation: bool whether or not to compute correlation as well as enrichment

	Returns:
	-------
		a dataframe[values_name x values_name] of % overlap of row values in col values in the consensus
		a dataframe[values_name x values_name] of correlation of values signal over the overlaps
	"""
	if correct:
		print("we will be correcting for fully similar lines/ columns by removing 1 on their last value")
	dat = bedfile[bedfile.columns[bedcol:]].values
	correlation = np.ones((dat.shape[1],dat.shape[1]))
	overlap = np.ones((dat.shape[1],dat.shape[1]))
	for i, col in enumerate(dat.T):
		#pdb.set_trace()
		overlapping = np.delete(dat,i,axis=1)[col!=0]
		col = col[col!=0]
		add=0
		for j, val in enumerate(overlapping.T):
			if j==i:
				add=1
			if docorrelation:
				if norm and not np.isnan(zscore(col[val!=0]).sum()) and not np.isnan(zscore(val[val!=0]).sum()) or not correct:
					correlation[i,j+add] = np.corrcoef(zscore(val[val!=0]),zscore(col)[val!=0])[0,1]
				else:
					tmp = np.corrcoef(val[val != 0], col[val != 0])[0, 1]
					if np.isnan(tmp) and correct:
						if len(col[val!=0]) == 0 or len(val[val!=0]) == 0:
						# one has no overlap
							correlation[i,j+add] =  0
						else:
							# one contains only the same value everywhere
							col[-1]-=max(0.01,abs(np.mean(col)))
							val[-1]-=max(0.01,abs(np.mean(val)))
							correlation[i,j+add] = np.corrcoef(val,col)[0,1]
			overlap[i,j+add]=len(val[val!=0])/len(col)
	if docorrelation:
		correlation = pd.DataFrame(data=correlation, index=bedfile.columns[bedcol:], columns=bedfile.columns[bedcol:])
		correlation[correlation.isna()] = 0
	overlap = pd.DataFrame(data=overlap, index=bedfile.columns[bedcol:], columns=bedfile.columns[bedcol:]).T
	return overlap, correlation if docorrelation else None


async def enrichment(bedfile: pd.DataFrame, bedcol=8, groups=None, correct=True, okpval=10**-3):
	"""
	compute pairwise enrichment and correlation for a set of peaks mappe to a conscensus 

	with each columns after the 7th one representing the signal of a given ChIP experiment
	over this conscensus. will present enrichment of row values in col values

	Args:
	----
		bedfile: df bed-like representing a conscensus set of peaks, and a set of values/foldchanges
			over it
		bedcol: int col where the bed information ends and signal information columns start
		correct: bool whether to correct for multiple hypothesis testing or not
		docorrelation: bool whether or not to compute correlation as well as enrichment
		okpval: float max pvalue over which to set the enrichment to 0

	Returns:
	-------
		a dataframe[values_name x values_name] of enrichment of row values in col values
		a dataframe[values_name x values_name] of correlation of values signal over the overlaps
	"""
	dat = bedfile[bedfile.columns[bedcol:]].values
	prob = dat.astype(bool).sum(0)/len(dat)
	enrichment = np.zeros((dat.shape[1] if groups is None else len(set(groups)), dat.shape[1]))
	pvals = np.zeros(
		(dat.shape[1] if groups is None else len(set(groups)), dat.shape[1]))
	if groups is not None:
		for i in set(groups):
			overlapping = dat[groups==i]
			for j,val in enumerate(overlapping.T):
				# enrichment of j in i
				e, p = fisher_exact([
					[len(val[val != 0]), len(val[val == 0])],
					[prob[j]*len(dat), (1-prob[j])*len(dat)]])
				enrichment[i, j] = np.log2(e)
				pvals[i, j] = p
	else:
		for i, col in enumerate(dat.T):
			overlapping = np.delete(dat,i,axis=1)[col!=0]
			col = col[col!=0]
			add=0
			for j, val in enumerate(overlapping.T):
				if j==i:
					add=1
					enrichment[i,i]=0
				e, p = fisher_exact([[len(val[val != 0]), len(val[val == 0])], [
					prob[j+add]*len(dat), (1-prob[j+add])*len(dat)]])
				enrichment[i, j+add] = np.log2(e)
				pvals[i, j+add] = p
		enrichment[i,i]=0
	enrichment = pd.DataFrame(data=enrichment, index=bedfile.columns[bedcol:] if groups is None else set(groups), columns=bedfile.columns[bedcol:]).T
	enrichment[enrichment==-np.inf] = -1000
	enrichment[enrichment.isna()] = 0
	enrichment[enrichment == np.inf] = 1000
	if correct:
		pvals = np.reshape(multipletests(pvals.ravel(),
										0.1, method="bonferroni")[1], pvals.shape)
	pvals = pd.DataFrame(
		data=pvals, index=bedfile.columns[bedcol:]  if groups is None else set(groups), columns=bedfile.columns[bedcol:]).T
	enrichment[pvals>okpval] = 0
	return enrichment, pvals


def findAdditionalCobindingSignal(conscensus, known=None, bigwigs=[], window=100):
	"""
	somewhat similar concept to computePeaksAt

	# get pysam data
	# ask for counts only at specific locus based on peak windows from mergedpeakset
	# append to an array
	# return array, normalized
	"""
	raise NotImplementedError()
	if known:
		print('getting '+ str(len(peaks.tf))+' peaks. Using the peaks values directly if \
				available and using the bigwigs otherwise.')
		res = known.values.astype(float)
	elif len(bigwigs)>0:
		print('getting '+str(len(bigwigs))+' bigwigs, no peaks passed. Will compute the cobinding values\
			across the conscensus for each bigwigs.')
		res = np.zeros((len(bigwigs), len(conscensus)), dtype=float)
	else:
		raise ValueError('you need to pass a list of path to bigwigs for each/some samples')
	for i, bw in enumerate(bigwigs):
		if known:
			found = False
			for j, val in enumerate(known.tf):
				if val in bw:
					if found:
						raise ValueError('found two or more matching tf for bigwig: '+str(bw))
					found = True
					i=j
					break
				if not found:
					print('no tf found in known for tf: '+bw)
					raise ValueError('you need to have an amount of known columns equal to your bigwigs')
		print('doing file ' + str(bw))
		bw = pyBigWig.open(bw)
		for k, val in conscensus.iterrows():
			if known:
				if res[i][k]!=0:
					continue
			start = max([val.start - window, 0])
			end = min(val.end + window, bw.chroms(str(val.chrom)))
			res[i][k] = bw.stats(str(val.chrom), start, end)[0]
	res = np.nan_to_num(res, 0)
	return conscensus.join(pd.Dataframe(data=(res.T / res.max(1)).T, columns=bigwigs if not known else known.columns))



def annotatePeaks():
	"""
	get H3k27AC peaks and compute Rose
	for each TF peak
	assign super enhancer if within it and create a super enhancer TFgroup
	for each peaks
	apply similar to rose and merge TF together. (create new TFgroup)
	for each TF groups say
		if within super enhancer
		if within high h3k27ac marks
		its nearest gene (distance to it, cis/trans)


	FIND A WAY TO FILTER TO ONLY PLACES WITH H3K27AC marks??

	TAD points where most contacts on one side happened on this side.
	specific distance. zone of most concentration of contacts for a peak region.
	"""
	# if way == 'Closest':

	# elif way == 'ClosestExpressed':

	# elif way == 'ActivityByContact':

	# else:
	#     raise ValueError('needs to be oneof Closest ClosestExpressed ActivityByContact')
	raise NotImplementedError()


# os.system()


def getCoLocalization():
	"""
	for each annotations (super enhancer & TFgroups)
	for each TF, find highest peak/meanCov. if above thresh add to localization
	"""
	raise NotImplementedError()


def refineGroupsWithHiC():
	"""
	given HiC data, for each loops (should be less than X Mb distance. create a Xkb zone around both sides
	and find + merge each TF/TFgroup at this location)
	"""
	raise NotImplementedError()


async def fullDiffPeak(bam1, bam2, control1, size=None, control2=None, scaling=None, directory='diffData/',
				 res_directory="diffPeaks/", isTF=False, compute_size=True, pairedend=True):
	"""
	will use macs3 to call differential peak binding from two bam files and their control

	one can also provide some spike in scaling information

	Args:
	-----
	bam1: str bamfilepath
	bam2: str bamfilepath
	control1: str control bam filepath (INPUT/IGG)
	control2: str control bam filepath (INPUT/IGG)
	scaling: tuple(float1,float2) of scaling of each bams 
	size: int the read length (extsize parameter)
	directory: str the directory where to save the produced data
	res_directory: str the directory where to save the results
	isTF: bool true if TF else false (will change the size)
	compute_size: bool whether to compute the extsize value with `macs3 predictd`
	pairedend: bool if paired end

	Returns:
	-------
		str the log of the macs3 bdgdiff command
	"""
	print("doing diff from " + bam1 + " and " + bam2)
	if scaling is not None:
		if max(scaling) > 1:
			raise ValueError("scalings need to be between 0-1")
	name1 = bam1.split('/')[-1].split('.')[0]
	name2 = bam2.split('/')[-1].split('.')[0]
	if size is None:
		if isTF:
			size = 147
		else:
			size = 200
	if compute_size:
		print('computing the fragment avg size')
		cmd = "macs3 predictd -i " + bam1
		ret = subprocess.run(cmd, capture_output=True, shell=True)
		size = re.findall("# predicted fragment length is (\d+)", str(ret.stderr))[0]
		print(size)
	else:
		print('using default|given size')
	pairedend = "BAMPE" if pairedend else "BAM"
	if control2 is None:
		control2 = control1
	cmd1 = "macs3 callpeak -B -t " + bam1 + " -c " + control1 + " --nomodel --extsize " + str(size) + " -n " + name1 + " --outdir " + directory + " -f " + pairedend
	cmd2 = "macs3 callpeak -B -t " + bam2 + " -c " + control2 + " --nomodel --extsize " + str(size) + " -n " + name2 + " --outdir " + directory + " -f " + pairedend
	print('computing the scaling values')
	ret = subprocess.run(cmd1, capture_output=True, shell=True)
	print(ret.stderr)
	scaling1a = int(re.findall(" after filtering in treatment: (\d+)", str(ret.stderr))[0])
	scaling1b = int(re.findall(" after filtering in control: (\d+)", str(ret.stderr))[0])
	scaling1 = scaling1a if scaling1a <= scaling1b else scaling1b
	ret = subprocess.run(cmd2, capture_output=True, shell=True)
	print(ret.stderr)
	scaling2a = int(re.findall(" after filtering in treatment: (\d+)", str(ret.stderr))[0])
	scaling2b = int(re.findall(" after filtering in control: (\d+)", str(ret.stderr))[0])
	scaling2 = scaling2a if scaling2a <= scaling2b else scaling2b
	if scaling is not None:
		scaling1 = int(scaling1/scaling[0])
		scaling2 = int(scaling2/scaling[1])
	print(scaling1, scaling2)
	return (await diffPeak(directory+name1+"_treat_pileup.bdg", directory+name2+"_treat_pileup.bdg",
		directory+name1+"_control_lambda.bdg", directory+name2+"_control_lambda.bdg",
		res_directory, scaling1, scaling2, size))


async def diffPeak(name1, name2, control1, control2, res_directory, scaling1, scaling2, size):
	"""
	calls MACS2 bdgdiff given some parameters

	Args:
	-----
		name1: str bamfilepath
		name2: str bamfilepath
		control1: str control bam filepath (INPUT/IGG)
		control2: str control bam filepath (INPUT/IGG)
		res_directory: str the directory where to save the results
		scaling1: float the scaling factor of the first bam
		scaling2: float the scaling factor of the second bam
		size: int the read length (extsize parameter)

	Returns:
	-------
		str the log of the macs3 bdgdiff command
	"""
	print("doing differential peak binding")
	cmd = "macs3 bdgdiff --t1 " + name1 + " --c1 "
	cmd += control1+" --t2 " + name2 +" --c2 " + control2
	cmd += " --d1 " + str(scaling1) + " --d2 " + str(scaling2) + " -g 60 "
	cmd += "-l " + str(size) + " --o-prefix " + name1.split('/')[-1].split('.')[0] + "_vs_"
	cmd += name2.split('/')[-1].split('.')[0] + " --outdir " + res_directory
	res = subprocess.run(cmd, capture_output=True, shell=True)
	return res


def AssignToClosestExpressed(bed,countFile,genelocFile):
	print("the bed file and genelocFile should use the same assembly")
	genelocFile = pd.read_csv(genelocFile,sep="\t", compression="", columns=['chrom','start','end',"name"])
	#for val in bed.iterrows():
	raise NotImplementedError()



async def makeSuperEnhancers(MACS2bed, bamFile, outdir, baiFile=None, rosePath=".",
	stitching_distance=None, TSS_EXCLUSION_ZONE_SIZE="2500", assembly="hg38",controlBam=None,controlBai=None):
	"""
	Calls super enhancer from H3K27ac with the ROSE algorithm

	Args:
	----
		MACS2bed: str he bed filepath to a MACS2 called set of peaks for an H3K27ac bam file
		bamFile: str bam filepath of the H3K27ac bam file
		outdir: str folderpath of the output. has to be an absolute path or a ~/path
		baiFile: str bam filepath of the H3K27ac bam file
		rosePath: str path to the ROSE folder containing the rose algorithm
		stitching_distance: int max stitching distance
		TSS_EXCLUSION_ZONE_SIZE: int size around TSS to excllude from SEs
		assembly: flag one of hg38, hg37...
		controlBam: str filepath to input bam file
		controlBai: str filepath to input bai file

	Returns:
	--------
		a bed-like dataframe of the superenhancers

	"""
	print("we are going to move your input files to "+rosePath)
	cmd = 'cp '+MACS2bed+' '+rosePath+MACS2bed.split('/')[-1]+'.bed'
	cmd += ' && mv '+bamFile+' '+rosePath
	baiFile = baiFile if baiFile else bamFile[:-1]+'i'
	cmd += ' && cp '+baiFile +' '+rosePath
	res = subprocess.run(cmd, capture_output=True, shell=True)
	if res.returncode != 0:
		raise SystemError('failed moving files: '+str(res))
	cmd = "cd "+rosePath+" && python ROSE_main.py -g "+assembly+" -i "+MACS2bed.split('/')[-1]+'.bed' + " -r " + bamFile.split('/')[-1] + " -o " + outdir
	if TSS_EXCLUSION_ZONE_SIZE:
		cmd+=" -t "+TSS_EXCLUSION_ZONE_SIZE
	if controlBam:
		os.system('mv '+controlBam+' '+rosePath)
		controlBai = controlBai if controlBai else controlBam[:-1]+'i'
		os.system('cp '+controlBai+' '+rosePath)
		cmd+=" -c "+controlBam.split('/')[-1]
	if stitching_distance:
		cmd+=" -s "+ stitching_distance

	res = subprocess.run(cmd, capture_output=True, shell=True)
	fail = False
	if res.returncode != 0:
		v = 'ROSE failed:' +str(res)
		fail = True
	print('finished.. moving them back to their original folder')
	cmd = 'rm '+rosePath+MACS2bed.split('/')[-1]+'.bed'
	cmd += ' && mv '+rosePath+bamFile.split('/')[-1]+' '+bamFile
	cmd+= ' && rm '+rosePath+baiFile.split('/')[-1]
	if controlBam:
		cmd += ' && mv '+rosePath+controlBam.split('/')[-1]+' '+controlBam
		cmd += ' && rm '+rosePath+controlBai.split('/')[-1]
	res = subprocess.run(cmd, capture_output=True, shell=True)
	if res.returncode != 0:
		raise SystemError('failed moving files: '+str(res))
	if fail:
		raise SystemError(v)
	print('worked')
	print(res)
	return ReadRoseSuperEnhancers(outdir ,bool(controlBam))



def runChromHMM(outdir, data, numstates=15, datatype='bed', folderPath=".", 
chromHMMFolderpath="~/ChromHMM/", assembly="hg38",control_bam_dir=None):
	"""
	runs the chromHMM algorithm

	Args:
	-----
		outdir str: an existing dir where the results should be saved
		data: df[cellname,markname,markbed|bam|bigwig, ?controlbed|bam|bigwig]
		numstates: number of states to use
		datatype: flag one of bed
		folderPath: str folder where to save chromHMM's work
		chromHMMFolderpath: str folderpath to chromHMM algorithm
		assembly: flag one of hg38, hg37 ...
		control_bam_dir: str directory where the control would be stored (if not given in the ddf)

	Returns:
	-------
		A dict of bed like dataframes containing the regions of the different states
	"""
	print("you need to have ChromHMM")
	chromHMM = "java -mx8000M -jar "+chromHMMFolderpath+"ChromHMM.jar "
	h.createFoldersFor(outdir+'binarized/')
	data.to_csv(outdir+"input_data.tsv", sep='\t',index=None,header=None)
	cmd = chromHMM
	if datatype=="bed":
		cmd+="BinarizeBed "
	elif datatype=="bigwig":
		cmd+="BinarizeSignal "
	elif datatype=="bam":
		cmd+="BinarizeBam "
	else:
		raise ValueError('you need to provide one of bam, bigwig, bed')
	cmd+= chromHMMFolderpath+"CHROMSIZES/"+assembly+".txt "+ folderPath+" "+outdir+"input_data.tsv "+outdir+"binarized"
	if control_bam_dir:
		cmd+=" -c "+control_bam_dir
	res1 = subprocess.run(cmd, capture_output=True, shell=True)
	print(res1)
	if res1.returncode!=0:
		raise ValueError(str(res1.stderr))
	cmd = chromHMM + "LearnModel -printposterior -noautoopen "
	if len(data)<10:
		cmd += '-init load -m '+chromHMMFolderpath+'model_15_coreMarks.txt '
	cmd += outdir+"binarized "+outdir+" "+str(numstates)+" "+assembly
	res2 = subprocess.run(cmd, capture_output=True, shell=True)
	print(res2)
	if res2.returncode!=0:
		raise ValueError(res2.stderr)
	ret = {}
	for v in set(data[0]):
		ret[v] = pd.read_csv(outdir+v+'_'+str(numstates)+'_dense.bed', sep='\t', header=None,
			skiprows=1).drop(columns=[4,5,6,7]).rename(columns=
			{0:'chrom',1:'start',2:'end',3:'state',8:"color"})
	return ret


def loadMEMEmotifs(file, tfsubset=[],motifspecies='HUMAN'):
	"""
	loads motif from the output file of MEME after running fimo.

	Args:
	----
		file: str location of the fimo gff|bed file
		tfsubset: list[str] subset of tf to load
		motifspecies: str name of the specied concerned (should be the same one as provided to MEME)

	Returns:
	--------
		df bed-like of each motifs across the genome
	"""
	if file.endswith('.gff'):
		print('converting to bed, you need to have "gfftobed" installed')
		cmd = 'gff2bed < '+file+' > '+file+'.bed'
		file = file+'.bed'
		res = subprocess.run(cmd,capture_output=True, shell=True)
		if res.returncode != 0:
			raise ValueError('issue with the command: ' + str(res.stderr))
		else:
			print(res.stdout.decode("utf-8"))
	## What are the motifs of our CRC members in ATACseq but not in our matrix
	merged_motif = pd.read_csv(file, sep='\t',skiprows=0,index_col=None, names=['pos',"fimo",
		"nucleotide_motif","relStart","relEnd","pval","strand",".","data"])
	merged_motif['tf']=[i[5:].split("_"+motifspecies)[0] for i in merged_motif.data]
	if tfsubset:
		merged_motif = merged_motif[merged_motif.tf.isin(tfsubset)]
	merged_motif['chrom'] = [i.split(':')[0][3:] for i in merged_motif.index]
	merged_motif['start'] = [i.split(':')[1].split('-')[0] for i in merged_motif.index]
	merged_motif['end'] = [i.split(':')[1].split('-')[1] for i in merged_motif.index]
	merged_motif = merged_motif.reset_index(drop=True)
	merged_motif = merged_motif.rename(columns={'pos':'relStart','fimo':'relEnd','nucleotide_motif':'pos',
		'relStart':'pval','relEnd':'strand','pval':'fimo','strand':'motif'})
	merged_motif['motif'] = [i.split('sequence=')[1].split(';')[0] for i in merged_motif.data]
	merged_motif['p_val'] = [i.split('pvalue=')[1].split(';')[0] for i in merged_motif.data]
	merged_motif['q_val'] = [i.split('qvalue=')[1].split(';')[0] for i in merged_motif.data]
	merged_motif = merged_motif.drop(columns=['pos','.','fimo','data'])
	merged_motif = merged_motif[merged_motif.columns[[6,7,8,0,1,3,2,9,10,5,4]]]
	merged_motif = merged_motif.sort_values(by=['chrom','start','end']).reset_index(drop=True)
	return merged_motif


def simpleMergeMotifs(motifs, window=0):
	"""
	aggregates the motifs if they overlap, into one motif file

	Args:
	----
		motifs: df bed-like of motifs locations
		window: int maxsize around motif for which to still merge

	Returns:
	--------
		df bedlike of merged motif:
		df bedlike of motifs that were not merged as they were different but still overlapping
	"""
	if type(motifs) is list:
		motifs = pd.concat(motifs)
	motifs = motifs.sort_values(by=['chrom', 'start'])
	toremove = []
	issues = []
	prevmotif = motifs.iloc[0]
	for i, (pos, motif) in enumerate(motifs.iloc[1:].iterrows()):
		print(str(i / len(motifs)), end="\r")
		if prevmotif['end'] + window > motif['start'] and prevmotif['chrom'] == motif['chrom']:
			# can be merged
			if motif['tf']!= prevmotif['tf'] or motif['motif'] != prevmotif['motif']:
				print('found different motifs overlapping')
				issues.extend([motif,prevmotif])
			else:
				toremove.append(pos)
		prevmotif = motif
	motifs = motifs.drop(index=toremove).reset_index(drop=True)
	issues = pd.concat(issues)
	return motifs, issues


def substractPeaksTo(peaks,loci, bp=50):
	"""
	removes all peaks that are not within a bp distance to a set of loci

	Args:
	----
		peaks: a bed file df with a chrom,start, end column at least
		loci: a df witth a chrom & loci column
		bp: the max allowed distance to the loci

	Returns:
	-------
		all the peaks that are within this distance
	"""
	i=0
	j=0
	keep=[]
	bp=50
	while j<len(peaks) and i<len(loci):
		h.showcount(j,len(peaks))
		if peaks.loc[j].chrom > loci.loc[i].chrom:
			i+=1
			continue
		if peaks.loc[j].chrom < loci.loc[i].chrom:
			j+=1
			continue
		if peaks.loc[j].start - bp > loci.loc[i].loci:
			i+=1
			continue
		if peaks.loc[j].end + bp< loci.loc[i].loci:
			j+=1
			continue
		if peaks.loc[j].end + bp >= loci.loc[i].loci and peaks.loc[j].start - bp <= loci.loc[i].loci:
			keep.append(j)
			j+=1
	return peaks.loc[set(keep)]

	
