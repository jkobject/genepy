import os
import pandas as pd
import pysam
import sys
import numpy as np
from JKBio import Helper as h
import re
from pybedtools import BedTool
import seaborn as sns
import pyBigWig
import matplotlib.pyplot as plt
import ipdb
import signal
from scipy.optimize import curve_fit,minimize
from sklearn.preprocessing import normalize
from scipy.stats import poisson
from scipy.special import factorial
import subprocess
from pandas.io.parsers import ParserError, EmptyDataError



size = {"GRCh37": 2864785220,
        "GRCh38": 2913022398}

cmaps = ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
         'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
         'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']

chroms = {'chr1','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
'chr2','chr20','chr21','chr22','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chrX',
 'chrY','1','10','11','12','13','14','15','16','17','18','19','2','20','21','22','3','4','5','6',
 '7','8','9','X','Y'}

def dropWeirdChromosomes(bedfile, keep=[], skip=0):
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

def extractPairedSingleEndFrom(folder, pattern='R1/R2', sep='-', namepos=2):
    """
    given a folder, find using a specific patterns, single end and paired end files
    return a list of single end files and a list of list of paired end files [[R1,R2]]
    """
    single = []
    paired = {}
    for val in os.listdir(folder):
        if 'R1' in val:
            name = val.split(sep)[namepos]
            paired[name] = {'R1': val}
        elif 'R2' in val:
            name = val.split(sep)[namepos]
            paired[name].update({'R2': val})
        else:
            single.append(val)

    return single, pd.DataFrame(paired)


def findReplicates(folder, sep='-', namings='-r([0-9])', namepos=2):
    """
    creates a dict of name and replicate files
    """
    rep = {}
    for val in os.listdir(folder):
        if val[-4:] == '.bam':
            match = re.search(namings, val)
            if match:
                number = match.groups()[0]
                name = val.split(sep)[namepos]
                if name in rep:
                    rep[name].append(val)
                else:
                    rep[name] = [val]

    return rep


def computeSingleEnd(singlend, folder="data/seqs/", numthreads=8, peaksFolder="peaks/",
                     ismapped=False, mappedFolder='mapped/', refFolder='data/reference/index'):
    """
    run the singleEnd pipeline
    for alignment etc, one can use pysam ready made implementation of samtools
    """
    for val in singlend:
        out1 = folder + mappedFolder + val.split('.')[0] + ".mapped.sam"
        if not ismapped:
            in1 = folder + val
            os.system("bowtie2 -x " + refFolder + " --threads " + str(numthreads) +
                      " -t -k 1 --very-sensitive -U " + in1 + " -S " + out1)
        out2 = folder + peaksFolder + val.split('.')[0]
        print(out1)
        os.system("macs2 callpeak -f SAM -t " + out1 + " --outdir " + out2)
        # it can take many TB so better delete


def computePairedEnd(pairedend, folder="", numthreads=8, peaksFolder="peaks/",
                     ismapped=False, mappedFolder='mapped/', refFolder='data/reference/index'):
    """
    # run the paired end pipeline
    """
    for key, val in pairedend.items():
        out1 = folder + mappedFolder + val[0].split('.')[0] + ".mapped.sam"
        in1 = folder + val[0]
        in2 = folder + val[1]
        os.system("bowtie2 -x " + refFolder + " --threads " + str(numthreads) + " -t -k 1 \
      --very-sensitive -1 " + in1 + " -2 " + in2 + " - S " + out1)
        out2 = folder + peaksFolder + val[0].split('.')[0]
        print(out1)
        changefrom = out1
        changeto = out1[:-4] + '.bam'
        os.system("samtools view -b " + changefrom + " -o " + changeto)
        os.system("macs2 callpeak --format 'BAMPE' --treatment " + changeto + " --outdir " + out2)
        # it can take many TB so better delete


def bigWigFrom(bams, folder="", numthreads=8, genome='GRCh37', scaling=None, verbose=1):
    """
    run the bigwig command line for a set of bam files in a folder
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


def mergeBams(rep):
    """
    uses samtools to merge a set of replicates considered into one file
    """
    for i, val in rep.items():
        out1 = i + '.merged.bam'
        for bam in val:
            in1 += ' ' + bam
        os.system("samtools merge " + out1 + in1)


def loadPeaks(peakfolder="data/peaks/", isMacs=True, CTFlist=[], skiprows=0):
    # for each data peak type file listed in a MACS2 way in a given MACS2 output folder, will merge them
    # all into one dataframe and output the dataframe
    bindings = pd.DataFrame()
    for folder in os.listdir(peakfolder):
        if isMacs:
            if any(tf in folder for tf in CTFlist) or not CTFlist:
                binding = pd.read_csv(peakfolder + folder + "/NA_peaks.narrowPeak", sep='\t', header=None) if\
                os.exists(peakfolder + folder + "/NA_peaks.narrowPeak") else peakfolder + folder + "/NA_peaks.broadPeak"
                binding['name'] = folder.replace('.narrowPeak', '').replace('.broadPeak','')
                bindings = bindings.append(binding)
        else:
            file = folder
            if file[-10:] in ["narrowPeak",".broadPeak"] and (any(tf in file for tf in CTFlist) or not CTFlist):
                binding = pd.read_csv(peakfolder + file, sep='\t', header=None, skiprows=skiprows)
                binding['name'] = file.replace('.narrowPeak', '').replace('.broadPeak','')
                bindings = bindings.append(binding)
    bindings = bindings.drop(5, 1).drop(4, 1)
    bindings = bindings.rename(columns={
        0: "chrom",
        1: 'start',
        2: 'end',
        3: 'peak_number',
        6: "foldchange",
        7: "-log10pvalue",
        8: "-log10qvalue",
        9: 'relative_summit_pos'})
    bindings = bindings.sort_values(by=["chrom", "start", "end"], axis=0)
    bindings.start = bindings.start.astype('int')
    bindings.end = bindings.end.astype('int')
    bindings['relative_summit_pos'] = bindings.relative_summit_pos.astype('int') if 'relative_summit_pos' in bindings.columns else bindings.end - bindings.start
    bindings.foldchange = bindings.foldchange.astype('float')
    bindings["-log10pvalue"] = bindings["-log10pvalue"].astype('float')
    bindings['-log10qvalue'] = bindings['-log10qvalue'].astype('float')
    return bindings.reset_index(drop=True)


def pysam_computePeaksAt(peaks, bams, folder='data/seqs/', window=1000, numpeaks=1000, numthreads=8):

    # get pysam data
    # ask for counts only at specific locus based on windows from center+-size from sorted MYC peaks
    # for each counts, do a rolling average (or a convolving of the data) with numpy
    # append to an array
    # return array, normalized
    loaded = {}
    res = {i: np.zeros((len(peaks), window * 2)) for i in bams}
    peaks = peaks.sort_values(by="foldchange", ascending=False).iloc[:numpeaks]
    peaks.chrom = peaks.chrom.astype(str)
    for val in bams:
        loaded.update({val: pysam.AlignmentFile(folder + val, 'rb', threads=numthreads)})
    for k, bam in loaded.items():
        for num, (i, val) in enumerate(peaks.iterrows()):
            print(int(num / len(peaks)), end='\r')
            center = int((val['start'] + val['end']) / 2)
            for pileupcolumn in bam.pileup(val['chrom'], start=center - window,
                                           stop=center + window, truncate=True):
                res[k][num][pileupcolumn.pos - (center - window)] = pileupcolumn.n
    fig, ax = plt.subplots(1, len(res))
    for i, (k, val) in enumerate(res.items()):
        sns.heatmap(val, ax=ax[i])
        ax[i].set_title(k.split('.')[0])
    fig.show()
    return res, fig


def bedtools_computePeaksAt(peaks, bams, folder='data/seqs/', window=1000, numpeaks=1000, numthreads=8):
    """
    get pysam data
    ask for counts only at specific locus based on windows from center+-size from sorted MYC peaks
    for each counts, do a rolling average (or a convolving of the data) with numpy
    append to an array
    return array, normalized
    """
    loaded = {}
    center = [int((val['start'] + val['end']) / 2) for k, val in peaks.iterrows()]
    peaks['start'] = [c - window for c in center]
    peaks['end'] = [c + window - 1 for c in center]
    peaks[peaks.columns[:3]].sort_values(by=['chrom', 'start']).to_csv('temp/peaks.bed', sep='\t', index=False, header=False)
    bedpeaks = BedTool('temp/peaks.bed')

    fig, ax = plt.subplots(1, len(bams))
    peakset = peaks["foldchange"].values.argsort()[::-1][:numpeaks]
    for i, val in enumerate(bams):
        coverage = BedTool(folder + val).intersect(bedpeaks).genome_coverage(bga=True, split=True)\
            .intersect(bedpeaks).to_dataframe(names=['chrom', 'start', 'end', 'coverage'])
        cov = np.zeros((len(peaks), window * 2), dtype=int)
        j = 0
        pdb.set_trace()
        for i, (k, val) in enumerate(peaks.iterrows()):
            print(i / len(peaks), end='\r')
            while coverage.iloc[j].start > val.start:
                j -= 1
            while coverage.iloc[j].start < val.end:
                cov[i][coverage.iloc[j].start - val.start:coverage.iloc[j].end - val.start] =\
                    coverage.iloc[j].coverage
                j += 1
        sns.heatmap(coverage, ax=ax[i])
        ax[i].set_title(val.split('.')[0])
    fig.show()
    return None, fig




def computePeaksAt(peaks, bigwigs, folder='', bigwignames=[], peaknames=[], window=1000, title='', numpeaks=4000, numthreads=8,
                   width=5, length=10,torecompute=False, name='temp/peaksat.png', refpoint="TSS", scale=None, 
                   sort=False, withDeeptools=True, onlyProfile=False, cluster=1):
    """
    get pysam data
    ask for counts only at specific locus based on windows from center+-size from sorted MYC peaks
    for each counts, do a rolling average (or a convolving of the data) with numpy
    append to an array
    return array, normalized
    """
    if withDeeptools:
        if isinstance(peaks, pd.DataFrame):
            peaks = 'peaks.bed '
            peaks.to_csv('peaks.bed', sep='\t', index=False, header=False)
        elif type(peaks) == list:
            pe = ''
            i=0
            for n, p in enumerate(peaks):
                if 20 < int(os.popen('wc -l ' + p).read().split(' ')[0]):
                    pe += p + ' '
                elif len(peaknames) > 0:
                    peaknames.pop(n-i)
                    i+=1
            peaks = pe
        elif type(peaks) == str:
            peaks += ' '
        else:
            raise ValueError(' we dont know this filetype')
        if type(bigwigs) is list:
            pe = ''
            for val in bigwigs:
                pe += folder + val + ' '
            bigwigs = pe
        else:
            bigwigs = folder + bigwigs + ' '
        h.createFoldersFor(name)
        cmd= ''
        if not os.path.exists(name) or torecompute:
            cmd += "computeMatrix reference-point -S "
            cmd += bigwigs
            cmd += " --referencePoint "+refpoint
            cmd += " --regionsFileName " + peaks
            cmd += " --missingDataAsZero"
            cmd += " --outFileName " + '.'.join(name.split('.')[:-1]) + ".gz"
            cmd += " --upstream " + str(window) + " --downstream " + str(window)
            cmd += " --numberOfProcessors " + str(numthreads) + ' && '
        cmd += "plotHeatmap " if not onlyProfile else 'plotProfile '
        cmd += "--matrixFile " + '.'.join(name.split('.')[:-1]) + ".gz"
        cmd += " --outFileName " + name
        cmd += " --refPointLabel "+ refpoint
        if cluster>1:
            cmd += " --perGroup --kmeans "+str(cluster)

        if len(peaknames) > 0:
            pe = ''
            for i in peaknames:
                pe += i + ' '
            peaknames = pe
            cmd += " --regionsLabel " + peaknames
        if len(bigwignames) > 0:
            pe = ''
            for i in bigwignames:
                pe += i + ' '
            bigwignames = pe
            cmd += " --samplesLabel " + bigwignames
        if title:
            cmd += " --plotTitle " + title
        data = subprocess.run(cmd, shell=True, capture_output=True)
        print(data)
    else:
        if 'relative_summit_pos' in peaks.columns:
            center = [int((val['start'] + val['relative_summit_pos'])) for k, val in peaks.iterrows()]
        else:
            center = [int((val['start'] + val['end']) / 2) for k, val in peaks.iterrows()]
        pd.set_option('mode.chained_assignment', None)
        peaks['start'] = [c - window for c in center]
        peaks['end'] = [c + window for c in center]
        fig, ax = plt.subplots(1, len(bigwigs), figsize=[width, length], title=title if title else 'Chip Heatmap')
        if sort:
            peaks = peaks.sort_values(by=["foldchange"], ascending=False)
        if numpeaks > len(peaks):
            numpeaks = len(peaks) - 1
        cov = {}
        maxs = []
        for num, bigwig in enumerate(bigwigs):
            bw = pyBigWig.open(folder + bigwig)
            co = np.zeros((numpeaks, window * 2), dtype=int)
            scale = scale[bigwig] if scale is dict else 1
            for i, (k, val) in enumerate(peaks.iloc[:numpeaks].iterrows()):
                try:
                    co[i] = np.nan_to_num(bw.values(str(val.chrom), val.start, val.end), 0)
                except RuntimeError as e:
                    print(str(val.chrom), val.start, val.end)
                    pass
            cov[bigwig] = co
            maxs.append(co.max())
        for num, bigwig in enumerate(bigwigs):
            sns.heatmap(cov[bigwig] * scale, ax=ax[num], vmax=max(maxs), yticklabels=[], cmap=cmaps[num],
                        cbar=True)
            ax[num].set_title(bigwig.split('.')[0])
        fig.subplots_adjust(wspace=0.1)
        fig.show()
        fig.savefig(name)
        return cov, fig


"""
def computeMeanCov(bigwigFolder, meanOnly=True, ref="GRCh38", averageFragSize=150, outputfolder='/data/coverage/'):
    meancov = {}
    for val in os.listdir(bigwigFolder):
        bw = pyBigWig.open(folder + bigwig)
        if bw:
            meancov[val.split('.')[0]] = bw.
    return meancov


def substractPeaks(peaks1, to):
    peaks1 = peaks1.sort_values(by=['chrom', 'start'])
    peaks2 = to.sort_values(by=['chrom', 'start'])
    j = 0
    i = 0
    newpeaks = pd.DataFrame(columns=peaks2.columns)
    while i < len(peaks2):
        if peaks1.iloc[j].chrom == peaks2.iloc[j].chrom:
            else
        intersection = helper.intersection([peaks2.iloc[i]['start'], peaks2.iloc[i]['end']],
                                           [peaks1.iloc[j]['start'], peaks1.iloc[j]['start']])
        if intersection:

            j += 1
        while peaks2.iloc[i]['end'] > peaks1.iloc[j]['start']:
            if peaks2.iloc[i]['end'] > peaks1.iloc[j]['end']:
                newpeaks.append(peaks2.iloc[i], ignore_index=True)
            else:
                newpeaks.append({'start': peaks2.iloc[i]['start'],
                                 'end': peaks1.iloc[j]['start'],
                                 'chromomsome': peaks1.iloc[j]['chrom']}, ignore_index=True)
"""


def simpleMergedPeaks(peaks, window=0, totpeaknumber=0, maxp=True):
    peaks = peaks.sort_values(by=['chrom', 'start'])
    tfs = list(set(peaks['name'].tolist()))
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
            foldchange.append(peak.get('foldchange', 1))
            log10pvalue.append(peak.get('-log10pvalue', 0))
            log10qvalue.append(peak.get('-log10qvalue', 0))
            if peak.get('foldchange', 1) > max(foldchange):
                relative_summit_pos = peak.get('relative_summit_pos', peaks['start'])
        else:
            # newpeak
            for k, val in tfmerged.items():
                val.append(0)
            peaknumber += 1
            merged_bed['chrom'].append(peak['chrom'])
            merged_bed['start'].append(peak['start'])
            merged_bed['end'].append(prev_end)
            merged_bed['peak_number'].append(peaknumber + totpeaknumber)
            merged_bed['foldchange'].append(np.mean(foldchange))
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
    merged_bed['foldchange'].append(np.mean(foldchange))
    merged_bed['-log10pvalue'].append(max(log10pvalue) if maxp else np.prod(log10pvalue))
    merged_bed['-log10qvalue'].append(max(log10qvalue) if maxp else np.prod(log10qvalue))
    merged_bed['relative_summit_pos'].append(relative_summit_pos)

    merged_bed = pd.DataFrame(merged_bed)
    tfmerged = pd.DataFrame(tfmerged)
    return pd.concat([merged_bed, tfmerged], axis=1, sort=False)


def findpeakpath(folder, peakname):
    res = None
    for val in os.listdir(folder):
        if str(peakname) in val:
            if res:
                raise ValueError('more than 1 bigwig file found')
            res= val
    if res:
        return res
    raise ValueError('no bigwig file found')


def mergeReplicatePeaks(peaks, bigwigfolder, markedasbad=None, window=100,
                        sampling=3000, mincov=4, doPlot=True, cov={}, minKL=8, use='max',
                        MINOVERLAP=0.3, mergewindow=1, lookeverywhere=True,only=''):
    """
    /!/ should only be passed peaks with at least one good replicate
    for each TFpeaksets,
    1. find the replicate that have the most peaks
    2. correlate peaks and get in highest correlation order with the replicate found in 1
    3. find overlap of both and get size of second replicate
    4. if small(er)-> use only to increase statistics
      1. if a lot of uncalled peaks in replicate 2 at replicate 1 peaks (flag for mergebam)
    5. if similar size -> get only intersect
      2. add to intersect, find uncalled peaks in both replicates which are called in the other
    6. repeat for all replicates
    -------------------------
    if full overlap of one of the peak replicate, only use the overlapped one to increase confidence on peak
    if >80% average non overlap,
      print warning and percentage of overlap

    if <20% average non overlap,
      take the overlap and increase confidence and avg logfold

    if one is <20%:
      if other <40% average non overlap,
        take the overlap and increase confidence and avg logfold
      else
        take

    gets the max cov at the genomic window and if above some threshold, accepts the peak.

    extend peak by X bp if no TSS
    remove TSS from peaks

    create a new data frame containing merged peak size, reassembled peak data (p value etc..) and
    a boolean value for presence of each TF listed in previous df
    ------------------------------------

    args:
    ----
    peaks: df[bed-like] all the peaks into the sameBam with a column containing the 'name'
    being the id of the sample, the 'replicate' number of this sample, the 'tf' chiped here
    bamfolder: str, foldername
    avgCov: dict(filename:int) a dict where for each bam filename is given an averageCoverage
    if use=='max':
        window:
        mincov:

    if use=='max':  


    returns:
    -------
    mergedpeaks: dict{df-peakslike}
    bamtomerge: [[bam1,bam2]]

    """
    print("/!/ should only be passed peaks with at least one good replicate")
    # for a df containing a set of peaks in bed format and an additional column of different TF
    tfs = list(set(peaks['tf']))
    totpeaknumber = 0
    mergedpeaksdict = {}
    remove = []
    tomergebam = []
    for tf in tfs:
        if only and tf!=only:
            continue
        cpeaks = peaks[peaks.tf==tf]
        if len(set(cpeaks['replicate'])) == 1:
            if cpeaks.name.tolist()[0] in markedasbad:
                print('the only replicate is considered bad!')
                print('wrong TF: '+tf)
                mergedpeaksdict.update({tf: cpeaks})
                remove.append(tf)
                continue
            print("we only have one replicate for " + tf + " .. pass")
            mergedpeaksdict.update({tf: cpeaks})
            continue
        print("merging " + tf + " peaks")
        merged = simpleMergedPeaks(cpeaks, window=mergewindow, maxp=False)
        merged_bed = merged[merged.columns[8:]]
        finalpeaks = merged[merged.columns[:8]]
        print('finish first overlaps lookup')
        # flag when  biggest is <1000 peaks
        if len(finalpeaks) < 1000:
            print('!TF has less than 1000 PEAKS!')
        # for each TF (replicates), compute number of peaks
        peakmatrix = merged_bed.values

        # compute overlap matrix (venn?)
        if peakmatrix.shape[1] < 7 and doPlot:
            presence = []
            for peakpres in peakmatrix.T:  # https://github.com/tctianchi/pyvenn
                presence.append(set([i for i, val in enumerate(peakpres) if val == 1]))
            h.venn(presence, [i+'_BAD' if i.split('-')[0] in markedasbad else i for i in merged_bed.columns], title=tf)
        else:
            print('too many replicates for Venn: '+str(peakmatrix.shape[1]))
        bigwigs = os.listdir(bigwigfolder)
        totpeak = np.sum(peakmatrix, 0)
        sort = np.argsort(totpeak)[::-1]
        print("found total peak for this replicate set: "+str(totpeak))
        foundgood=False
        for ib,sb in enumerate(sort):
            if merged_bed.columns[sb].split('-')[0] not in markedasbad:
                foundgood=True
                break
        if not foundgood:
            print('no peaks were good enough quality')
            print('wrong TF: '+tf)
            remove.append(tf)
        biggest_ind = sort[ib]
        peakmatrix = peakmatrix.T
        biggest = merged_bed.columns[biggest_ind]
        print('main rep is: '+str(biggest))
        tot = peakmatrix[biggest_ind].copy()
        # starts with highest similarity and go descending
        j = 0
        recovered = 0
        additionalpeaksinbig = np.array([])
        for i, val in enumerate(sort):
            if i==ib:
                continue
            j+=1
            # if avg non overlap > 60%, and first, and none small flag TF as unreliable.
            overlap = len(presence[biggest_ind] & presence[val]) / len(presence[biggest_ind])
            peakname = merged_bed.columns[val]
            print(peakname)
            print('overlap: ' + str(overlap*100)+"%")
            if overlap < MINOVERLAP:
                smallsupport = len(presence[biggest_ind] & presence[val]) / len(presence[val])
                print('not enough overlap')
                if smallsupport < MINOVERLAP:
                    # if the secondary does not have itself the required support
                    if j == 1 and merged_bed.columns[val].split('-')[0] not in markedasbad:
                        print("Wrong TF: "+tf)
                        remove.append(tf)
                        break
                    # if not first, throw the other replicate and continue
                    print("not using this replicate from the peakmatrix")
                    continue
            if lookeverywhere:
                tolookfor = peakmatrix[val] == 0
            else:
                tolookfor = np.logical_and(peakmatrix[biggest_ind], peakmatrix[val] == 0)
            # ones that we have in the Primary but not in the secondary
            additionalpeaksinsec = findAdditionalPeaks(finalpeaks, tolookfor, bigwigfolder + findpeakpath(bigwigfolder, peakname), sampling=sampling, mincov=mincov, window=window, minKL=minKL, use=use)
            # for testing purposes mainly
            finalpeaks[additionalpeaksinsec].to_csv('additionalpeaksinsec_mp'+merged_bed.columns[val]+'.bed',sep='\t',index=None,header=False)
            peakmatrix[val] = np.logical_or(peakmatrix[val], additionalpeaksinsec)
            overlap = np.sum(np.logical_and(peakmatrix[val],peakmatrix[biggest_ind]))/np.sum(peakmatrix[biggest_ind])
            if overlap < MINOVERLAP:
                newsmalloverlap = np.sum(np.logical_and(peakmatrix[val],peakmatrix[biggest_ind]))/np.sum(peakmatrix[val])
                print("we did not had enough initial overlap.")
                if newsmalloverlap < MINOVERLAP:
                    if merged_bed.columns[val].split('-')[0] in markedasbad:
                        print('replicate ' + merged_bed.columns[val] + ' was too bad and had not enough overlap')
                        continue
                    elif h.askif("we have two good quality peaks that don't merge well at all: "+merged_bed.columns[val] +\
                            " and " +merged_bed.columns[biggest_ind]+ " can the first one be removed?:\n\
                            overlap: "+str(overlap*100)+'%\nsmalloverlap: '+str(smalloverlap*100)+'%\nnew smalloverlap: '+str(newsmalloverlap*100)+"%"):
                        continue
                    else:
                        print("enough from small overlaps")
            print('enough overlap')
            recovered += np.sum(additionalpeaksinsec)
            if merged_bed.columns[val].split('-')[0] not in markedasbad:
                tot += peakmatrix[val].astype(int) 
            # ones that we have in the Primary but not in the secondary
            if not lookeverywhere or len(additionalpeaksinbig)==0:
                tolookfor = peakmatrix[biggest_ind] == 0 if lookeverywhere else np.logical_and(peakmatrix[biggest_ind]==0, peakmatrix[val])
                additionalpeaksinbig = findAdditionalPeaks(finalpeaks, tolookfor, bigwigfolder + findpeakpath(bigwigfolder, biggest), sampling=sampling, mincov=mincov, window=window, minKL=minKL, use=use)
                peakmatrix[biggest_ind] = np.logical_or(peakmatrix[biggest_ind], additionalpeaksinbig)        
                tot +=additionalpeaksinbig.astype(int)
                recovered += np.sum(additionalpeaksinbig)
            print('we have recovered ' + str(recovered)+' peaks, equal to '+ str(100*recovered/np.sum(peakmatrix[biggest_ind]))+\
                '% of the peaks in main replicate')
            if overlap < (MINOVERLAP+0.2)/1.2:
                # we recompute to see if the overlap changed
                newoverlap = np.sum(np.logical_and(peakmatrix[val],peakmatrix[biggest_ind]))/np.sum(peakmatrix[biggest_ind])
                smalloverlap = np.sum(np.logical_and(peakmatrix[val],peakmatrix[biggest_ind]))/np.sum(peakmatrix[val])
                if newoverlap < (MINOVERLAP+0.2)/1.2:
                    if smalloverlap < (2+MINOVERLAP)/3:
                        print("not enough overlap to advice to merge the bams.\noldnew overlap: "+str(overlap*100)+'%\n\
                            new overlap: '+str(newoverlap*100)+"%")
                        continue
                    else:
                        print('enough from small overlap to advice to merge the peaks')
            tomergebam.append([biggest, peakname])
        if len(peakmatrix.shape) > 1 and doPlot and tf not in remove:
            if peakmatrix.shape[0] < 7:
                presence = []
                for peakpres in peakmatrix:  # https://github.com/tctianchi/pyvenn
                    presence.append(set([i for i, val in enumerate(peakpres) if val == 1]))
                h.venn(presence, [i+'_BAD' if i in markedasbad.split('-')[0] else i for i in merged_bed.columns],title=tf+'_recovered')
            else:
                print('too many replicates for Venn')    
        if tf not in remove:
            finalpeaks = finalpeaks[peakmatrix[biggest_ind].astype(bool)] if len(peakmatrix.shape) < 3 else np.logical_or(tot>1,peakmatrix[biggest_ind])
            finalpeaks['name'] = biggest
            finalpeaks['tf'] = tf
            mergedpeaksdict.update({tf: finalpeaks})
            print(tf,len(finalpeaks),set(finalpeaks.tf))
    mergedpeak = pd.concat([peaks for _, peaks in mergedpeaksdict.items()]).reset_index(drop=True)
    return mergedpeak, tomergebam, remove


def findAdditionalPeaks(peaks, tolookfor, filepath, sampling=1000, mincov=4,
                        window=100, cov={}, minKL=8, use='max'):
    """
    findAdditionalPeaks: for all peaks in A and/or B find in coverage file if zone has relative cov
    of more than thresh then add to peak
    if B is small and > 20% of peaks in A are found back, increase confidence and
    flag for mergeBams
    if < 20% don't flag for merge bam
    f B is big and now mean non overlap < 40%, take union and flag for mergeBam else, throw B.

    returns:
    -------
    peaks: np.array(bool) for each peaks in peakset, returns a binary
    """
    # def poisson(k, lamb, scale): return scale * (lamb**k / factorial(k)) * np.exp(-lamb)

    def KLpoisson(lamb1, lamb2): return lamb1 * np.log(lamb1 / lamb2) + lamb2 - lamb1
    def poisson(k, lamb): return (lamb**k/factorial(k)) * np.exp(-lamb)
    def negLogLikelihood(params, data): return - np.sum(np.log(poisson(data, params[0])))
    def poissonFit(data): return float(minimize(negLogLikelihood,x0=np.ones(1),args=(data,),method='Powell').x)
    bw = pyBigWig.open(filepath)
    res = np.zeros(len(peaks))
    prevchrom = ''
    lamb = {}
    cov = {}
    for i, has in enumerate(tolookfor):
        if has:
            val = peaks.iloc[i]
            if val.chrom not in chroms:
                continue
            if val.chrom != prevchrom:
                if val.chrom not in cov:
                    cov[val.chrom] = bw.stats(str(val.chrom))[0]
                    prevchrom = val.chrom
                    if use == 'poisson':
                        samples = np.zeros(window * sampling)
                        sam = np.random.rand(sampling)
                        sam = sam * (bw.chroms(str(val.chrom))-window)
                        for j, sample in enumerate(sam.astype(int)):
                            samples[j*window:(j + 1)*window] = np.nan_to_num(bw.values(str(val.chrom), sample, sample + window), 0)
                        scale = np.unique(samples)[1]
                        samples = (samples/scale).astype(int)
                        lamb[val.chrom] = (poissonFit(samples),scale)

            start = max([val.start - window, 0])
            end = min(val.end + window, bw.chroms(str(val.chrom)))
            zone = np.nan_to_num(bw.values(str(val.chrom), start, end), 0)
            if use == 'max':
                if max(zone) / cov[val.chrom] > mincov*1.5 or sum(zone) / (cov[val.chrom] * (end - start)) > mincov:
                    res[i] = 1
            elif use == 'poisson':
                #TODO: compute foldchange and -log10pvalue
                zone = (zone/lamb[val.chrom][1]).astype(int)
                la = poissonFit(zone)
                kl=KLpoisson(la, lamb[val.chrom][0])
                if kl > minKL:
                    res[i] = 1

    return res.astype(bool)


def findAdditionalCobindingSignal(conscensus, known=None, bigwigs=[], window=100):
    """
    somewhat similar concept to computePeaksAt

    # get pysam data
    # ask for counts only at specific locus based on peak windows from mergedpeakset
    # append to an array
    # return array, normalized
    """
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


# os.system()


def getCoLocalization():
    """
    for each annotations (super enhancer & TFgroups)
    for each TF, find highest peak/meanCov. if above thresh add to localization
    """


def refineGroupsWithHiC():
    """
    given HiC data, for each loops (should be less than X Mb distance. create a Xkb zone around both sides
    and find + merge each TF/TFgroup at this location)
    """


def getPeaksOverlap(peaks, isMerged=False, correlationMatrix=None, countMatrix=None, extend=1000,
                    onlyOn=[], quality={}):
    """
    generates a venn diagram of overlap between a replicate (quality anotation can be added)

    args:
    ----
      peakset: df of Bed Peak format
      onlyOn: list[string] of TF name
      quality: dict[string: {1,2,3}] a dict with a value between 1: good, 2: bad, 3: throw for each TF
    """


# def assignGene(peaks, bedFolder):


def fullDiffPeak(bam1, bam2, control1, control2=None, scaling=None, directory='diffData/',
                 res_directory="diffPeaks/", isTF=False, compute_size=True, pairedend=True):
    """
    will use macs2 to call differential peak binding

    Args:
    -----
    bam1
    bam2()
    control1
    control2
    scaling
    """
    print("doing diff from " + bam1 + " and " + bam2)
    name1 = bam1.split('/')[-1].split('.')[0]
    name2 = bam2.split('/')[-1].split('.')[0]
    if isTF:
        size = 147
    else:
        size = 200
    if compute_size:
        print('computing the fragment avg size')
        cmd = "macs2 predictd -i " + bam1
        ret = subprocess.run(cmd, capture_output=True, shell=True)
        size = re.findall("# predicted fragment length is (\d+)", str(ret.stderr))[0]
        print(size)
    else:
        print('using default size')
    pairedend = "BAMPE" if pairedend else "BAM"
    if control2 is None:
        control2 = control1
    cmd1 = "macs2 callpeak -B -t " + bam1 + " -c " + control1 + " --nomodel --extsize " + size + " -n " + name1 + " --outdir " + directory + " -f " + pairedend
    cmd2 = "macs2 callpeak -B -t " + bam2 + " -c " + control2 + " --nomodel --extsize " + size + " -n " + name2 + " --outdir " + directory + " -f " + pairedend
    if scaling is None:
        print('computing the scaling values')
        ret = subprocess.run(cmd1, capture_output=True, shell=True)
        print(ret.stderr)
        scaling1a = re.findall("tags after filtering in treatment: (\d+)", str(ret.stderr))[0]
        scaling1b = re.findall("tags after filtering in control: (\d+)", str(ret.stderr))[0]
        scaling1 = scaling1a if scaling1a <= scaling1b else scaling1b
        ret = subprocess.run(cmd2, capture_output=True, shell=True)
        print(ret.stderr)
        scaling2a = re.findall("tags after filtering in treatment: (\d+)", str(ret.stderr))[0]
        scaling2b = re.findall("tags after filtering in control: (\d+)", str(ret.stderr))[0]
        scaling2 = scaling2a if scaling2a <= scaling2b else scaling2b
    else:
        res = subprocess.run(cmd1, capture_output=True, shell=True)
        scaling1 = scaling[0]
        res = subprocess.run(cmd2, capture_output=True, shell=True)
        scaling2 = scaling[1]
    print(scaling1, scaling2)
    diffPeak(name1, name2, res_directory, directory, scaling1, scaling2, size)


def diffPeak(name1, name2, res_directory, directory, scaling1, scaling2, size):
    print("doing differential peak binding")
    cmd = "macs2 bdgdiff --t1 " + directory + name1 + "_treat_pileup.bdg --c1 "
    cmd += directory + name1 + "_control_lambda.bdg --t2 " + directory + name2
    cmd += "_treat_pileup.bdg --c2 " + directory + name2 + "_control_lambda.bdg "
    cmd += "--d1 " + str(scaling1) + " --d2 " + str(scaling2) + " -g 60 "
    cmd += "-l " + str(size) + " --o-prefix " + name1 + "_vs_" + name2 + " --outdir " + res_directory
    res = subprocess.run(cmd, capture_output=True, shell=True)
    return res