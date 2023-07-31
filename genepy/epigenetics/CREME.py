import os
import pandas as pd
import numpy as np
from genepy.utils import helper as h
from genepy.utils import plot
from genepy.epigenetics.chipseq import *
import seaborn as sns
import pyBigWig
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.special import factorial
import warnings
import itertools


def findpeakpath(folder, proteiname):
    """
    given a folder of bigwigs and a protein name, finds the right bigwig
    """
    res = None
    for val in os.listdir(folder):
        if str(proteiname) in val:
            if res:
                raise ValueError('more than 1 bigwig file found')
            res = val
    if res:
        return res
    raise ValueError('no bigwig file found')


def findBestPeak(presence):
    """
    given a list of -sets of peak locations for each replicate- will return the best replicate given a simple metric
    """
    tot = []
    for ind, el in enumerate(presence):
        val = len(el)
        pres = [x for j, x in enumerate(presence) if j != ind]
        for jnd in range(1, len(pres)+1):
            for comb in itertools.combinations(pres, jnd):
                ov = el
                for knd in range(jnd):
                    ov = ov & comb[knd]
                val += len(ov)*(jnd+1)
        tot.append(val)
    return np.argsort(tot)[::-1]


def mergeReplicatePeaks(peaks, bigwigfolder, markedasbad=None, window=100,
                                                sampling=3000, mincov=4, doPlot=True, cov={}, minKL=8, use='max',
                                                MINOVERLAP=0.3, lookeverywhere=True, only='', saveloc=''):
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
        a the value for presence of each TF listed in previous df
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
    def col_nan_scatter(x, y, **kwargs):
        df = pd.DataFrame({'x': x[:], 'y': y[:]})
        df = df[df.sum(0) != 0]
        x = df['x']
        y = df['y']
        plt.gca()
        plt.scatter(x, y)

    def col_nan_kde_histo(x, **kwargs):
        df = pd.DataFrame({'x': x[:]})
        df = df[df['x'] != 0]
        x = df['x']
        plt.gca()
        sns.kdeplot(x)
    print("/!/ should only be passed peaks with at least one good replicate")
    # for a df containing a set of peaks in bed format and an additional column of different TF
    tfs = list(set(peaks['tf']))
    totpeaknumber = 0
    mergedpeaksdict = {}
    remove = []
    tomergebam = []
    ratiosofunique = {}
    h.createFoldersFor(saveloc)
    f = open(saveloc+'results.txt', 'w')
    warnings.simplefilter("ignore")
    for tf in tfs:
        if only and tf != only:
            continue
        cpeaks = peaks[peaks.tf == tf]
        print('_____________________________________________________')
        f.write('_____________________________________________________' + '\n')
        if len(set(cpeaks['replicate'])) == 1:
            if cpeaks.name.tolist()[0] in markedasbad:
                print('the only replicate is considered bad!')
                f.write('the only replicate is considered bad!'+"\n")
                print('wrong TF: '+tf)
                f.write('wrong TF: '+tf+"\n")
                mergedpeaksdict.update({tf: cpeaks})
                remove.append(tf)
                continue
            print("we only have one replicate for " + tf + " .. pass")
            f.write("we only have one replicate for " + tf + " .. pass"+"\n")
            mergedpeaksdict.update({tf: cpeaks})
            continue
        print("merging " + tf + " peaks")
        f.write("merging " + tf + " peaks"+"\n")
        merged = simpleMergePeaks(cpeaks, window=window, maxp=False)
        merged_bed = merged[merged.columns[8:]]
        finalpeaks = merged[merged.columns[:8]]
        print('--> finish first overlaps lookup')
        f.write('--> finish first overlaps lookup'+"\n")
        # flag when  biggest is <1000 peaks
        if len(finalpeaks) < 1000:
            print('!TF has less than 1000 PEAKS!')
            f.write('!TF has less than 1000 PEAKS!'+"\n")
        # for each TF (replicates), compute number of peaks
        peakmatrix = merged_bed.values.astype(bool)

        presence = []
        for peakpres in peakmatrix.T:  # https://github.com/tctianchi/pyvenn
            presence.append(set([i for i, val in enumerate(peakpres) if val == 1]))
        # compute overlap matrix (venn?)
        if peakmatrix.shape[1] < 7 and doPlot:
            plot.venn(presence, [i+'_BAD' if i.split('-')[0]
                                                in markedasbad else i for i in merged_bed.columns], title=tf+"_before_venn", folder=saveloc)
            plt.show()
        else:
            print('too many replicates for Venn: '+str(peakmatrix.shape[1]))
            f.write('too many replicates for Venn: '+str(peakmatrix.shape[1])+"\n")
        if doPlot:
            fig = sns.pairplot(merged_bed, corner=True, diag_kind="kde",
                                                 kind="reg", plot_kws={"scatter_kws": {"alpha": .05}})
            #fig = fig.map_upper(col_nan_scatter)
            #fig = fig.map_upper(col_nan_kde_histo)
            plt.suptitle("correlation of peaks in each replicate", y=1.08)
            if saveloc:
                fig.savefig(saveloc+tf+"_before_pairplot.pdf")
            plt.show()
            for i, val in enumerate(merged_bed):
                unique_inval = np.logical_and(
                    np.delete(peakmatrix, i, axis=1).sum(1).astype(bool) == 0, peakmatrix[:, i])
                sns.kdeplot(merged_bed[val][unique_inval], legend=True).set(xlim=(0, None))
            plt.title("distribution of unique peaks in each replicate")
            if saveloc:
                plt.savefig(saveloc+tf+"_before_unique_kdeplot.pdf")
            plt.show()

        bigwigs = os.listdir(bigwigfolder)

        foundgood = False
        sort = findBestPeak(presence)
        for ib, sb in enumerate(sort):
            if merged_bed.columns[sb].split('-')[0] not in markedasbad:
                foundgood = True
                break
        if not foundgood:
            print('no peaks were good enough quality')
            f.write('no peaks were good enough quality'+"\n")
            print('bad TF: '+tf)
            f.write('bad TF: '+tf+"\n")
            remove.append(tf)
            ib = 0
        # distplot
        # correlation plot

        biggest_ind = sort[ib]
        peakmatrix = peakmatrix.T
        biggest = merged_bed.columns[biggest_ind]
        print('-> main rep is: '+str(biggest))
        f.write('-> main rep is: '+str(biggest)+'\n')
        tot = peakmatrix[biggest_ind].copy().astype(int)
        # starts with highest similarity and go descending
        j = 0
        recovered = 0
        additionalpeaksinbig = np.array([])
        for i, val in enumerate(sort):
            if i == ib:
                continue
            j += 1
            # if avg non overlap > 60%, and first, and none small flag TF as unreliable.
            overlap = len(presence[biggest_ind] & presence[val]
                                        ) / len(presence[biggest_ind])
            peakname = merged_bed.columns[val]
            print('- '+peakname)
            f.write('- '+peakname+'\n')
            print('  overlap: ' + str(overlap*100)+"%")
            f.write('  overlap: ' + str(overlap*100)+"%"+'\n')
            if overlap < MINOVERLAP:
                smallsupport = len(presence[biggest_ind] &
                                                     presence[val]) / len(presence[val])
                print(' --> not enough overlap')
                f.write(' --> not enough overlap'+'\n')
                if smallsupport < MINOVERLAP:
                    # if the secondary does not have itself the required support
                    if j == 1 and merged_bed.columns[val].split('-')[0] not in markedasbad:
                        print("  Wrong TF: "+tf)
                        f.write("  Wrong TF: "+tf+'\n')
                        remove.append(tf)
                        break
                    # if not first, throw the other replicate and continue
                    print("  not using this replicate from the peakmatrix")
                    f.write("  not using this replicate from the peakmatrix"+'\n')
                    continue
            if lookeverywhere:
                tolookfor = peakmatrix[val] == 0
            else:
                tolookfor = np.logical_and(peakmatrix[biggest_ind], peakmatrix[val] == 0)
            # ones that we have in the Primary but not in the secondary
            additionalpeaksinsec = findAdditionalPeaks(finalpeaks, tolookfor, bigwigfolder + findpeakpath(
                bigwigfolder, peakname), sampling=sampling, mincov=mincov, window=window, minKL=minKL, use=use)
            if len(additionalpeaksinsec[additionalpeaksinsec > 0]) > 0:
                sns.kdeplot(additionalpeaksinsec[additionalpeaksinsec > 0],
                                        label=peakname, legend=True).set(xlim=(0, None))
                print('  min,max from newly found peaks: ' +
                            str((additionalpeaksinsec[additionalpeaksinsec > 0].min(), additionalpeaksinsec[additionalpeaksinsec > 0].max())))
                f.write('  min,max from newly found peaks: '+str((additionalpeaksinsec[additionalpeaksinsec > 0].min(
                ), additionalpeaksinsec[additionalpeaksinsec > 0].max()))+'\n')
            # for testing purposes mainly
            finalpeaks[additionalpeaksinsec.astype(bool)].to_csv(
                'additionalpeaksinsec_mp'+merged_bed.columns[val]+'.bed', sep='\t', index=None, header=False)
            peakmatrix[val] = np.logical_or(
                peakmatrix[val], additionalpeaksinsec.astype(bool))
            overlap = np.sum(np.logical_and(
                peakmatrix[val], peakmatrix[biggest_ind]))/np.sum(peakmatrix[biggest_ind])
            if overlap < MINOVERLAP:
                newsmalloverlap = np.sum(np.logical_and(
                    peakmatrix[val], peakmatrix[biggest_ind]))/np.sum(peakmatrix[val])
                print("  we did not had enough initial overlap.")
                f.write("  we did not had enough initial overlap."+'\n')
                if newsmalloverlap < MINOVERLAP:
                    if merged_bed.columns[val].split('-')[0] in markedasbad:
                        print('  replicate ' +
                                    merged_bed.columns[val] + ' was too bad and had not enough overlap')
                        f.write('  replicate ' +
                                        merged_bed.columns[val] + ' was too bad and had not enough overlap'+'\n')
                        continue
                    elif h.askif("we have two good quality peaks that don't merge well at all: "+merged_bed.columns[val] +
                                                                                        " and " + merged_bed.columns[biggest_ind] + " can the first one be removed?:\n  \
                            overlap: "+str(overlap*100)+'%\n  smalloverlap: '+str(smalloverlap*100)+'%\n  new smalloverlap: '+str(newsmalloverlap*100)+"%"):
                        continue
                    else:
                        print("  enough from small overlaps")
                        f.write("  enough from small overlaps"+'\n')
            print(' --> enough overlap')
            f.write(' --> enough overlap'+'\n')
            recovered += np.sum(additionalpeaksinsec.astype(bool))
            if merged_bed.columns[val].split('-')[0] not in markedasbad:
                tot += peakmatrix[val].astype(int)
            # ones that we have in the Primary but not in the secondary
            if not lookeverywhere or len(additionalpeaksinbig) == 0:
                tolookfor = peakmatrix[biggest_ind] == 0 if lookeverywhere else np.logical_and(
                    peakmatrix[biggest_ind] == 0, peakmatrix[val])
                additionalpeaksinbig = findAdditionalPeaks(finalpeaks, tolookfor, bigwigfolder + findpeakpath(
                    bigwigfolder, biggest), sampling=sampling, mincov=mincov, window=window, minKL=minKL, use=use)
                if len(additionalpeaksinbig[additionalpeaksinbig > 0]) > 0:
                    sns.kdeplot(additionalpeaksinbig[additionalpeaksinbig > 0],
                                            label=biggest, legend=True).set(xlim=(0, None))
                    print('  min,max from newly found peaks: ' +
                                str((additionalpeaksinbig[additionalpeaksinbig > 0].min(), additionalpeaksinbig[additionalpeaksinbig > 0].max())))
                    f.write('  min,max from newly found peaks: '+str((additionalpeaksinbig[additionalpeaksinbig > 0].min(
                    ), additionalpeaksinbig[additionalpeaksinbig > 0].max()))+'\n')

                peakmatrix[biggest_ind] = np.logical_or(
                    peakmatrix[biggest_ind], additionalpeaksinbig)
                tot += additionalpeaksinbig.astype(bool).astype(int)
                recovered += np.sum(additionalpeaksinbig.astype(bool))
            print('  we have recovered ' + str(recovered)+' peaks, equal to ' + str(100*recovered/np.sum(peakmatrix[biggest_ind])) +
                                                        '% of the peaks in main replicate')
            f.write('  we have recovered ' + str(recovered)+' peaks, equal to ' + str(100*recovered/np.sum(peakmatrix[biggest_ind])) +
                                                        '% of the peaks in main replicate'+'\n')
            if overlap < (MINOVERLAP+0.2)/1.2:
                # we recompute to see if the overlap changed
                newoverlap = np.sum(np.logical_and(
                    peakmatrix[val], peakmatrix[biggest_ind]))/np.sum(peakmatrix[biggest_ind])
                smalloverlap = np.sum(np.logical_and(
                    peakmatrix[val], peakmatrix[biggest_ind]))/np.sum(peakmatrix[val])
                if newoverlap < (MINOVERLAP+0.2)/1.2:
                    if smalloverlap < (2+MINOVERLAP)/3:
                        print("  not enough overlap to advice to merge the bams.\n  oldnew overlap: "+str(overlap*100)+'%\n  \
                            new overlap: '+str(newoverlap*100)+"%")
                        f.write("  not enough overlap to advice to merge the bams.\n  oldnew overlap: "+str(overlap*100)+'%\n  \
                            new overlap: '+str(newoverlap*100)+"%"+'\n')
                        continue
                    else:
                        print('  enough from small overlap to advice to merge the peaks')
                        f.write('  enough from small overlap to advice to merge the peaks'+'\n')
            tomergebam.append([biggest, peakname])
            #the quality is good enough in the end we can pop from the list if it exists
            if tf in remove:
                remove.remove(tf)
        plt.title('distribution of new found peaks')
        if saveloc:
            plt.savefig(saveloc+tf+"_new_found_peaks_kdeplot.pdf")
        plt.show()
        # new distplot
        # new correlation plot
        ratiosofunique[tf] = len(np.argwhere(
            peakmatrix.sum(0) == 1))/peakmatrix.shape[1]
        if doPlot:
            sns.pairplot(merged_bed, corner=True, diag_kind="kde",
                                     kind="reg", plot_kws={"scatter_kws": {"alpha": .05}})
            #fig = fig.map_upper(col_nan_scatter)
            #fig = fig.map_upper(col_nan_kde_histo)
            plt.suptitle("correlation and distribution of peaks after recovery", y=1.08)
            if saveloc:
                plt.savefig(saveloc+tf+"_after_pairplot.pdf")
            plt.show()
            for i, val in enumerate(merged_bed):
                unique_inval = np.logical_and(
                    np.delete(peakmatrix, i, axis=0).sum(0).astype(bool) == 0, peakmatrix[i])
                sns.kdeplot(merged_bed[val][unique_inval], legend=True).set(xlim=(0, None))
            plt.title("distribution of unique peaks in each replicate after recovery")
            if saveloc:
                plt.savefig(saveloc+tf+"_after_unique_kdeplot.pdf")
            plt.show()
        if len(peakmatrix.shape) > 1 and doPlot:
            if peakmatrix.shape[0] < 7:
                presence = []
                for peakpres in peakmatrix:  # https://github.com/tctianchi/pyvenn
                    presence.append(set([i for i, val in enumerate(peakpres) if val == 1]))
                title = tf + '_recovered (TOREMOVE)' if tf in remove else tf+'_recovered'
                plot.venn(presence, [i+'_BAD' if i.split('-')[0]
                                                                        in markedasbad else i for i in merged_bed.columns], title=title, folder=saveloc)
                plt.show()
            else:
                print('too many replicates for Venn')
                f.write('(too many replicates for Venn)'+'\n')
            finalpeaks = finalpeaks[np.logical_or(tot > 1, peakmatrix[biggest_ind])]
        finalpeaks['name'] = biggest
        finalpeaks['tf'] = tf
        mergedpeaksdict.update({tf: finalpeaks})
        print(str((tf, len(finalpeaks))))
        f.write(str((tf, len(finalpeaks)))+'\n')
    mergedpeak = pd.concat(
        [peaks for _, peaks in mergedpeaksdict.items()]).reset_index(drop=True)
    if doPlot:
        df = pd.DataFrame(data=ratiosofunique, index=['percentage of unique'])
        df['proteins'] = df.index
        fig = sns.barplot(data=df)
        plt.xticks(rotation=60, ha='right')
        plt.title("ratios of unique in replicates across experiments")
        if saveloc:
            plt.savefig(saveloc+"All_ratios_unique.pdf")
        plt.show()
    f.close()
    mergedpeak['name'] = mergedpeak.tf
    return mergedpeak, tomergebam, remove, ratiosofunique


def findAdditionalPeaks(peaks, tolookfor, filepath, sampling=1000, mincov=4,
                                                window=100, cov={}, minKL=8, use='max'):

    """
    findAdditionalPeaks: for all peaks in A and/or B find in coverage file if zone has relative cov
    of more than thresh then add to peak
    if B is small and > 20% of peaks in A are found back, increase confidence and
    flag for mergeBams
    if < 20% don't flag for merge bam
    f B is big and now mean non overlap < 40%, take union and flag for mergeBam else, throw B.

    Args:
    -----
        peaks
        tolookfor
        filepath
        sampling
        mincov
        window
        cov
        minKL
        use
    returns:
    -------
        np.array(bool) for each peaks in peakset, returns a binary
    """
    # def poisson(k, lamb, scale): return scale * (lamb**k / factorial(k)) * np.exp(-lamb)

    def KLpoisson(lamb1, lamb2): return lamb1 * \
                        np.log(lamb1 / lamb2) + lamb2 - lamb1

    def poisson(k, lamb): return (lamb**k/factorial(k)) * np.exp(-lamb)

    def negLogLikelihood(params, data): return - \
                        np.sum(np.log(poisson(data, params[0])))

    def poissonFit(data): return float(
        minimize(negLogLikelihood, x0=np.ones(1), args=(data,), method='Powell').x)
    bw = pyBigWig.open(filepath)
    res = np.zeros(len(peaks))
    prevchrom = ''
    lamb = {}
    cov = {}
    #ignore by message
    warnings.filterwarnings("ignore", message="encountered in")
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
                        #TODO: compute on INPUT file instead
                        samples = np.zeros(window * sampling)
                        sam = np.random.rand(sampling)
                        sam = sam * (bw.chroms(str(val.chrom))-window)
                        for j, sample in enumerate(sam.astype(int)):
                            samples[j*window:(j + 1)*window] = np.nan_to_num(
                                bw.values(str(val.chrom), sample, sample + window), 0)
                        scale = np.unique(samples)[1]
                        samples = (samples/scale).astype(int)
                        lamb[val.chrom] = (poissonFit(samples), scale)

            start = max([val.start - window, 0])
            end = min(val.end + window, bw.chroms(str(val.chrom)))
            zone = np.nan_to_num(bw.values(str(val.chrom), start, end), 0)
            if use == 'max':
                if max(zone) / cov[val.chrom] > mincov*1.5 or sum(zone) / (cov[val.chrom] * (end - start)) > mincov:
                    res[i] = max(zone) / cov[val.chrom]
            elif use == 'poisson':
                #TODO: compute -log10pvalue
                la = poissonFit((zone/lamb[val.chrom][1]).astype(int))
                kl = KLpoisson(la, lamb[val.chrom][0])
                if kl > minKL:
                    res[i] = max(zone) / cov[val.chrom]  # foldchange from macs3

    return res
