# Utils

1. helper functions to save data, generate random strings, run tasks in parallel etc.
2. a set of plotting tools based on [matplotlib]() and [bokeh]() 

## Contains:

_in ./helper.py_

- fileToList: convert a txt with a list of values to a python list
- listToFile: converts a list of values to a txt
- dictToFile: converts a dict to a json
- fileToDict: converts a json to a dict
- batchMove: move a lot of file in batch (can provide different locations)
- batchRename: rename a bunch of files in batch
- createFoldersFor: makes the required folders for a given filepath
- grouped: to use in a forloop to group values in batch
- overlap: given two tuples, returns the overlap
- union: given two tuples, returns the union
- nans: gets nans from a panda df
- randomString: generate a random string for namings
- parrun: runs list of commands in parallel
- askif: ask the user a questions and returns the y/n answer
- inttodate: converts an int to a string date.
- datetoint: converts a date to an int
- showcount: pretty print of i/size%, to put in a for loop
- combin: outputs the number of comabination of n object taken k at a time
- dups: shows the duplicates in a list
- makeCombinations: produces probability of X event happening at the same time. wil compute it given binomial probabilities of each event occuring and the number of trials
- closest: returns the index of the value closest to K in a lst
- compareDfs: compares df1 to df2. Shows col difference, index difference, nans & 0s differences

_in ./plot.py_

- scatter: makes a hoverable/zoomable bokeh scatter plot
- bigScatter: 
- CNV_Map: makes a hoverable Copy Number plot using bokeh
- volcano: makes a searchable volcano plot for a Differential expression experiment using bokeh
- correlationMatrix: makes a hoverable bokeh correlation matrix, works with annotations, pvalues, clusters, etc.
- venn: makes a venn diagram from a list of sets
- mergeImages: merge multiple pngs/pdfs together into one
- addTextToImage: adds a text in an image to a specific location
- SOMPlot: a tool that uses simpSOM's package output (which produces self organizing maps), to plot its output in an interactive fashion.

## other necessary tools

_I am not creating anything that overlaps with that/ I am using these tools_

- os (python)
- subprocess (python)
- sns (python)
- bokeh (python)
