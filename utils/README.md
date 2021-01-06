# Utils

1. a set of plotting tools based on [matplotlib]() and [bokeh]() 
2. helper functions to save data, generate random strings, run tasks in parallel etc.. 

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


_in ./plot.py_

- scatter: makes a hoverable/zoomable bokeh scatter plot
- CNV_Map: makes a hoverable Copy Number plot using bokeh
- volcano: makes a searchable volcano plot for a Differential expression experiment using bokeh
- plotCorrelationMatrix: makes a hoverable bokeh correlation matrix
- venn: makes a venn diagram from a list of sets
- mergeImages: merge multiple pngs/pdfs together into one
- addTextToImage adds a text in an image to a specific location


## other necessary tools

_I am not creating anything that overlaps with that/ I am using these tools_

- os (python)
- subprocess (python)
- sns (python)
- bokeh (python)