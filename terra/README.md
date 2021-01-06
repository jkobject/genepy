# terra

a file containing a set of functions that uses [dalmatian](github.com/broadinstitute/dalmatian) to interact with the [GCP](https://cloud.google.com/storage/docs/gsutil) powered genomics HPC platform: [Terra](www.terra.bio). 
They contain a list of additional functions to do more than what is available in dalmatian

#### Available functions:

- createManySubmissions: 
- waitForSubmission: 
- uploadFromFolder: 
- updateAllSampleSet: 
- addToSampleSet: 
- addToPairSet: 
- changeGSlocation: 
- renametsvs: 
- findBackErasedDuplicaBamteFromTerraBucket: 
- shareTerraBams: 
- saveConfigs: 
- cleanWorkspace: 
- saveOmicsOutput: 
- shareCCLEbams: 
- changeToBucket:

## highly recommanded

- firecloud-dalmatian (python)
- gsutil
- nextflow (better than terra)