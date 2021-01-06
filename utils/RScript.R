MofaRun <- function(valueList){
	library(reticulate)
	MOFAobject <- createMOFAobject(valueList)
	DataOptions <- getDefaultDataOptions()
	ModelOptions <- getDefaultModelOptions(MOFAobject)
	TrainOptions <- getDefaultTrainOptions()
	ModelOptions$numFactors <- 200
	TrainOptions$DropFactorThreshold <- 0.02
	MOFAobject <- prepareMOFA(
	  MOFAobject, 
	  DataOptions = DataOptions,
	  ModelOptions = ModelOptions,
	  TrainOptions = TrainOptions
	)
	MOFAobject <- runMOFA(MOFAobject)
	return(MOFAobject)
}

# Function that converts the segmented data to be continuous (so can plot chromosomes in 1, 2, 3, 4... order)
generate_chromosome_cutoffs_list <- function(cyto_band_file="data/hg38_cytoband.gz") {
  # Have to edit the chr values to 
  chr_bp_cutoffs <- read_tsv(cyto_band_file, col_names = F)
  cutoffs <- chr_bp_cutoffs %>% 
    group_by(X1) %>% 
    dplyr::summarize(pos=max(X3)) %>%
    mutate(X1=gsub('chr', '', X1)) %$% 
    setNames(pos, ifelse(X1 %in% seq(1,21), paste0('chr', as.integer(X1) + 1), ifelse(X1==22, 'chrX', ifelse(X1=='X', 'chrY', 'chrZ'))))
  
  cutoffs_final <- cutoffs[paste0('chr',c(seq(2, 22), 'X', 'Y'))] %>% cumsum()
  cutoffs_final['chr1'] = 0
  
  return(cutoffs_final)
}