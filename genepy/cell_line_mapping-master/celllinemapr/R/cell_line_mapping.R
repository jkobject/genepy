find.non.unique <- function(x) {
    b <- table(x) > 1
    names(b)[b]
}

make.mapping <- function(full.map, input.type, output.type, id.subset, check.unique.mapping) {
    x <- as.character(full.map[[output.type]])
    n <- as.character(full.map[[input.type]])

    mask <- n %in% id.subset

    # collapse duplicate rows
    df <- data.frame(x=x[mask], n=n[mask], stringsAsFactors =F)
    df <- unique(df)

    x <- df$x
    n <- df$n

    if(check.unique.mapping) {
        non.unique.inputs <- find.non.unique(n)
        non.unique.outputs <- find.non.unique(x)
        if(length(non.unique.inputs) > 0) {
            stop(paste0("The following had nonunique values: ", paste0(non.unique.inputs, collapse=", ")))
        }

        if(length(non.unique.outputs) > 0) {
            stop(paste0("The following had nonunique values: ", paste0(non.unique.outputs, collapse=", ")))
        }
    }

    names(x) <- n
    x
}

read.mapping <- (function() {
    # save value from previous read to avoid fetching from url every time
    cell.line.mapping.cache <- NULL
    mapping.url <- getOption("celllinemapr.url", "../naming.csv")
    cache.path <- getOption("celllinemapr.cache.path", "../.celllinemapr.Rds")

    function(force=F) {
      if(is.null(cell.line.mapping.cache) || force) {
            mapping <- try(read.csv("../naming.csv"))
            if(class(mapping) == "try-error") {
                # if we got an error, then warn user that this failed and try loading from cache file.
                warning(paste0("Could not fetch mapping from ", mapping.url, ", attempting to read most recent cached mapping from ", cache.path))
                mapping <- readRDS(cache.path)
            } else {
              stopifnot(is.data.frame(mapping))
              saveRDS(mapping, file=cache.path)              
            }
            stopifnot(is.data.frame(mapping))
            cell.line.mapping.cache <<- mapping
        }
      cell.line.mapping.cache
    }
})()


name.mapper <- function(input.type, input.names, output.type, ignore.problems, check.unique.mapping, read.mapping.fn) {
    full.mapping <- read.mapping.fn()
    mapping <- make.mapping(full.mapping, input.type, output.type, input.names, check.unique.mapping)
    result <- mapping[input.names]
    if(!ignore.problems) {
        bad.names <- input.names[is.na(result)]
        if(length(bad.names) > 5) {
            bad.names <- c(bad.names[1:5], "...")
        }
        if(length(bad.names) > 0) {
            stop(paste0("Could not find cell lines (searching by ", input.type, ") for ", paste(bad.names, collapse=", ")))
        }
    }
    result
}

# returns a function to get cell line mapping. Returns the default function if 
pick.mapping.fn <- function(mapping) {
    if(!is.null(mapping)) {
        stopifnot(is.data.frame(mapping))
        return (function() {
            return (mapping)
        } )
    } else {
        return(read.mapping)
    }
}

#' Map cell line Broad ID (aka ArxSpan IDs) to the latest CCLE names
#'
#' @param arxspan.ids A vector of arxspan ids. These are always of the form "ACH-XXXXXX"
#' @param ignore.problems if not set to True, any unknown cell lines will result in an error being thrown. If you set to True, then you'll get NA for unknown lines instead.
#' @param check.unique.mapping if set, will throw an error if it discovers two different IDs which map to the same CCLE name (which could cause issues downstream)
#' @param mapping if set, will use this dataframe for the mapping instead of fetching the latest
#' @examples
#' ccle_names <- arxspan.to.ccle(c('ACH-000007', 'ACH-000008'))
#' @export arxspan.to.ccle
arxspan.to.ccle <- function(arxspan.ids, ignore.problems=F, check.unique.mapping=T, mapping=NULL) {
    name.mapper('broad_id', arxspan.ids, 'canonical_ccle_name', ignore.problems, check.unique.mapping, pick.mapping.fn(mapping))
}

#' Map ccle names to Broad ID (aka ArxSpan IDs)
#'
#' @param ccle.names A vector of CCLE names
#' @param ignore.problems if not set to True, any unknown cell lines will result in an error being thrown. If you set to True, then you'll get NA for unknown lines instead.
#' @param check.unique.mapping if set, will throw an error if it discovers two different IDs which map to the same arxspan id (which could cause issues downstream)
#' @param mapping if set, will use this dataframe for the mapping instead of fetching the latest
#' @examples
#' broad_ids <- ccle.to.arxspan(c('HS294T_SKIN','NCIH1581_LUNG'))
#' @export ccle.to.arxspan
ccle.to.arxspan <- function(ccle.names, ignore.problems=F, check.unique.mapping=T, mapping=NULL) {
    name.mapper('ccle_name', ccle.names, 'broad_id', ignore.problems, check.unique.mapping, pick.mapping.fn(mapping))
}

#' Map any ccle names to the current/latest ccle names. Useful for updating old names and correcting lines which have been renamed.
#'
#' @param ccle.names A vector of CCLE names
#' @param ignore.problems if not set to True, any unknown cell lines will result in an error being thrown. If you set to True, then you'll get NA for unknown lines instead.
#' @param check.unique.mapping if set, will throw an error if it discovers two different IDs which map to the same ccle name (which could cause issues downstream)
#' @param mapping if set, will use this dataframe for the mapping instead of fetching the latest
#' @examples
#' ccle_names <- ccle.to.latest('HEL9217_2013_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE')
#' @export ccle.to.latest
ccle.to.latest <- function(arxspan.ids, ignore.problems=F, check.unique.mapping=T, mapping=NULL) {
    name.mapper('ccle_name', arxspan.ids, 'canonical_ccle_name', ignore.problems, check.unique.mapping, pick.mapping.fn(mapping))
}

