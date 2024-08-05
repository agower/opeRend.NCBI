#' Permissions for GEO Operations
#'
#' This object contains permission settings for adding GSE, GPL, GSM,
#' affymetrixCEL, and affymetrixCELSet Entities, as well as associated work files.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{geo}{A character vector with permissions "R", "U", "PR", "PU".}
#'   \item{_other}{A character vector with permission "R".}
#' }
#' @examples
#' opeRend.NCBI::geoPermissions
#' @seealso \code{\link[opeRend]{operendPermissions}}
#' @import opeRend
#' @export
geoPermissions <- opeRend::operendPermissions(
  "geo" = c("R", "U", "PR", "PU"),
  "_other" = c("R")
)

#' @name removeNull
#' @title Remove \code{NULL} values from a list
#' @param list_obj A list.
#' @returns A list with all \code{NULL} values removed.
removeNull <- function(list_obj) {
  Filter(function(x) length(x) != 0, list_obj)
}

#' @name addGEO
#' @title Add GEO Entities to Operend
#' @description
#' Wrapper for adding GEO entities.
#' @param class 
#' A character string specifying the class of the Entity record
#' @param variables
#' A list specifying variables of the Entity record to be added.
#' \code{NULL} values will be removed from the list prior to adding to Operend.
#' @returns
#' If the operation is successful, an \code{\linkS4class{operendEntity}} object.
#' @seealso \code{\link[opeRend]{addEntity}}
#' @author Dylan L. Tamayo \email{dltamayo@@bu.edu}
#' @import opeRend
addGEO <- function(class, variables) {
  opeRend::addEntity(
    class = class,
    variables = removeNull(variables),
    permissions = opeRend.NCBI::geoPermissions
  )
}

#' @name addGPL
#' @title Add GEO Platform to Operend
#' @description
#' This function queries Operend for the accession number of a GEO Platform; if 
#' the GPL has previously been uploaded, retrieves the \code{Entity} record,
#' otherwise extracts and uploads metadata from a GEO Platform as a new Operend
#' \code{Entity} record.
#' @param GPL 
#' A GEOquery GEO Platform entity.
#' @returns
#' If the operation is successful, an \code{\linkS4class{operendEntity}} object of the GEO Platform.
#' @examples
#' GSE <- getGEO("GSE994", GSEMatrix = FALSE)
#' lapply(GPLList(GSE), addGPL)
#' @author Dylan L. Tamayo \email{dltamayo@@bu.edu}
#' @import GEOquery
#' @import opeRend
#' @export
addGPL <- function(GPL) {
  # Retrieve GPL data
  metadata <- GEOquery::Meta(GPL)
  geo_accession <- metadata$geo_accession
  
  # Query if GPL accession is present in Operend
  query <- opeRend::listEntities("GEOPlatform", variables = list(geo_accession = geo_accession))
  
  # If GPL accession is not present in Operend, add GPL entity
  if(length(query)==0) {
    # Add to Operend
    gplList <- list(
      geo_accession = geo_accession,
      manufacturer = metadata$manufacturer,
      organism = metadata$organism,
      taxid = as.numeric(metadata$taxid),
      technology = metadata$technology,
      title = metadata$title
    )
    
    gplEntity <- addGEO(
      class = 'GEOPlatform',
      variables = gplList
    )
  } else {
    # If GPL accession is already present in Operend, return the oldest GPL entity
    gplEntity <- query[[length(query)]]
    cat(c('GEOPlatform record', objectId(gplEntity), 'retrieved.\n'))
  }
  
  # Return GEO Platform Entity
  return(gplEntity)
}

#' @name addAffymetrixCEL
#' @title Add affymetrixCEL file to Operend
#' @description
#' This function extracts and uploads the metadata and workfile of an affymetrixCEL
#' file as an Operend \code{Entity} record and \code{WorkFile}.
#' @param ftpUrl
#' Character vector of length 1 specifying the ftp link of the CEL file of a GEO Sample.
#' @returns
#' If the operation is successful, an \code{\linkS4class{operendEntity}} object of the affymetrixCEL.
#' @examples
#' addAffymetrixCEL(
#'   "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM15nnn/GSM15684/suppl/GSM15684.CEL.gz"
#' )
#' @author Dylan L. Tamayo \email{dltamayo@@bu.edu}
#' @import GEOquery
#' @import opeRend
#' @importFrom affyio read.celfile.header
#' @export
addAffymetrixCEL <- function(ftpUrl) {
  # Define a temporary file path to download the file. Clear file upon function exit
  filenames <- basename(ftpUrl)
  celFile <- file.path(tempdir(), filenames)
  on.exit(unlink(celFile))
  
  # Download the file
  download.file(ftpUrl, destfile=celFile)
  
  # Read the .CEL.gz file using affyio
  metadata <- processCELHeader(celFile)
  
  # Add CEL work file
  workFile <- opeRend::addWorkFile(celFile, fileType = ".CEL.gz",
                                   storage = "geo", permissions = opeRend.NCBI::geoPermissions)
  
  # Add to Operend
  affymetrixCELList <- list(
    workFile = opeRend::objectId(workFile),
    filenames = filenames,
    
    cdfName = metadata$cdfName,
    dimensions = as.numeric(metadata$dimensions),
    Datmetadata = metadata$Datmetadata,
    ScanDate = metadata$ScanDate
  )
  
  affymetrixCELEntity <- addGEO(
    class = 'AffymetrixCEL',
    variables = affymetrixCELList
  )
  
  # Return affymetrixCEL Entity
  return(affymetrixCELEntity)
}

#' @name addGSM
#' @title Add GEO Sample to Operend
#' @description
#' This function extracts and uploads metadata from a GEO Sample, including associated
#' affymetrixCEL files, as an Operend \code{Entity} record.
#' @param GSM
#' A GEOquery GEO Sample entity.
#' @returns
#' If the operation is successful, an \code{\linkS4class{operendEntity}} object of the GEO Sample.
#' @examples
#' GSE <- getGEO("GSE994", GSEMatrix = FALSE)
#' lapply(GSMList(GSE)[1:2], addGPL)
#' @author Dylan L. Tamayo \email{dltamayo@@bu.edu}
#' @import GEOquery
#' @import opeRend
#' @export
addGSM <- function(GSM) {
  # Retrieve ftp link from GSM metadata
  metadata <- GEOquery::Meta(GSM)
  supplementary_file <- metadata$supplementary_file
  
  # Download CEL file from ftp link, then retrieve UUID of object after adding.
  # If no CEL file is specified in metadata, set value to NULL
  if (is.null(supplementary_file)) {
    affymetrixCEL <- NULL
    affymetrixCELUuid <- NULL
  } else {
    affymetrixCEL <- addAffymetrixCEL(supplementary_file)
    affymetrixCELUuid <- opeRend::objectId(affymetrixCEL)
  }
  
  # Add to Operend
  gsmList <- list(
    affymetrixCEL = affymetrixCELUuid,
    
    channel_count = as.numeric(metadata$channel_count),
    characteristics_ch1 = metadata$characteristics_ch1,
    characteristics_ch2 = metadata$characteristics_ch2,
    data_processing = metadata$data_processing,
    description = metadata$description,
    geo_accession = metadata$geo_accession,
    molecule_ch1 = metadata$molecule_ch1,
    molecule_ch2 = metadata$molecule_ch2,
    platform_id = metadata$platform_id,
    series_id = metadata$series_id,
    source_name_ch1 = metadata$source_name_ch1,
    source_name_ch2 = metadata$source_name_ch2,
    title = metadata$title,
    type = metadata$type
  )
  
  gsmEntity <- addGEO(
    class = 'GEOSample',
    variables = gsmList
  )
  
  # Return affymetrixCEL and GSM Entity records. If GSM does not contain
  # a supplementary file, affymetrixCEL value will be NULL.
  return(list(affymetrixCEL = affymetrixCEL, geoSample = gsmEntity))
}

#' @name addGSE
#' @title Add GEO Series to Operend
#' @description
#' This function extracts and uploads metadata from a GEO Series, including associated
#' Platforms, Samples, and affymetrixCEL files, as an Operend \code{Entity} record.
#' @param GSE
#' A character vector of length 1 specifying a GEO Series accession number.
#' @returns
#' If the operation is successful, an \code{\linkS4class{operendEntity}} object of the GEO Series.
#' @examples
#' addGSE("GSE994")
#' @author Dylan L. Tamayo \email{dltamayo@@bu.edu}
#' @import GEOquery
#' @import opeRend
#' @export
addGSE <- function(GSE) {
  # Retrieve GSE data
  cat("Retrieving GEOSeries:\n")
  gseObj <- GEOquery::getGEO(GSE, GSEMatrix = FALSE)
  
  # Retrieve and add GPL data
  cat("Retrieving GPL data:\n")
  gplList <- GEOquery::GPLList(gseObj)
  
  # Add GPLs entities; unlist GPL IDs into vector
  cat("Uploading GPLs:\n")
  geoPlatforms <- unlist(lapply(gplList, function(x) opeRend::objectId(addGPL(x))))
  
  # Retrieve and add GSM data and corresponding affymetrixCEL file
  cat("Retrieving GSM data:\n")
  gsmList <- GEOquery::GSMList(gseObj)[1:2]                 #truncate during testing
  
  cat("Uploading GSM and affyCEL:\n")
  affyGsmIds <- lapply(gsmList, addGSM)
  
  # Separate affymetrixCEL and GSM Entity IDs; unlist into vectors
  affyIds <- removeNull(lapply(affyGsmIds, function(x) x$affymetrixCEL))
  affymetrixCELs <- unlist(lapply(affyIds, objectId))
  geoSamples <- unlist(lapply(affyGsmIds, function(x) opeRend::objectId(x$geoSample)))
  
  # Retrieve GSE metadata
  metadata <- GEOquery::Meta(gseObj)
  
  # Add affymetrixCELSet to Operend
  cat("Adding CEL set:\n")
  affymetrixCELSetList <- list(
    affymetrixCELs = affymetrixCELs,
    
    description = metadata$title,
    name = metadata$geo_accession
  )
  
  affymetrixCELSetEntity <- addGEO(
    class = 'AffymetrixCELSet',
    variables = affymetrixCELSetList
  )
  
  # Add GSE to Operend
  cat("Adding GSE:\n")
  gseList <- list(
    affymetrixCELSet = opeRend::objectId(affymetrixCELSetEntity),
    geoSamples = geoSamples,
    geoPlatforms = geoPlatforms,
    
    geo_accession = metadata$geo_accession,
    platform_id = metadata$platform_id,
    pubmed_id = as.numeric(metadata$pubmed_id),
    summary = metadata$summary,
    title = metadata$title,
    type = metadata$type
  )
  
  gseEntity <- addGEO(
    class = 'GEOSeries',
    variables = gseList
  )
  
  # List all files in the cache directory
  # cached_files <- list.files(tempdir(), full.names = TRUE)
  
  # Remove all cached files
  # unlink(cached_files, recursive = TRUE)
  
  # Return GEO Series Entity
  return(gseEntity)
}
