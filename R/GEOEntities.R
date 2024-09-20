#' @name addGPL
#' @title Add GEO Platform to Operend
#' @description
#' This function queries Operend for the accession number of a GEO Platform. If 
#' the GPL has previously been uploaded, retrieves the \code{Entity} record and
#' checks for any updates. If the GPL has not been uploaded, extracts and
#' uploads metadata from a GEO Platform as a new Operend \code{Entity} record.
#' @param GPL 
#' A character vector of length 1 specifying a GEO Platform accession number, OR
#' a GEOquery GEO Platform object.
#' @returns
#' If the operation is successful, an \code{\linkS4class{operendEntity}} object of the GEO Platform.
#' @examples
#' GPL <- getGEO("GPL96", GSEMatrix = FALSE)
#' addGPL(GPL)
#' @author Dylan L. Tamayo \email{dltamayo@@bu.edu}
#' @import GEOquery
#' @import opeRend
#' @export
addGPL <- function(GPL) {
  # Retrieve GPL data
  metadata <- GEOquery::Meta(retrieveGEOquery(GPL, "GPL"))
  geo_accession <- metadata$geo_accession
  last_update_date <- processDate(metadata$last_update_date)
  
  # Query if GPL accession is present in Operend
  query <- queryOperend("GEOPlatform", list(geo_accession = geo_accession))
  
  # Stage metadata
  gplList <- list(
    geo_accession = geo_accession,
    last_update_date = last_update_date,
    manufacturer = metadata$manufacturer,
    organism = metadata$organism,
    taxid = as.numeric(metadata$taxid),
    technology = metadata$technology,
    title = metadata$title
  )
  
  # If GPL accession is not present in Operend, add GPL entity
  if(is.null(query)) {
    gplEntity <- addGEO(
      class = "GEOPlatform",
      variables = gplList
    )
    
    return(gplEntity)
  }
    
  # If GPL accession is already present in Operend, return the oldest GPL entity
  cat(c(geo_accession, "already in Operend, retrieving record:\n"))
  cat(c("GEOPlatform record", opeRend::objectId(query), "retrieved.\n"))
    
  # Check if GPL entity requires updating
  if(needsUpdate(query, last_update_date)) {
    cat("Updating GEOPlatform entity:\n")
    gplEntity <- updateGEO(
      id = opeRend::objectId(query),
      variables = gplList
    )
    
    return(gplEntity)
  }
  
  # If no update is needed, return the existing entry
  return(query)
}

#' @name addAffymetrixCEL
#' @title Add affymetrixCEL file to Operend
#' @description
#' This function extracts and uploads the metadata and workfile of an affymetrixCEL
#' file as an Operend \code{Entity} record and \code{WorkFile}.
#' @param ftpUrl
#' Character vector of length 1 specifying the ftp link of the CEL file of a GEO Sample.
#' @param existingEntity
#' An optional parameter for an AffymetrixCEL Operend Entity, which, if provided,
#' will be compared with the ftpUrl CEL file to determine if the existingEntity
#' requires updating.
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
#' @importFrom utils download.file
#' @export
addAffymetrixCEL <- function(ftpUrl, existingEntity = NULL) {
  # Define a temporary file path to download the file. Clear file upon function exit
  filenames <- basename(ftpUrl)
  celFile <- file.path(tempdir(), filenames)
  on.exit(unlink(celFile))
  download.file(ftpUrl, destfile = celFile)
  
  # Read the .CEL.gz file using affyio
  metadata <- processCELHeader(celFile)
  
  # Stage metadata; workFile will be added later
  affymetrixCELList <- list(
    workFile = NULL,
    filenames = filenames,
    
    cdfName = metadata$cdfName,
    dimensions = as.numeric(metadata$dimensions),
    Datmetadata = metadata$Datmetadata,
    ScanDate = metadata$ScanDate
  )
  
  # If no existing entity provided, create new entity
  if(is.null(existingEntity)){
    workFile <- opeRend::addWorkFile(celFile, fileType = ".CEL.gz",
                                     storage = "geo", permissions = geoPermissions)
    
    affymetrixCELList$workFile <- opeRend::objectId(workFile)
    
    affymetrixCELEntity <- addGEO(
      class = "AffymetrixCEL",
      variables = affymetrixCELList
    )
    
    return(affymetrixCELEntity)
  }
  
  # Attempt to update existingEntity provided
  if (tools::md5sum(celFile) != getWorkFileProperties(existingEntity$workFile)@hash) {
    workFile <- opeRend::addWorkFile(celFile, fileType = ".CEL.gz",
                                     storage = "geo", permissions = geoPermissions)
    
    affymetrixCELList$workFile <- opeRend::objectId(workFile)
    
    affymetrixCELEntity <- updateGEO(
      id = objectId(existingEntity),
      variables = affymetrixCELList
    )
    
    # Trash old workFile
    updateWorkFileProperties(existingEntity$workFile, isTrashed = TRUE)
    
    return(affymetrixCELEntity)
  }
  
  # If existing entity does not require updating, return existing entity
  return(existingEntity)
}

#' @name addGSM
#' @title Add GEO Sample to Operend
#' @description
#' This function queries Operend for the accession number of a GEO Sample. If 
#' the GSM has previously been uploaded, retrieves the \code{Entity} record and
#' checks for any updates. If the GSM has not been uploaded, extracts and
#' uploads metadata from a GEO Sample and associated affymetrixCEL files as a
#' new Operend \code{Entity} record.
#' @param GSM
#' A character vector of length 1 specifying a GEO Sample accession number, OR
#' a GEOquery GEO Sample object.
#' @returns
#' If the operation is successful, an \code{\linkS4class{operendEntity}} object of the GEO Sample.
#' @examples
#' GSE <- GEOquery::getGEO("GSE994", GSEMatrix = FALSE)
#' lapply(GSMList(GSE)[1:2], addGPL)
#' @author Dylan L. Tamayo \email{dltamayo@@bu.edu}
#' @import GEOquery
#' @import opeRend
#' @export
addGSM <- function(GSM) {
  # Retrieve GSM data
  metadata <- GEOquery::Meta(retrieveGEOquery(GSM, "GSM"))
  geo_accession <- metadata$geo_accession
  last_update_date = processDate(metadata$last_update_date)
  ftpUrl <- getCELurl(metadata$supplementary_file)
  
  # Query if GSM accession is present in Operend
  query <- queryOperend("GEOSample", list(geo_accession = geo_accession))
  
  # Stage metadata. Add affymetrixCEL ID later
  gsmList <- list(
    affymetrixCEL = NULL,
    
    channel_count = as.numeric(metadata$channel_count),
    characteristics_ch1 = metadata$characteristics_ch1,
    characteristics_ch2 = metadata$characteristics_ch2,
    data_processing = metadata$data_processing,
    description = metadata$description,
    geo_accession = geo_accession,
    last_update_date = last_update_date,
    molecule_ch1 = metadata$molecule_ch1,
    molecule_ch2 = metadata$molecule_ch2,
    platform_id = metadata$platform_id,
    series_id = metadata$series_id,
    source_name_ch1 = metadata$source_name_ch1,
    source_name_ch2 = metadata$source_name_ch2,
    title = metadata$title,
    type = metadata$type
  )
  
  # If GPL accession is not present in Operend, add GPL entity
  if(is.null(query)) {
    # If CEL file specified in metadata, download CEL file,
    # then add to Operend and retrieve ID.
    if (!is.null(ftpUrl)) {
      gsmList$affymetrixCEL <- opeRend::objectId(addAffymetrixCEL(ftpUrl))
    }
    
    gsmEntity <- addGEO(
      class = "GEOSample",
      variables = gsmList
    )
    
    return(gsmEntity)
  }
  
  # If GSM accession is already present in Operend, return the oldest GPL entity
  cat("GSM accession number already in Operend, retrieving record:\n")
  cat(c("GEOSample record", opeRend::objectId(query), "retrieved.\n"))
  
  # Check if GSM entity requires updating
  if(needsUpdate(query, last_update_date)) {
    cat("Updating GEOSample entity:\n")
    
    # If supplementary file present, check if affymetrixCEL requires updating
    if (!is.null(ftpUrl)) {
      existingEntity <- getEntity(query$affymetrixCEL)
      gsmList$affymetrixCEL <- opeRend::objectId(addAffymetrixCEL(ftpUrl, existingEntity))
    }

    gsmEntity <- updateGEO(
      id = opeRend::objectId(query),
      variables = gsmList
    )
    
    return(gsmEntity)
  }
  
  # If no update is needed, return the existing entry
  return(query)
}

#' @name addGSE
#' @title Add GEO Series to Operend
#' @description
#' This function queries Operend for the accession number of a GEO Series. If 
#' the GSE has previously been uploaded, retrieves the \code{Entity} record and
#' checks for any updates. If the GSE has not been uploaded, extracts and
#' uploads metadata from a GEO Series and associated Platforms, Samples, and 
#' affymetrixCEL files as a new Operend \code{Entity} record.
#' @param GSE
#' A character vector of length 1 specifying a GEO Series accession number OR
#' a GEOquery GEO Series object.
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
  gseObj <- retrieveGEOquery(GSE, "GSE")
  
  # Remove all cached files after completion
  on.exit({
    cached_files <- list.files(tempdir(), full.names = TRUE)
    unlink(cached_files, recursive = TRUE)
  })
  
  metadata <- GEOquery::Meta(gseObj)
  geo_accession <- metadata$geo_accession
  last_update_date <- processDate(metadata$last_update_date)
  
  # Query if GSE accession is present in Operend
  query <- queryOperend("GEOSeries", list(geo_accession = geo_accession))
  
  # Helper function to process AffymetrixCELSets.
  # Returns NULL if no GSMs in the GSE contain AffymetrixCELs
  addAffyCELSets <- function(geoSamples, metadata) {
    GSMsplitByGPL <- split(geoSamples, sapply(geoSamples, function(x) x$platform_id))
    
    CELSetArray <- lapply(metadata$platform_id, function(GPL){
      # Retrieve affymetrixCELs corresponding to each GPL
      affys <- lapply(GSMsplitByGPL[[GPL]], function(GSM){
        GSM$affymetrixCEL
      })
      affymetrixCELs <- unlist(removeNull(affys))
      
      # Early return if there are no affymetrixCELs to upload as a set
      if(is.null(affymetrixCELs)) {
        return(NULL)
      }
      
      cat("Uploading associated affymetrixCEL sets:\n")
      affymetrixCELSetList <- list(
        affymetrixCELs = affymetrixCELs,
        
        description = metadata$title,
        name = metadata$geo_accession
      )
      
      affymetrixCELSetEntity <- addGEO(
        class = "AffymetrixCELSet",
        variables = affymetrixCELSetList
      )
      
      return(affymetrixCELSetEntity)
    })
    
    return(removeNull(CELSetArray))
  }
  
  # Helper function to process and add GPL and GSM data
  stageGSEList <- function(gseObj, metadata) {
    gseList <- list(
      affymetrixCELSet = NULL,
      geoSamples = NULL,
      geoPlatforms = NULL,
      
      geo_accession = metadata$geo_accession,
      last_update_date = last_update_date,
      platform_id = metadata$platform_id,
      pubmed_id = as.numeric(metadata$pubmed_id),
      relation = metadata$relation,
      summary = metadata$summary,
      title = metadata$title,
      type = metadata$type
    )
    
    # Retrieve and add GPL data
    cat("Retrieving GPL data:\n")
    gplList <- GEOquery::GPLList(gseObj)
    
    cat("Uploading GPLs:\n")
    geoPlatforms <- lapply(gplList, addGPL)
    gseList$geoPlatforms <- unlist(lapply(geoPlatforms, opeRend::objectId))
    
    # Retrieve and add GSM data
    cat("Retrieving GSM data:\n")
    gsmList <- GEOquery::GSMList(gseObj) #[1:2]                 #truncate during testing
    
    cat("Uploading GSMs:\n")
    geoSamples <- lapply(gsmList, addGSM)
    gseList$geoSamples <- unlist(lapply(geoSamples, opeRend::objectId))
    
    # If GSE is a SuperSeries, add SubSeries and skip adding CEL sets
    subseries <- getSubseries(metadata$relation)
    if(!is.null(subseries)) {
      cat("Uploading SubSeries:\n")
      lapply(subseries, addGSE)
      return(gseList)
    }
    
    # Add associated CEL sets, if present
    affymetrixCELSets <- addAffyCELSets(geoSamples, metadata)
    if(!is.null(affymetrixCELSets)) {
      gseList$affymetrixCELSet <- unlist(lapply(affymetrixCELSets, opeRend::objectId))
    }
    
    return(gseList)
  }
  
  # If GSE accession is not present in Operend, add GSE entity
  if(is.null(query)) {
    cat("Adding GSE:\n")
    gseList <- stageGSEList(gseObj, metadata)
    
    gseEntity <- addGEO(
      class = "GEOSeries",
      variables = gseList
    )
    
    return(gseEntity)
  }
  
  # If GSE accession is already present in Operend, return the oldest GSE entity
  cat("GSE accession number already in Operend, retrieving record:\n")
  cat(c("GEOSeries record", opeRend::objectId(query), "retrieved.\n"))
  
  # Check if GSE entity requires updating
  if(needsUpdate(query, last_update_date)) {
    cat(c("Updating GEOSeries record", opeRend::objectId(query), ".\n"))
    gseList <- stageGSEList(gseObj, metadata)
    
    # Delete old affymetrixCELSets after setting variable to NULL
    old_affymetrixCELSet <- query$affymetrixCELSet
    opeRend::updateEntity(id = opeRend::objectId(query), variables = list(affymetrixCELSet = NULL))
    lapply(old_affymetrixCELSet, opeRend::deleteEntity)
    
    gseEntity <- updateGEO(
      id = opeRend::objectId(query),
      variables = gseList
    )
    
    return(gseEntity)
  }
  
  # If no update is needed, return the existing entry
  return(query)
}
