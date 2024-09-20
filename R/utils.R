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
#' geoPermissions
#' @seealso \code{\link[opeRend]{operendPermissions}}
#' @import opeRend
#' @export
geoPermissions <- opeRend::operendPermissions(
  "geo" = c("R", "U", "PR", "PU"),
  "_other" = c("R")
)

#' @name processDate
#' @title Format dates for Operend date variables.
#' @param d A character vector of length 1 representing a date.
#' @returns A character vector of length 1 representing a formatted operendDate. 
#' @examples 
#' processDate("Feb 19 2004")
#' @importFrom methods as
processDate <- function(d) {
  as(as.Date(d, format="%b %d %Y"), "operendDate")
}

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
    permissions = geoPermissions
  )
}

#' @name updateGEO
#' @title Update Operend GEO Entities
#' @description
#' Wrapper for updating GEO entities.
#' @param id 
#' A character string specifying the id of the Entity record
#' @param variables
#' A list specifying variables of the Entity record to be added.
#' \code{NULL} values will be removed from the list prior to adding to Operend.
#' @returns
#' If the operation is successful, an \code{\linkS4class{operendEntity}} object.
#' @seealso \code{\link[opeRend]{addEntity}}
#' @author Dylan L. Tamayo \email{dltamayo@@bu.edu}
#' @import opeRend
updateGEO <- function(id, variables) {
  opeRend::updateEntity(
    id = id,
    variables = removeNull(variables)
  )
}

#' @name needsUpdate
#' @title Logic for comparing last updated dates.
#' @description
#' Helper function that returns TRUE if two dates are not equal (therefore
#' requiring an update), or FALSE if two dates are equal (therefore indicating
#' that the provided Operend entity is up to date with the NCBI entry).
#' @param entity
#' An Operend GEO entity containing a last_update_date variable.
#' @param last_update_date 
#' A character vector of length 1 representing an operendDate.
#' @returns A character vector of length 1 representing a formatted operendDate. 
needsUpdate <- function(entity, last_update_date) {
  # If Operend entry date is NULL, it needs an update
  if (is.null(entity$last_update_date) && !is.null(last_update_date)) {
    return(TRUE)
  }
  
  # If GEO entry is NULL but Operend entry is not NULL, assume update is needed
  if (!is.null(entity$last_update_date) && is.null(last_update_date)) {
    warning("Warning: GEO entry update date is NULL, updating Operend entry:\n")
    return(TRUE)
  }
  
  # If both entries are NULL, assume no update needed
  if (is.null(entity$last_update_date) && is.null(last_update_date)) {
    warning("Warning: Update dates for both Operend and GEO entries are NULL, skipping update:\n")
    return(FALSE)
  }
  
  # TRUE if entity is not up to date, FALSE if up to date
  entity$last_update_date != last_update_date
}

#' @name getCELurl
#' @title Return first matching CEL file ftp link
#' @description
#' Helper function that checks character vector for an ftp link to a CEL file.
#' @param supplementaryFile
#' A character vector or \code{NULL}.
#' @returns A character vector of length 1 representing a ftp link to a CEL file,
#' OR \code{NULL}.
getCELurl <- function(supplementaryFile) {
  # Check if supplementaryFile is not NULL
  if (!is.null(supplementaryFile)) {
    # Find matches for "ftp://(...).CEL.gz" (case insensitive)
    matches <- supplementaryFile[grepl("ftp://(.+)\\.CEL\\.gz$", supplementaryFile, ignore.case = TRUE)]
    
    # If there's at least one match, return the first one
    if (length(matches) > 0) {
      return(matches[1])
    }
  }
  
  # If supplementaryFile is NULL or no matches are found, return NULL
  return(NULL)
}

#' @name getSubseries
#' @title Return GSE accession numbers of subseries
#' @description
#' Helper function that extracts
#' @param relationVector
#' A character vector corresponding to the relation metadata of a GSE
#' @returns A character vector representing matches for GSE subseries, OR \code{NULL}.
getSubseries <- function(relationVector) {
  subseries <- relationVector[grepl("SuperSeries of:", relationVector)]
  subseriesGSEs <- unlist(regmatches(subseries, gregexpr("GSE[0-9]+", subseries)))
  
  # Return NULL if no matches are found
  if (length(subseriesGSEs) == 0) {
    return(NULL)
  }
  
  return(subseriesGSEs)
}

# Wrapper function that returns either the submitted GEOquery object or retrieves GEOquery object based on accession number
#' @name retrieveGEOquery
#' @title Retrieve GEOquery object based on th
#' @param GEO A character vector of length 1 specifying a GEO accession number,
#' or a GEOquery GEO object.
#' @param GEOqueryClass A character vector of length 1 specifying the class of
#' the GEOquery GEO object to be returned.
#' @importFrom methods is
#' @returns A GEOquery GEO object.
retrieveGEOquery <- function(GEO, GEOqueryClass) {
  # Check if the input is already of the desired class
  if (is(GEO, GEOqueryClass)) {
    return(GEO)
  }
  
  # Check for valid inputs
  if (!(is.character(GEO) && length(GEO) == 1)) {
    stop("Argument 'GEO' must be a character vector of length 1") 
  }
  
  if (!grepl(paste0("^", GEOqueryClass), GEO)) {
    stop(paste("Invalid GEO accession number. Expected prefix:", GEOqueryClass))
  }
  
  # Retrieve the GEOquery object if all checks pass
  GEOquery::getGEO(GEO, GSEMatrix = FALSE)
}

#' @name queryOperend
#' @title Query Operend for GEO accession number
#' @description
#' Query Operend for the GEO accession number of a GEO entity, returning the
#' least recently added entry.
#' @param class A character vector of length 1 specifying the Operend \code{EntityClass}
#' of a GEO entity.
#' @param variables A named list indicating variables in the Operend Entity to search for.
#' @returns An \code{\linkS4class{operendEntity}} object, or \code{NULL}.
queryOperend <- function(class, variables) {
  # Query if accession number for class is present in Operend
  query <- opeRend::listEntities(class, variables)
  queryLength <- length(query)
  
  # Return NULL if no entity available
  if(queryLength == 0) {
    return(NULL)
  }
  
  if(queryLength > 1) {
    warning("More than one entry found, returning least recent entity:\n")
  }
  
  # Return least recent Operend entity
  return(query[[queryLength]])
}  
