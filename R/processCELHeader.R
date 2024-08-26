#' @importFrom affyio read.celfile.header

#' @rdname processCELHeader
#' @name Process CEL header
#' @title Process the header of a CEL file
#' @description
#' This function extracts and processes the header of a CEL file.
#' @param filename
#' A character vector of length 1 containing the full path to a CEL file
#' @return
#' This function returns a named list containing the following elements:
#' \describe{
#'   \item{cdfName}{
#'     A character string containing the \code{cdfName} field of the header
#'   }
#'   \item{dimensions}{
#'     A two-element, unnamed integer vector containing the
#'     \code{CEL dimensions} field of the header
#'   }
#'   \item{DatHeader}{
#'     A character string containing the \code{DatHeader} field of the header
#'   }
#'   \item{ScanDate}{
#'     A character string containing the \code{ScanDate} field of the header,
#'     with extraneous characters removed
#'   }
#' }
#' @author Adam C. Gower \email{agower@@bu.edu}

processCELHeader <- function (filename)
{
  # If the CEL file header can be read, the 'ScanDate' field of the file
  # header is returned; otherwise, the function terminates with an error
  # message.
  header <- try(affyio::read.celfile.header(filename, info="full"), silent=TRUE)
  if (!inherits(header, "try-error")) {
    # Extract selected fields from the CEL file header to a named list
    result <- list(
      cdfName = header$cdfName,
      dimensions = unname(header[["CEL dimensions"]]),
      DatHeader = header$DatHeader,
      # ScanDate entries have been observed to take the following forms
      # (sometimes with leading spaces):
      # "mm/dd/yy hh:mm:ss" (most common), "yyyy-mm-ddThh:mm:ssZ",
      # "yyyy-mm-ddThh:mm:ss", "mm/dd/ hh:mm:ss", as.character(NA)
      ScanDate = sub(
        # Replace any instances of '/ ' or 'T' with a space
        "/ |T", " ",
        # Remove any leading spaces or trailing 'Z' characters
        sub("^[ ]*([^Z]+)Z?", "\\1", header$ScanDate)
      )
    )
  } else {
    stop("A CEL file header could not be read from file ", sQuote(filename))
  }
  result
}
