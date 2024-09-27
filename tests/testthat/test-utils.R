library(testthat)
library(opeRend.NCBI)

test_that("processDate handles various cases", {
  # Standard date format (valid input)
  date1 <- "Jan 01 2000"
  result <- processDate(date1)
  expected <- as(as.Date("2000-01-01"), "operendDate")
  expect_equal(result, expected)
  
  # Another valid date input (February 29th on a leap year)
  date2 <- "Feb 29 2020"
  result <- processDate(date2)
  expected <- as(as.Date("2020-02-29"), "operendDate")
  expect_equal(result, expected)
  
  # Invalid date format (missing day)
  date3 <- "Feb 2020"
  expect_error(processDate(date3))
  
  # Invalid date (February 30th)
  date4 <- "Feb 30 2020"
  expect_error(processDate(date4))

  # Empty string as input (invalid input)
  date5 <- ""
  expect_error(processDate(date5))

  # NULL input (edge case handling)
  expect_error(processDate(NULL))
  
  # Multiple string input
  expect_error(processDate(c(date1, date1)))

  # Edge case - handling very old date
  date6 <- "Jan 01 1900"
  result <- processDate(date6)
  expected <- as(as.Date("1900-01-01"), "operendDate")
  expect_equal(result, expected)

  # Edge case - handling far future date
  date7 <- "Dec 31 2999"
  result <- processDate(date7)
  expected <- as(as.Date("2999-12-31"), "operendDate")
  expect_equal(result, expected)
})
  
test_that("removeNull handles various cases", {
  list1 <- list(a = NULL, b = character(0), c = numeric(0), d = "test")
  list2 <- list(NULL)
  list3 <- list(a = NULL, b = character(0), c = c(numeric(0), NULL), d = c("test", NULL))
  list4 <- list(a = "test1", b = list(b1 = "test2", b2 = "test3"))
  list5 <- list(a = c("", ""), b = "non-empty", c = character(0))
  list6 <- list(a = numeric(0), b = NULL, c = "", d = 5)
  
  # List with NULL, empty character, empty numeric, and valid string
  result <- removeNull(list1)
  expect_equal(result, list(d = "test"))
  
  # List containing only NULL should return an empty list
  result <- removeNull(list2)
  expect_equal(result, list())
  
  # Nested list with NULLs and empty values should return only the valid string
  result <- removeNull(list3)
  expect_equal(result, list(d = "test"))
  
  # Nested list with valid strings should remain unchanged
  result <- removeNull(list4)
  expect_equal(result, list4)
  
  # List with an empty string and valid string, should only return the valid string
  result <- removeNull(list5)
  expect_equal(result, list(a = c("", ""), b = "non-empty"))
  
  # List with mixed data types: numeric(0), NULL, empty string, and valid numeric
  result <- removeNull(list6)
  expect_equal(result, list(c = "", d = 5))
})

test_that("needsUpdate handles various cases", {
  entity1 <- list(last_update_date = processDate("March 09 2021"))
  entity2 <- list(last_update_date = NULL)
  date1 <- processDate("Mar 9 2021")
  date2 <- processDate("January 01 1900")
  date3 <- processDate("January 01 2300")
  
  # Operend and GEO dates are the same (no update needed)
  result <- needsUpdate(entity1, date1)
  expect_equal(result, FALSE)
  
  # GEO date is earlier than Operend date (update needed)
  result <- needsUpdate(entity1, date2)
  expect_equal(result, TRUE)
  
  # GEO date is later than Operend date (update needed)
  result <- needsUpdate(entity1, date3)
  expect_equal(result, TRUE)
  
  # Operend date is NULL, GEO date is present (update needed)
  result <- needsUpdate(entity2, date1)
  expect_equal(result, TRUE)
  
  # GEO date is NULL, Operend date is present (warning, update needed)
  warning1 <- "updating Operend entry:\n"
  expect_warning(needsUpdate(entity1, NULL), warning1)
  
  # Both GEO and Operend dates are NULL (warning, no update needed)
  warning2 <- "skipping update:\n"
  expect_warning(needsUpdate(entity2, NULL), warning2)
})

test_that("getCELurl handles various cases", {
  ftpLink1 <- "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM15nnn/GSM15684/suppl/GSM15684.CEL.gz"
  ftpLink2 <- "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM15nnn/GSM15684/suppl/GSM15685.CEL.gz"
  nonCelLink <- "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM15nnn/GSM15684/suppl/GSM15684.txt.gz"
  
  # Correct matching
  result <- getCELurl(ftpLink1)
  expect_equal(result, ftpLink1)
  
  # Multiple valid links, expect first
  result <- getCELurl(c(ftpLink2, ftpLink1))
  expect_equal(result, ftpLink2)
  
  # No matching CEL file, return NULL
  result <- getCELurl(nonCelLink)
  expect_equal(result, NULL)
  
  # Input contains NULL
  result <- getCELurl(c(NULL, ftpLink1, ftpLink2))
  expect_equal(result, ftpLink1)
  
  # NULL input
  result <- getCELurl(NULL)
  expect_equal(result, NULL)
  
  # Empty input should return NULL
  result <- getCELurl(character(0))
  expect_equal(result, NULL)
  
  # Mixed with a non-CEL link
  result <- getCELurl(c(nonCelLink, ftpLink1))
  expect_equal(result, ftpLink1)
})

test_that("getSubseries handles various cases", {
  relationVector1 <- c("SubSeries of: GSE210661", "BioProject: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA864867")
  relationVector2 <- c("SuperSeries of: GSE210271", "SuperSeries of: GSE210272", "BioProject: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA866709")
  relationVector3 <- c("SuperSeries of: GSM210271", "BioProject: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA866709")
  
  # Case where there are no "SuperSeries of" matches
  result <- getSubseries(relationVector1)
  expect_equal(result, NULL)
  
  # Case where there are valid "SuperSeries of" matches
  result <- getSubseries(relationVector2)
  expect_equal(result, c("GSE210271", "GSE210272"))

  # Case where input is empty
  result <- getSubseries(character(0))
  expect_equal(result, NULL)

  # Case where input is NULL
  result <- getSubseries(NULL)
  expect_equal(result, NULL)

  # Case with incorrect formatting (no GSE)
  result <- getSubseries(relationVector3)
  expect_equal(result, NULL)
  
})

test_that("retrieveGEOquery handles various cases", {
  GPLObject <- GEOquery::getGEO("GPL96", GSEMatrix = FALSE)
  
  # Check if passing an existing GEO object returns the same object without re-retrieving it
  result <- retrieveGEOquery(GPLObject, "GPL")
  expect_equal(result, GPLObject)
  
  # Check if an error is thrown when passing a GEO object with the wrong class
  expect_error(retrieveGEOquery(GPLObject, "GSM"))
  
  
  # Check if a valid GEO character string returns the expected GEO object
  GPLCharacter <- "GPL96"
  result <- retrieveGEOquery(GPLCharacter, "GPL")
  expect_equal(result, GPLObject)

  # Check if an error is thrown when "GEO" is not a character string (e.g., passing an integer)
  expect_error(retrieveGEOquery(1, "GPL"))
  
  # Check if an error is thrown when "GEO" is a character vector with length greater than 1
  expect_error(retrieveGEOquery(c("GPL96", "GPL96"), "GPL"))
  
  # Check if an error is thrown when "GEO" specifies a different class
  expect_error(retrieveGEOquery(GPLCharacter, "GSE"))
})

test_that("queryOperend handles various cases", {
  expected <- opeRend::listEntities("GEOPlatform", list(geo_accession = "GPL96"))[[1]]
  
  # Check querying using characters, OperendDates, and numerics
  result <- queryOperend("GEOPlatform", list(geo_accession = "GPL96"))
  expect_equal(result, expected)
  
  result <- queryOperend("GEOPlatform", list(last_update_date = processDate("Mar 9 2021")))
  expect_equal(result, expected)
  
  result <- queryOperend("GEOPlatform", list(taxid = 9606))
  expect_equal(result, expected)
  
  # Check if queries without results return NULL
  result <- queryOperend("GEOPlatform", list(geo_accession = "invalid"))
  expect_equal(result, NULL)
  
  # Check if invalid class names cause error
  expect_error(queryOperend("invalidClassName"))
})
