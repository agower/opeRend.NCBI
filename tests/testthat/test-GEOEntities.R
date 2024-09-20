library(testthat)
library(opeRend)
library(opeRend.NCBI)

test_that("addGPL uploads and updates properly", {
  test <- "GPL4"
  
  # Find another GPL to test if already in Operend
  if(!is.null(queryOperend("GEOPlatform", list(geo_accession = test)))) {
    stop("Test GPL already in Operend, change to another GPL\n")
  }
  
  # Check if GPL is added to Operend properly
  result <- addGPL(test)
  testId <- opeRend::objectId(result)
  metadata <- GEOquery::Meta(GEOquery::getGEO(test, GSEMatrix = FALSE))
  expected <- removeNull(list(
    geo_accession = metadata$geo_accession,
    last_update_date = processDate(metadata$last_update_date),
    manufacturer = metadata$manufacturer,
    organism = metadata$organism,
    taxid = as.numeric(metadata$taxid),
    technology = metadata$technology,
    title = metadata$title
  ))
  
  expect_equal(result@.Data, expected, ignore_attr = TRUE)
  
  # Adding same GPL does not result in additional entries
  addGPL(test)
  expect_equal(length(opeRend::listEntities("GEOPlatform", list(geo_accession = test))), 1)
  
  # Test updating variables
  expected2 <- expected
  expected2[["last_update_date"]] <- processDate("Jan 01 1990")
  expected2[["organism"]] <- "test"
  
  result2 <- opeRend::updateEntity(testId, variables = expected2)
  expect_equal(result2@.Data, expected2, ignore_attr = TRUE)
  
  # Test updating Operend GPL entry
  result3 <- addGPL(test)
  expect_equal(result3@.Data, expected, ignore_attr = TRUE)
  
  # Clean up after testing
  opeRend::deleteEntity(testId)
})
