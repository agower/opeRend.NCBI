library(testthat)
library(opeRend)
library(opeRend.NCBI)
library(GEOquery)

test_that("addGPL uploads and updates properly", {
  test <- "GPL4"
  
  # Find another GPL to test if already in Operend
  if(!is.null(queryOperend("GEOPlatform", list(geo_accession = test)))) {
    stop("Test GPL already in Operend, change to another GPL\n")
  }
  
  # Check if GPL is added to Operend properly
  result <- addGPL(test)
  testId <- opeRend::objectId(result)
  
  # Clean up after testing
  on.exit(opeRend::deleteEntity(testId))
  
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
  expected2[["title"]] <- "test"
  
  result2 <- opeRend::updateEntity(testId, variables = expected2)
  expect_equal(result2@.Data, expected2, ignore_attr = TRUE)
  
  # Test updating Operend GPL entry
  result3 <- addGPL(test)
  expect_equal(result3@.Data, expected, ignore_attr = TRUE)
})

test_that("addGSM uploads and updates properly", {
  test <- "GSM1"
  
  # Find another GPL to test if already in Operend
  if(!is.null(queryOperend("GEOSample", list(geo_accession = test)))) {
    stop("Test GSM already in Operend, change to another GSM\n")
  }
  
  # Check if GSM is added to Operend properly
  result <- addGSM(test)
  testId <- opeRend::objectId(result)
  
  # Clean up after testing
  on.exit(opeRend::deleteEntity(testId))
  
  metadata <- GEOquery::Meta(GEOquery::getGEO(test, GSEMatrix = FALSE))
  expected <- removeNull(list(
    affymetrixCEL = NULL,
    
    channel_count = as.numeric(metadata$channel_count),
    characteristics_ch1 = metadata$characteristics_ch1,
    characteristics_ch2 = metadata$characteristics_ch2,
    data_processing = metadata$data_processing,
    description = metadata$description,
    geo_accession = metadata$geo_accession,
    last_update_date = processDate(metadata$last_update_date),
    molecule_ch1 = metadata$molecule_ch1,
    molecule_ch2 = metadata$molecule_ch2,
    platform_id = metadata$platform_id,
    series_id = metadata$series_id,
    source_name_ch1 = metadata$source_name_ch1,
    source_name_ch2 = metadata$source_name_ch2,
    title = metadata$title,
    type = metadata$type
  ))
  expect_equal(result@.Data, expected, ignore_attr = TRUE)
  
  # Adding same GSM does not result in additional entries
  addGSM(test)
  expect_equal(length(opeRend::listEntities("GEOSample", list(geo_accession = test))), 1)
  
  # Test updating variables
  expected2 <- expected
  expected2[["last_update_date"]] <- processDate("Jan 01 1990")
  expected2[["description"]] <- "test"
  
  result2 <- opeRend::updateEntity(testId, variables = expected2)
  expect_equal(result2@.Data, expected2, ignore_attr = TRUE)
  
  # Test updating Operend GSM entry
  result3 <- addGSM(test)
  expect_equal(result3@.Data, expected, ignore_attr = TRUE)
})

test_that("addGSE uploads and updates properly", {
  # Setting up GPL and GSMs
  gpl <- opeRend::objectId(addGPL("GPL9854"))
  gsm1 <- opeRend::objectId(addGSM("GSM764595"))
  gsm2 <- opeRend::objectId(addGSM("GSM764596"))
  
  test <- "GSE30827"
  
  # Find another GPL to test if already in Operend
  if(!is.null(queryOperend("GEOSeries", list(geo_accession = test)))) {
    stop("Test GSE already in Operend, change to another GSE\n")
  }
  
  # Check if GSM is added to Operend properly
  result <- addGSE(test)
  testId <- opeRend::objectId(result)
  gsmIds <- result$geoSamples
  
  # Clean up after testing
  on.exit({
    lapply(c(testId, gsm1, gsm2, gpl), opeRend::deleteEntity)
  })
  
  metadata <- GEOquery::Meta(GEOquery::getGEO(test, GSEMatrix = FALSE))
  expected <- removeNull(list(
    affymetrixCELSet = NULL,
    geoPlatforms = gpl,
    geoSamples = c(gsm1, gsm2),
    
    geo_accession = metadata$geo_accession,
    last_update_date = processDate(metadata$last_update_date),
    platform_id = metadata$platform_id,
    pubmed_id = as.numeric(metadata$pubmed_id),
    relation = metadata$relation,
    summary = metadata$summary,
    title = metadata$title,
    type = metadata$type
  ))
  # print(expected)
  expect_equal(result@.Data, expected, ignore_attr = TRUE)
  
  # Adding same GSM does not result in additional entries
  # addGSE(test)
  # expect_equal(length(opeRend::listEntities("GEOSeries", list(geo_accession = test))), 1)
  
  # Test updating variables
  expected2 <- expected
  expected2[["last_update_date"]] <- processDate("Jan 01 1990")
  # expected2[["summary"]] <- "test"

  result2 <- opeRend::updateEntity(testId, variables = expected2)
  expect_equal(result2@.Data, expected2, ignore_attr = TRUE)

  # Test updating Operend GPL entry

  result3 <- addGSE(test)
  expect_equal(result3@.Data, expected, ignore_attr = TRUE)
})

test_that("addAffymetrixCEL uploads and updates properly", {
  ftpUrl <- "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM603nnn/GSM603984/suppl/GSM603984_CR379.CEL.gz"
  result <- addAffymetrixCEL(ftpUrl)
  on.exit({
    # updateWorkFileProperties(result$workFile, isTrashed = TRUE)
    deleteEntity(objectId(result))
  })
  
  # Check if affymetrixCEL is added to Operend properly
  expected <- list(
    ScanDate = "11/04/08 11:43:22",
    cdfName = "HG-U133A",
    dimensions = c(712, 712),
    filenames = "GSM603984_CR379.CEL.gz",
    workFile = result$workFile
  )
  expect_equal(result@.Data, expected, ignore_attr = TRUE)
  
  # Check if hashes match between uploaded and original CEL file
  ftpUrl <- "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM603nnn/GSM603984/suppl/GSM603984_CR379.CEL.gz"
  filenames <- basename(ftpUrl)
  celFile <- file.path(tempdir(), filenames)
  on.exit(unlink(celFile), add = TRUE)
  download.file(ftpUrl, destfile = celFile)
  result2 <- tools::md5sum(celFile)

  workFile1 <- getWorkFileProperties(result$workFile)
  expect_equal(result2, workFile1@hash, ignore_attr = TRUE)
  
  # Check that comparing CEL file with uploaded entity returns uploaded entity
  result3 <- addAffymetrixCEL(ftpUrl, result)
  expect_equal(result3@.Data, expected, ignore_attr = TRUE)
  
  # Supplying an existing entity and a CEL file with a different hash modifies the existing entity
  ftpUrl2 <- "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM603nnn/GSM603985/suppl/GSM603985_CR446B-1.CEL.gz"
  result4 <- addAffymetrixCEL(ftpUrl2, result)
  on.exit(
    updateWorkFileProperties(result4$workFile, isTrashed = TRUE),
    add = TRUE
  )
  expected <- list(
    ScanDate = "11/04/08 11:49:02",
    cdfName = "HG-U133A",
    dimensions = c(712, 712),
    filenames = "GSM603985_CR446B-1.CEL.gz",
    workFile = result4$workFile
  )
  expect_equal(result4@.Data, expected, ignore_attr = TRUE)
  
  # Previous workfile is trashed
  workFile1 <- getWorkFileProperties(result$workFile)
  expect_equal(workFile1@isTrashed, TRUE)
  
  # New workFile has proper hash
  filenames2 <- basename(ftpUrl2)
  celFile2 <- file.path(tempdir(), filenames2)
  on.exit(unlink(celFile2), add = TRUE)
  download.file(ftpUrl2, destfile = celFile2)
  result5 <- tools::md5sum(celFile2)
  
  workFile2 <- getWorkFileProperties(result4$workFile)
  expect_equal(result5, workFile2@hash, ignore_attr = TRUE)
})

test_that("addGSM with CEL files uploads and updates properly", {
  test <- "GSM603984"
  
  # Find another GPL to test if already in Operend
  if(!is.null(queryOperend("GEOSample", list(geo_accession = test)))) {
    stop("Test GSM already in Operend, change to another GSM\n")
  }
  
  # Check if GSM is added to Operend properly
  result <- addGSM(test)
  testId <- opeRend::objectId(result)
  affy1 <- result$affymetrixCEL
  workFile1 <- getEntity(affy1)$workFile
  
  # Clean up after testing
  on.exit({
    opeRend::deleteEntity(testId)
    # updateWorkFileProperties(opeRend::getEntity(affy1)$workFile, isTrashed = TRUE)
    # opeRend::deleteEntity(affy1)
  })
  
  metadata <- GEOquery::Meta(GEOquery::getGEO(test, GSEMatrix = FALSE))
  expected <- removeNull(list(
    affymetrixCEL = result$affymetrixCEL,
    
    channel_count = as.numeric(metadata$channel_count),
    characteristics_ch1 = metadata$characteristics_ch1,
    characteristics_ch2 = metadata$characteristics_ch2,
    data_processing = metadata$data_processing,
    description = metadata$description,
    geo_accession = metadata$geo_accession,
    last_update_date = processDate(metadata$last_update_date),
    molecule_ch1 = metadata$molecule_ch1,
    molecule_ch2 = metadata$molecule_ch2,
    platform_id = metadata$platform_id,
    series_id = metadata$series_id,
    source_name_ch1 = metadata$source_name_ch1,
    source_name_ch2 = metadata$source_name_ch2,
    title = metadata$title,
    type = metadata$type
  ))
  expect_equal(result@.Data, expected, ignore_attr = TRUE)
  
  # Check if affymetrixCEL is added to Operend properly
  result2 <- getEntity(affy1)
  expected2 <- list(
    ScanDate = "11/04/08 11:43:22",
    cdfName = "HG-U133A",
    dimensions = c(712, 712),
    filenames = "GSM603984_CR379.CEL.gz",
    workFile = result2$workFile
  )
  expect_equal(result2@.Data, expected2, ignore_attr = TRUE)
  
  # Changes to metadata only does not affect affymetrixCEL
  expected3 <- expected
  expected3[["last_update_date"]] <- processDate("Jan 01 1990")
  expected3[["description"]] <- "test"
  
  opeRend::updateEntity(testId, variables = expected3)
  result3 <- addGSM(test)
  expect_equal(result3@.Data, expected, ignore_attr = TRUE)
  
  # Removing CEL file in GEO removes CEL file
  testGEO <- GEOquery::getGEO(test, GSEMatrix = FALSE)
  testGEO@header$supplementary_file <- NULL
  
  expected4 <- expected
  expected4$affymetrixCEL <- NULL
  expected4 <- removeNull(expected4)
  
  result4 <- addGSM(testGEO)
  
  # affymetrixCEL has been deleted
  expect_equal(result4@.Data, expected4, ignore_attr = TRUE)
  expect_error(opeRend::getEntity(affy1), "No visible entity with id")
  
  # Previous workfile is trashed
  workFile <- getWorkFileProperties(workFile1)
  expect_equal(workFile@isTrashed, TRUE)
  
  # Adding CEL file in GEO adds CEL file
  testGEO@header$supplementary_file <- "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM603nnn/GSM603985/suppl/GSM603985_CR446B-1.CEL.gz"
  result5 <- addGSM(testGEO)
  affy2 <- result5$affymetrixCEL
  workFile2 <- getEntity(affy2)$workFile
  
  on.exit({
    # opeRend::deleteEntity(testId)
    # updateWorkFileProperties(workFile2, isTrashed = TRUE)
    opeRend::deleteEntity(affy2)
  }, add = TRUE)
  
  
  expected5 <- expected
  expected5$affymetrixCEL <- affy2
  expect_equal(result5@.Data, expected5, ignore_attr = TRUE)
  
  # Changing CEL file in GEO changes CEL file
  result6 <- addGSM("GSM603984")
  workFile3 <- getEntity(affy2)$workFile
  
  on.exit({
    # opeRend::deleteEntity(testId)
    updateWorkFileProperties(workFile3, isTrashed = TRUE)
    # opeRend::deleteEntity(affy2)
  }, add = TRUE)
  
  expected6 <- expected5
  expect_equal(result6@.Data, expected6, ignore_attr = TRUE)
  
  # Check if hashes match between uploaded and original CEL file
  ftpUrl <- "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM603nnn/GSM603984/suppl/GSM603984_CR379.CEL.gz"
  filenames <- basename(ftpUrl)
  celFile <- file.path(tempdir(), filenames)
  on.exit(unlink(celFile), add = TRUE)
  download.file(ftpUrl, destfile = celFile)
  result7 <- tools::md5sum(celFile)
  
  # Hash from initial upload
  expect_equal(result7, getWorkFileProperties(workFile1)@hash, ignore_attr = TRUE)
  
  # New workFile has proper hash
  ftpUrl2 <- "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM603nnn/GSM603985/suppl/GSM603985_CR446B-1.CEL.gz"
  filenames2 <- basename(ftpUrl2)
  celFile2 <- file.path(tempdir(), filenames2)
  on.exit(unlink(celFile2), add = TRUE)
  download.file(ftpUrl2, destfile = celFile2)
  result8 <- tools::md5sum(celFile2)
  
  # Hash from newly added upload
  expect_equal(result8, getWorkFileProperties(workFile2)@hash, ignore_attr = TRUE)
  
  # Hash after supplementary file link change
  expect_equal(result7, getWorkFileProperties(workFile3)@hash, ignore_attr = TRUE)
})
