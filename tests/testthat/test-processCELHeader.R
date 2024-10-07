library(testthat)
library(opeRend.NCBI)
library(affyio)

test_that("processCELHeader handles various cases", {
  # Check if CEL files are successfully read
  ftpUrl <- "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM603nnn/GSM603984/suppl/GSM603984_CR379.CEL.gz"
  filenames <- basename(ftpUrl)
  celFile <- file.path(tempdir(), filenames)
  on.exit(unlink(celFile))
  download.file(ftpUrl, destfile = celFile)
  result <- processCELHeader(celFile)

  expected <- list(
    cdfName = "HG-U133A",
    dimensions = c(712, 712),
    DatHeader = paste(
      c("[0..65534]  CR379:CLS=5464 RWS=5464 XIN=2  YIN=2  VE=30        2.0",
        " 11/04/08 11:43:22 50201191  M10   \024  \024 HG-U133A.1sq \024  ",
        "\024  \024  \024  \024 570 \024 25348.910156 \024 3.500000 \024 ",
        "2.5000 \024 6"
      ),
      collapse = ""
    ),
    ScanDate = "11/04/08 11:43:22"
  )
  expect_equal(result, expected)
  
  expect_error(processCELHeader(celFile2))
  
  # Check if non-CEL files cause error
  nonCelLink <- "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM227nnn/GSM227858/suppl/GSM227858.EXP.gz"
  filenames <- basename(nonCelLink)
  celFile2 <- file.path(tempdir(), filenames)
  on.exit(unlink(celFile2))
  download.file(nonCelLink, destfile = celFile2)
  expect_error(processCELHeader(celFile2))
  
  # Passing string instead of file as argument causes error
  expect_error(processCELHeader(ftpUrl))
})




