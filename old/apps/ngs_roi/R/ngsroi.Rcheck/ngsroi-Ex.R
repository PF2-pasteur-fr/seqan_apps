pkgname <- "ngsroi"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('ngsroi')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("ngsroi-package")
### * ngsroi-package

flush(stderr()); flush(stdout())

### Name: ngsroi-package
### Title: NGS Regions of Interest Analysis
### Aliases: ngsroi-package ngsroi
### Keywords: package

### ** Examples

library(ngsroi)

# Load ROI file into data.frame.
#roi = readROI("dmel.bowtie.sam.roi")

# Compute some metrics.
#roi$min = as.numeric(lapply(roi$counts ,min))
#roi$median =  as.numeric(lapply(roi$counts ,median))
#roi$mean =  as.numeric(lapply(roi$counts ,mean))
#roi$quantile75 =  as.numeric(lapply(roi$counts ,quantile, probs=0.75))
#roi$quantile95 =  as.numeric(lapply(roi$counts ,quantile, probs=0.95))

# Write data.frame into ROI file again.
#writeROI(roi, "dmel.bowtie.sam.trans.roi");



cleanEx()
nameEx("readROI")
### * readROI

flush(stderr()); flush(stdout())

### Name: readROI
### Title: Read ROI file.
### Aliases: readROI
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (fpname) 
{
    fi <- function(x, i) {
        x[i]
    }
    fni <- function(x, i) {
        as.numeric(x[[i]])
    }
    getVec <- function(x) {
        mylen = as.numeric(x[5])
        veclen = length(x)
        as.numeric(unlist(strsplit(x[veclen], ",")))
    }
    getVals <- function(y, x, columnNames) {
        veclen = length(columnNames)
        t1 = veclen - 1
        for (colN in c(7:t1)) {
            y[, columnNames[colN]] = unlist(lapply(x, fni, colN))
        }
        return(y)
    }
    con = gzfile(fpname)
    rLines = readLines(con)
    close(con)
    values = rLines[substr(rLines, 1, 1) != "#"]
    values = strsplit(values, "\t")
    if (length(rLines[substr(rLines, 1, 2) == "##"]) == 0) {
        columnNames = unlist(list("##ref", "begin_pos", "end_pos", 
            "region_name", "length", "strand", "max_count", "cg_content", 
            "counts"))
    }
    else {
        columnNames = unlist(strsplit(rLines[substr(rLines, 1, 
            2) == "##"], "\t"))
    }
    df = data.frame(ref = unlist(lapply(values, fi, 1)), begin_pos = as.integer(unlist(lapply(values, 
        fni, 2))), end_pos = as.integer(unlist(lapply(values, 
        fni, 3))), region_name = unlist(lapply(values, fi, 4)), 
        length = as.integer(unlist(lapply(values, fni, 5))), 
        strand = unlist(lapply(values, fi, 6)))
    df = getVals(df, values, columnNames)
    df$counts = lapply(values, getVec)
    df$counts = lapply(df$counts, unlist)
    roiNames = names(df)
    return(df)
  }



cleanEx()
nameEx("writeROI")
### * writeROI

flush(stderr()); flush(stdout())

### Name: writeROI
### Title: Write ROI file.
### Aliases: writeROI
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (roi, fpname) 
{
    fpConn <- file(fpname, "w")
    colCount = 1
    write("# ROI written from R", fpConn, append = F)
    colNames = names(roi)
    colNames = colNames[!colNames == "counts"]
    outStr = ""
    outStr = paste(c(outStr, "##", colNames[1]), collapse = "")
    for (colIds in c(2:length(colNames))) {
        outStr = paste(c(outStr, colNames[colIds]), collapse = "\t")
    }
    outStr = paste(c(outStr, "counts"), collapse = "\t")
    write(outStr, file = fpConn, append = TRUE)
    apply(roi, 1, writeRoiLine, fpConn)
    close(fpConn)
  }



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
