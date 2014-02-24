#! /usr/bin/Rscript

args <- commandArgs(TRUE)
inFile=args[1]
outFile=args[2]
source('/pasteur/projets/specific/PF2_ngs/protected/scripts/Rnw/rfunctions/readROI.R')
source('/pasteur/projets/specific/PF2_ngs/protected/scripts/Rnw/rfunctions/writeROI.R')
print(inFile)
print(outFile)

roi=readROI(inFile)

data=roi[,8:length(roi)-1]
pca=princomp(data, na.action=na.exclude)
prj=predict(pca, data)
roi$pca1=prj[,1]
roi$pca2=prj[,2]

writeROI(roi, outFile)

q(status=0)

