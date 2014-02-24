#! /usr/bin/Rscript
library('kohonen')

args <- commandArgs(TRUE)
inFile=args[1]
outFile=args[2]
source('/pasteur/projets/specific/PF2_ngs/protected/scripts/Rnw/rfunctions/readROI.R')
source('/pasteur/projets/specific/PF2_ngs/protected/scripts/Rnw/rfunctions/writeROI.R')
print(inFile)
print(outFile)

roi=readROI(inFile)

dataFound=roi$max_count>0
data=roi[dataFound,7:(length(roi)-1)]
roi$kohonenClass=-1
if(dim(data)[[1]]<50){
   roi$kohonenClass=-1
}else{
   data.sc = scale(data)
   data.som = som(data=data.sc, rlen=10000, grid = somgrid(5,5,"hexagonal"))
   roi[dataFound,]$kohonenClass=data.som$unit.classif

   png(file=paste(outFile,".png", sep=""))
   par(mfrow = c(2, 1))
   plot(data.som, main = paste(inFile, "data"))
   plot(data.som, type="changes", main="Learning curve")
   dev.off()
   save(data.som, file=paste(outFile,".rData",sep=""))
}
writeROI(roi, outFile)

q(status=0)

