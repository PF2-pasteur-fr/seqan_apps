
for fp in tophat/*/accepted_hits.bam ; do
   tmp=${fp%/accepted_hits.bam}
   project=${tmp#tophat/*}
   if [ ! -d rois/$project ] ;then
      mkdir rois/$project
   fi
   echo "bam2roi for $project"
   echo time /pasteur/projets/specific/PF2_ngs/protected/programs/seqan/seqan-trunk-build/Release/bin/bam2roi \
      -if tophat/$project/accepted_hits.bam \
      -of rois/$project/accepted_hits.roi \
      -ss -ls > qsubs/$project.qsub
   echo gzip -f9 rois/$project/accepted_hits.roi  >> qsubs/$project.qsub
   echo time /pasteur/projets/specific/PF2_ngs/protected/scripts/roiTransforms.sh rois/$project/accepted_hits.roi.gz  Rtransforms/${project}.Rtransforms.roi  >> qsubs/$project.qsub
   echo gzip -f9 Rtransforms/${project}.Rtransforms.roi >> qsubs/$project.qsub
   echo /pasteur/projets/specific/PF2_ngs/protected/scripts/roiPCA.sh Rtransforms/${project}.Rtransforms.roi.gz Rtransforms/${project}.pca.roi >> qsubs/$project.qsub
   echo gzip -9f Rtransforms/${project}.pca.roi >> qsubs/$project.qsub
done


