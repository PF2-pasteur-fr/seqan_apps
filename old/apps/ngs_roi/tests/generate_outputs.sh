#!/bin/sh
#
# Output generation for roi_intersect

INTERSECT="../../../../../../seqan-pasteur-build/debug/bin/roi_intersect"

# ----------------------------------------------------------------------------
# Intersection with BED or GFF/GTF in BED style.
# ----------------------------------------------------------------------------

for mode in intersection projection union difference; do
    for format in bed gff gtf; do
        echo "${INTERSECT} --mode ${mode} -ir small.roi -if small.${format} -or out_small_${format}_m${mode}.roi >out_small_${format}_m${mode}.stdout 2>out_small_${format}_m${mode}.stderr"
        ${INTERSECT} --mode ${mode} -ir small.roi -if small.${format} -or out_small_${format}_m${mode}.roi >out_small_${format}_m${mode}.stdout 2>out_small_${format}_m${mode}.stderr
        echo "${INTERSECT} --strand-specific --mode ${mode} -ir small.roi -if small.${format} -or out_small_${format}_m${mode}_ss.roi >out_small_${format}_m${mode}_ss.stdout 2>out_small_${format}_m${mode}_ss.stderr"
        ${INTERSECT} --strand-specific --mode ${mode} -ir small.roi -if small.${format} -or out_small_${format}_m${mode}_ss.roi >out_small_${format}_m${mode}_ss.stdout 2>out_small_${format}_m${mode}_ss.stderr
    done
done
