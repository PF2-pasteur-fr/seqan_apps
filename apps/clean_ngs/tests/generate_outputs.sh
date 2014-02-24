#!/bin/sh
#
# Output generation for the adaptor removal tool.

REMOVER="../../../../../../seqan-pasteur-build/release/bin/clean_ngs"

# ----------------------------------------------------------------------------
# Single-End Adapter Removal
# ----------------------------------------------------------------------------

# Default settings.
${REMOVER} -adf adapters.txt -if1 reads_1.fq -of1 reads_1_out_se.fq -rf1 reads_1_out_se_rej.fq >reads_1_out_se_stdout.txt 2>reads_1_out_se_stderr.txt

# Setting minimal length.
${REMOVER} -adf adapters.txt -if1 reads_1.fq -of1 reads_1_out_se_l90.fq -rf1 reads_1_out_se_l90_rej.fq >reads_1_out_se_l90_stdout.txt 2>reads_1_out_se_l90_stderr.txt -l 90

# ----------------------------------------------------------------------------
# Paired-End Adaptor Removal
# ----------------------------------------------------------------------------

${REMOVER} -adf adapters.txt -if1 reads_1.fq -if2 reads_2.fq -of1 reads_1_out_pe.fq -of2 reads_2_out_pe.fq -rf1 reads_1_out_pe_rej.fq -rf2 reads_2_out_pe_rej.fq >reads_1_out_pe_stdout.txt 2>reads_1_out_pe_stderr.txt
