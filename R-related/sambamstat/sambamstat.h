// ==========================================================================
//                                 sambamstat
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Your Name <your.email@example.net>
// ==========================================================================
// Version 0.1 = NH record
// Version 0.7 = replace BamSamheader instead of appending to it

#ifndef SANDBOX_JAGLA_APPS_SAMBAMSTAT_SAMBAMSTAT_H_
#define SANDBOX_JAGLA_APPS_SAMBAMSTAT_SAMBAMSTAT_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <sstream>
#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/stream/stream_bgzf.h>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <seqan/bam_io/read_bam.h>
#include <seqan/bam_io/write_bam.h>
#include <seqan/bam_io/write_sam.h>


#include <seqan/misc/misc_cmdparser.h>
#include <seqan/file.h>      // For printing SeqAn Strings.
#include <seqan/misc/misc_cmdparser.h>

//#include "dataanalysis.h"

#define VERSION "0.14"
#define PROGNAME_ID "sambamstats-ID"
#define PROGNAME_POS "sambamstats-POS"

//#if SEQAN_HAS_ZLIB

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================
enum Format {
	FORMAT_AUTO, FORMAT_SAM, FORMAT_BAM
};

enum SortOrder {
	SORTORDER_NA, SORTORDER_ID, SORTORDER_POS
};

struct ROI {
	int min; // start position
	int max; // end position
	int rId;
	int roiID;
	String<unsigned> count;
	unsigned len;
	unsigned countMax;
	ROI() {
		min = INT_MAX;
		max = 0;
		rId = -2;
		len = 0;
		roiID = 0;
		countMax = 0;
	}
};

struct Stats {
	__uint64 numRecords;
	__uint64 alignedRecords;
	String<double> avrgQuality;
	String<unsigned> nhHisto;
	String<unsigned> intronLHisto;
	String<unsigned> exonLHisto;
	String<unsigned> intronCHisto;
	String<unsigned> lengthsHisto;
	String<unsigned> cigarHisto;
	String<unsigned> flagsHisto;
	String<unsigned> rIDs;
	String<unsigned> insertHisto;
	String<unsigned> deletionHisto;
	String<unsigned> mismatchHisto;
	String<unsigned> editDistanceHisto;
	String<unsigned> roicountMax;
	String<unsigned> roiLenHisto;
	Stats() :
			numRecords(0), alignedRecords(0) {
	}
};

struct Options {
	bool showHelp;
	bool showVersion;
	unsigned int verbosity;
	CharString inFileName;
	CharString outFileName;
	CharString refFileName;
	CharString roiFileName;
	CharString tmpDir;
//	bool statistics;
	Format inFormat;
	SortOrder sortOrder;
	bool forceRedo;
	CharString cmd;
	unsigned maxCoverage;
	bool strandSpecific;

	Options() {
		// Set defaults.
		showHelp = false;
		showVersion = false;
		verbosity = 0;
		//statistics = false;
		forceRedo = false;
		inFileName = "";
		outFileName = "";
		refFileName = "";
		roiFileName = "";
		maxCoverage = 100000000;
		//statistics = false;
		inFormat = FORMAT_AUTO;
		sortOrder = SORTORDER_NA;
		tmpDir = "";
		strandSpecific = false;
	}
};

double programStartTime=0.0;
int roiCount = 0;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

void setupCommandLineParser(CommandLineParser & parser,
		Options const & options) {
	addVersionLine(parser, VERSION);

	addTitleLine(parser, "**********************");
	addTitleLine(parser, "* sambamstat *");
	addTitleLine(parser, "**********************");
	addTitleLine(parser, "");
	addTitleLine(parser, "(c) 2012 by Bernd Jagla <Bernd.Jagla@pasteur.fr>");
	addTitleLine(parser,
			"this program makes use of the env variable TMPDIR [/tmp]");

	addUsageLine(parser, "[OPTIONS] -if inputFile -of outputFile");

	addSection(parser, "Main Options");
	addOption(parser,
			CommandLineOption("if", "input-file", "Input file.",
					OptionType::String | OptionType::Label,
					options.inFileName));
	addOption(parser,
			CommandLineOption("ref", "ref_file", "Reference Genome fasta",
					OptionType::String | OptionType::Label,
					options.refFileName));
	addOption(parser,
			CommandLineOption("S", "input-sam",
					"Input file is SAM (default: auto).", OptionType::Bool,
					options.inFormat == FORMAT_SAM));
	addOption(parser,
			CommandLineOption("B", "input-bam",
					"Input file is BAM (default: auto).", OptionType::Bool,
					options.inFormat == FORMAT_BAM));

	addOption(parser,
			CommandLineOption("sP", "sortOrderPos", "assume sorted by position",
					OptionType::Bool, options.sortOrder == SORTORDER_POS));
	addOption(parser,
			CommandLineOption("sI", "sortOrderID", "assume sorted by ID",
					OptionType::Bool, options.sortOrder == SORTORDER_ID));

	addOption(parser,
			CommandLineOption("of", "output-file", "Output file.",
					OptionType::String | OptionType::Label,
					options.outFileName));

	addOption(parser,
			CommandLineOption("rf", "roi-output-file",
					"Output file for region of interest.",
					OptionType::String | OptionType::Label,
					options.roiFileName));

	addSection(parser, "optional parameters:");
//	addOption(parser,
//			CommandLineOption("s", "stats", "Perform stats", OptionType::Bool,
//					options.statistics));
	addOption(parser,
			CommandLineOption("ss", "strandspecific",
					"calculate strandspecific stats (only position sorted)",
					OptionType::Bool));
	addOption(parser,
			CommandLineOption("f", "forceRedo", "enforce recalculation",
					OptionType::Bool, options.forceRedo));
	addOption(parser,
			CommandLineOption("T", "tempdirectory",
					"use DIR for temporaries, not /tmp; multiple options specify multiple directories",
					OptionType::String, options.tmpDir));
	addOption(parser,
			CommandLineOption("mC", "maxCoverage",
					"maximum number of reads per ROI [100,000,000]",
					OptionType::Integer, options.maxCoverage));
	addSection(parser, "General Options");
	addOption(parser,
			CommandLineOption("vv", "very-verbose", "Very verbose output.",
					OptionType::Bool));
	addOption(parser,
			CommandLineOption("v", "verbose", "verbose output.",
					OptionType::Bool));

	requiredArguments(parser, 0);
}

int parseCommandLineAndCheck(Options & options, CommandLineParser & parser,
		int argc, char const ** argv) {
	bool stop = !parse(parser, argc, argv);
	if (stop)
		return 1;
	if (isSetLong(parser, "help")) {
		options.showHelp = true;
		return 0;
	}
	if (isSetLong(parser, "version")) {
		options.showVersion = true;
		return 0;
	}

	getOptionValueLong(parser, "input-file", options.inFileName);
	getOptionValueLong(parser, "output-file", options.outFileName);
	getOptionValueLong(parser, "ref_file", options.refFileName);
	getOptionValueLong(parser, "roi-output-file", options.roiFileName);
	getOptionValueLong(parser, "maxCoverage", options.maxCoverage);
	getOptionValueLong(parser, "tempdirectory", options.tmpDir);
	if (isSetLong(parser, "verbose"))
		options.verbosity = 1;
	if (isSetLong(parser, "very-verbose"))
		options.verbosity = 2;
//	if (isSetLong(parser, "stats"))
//		options.statistics = true;
	if (isSetLong(parser, "sortOrderID"))
		options.sortOrder = SORTORDER_ID;
	if (isSetLong(parser, "sortOrderPos"))
		options.sortOrder = SORTORDER_POS;
	if (isSetLong(parser, "forceRedo"))
		options.forceRedo = true;
	if (isSetLong(parser, "strandspecific"))
		options.strandSpecific = true;
	for (int i = 0; i < argc; i++) {
		options.cmd += argv[i];
		options.cmd += " ";
	}

	return 0;
}

// record is a BAM record of a possibly split mapped read, sets tags "il" to number of deleted characters from
// CIGAR string, set "cl" to number of covered bases in genome, including intron.

void intronLength(BamAlignmentRecord & record, Stats & stats,
		BamTagsDict & tagDict) {

	String<CigarElement<> > cigar = record.cigar;
	if (length(cigar) > 0) {
		int anzN = 0;
		unsigned coverageLength = 0;
		for (unsigned long i = 0; i < length(cigar); i++) {
			unsigned cLen = length(stats.cigarHisto);
			char op = ((CigarElement<> ) cigar[i]).operation;
			resize(stats.cigarHisto, std::max(cLen, (unsigned) op + 1), 0);
			stats.cigarHisto[((CigarElement<> ) cigar[i]).operation]++;
			coverageLength += ((CigarElement<> ) cigar[i]).count;
			if (((CigarElement<> ) cigar[i]).operation != 'M') {
				unsigned len = length(stats.intronLHisto);
				resize(stats.intronLHisto,
						std::max(len, ((CigarElement<> ) cigar[i]).count + 1),
						0);
				anzN += ((CigarElement<> ) cigar[i]).count;
				stats.intronLHisto[((CigarElement<> ) cigar[i]).count]++;
			} else {
				unsigned len = length(stats.exonLHisto);
				resize(stats.exonLHisto,
						std::max(len, ((CigarElement<> ) cigar[i]).count + 1),
						0);
				stats.exonLHisto[((CigarElement<> ) cigar[i]).count]++;
			}
		}
		setTagValue(tagDict, "il", anzN, 'i'); //intron length
		//setTagValue(tagDict, "nl", record.tLen - anzN, 'i'); //inner length target length - introns for this read
		setTagValue(tagDict, "cl", coverageLength, 'i');
	}
}

// Set tag "ic" to the number of introns in a record.

void intronCount(BamAlignmentRecord & record, Stats & stats,
		BamTagsDict & tagDict) {

	String<CigarElement<> > cigar = record.cigar;
	if (length(cigar) > 0) {
		unsigned anzN = 0;
		for (unsigned long i = 0; i < length(cigar); i++) {
			if (((CigarElement<> ) cigar[i]).operation != 'M')
				anzN++;
		}
		setTagValue(tagDict, "ic", anzN, 'i');
		unsigned len = length(stats.intronCHisto);
		resize(stats.intronCHisto, std::max(len, anzN + 1), 0);
		stats.intronCHisto[anzN]++;
	}
}

// Compute average quality value of quality string.

void averageQV(BamAlignmentRecord & record, Stats & stats,
		BamTagsDict & tagDict) {

	CharString scores = record.qual;

	if (length(scores) > 0) {
		unsigned len = length(stats.avrgQuality);
		resize(stats.avrgQuality,
				std::max(len, unsigned(length(record.seq)) + 1), 0);
		unsigned long sum = 0;
		for (unsigned long i = 0; i < length(scores); i++) {
			sum += scores[i];
			stats.avrgQuality[i + 1] += scores[i];
		}
		setTagValue(tagDict, "aq", float(int(sum)) / float(length(scores)),
				'f');
	}
}

//TODO rename is not working for me....
int myrename(char * oldfp, char * newfp) {
	char cmd[2000];
#if defined(_WIN32)
	sprintf(cmd, "move %s %s", oldfp, newfp);
#else
	sprintf(cmd, "mv %s %s", oldfp, newfp);
#endif
	return system(cmd);
}

// Update the length histogram in stats given the record.

void lengthHist(Stats& stats, BamAlignmentRecord& record) {
	//Lengths histogram
	unsigned len = length(stats.lengthsHisto);
	unsigned seqLen = length(record.seq);
	resize(stats.lengthsHisto, std::max(len, seqLen + 1), 0);
	for (unsigned i = 1; i <= seqLen; i++) {
		stats.lengthsHisto[i] += 1;
	}
}

// Compute alignment statistics from the record and the contigs.  This does currently not work.

void fastaStats(const Options& options, BamAlignmentRecord& record,
		unsigned maxLen, Stats& stats, StringSet<Dna5String> faSeqs,
		Align<Dna5String> align, String<__uint64> qualSum, __int64& reads)
		{
			//stats with fasta reference genome...
		if (options.refFileName != "") {
			if (length(record.seq) > maxLen) {
				maxLen = length(record.seq);
				resize(stats.insertHisto, maxLen + 1, 0);
				resize(stats.deletionHisto, maxLen + 1, 0);
				resize(stats.mismatchHisto, maxLen + 1, 0);
				resize(stats.editDistanceHisto, maxLen + 1, 0);
			}

			if (record.rId != BamAlignmentRecord::INVALID_REFID && record.rId < length(faSeqs)) {
				//std::cerr << header.sequenceInfos[record.rId] << "\t"  << record.seq << std::endl;
				bamRecordToAlignment(align, faSeqs[record.rId], record);
				if (options.verbosity >= 2)
				std::cerr << align << std::endl;

				typedef Align<Dna5String> TAlign;
				typedef Row<TAlign>::Type TRow;
				typedef Iterator<TRow>::Type TRowIter;
				unsigned editDistance = 0;
				unsigned posReadFwd = 0;
				for (TRowIter it0 = begin(row(align, 0)), it1 = begin(row(align, 1));!atEnd(it0) && !atEnd(it1);goNext(it0), goNext(it1)) {
					unsigned posRead = posReadFwd;
					// is read aligned to reverse strand?
					if ((record.flag & 0x10) != 0)
					posRead = (length(record.seq) - 1) - posReadFwd;

					if (posRead > maxLen) {
						maxLen = length(record.seq);
						resize(stats.insertHisto, maxLen + 1, 0);
						resize(stats.deletionHisto, maxLen + 1, 0);
						resize(stats.mismatchHisto, maxLen + 1, 0);
					}
					if (isGap(it0) && isGap(it1))
					continue;

					if (isGap(it0) || isGap(it1)) {
						if (isGap(it0)) {
							stats.insertHisto[posRead] += 1;
							posReadFwd += 1;
						} else {
							stats.deletionHisto[posRead] += 1;
						}
						editDistance += 1;
						continue;
					}

					if (*it0 != *it1) {
						stats.mismatchHisto[posRead] += 1;
						editDistance += 1;
					}
					posReadFwd += 1;
				}

				resize(qualSum, std::max(length(qualSum), length(record.qual)), 0);
				if ((record.flag & 0x10) == 0) {
					// read aligns to forward strand
					for (unsigned i = 0;i < length(record.qual);++i)
					qualSum[i] += record.qual[i] - '!';
				} else {
					// read aligns to reverse strand
					for (unsigned i = 0;i < length(record.qual);++i)
					qualSum[(length(record.qual) - 1) - i] += record.qual[i] - '!';
				}
				++reads;
				if (options.verbosity >= 2)
				std::cerr << "edit distance: " << editDistance << std::endl;

				unsigned len = length(stats.editDistanceHisto);
				resize(stats.editDistanceHisto, std::max(len, editDistance + 1), 0);
				stats.editDistanceHisto[editDistance] += 1;
			}
		}

	}

// Create temporary BAM file.

void createTmpFile(char * tmpFPnam, Options const & options, FILE ** sFile,
		Stream<Bgzf> & bFile) {
	char * tmpDir = NULL;
	char defaultTMPDIR[4095] = "/tmp";
	tmpDir = getenv("TMPDIR");
	if (tmpDir == NULL) {
		tmpDir = defaultTMPDIR;
	}
	if (options.tmpDir != "") {
		tmpDir = toCString(options.tmpDir);
		std::cerr << "tmpdir  " << tmpDir << " is set" << std::endl;
	}
	std::cerr << "this is options.tempDir " << options.tmpDir << " is set"
			<< std::endl;
	tmpFPnam = strcpy(tmpFPnam, tmpDir);
	tmpFPnam = strcat(tmpFPnam, "/samBamStatXXXXXX");
#if defined(_WIN32)
	tmpFPnam = strcat(toCString(options.outFileName), ".tmp");
#else
	if (mkstemp(tmpFPnam) == 0) {
		std::cerr << "Cannot create temp file" << std::endl;
		//the tmp file has not been created yet so we can exit
		exit(-1);
	}
#endif
	std::cerr << "Using " << tmpFPnam << " as temp file" << std::endl;
	if (options.verbosity >= 2)
		printf("Tempname #1: %s\n", tmpFPnam);
	if (options.inFormat == FORMAT_SAM) {
		*sFile = fopen(tmpFPnam, "w");
		if (*sFile == NULL) {
			std::cerr << "Could not create temp file" << std::endl;
			//the tmp file has not been created yet so we can exit
			exit(-1);
		}
	} else {
		if (!open(bFile, tmpFPnam, "w")) {
			std::cerr << "Could not create temp file" << std::endl;
			//the tmp file has not been created yet so we can exit
			exit(-1);
		}
	}

}

void printStatsIDsorted(Stats stats, BamHeader header, const Options& options,
		__int64 reads) {
	// Print stats
	std::cout << ">>Basic Statistics" << std::endl;
	std::cout << "#Measure\tValue" << std::endl;
	std::cout << "time\t" << sysTime() - programStartTime << "\n";
	std::cout << "num records\t" << stats.numRecords << std::endl;
	std::cout << "aligned records\t" << stats.alignedRecords << std::endl;
	std::cout << "aligned record %\t"
			<< 100.0 * stats.alignedRecords / stats.numRecords << std::endl;
	std::cout << ">>END_MODULE" << std::endl;

	std::cout << ">>NH Histogram" << std::endl;
	std::cout << "#NH value\tcount" << std::endl;
	for (unsigned i = 1; i < length(stats.nhHisto); ++i) {
		if (stats.nhHisto[i] > 0) {
			std::cout << i << '\t';
			std::cout << stats.nhHisto[i] << std::endl;
		}
	}
	std::cout << ">>END_MODULE" << std::endl;

	std::cout << ">>Flags Histogram" << std::endl;
	std::cout << "#Flag\tcount" << std::endl;
	for (unsigned i = 0; i < length(stats.flagsHisto); ++i) {
		if (stats.flagsHisto[i] > 0) {
			std::cout << i << '\t';
			std::cout << stats.flagsHisto[i] << std::endl;
		}
	}
	std::cout << ">>END_MODULE" << std::endl;

	std::cout << ">>Cigar Histogram" << std::endl;
	std::cout << "#Cigar character\tcount" << std::endl;
	for (unsigned i = 1; i < length(stats.cigarHisto); ++i) {
		if (stats.cigarHisto[i] > 0) {
			std::cout << (char) (i) << '\t';
			std::cout << stats.cigarHisto[i] << std::endl;
		}
	}
	std::cout << ">>END_MODULE" << std::endl;

	std::cout << ">>Chromosome coverage" << std::endl;
	std::cout << "#Chromosome\tChr length\treads\tPercentage" << std::endl;
	std::cout << "*\t" << stats.rIDs[0] << "\t0\t0" << std::endl;
	for (unsigned i = 1; i < length(stats.rIDs); i++) {
		std::cout << header.sequenceInfos[i - 1].i1 << "\t"
				<< header.sequenceInfos[i - 1].i2 << "\t" << stats.rIDs[i];
		std::cout << "\t"
				<< 100.0 * double(stats.rIDs[i])
						/ double(header.sequenceInfos[i - 1].i2) << std::endl;
	}
	std::cout << ">>END_MODULE" << std::endl;

	std::cout << ">>Intron count Distribution" << std::endl;
	std::cout << "#Number of introns\tcount" << std::endl;
	for (unsigned i = 0; i < length(stats.intronCHisto); ++i) {
		if (stats.intronCHisto[i] > 0) {
			std::cout << i << '\t';
			std::cout << stats.intronCHisto[i] << std::endl;
		}
	}
	std::cout << ">>END_MODULE" << std::endl;

	std::cout << ">>Intron Lengths Distribution" << std::endl;
	std::cout << "#Intron length\tcount" << std::endl;
	for (unsigned i = 0; i < length(stats.intronLHisto); ++i) {
		if (stats.intronLHisto[i] > 0) {
			std::cout << i << '\t';
			std::cout << stats.intronLHisto[i] << std::endl;
		}
	}
	std::cout << ">>END_MODULE" << std::endl;

	std::cout << ">>Exon Lengths Distribution" << std::endl;
	std::cout << "#Exon length\tcount" << std::endl;
	for (unsigned i = 0; i < length(stats.exonLHisto); ++i) {
		if (stats.exonLHisto[i] > 0) {
			std::cout << i << '\t';
			std::cout << stats.exonLHisto[i] << std::endl;
		}
	}
	std::cout << ">>END_MODULE" << std::endl;

	std::cout << ">>Per position stats" << std::endl;
	std::cout << "#position\t";
	if (options.refFileName != "") {
		std::cout
				<< "mismatches\tinsertions\tdeletions\terror prob\tquality-based error prob\t";
	}
	std::cout << "quality" << std::endl;
	for (unsigned i = 1; i < length(stats.lengthsHisto); ++i) {
		std::cout << i << '\t';
		if (options.refFileName != "") {
			std::cout << stats.mismatchHisto[i] << '\t';
			std::cout << stats.insertHisto[i] << '\t';
			std::cout << stats.deletionHisto[i] << '\t';
			std::cout
					<< (stats.mismatchHisto[i] + stats.insertHisto[i]
							+ stats.deletionHisto[i]) / (double) (reads)
					<< '\t';
		}
		double e = stats.avrgQuality[i] / stats.lengthsHisto[i];
		std::cout << e << std::endl;
	}
	std::cout << ">>END_MODULE" << std::endl;

	if (options.refFileName != "") {
		std::cout << ">>Edit distance" << std::endl;
		std::cout << "#distance\tnumber of reads" << std::endl;
		for (unsigned i = 0; i < length(stats.editDistanceHisto); ++i) {
			if (stats.editDistanceHisto[i] > 0)
				std::cout << i << '\t' << stats.editDistanceHisto[i]
						<< std::endl;
		}
		std::cout << ">>END_MODULE" << std::endl;
	}

}

void printStatsPosSorted(Stats stats, BamHeader header, const Options& options,
		__int64 reads) {
	// Print stats
	std::cout << ">>Basic Statistics" << std::endl;
	std::cout << "#Measure\tValue" << std::endl;
	std::cout << "time\t" << sysTime() - programStartTime << "\n";
	std::cout << "num records\t" << stats.numRecords << std::endl;
	std::cout << "aligned records\t" << stats.alignedRecords << std::endl;
	std::cout << "aligned record %\t"
			<< 100.0 * stats.alignedRecords / stats.numRecords << std::endl;
	std::cout << ">>END_MODULE" << std::endl;
	std::cout << ">>ROI length Histogram" << std::endl;
	std::cout << "#roi length\tcount" << std::endl;
	for (unsigned i = 0; i < length(stats.roiLenHisto); ++i) {
		if (stats.roiLenHisto[i] > 0) {
			std::cout << i << '\t';
			std::cout << stats.roiLenHisto[i] << std::endl;
		}
	}
	std::cout << ">>END_MODULE" << std::endl;

	std::cout << ">>ROI max count Histogram" << std::endl;
	std::cout << "#maxCoverage\tcount" << std::endl;
	for (unsigned i = 0; i < length(stats.roicountMax); ++i) {
		if (stats.roicountMax[i] > 0) {
			std::cout << i << '\t';
			std::cout << stats.roicountMax[i] << std::endl;
		}
	}
	std::cout << ">>END_MODULE" << std::endl;

	std::cout << ">>Flags Histogram" << std::endl;
	std::cout << "#Flag\tcount" << std::endl;
	for (unsigned i = 0; i < length(stats.flagsHisto); ++i) {
		if (stats.flagsHisto[i] > 0) {
			std::cout << i << '\t';
			std::cout << stats.flagsHisto[i] << std::endl;
		}
	}
	std::cout << ">>END_MODULE" << std::endl;

	std::cout << ">>Cigar Histogram" << std::endl;
	std::cout << "#Cigar character\tcount" << std::endl;
	for (unsigned i = 1; i < length(stats.cigarHisto); ++i) {
		if (stats.cigarHisto[i] > 0) {
			std::cout << (char) (i) << '\t';
			std::cout << stats.cigarHisto[i] << std::endl;
		}
	}
	std::cout << ">>END_MODULE" << std::endl;

	std::cout << ">>Chromosome coverage" << std::endl;
	std::cout << "#Chromosome\tChr length\treads\tPercentage" << std::endl;
	std::cout << "*\t" << stats.rIDs[0] << "\t0\t0" << std::endl;
	for (unsigned i = 1; i < length(stats.rIDs); i++) {
		std::cout << header.sequenceInfos[i - 1].i1 << "\t"
				<< header.sequenceInfos[i - 1].i2 << "\t" << stats.rIDs[i];
		std::cout << "\t"
				<< 100.0 * double(stats.rIDs[i])
						/ double(header.sequenceInfos[i - 1].i2) << std::endl;
	}
	std::cout << ">>END_MODULE" << std::endl;

	std::cout << ">>Intron count Distribution" << std::endl;
	std::cout << "#Number of introns\tcount" << std::endl;
	for (unsigned i = 0; i < length(stats.intronCHisto); ++i) {
		if (stats.intronCHisto[i] > 0) {
			std::cout << i << '\t';
			std::cout << stats.intronCHisto[i] << std::endl;
		}
	}
	std::cout << ">>END_MODULE" << std::endl;

	std::cout << ">>Intron Lengths Distribution" << std::endl;
	std::cout << "#Intron length\tcount" << std::endl;
	for (unsigned i = 0; i < length(stats.intronLHisto); ++i) {
		if (stats.intronLHisto[i] > 0) {
			std::cout << i << '\t';
			std::cout << stats.intronLHisto[i] << std::endl;
		}
	}
	std::cout << ">>END_MODULE" << std::endl;

	std::cout << ">>Exon Lengths Distribution" << std::endl;
	std::cout << "#Exon length\tcount" << std::endl;
	for (unsigned i = 0; i < length(stats.exonLHisto); ++i) {
		if (stats.exonLHisto[i] > 0) {
			std::cout << i << '\t';
			std::cout << stats.exonLHisto[i] << std::endl;
		}
	}
	std::cout << ">>END_MODULE" << std::endl;

	std::cout << ">>Per position stats" << std::endl;
	std::cout << "#position\t";
	if (options.refFileName != "") {
		std::cout
				<< "mismatches\tinsertions\tdeletions\terror prob\tquality-based error prob\t";
	}
	std::cout << "quality" << std::endl;
	for (unsigned i = 1; i < length(stats.lengthsHisto); ++i) {
		std::cout << i << '\t';
		if (options.refFileName != "") {
			std::cout << stats.mismatchHisto[i] << '\t';
			std::cout << stats.insertHisto[i] << '\t';
			std::cout << stats.deletionHisto[i] << '\t';
			std::cout
					<< (stats.mismatchHisto[i] + stats.insertHisto[i]
							+ stats.deletionHisto[i]) / (double) (reads)
					<< '\t';
		}
		double e = stats.avrgQuality[i] / stats.lengthsHisto[i];
		std::cout << e << std::endl;
	}
	std::cout << ">>END_MODULE" << std::endl;

	if (options.refFileName != "") {
		std::cout << "#Edit distance" << std::endl;
		std::cout << ">>distance\tnumber of reads" << std::endl;
		for (unsigned i = 0; i < length(stats.editDistanceHisto); ++i) {
			if (stats.editDistanceHisto[i] > 0)
				std::cout << i << '\t' << stats.editDistanceHisto[i]
						<< std::endl;
		}
		std::cout << ">>END_MODULE" << std::endl;
	}

}
// Why did we do this???
// the only thing that comes to mind is that I believe I had once problem with writing to /dev/null
template<class T> const T& notZero(const T& a, const T& b) {
	if (a != 0)
		return a;
	return b;
}

/*
 * Analyze read that are sorted by ID
 * - check if alignment where unique
 * - can we check if there are biases on the flow cell, i.e. analyze the chip position in correlation with the number of reads
 *   or their quality value
 * - check NM/NH values
 * - analyze cigar string (intron length distribution, length as field
 * - count per chromosome unique/non-unique
 * - number of mutations
 * - quality string analysis (min, mean)
 * - quality over postion
 *
 *
 */
template<typename TStreamOrReader, typename TSeqString, typename TSpec,
		typename TFormat>
int analyze_idSorted(TStreamOrReader & reader,
		StringSet<TSeqString, TSpec> & seqs, Options const & options,
		TFormat const & tag, BamHeader header,
		BamIOContext<StringSet<CharString> > & context) {
	//StringSet<CharString> refNames;
	//NameStoreCache<StringSet<CharString> > refNamesCache(refNames);
	// context(refNames, refNamesCache);
	double programStartTime = sysTime();
	String<__uint64> qualSum;
	Align<Dna5String> align;
	__int64 reads = 0;
	//
	int writeFail = 0;
	Stats stats;
	resize(stats.flagsHisto, 65536, 0);

	// Create temporary file to ensure that the file is ordered by ID and not overwrite output file
	// if file is not properly sorted the header section will not be written correctly

	char tmpFPnam[2000];
	FILE *sFile;
	Stream<Bgzf> bFile;
	createTmpFile(tmpFPnam, options, &sFile, bFile);

	// If reference is given load it...
	StringSet<CharString> seqIds;
	//StringSet<TIdString, TIdSpec> seqIds;
	StringSet<Dna5String> faSeqs;
	if (options.refFileName != "") {
		std::cerr << "Reading references from " << options.refFileName
				<< std::endl;
		SequenceStream seqIO(toCString(options.refFileName),
				SequenceStream::READ, SequenceStream::FASTA);
		int read2out = readAll(seqIds, faSeqs, seqIO);
		if (read2out != 0) {
			std::cerr << "Could not read reference from " << options.refFileName
					<< " because: " << read2out << std::endl;
			return 1;
		}
		if (options.verbosity >= 2) {
			std::cerr << "Finished reading references from "
					<< options.refFileName << std::endl;
			std::cerr << " this took: " << sysTime() - programStartTime << "\n";
		}
		std::cerr << " read FASTA sequence " << length(faSeqs) << "\n";
	}

	// check if already run
	// Since there is a bug in IGV that doesn't allow for multiple versions of a given program to be exectued we are removing all old records.
	BamHeader newHeader;
	newHeader.sequenceInfos = header.sequenceInfos;
	unsigned myIdx = 0;
	for (unsigned i = 0; i < length(header.records); i++) {
		String<char> porgID = "ID";
		if (findTagKey(myIdx, porgID, header.records[i])) {
			String<char> tagValue;
			getTagValue(tagValue, myIdx, header.records[i]);
			if (tagValue == PROGNAME_ID) {
				if (options.verbosity >= 1)
					std::cerr << "already run " << std::endl;
				if (options.forceRedo) {
					erase(header.records, myIdx);
				} else {
					std::cerr
							<< "Program has already been run, try using -f to force execution "
							<< "\n";
					exit(-1);
				}
			} else {
				appendValue(newHeader.records, header.records[i]);
			}
		} else {
			appendValue(newHeader.records, header.records[i]);
		}

	}

	//new bam header line
	typedef typename BamHeaderRecord::TTag TTag;
	BamHeaderRecord progRecord;
	progRecord.type = BAM_HEADER_PROGRAM;
	appendValue(progRecord.tags, TTag("ID", PROGNAME_ID));
	appendValue(progRecord.tags, TTag("VN", VERSION));
	appendValue(progRecord.tags, TTag("CL", toCString(options.cmd)));
	appendValue(newHeader.records, progRecord);
	header = newHeader;

	if (options.inFormat == FORMAT_SAM) {
		writeFail = notZero(writeFail, write2(sFile, newHeader, context, tag));
		//fflush(sFile);
	} else {
		writeFail = notZero(writeFail, write2(bFile, newHeader, context, tag));
		//fflush(bFile);
	}

	// Some sanity checks:
	if (options.refFileName != "") {
		if (length(faSeqs) < length(header.sequenceInfos)) {
			std::cerr
					<< "There are more sequences in the SAM/BAM file than there are in the fasta file"
					<< std::endl;
			exit(-1);
		}
		if (options.verbosity >= 2)
			std::cerr << "sequence lengths for header "
					<< length(header.sequenceInfos) << " and fasta "
					<< length(faSeqs) << std::endl;

		for (unsigned i = 0; i < length(header.sequenceInfos); i++) {
			if (unsigned(header.sequenceInfos[i].i2) != length(faSeqs[i])) {
				std::cerr << "sequence lengths for "
						<< header.sequenceInfos[i].i1 << " and " << seqIds[i]
						<< " don't match (" << header.sequenceInfos[i].i2
						<< " : " << length(faSeqs[i]) << ")" << std::endl;
				exit(-1);
			}
			if (options.verbosity >= 2)
				std::cerr << "sequence lengths for "
						<< header.sequenceInfos[i].i1 << " and " << seqIds[i]
						<< " do match (" << header.sequenceInfos[i].i2 << " : "
						<< length(faSeqs[i]) << ")" << std::endl;
		}
		if (options.verbosity >= 2)
			std::cerr << "fasta file and header seem to correspond"
					<< std::endl;
	}

	resize(stats.rIDs, length(header.sequenceInfos) + 1, 0);
	StringSet<BamAlignmentRecord> idSet;
	BamAlignmentRecord oldRec;
	oldRec.qName = "";
	BamAlignmentRecord record;

	// Read alignments.
	if (options.verbosity >= 2)
		std::cerr << "Reading alignments" << std::endl;
	unsigned maxLen = 1;
	unsigned readCounter = 0;
	bool foundMultiHit = false;
	while (!atEnd(reader)) {
		readCounter++;
		if (!(readCounter % 100000)) {
			if (readCounter == 100000) {
				std::cerr
						<< "N(seq)\ttime[sec]\tRSS[kB]\tShared Memory[kB]\tPrivate Memory[kB]\n";
			}
			std::cerr << readCounter << "\t" << sysTime() - programStartTime
					<< "\t";
#if defined(_WIN32)
#else
			int tSize = 0, resident = 0, share = 0;
			std::ifstream buffer("/proc/self/statm");
			buffer >> tSize >> resident >> share;
			buffer.close();
			long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
			double rss = resident * page_size_kb;
			std::cerr << rss << "\t";
			double shared_mem = share * page_size_kb;
			std::cerr << shared_mem << "\t";
			std::cerr << rss - shared_mem;
#endif
			std::cerr << "\n";
		}
		// Read alignment record.
		if (readRecord(record, context, reader, tag) != 0) {
			std::cerr << "Could not read alignment!" << std::endl;
			return 1;
		}
		stats.numRecords++;
		stats.alignedRecords += !hasFlagUnmapped(record);
		stats.flagsHisto[record.flag]++;

		BamTagsDict tagDict(record.tags);

		//Lengths histogram
		lengthHist(stats, record);

		// count reference ids
		{
			if (record.rId == BamAlignmentRecord::INVALID_REFID) {
				stats.rIDs[0]++;
			} else {
				stats.rIDs[record.rId + 1]++;
			}
		}

		//stats with fasta reference genome...
		fastaStats(options, record, maxLen, stats, faSeqs, align, qualSum,
				reads);

		/*
		 * we basically can only handle one group at a time any other property has to be calculated for each record
		 */
		// MD,NM distance to reference (should already be calculated... We would need the reference genome for this and compare the it
		// therefore we assume that this has already been correctly calculated.
		// H0, H1, H2 this would be cool to have, but it is too expensive to calculate
		// ni number of introns
		intronCount(record, stats, tagDict);

		// il total intron length
		intronLength(record, stats, tagDict);

		// av average quality value
		averageQV(record, stats, tagDict);

		if (options.verbosity >= 2)
			write2(std::cout, record, context, Sam());
		/*
		 * IH tag: Query hit index, indicating the alignment record is the i-th one stored
		 * NH tag: Number of reported alignments that contains the query in the current record
		 */
		//Only occurs the first time
		if (oldRec.qName == "") {
			oldRec = record;
			appendValue(idSet, record);
		} else {
			// same query Name => collect in set
			if (oldRec.qName == record.qName) {
				foundMultiHit = true;
				appendValue(idSet, record);
			} else {
				// different Name => write out old...
				if (oldRec.rId == BamAlignmentRecord::INVALID_REFID) {
					BamTagsDict tagDictid(oldRec.tags);
					setTagValue(tagDictid, "NH", 0, 'i');
					if (options.inFormat == FORMAT_SAM) {
						writeFail = notZero(writeFail,
								write2(sFile, oldRec, context, tag));
					} else {
						writeFail = notZero(writeFail,
								write2(bFile, oldRec, context, tag));
					}
				} else {
					unsigned nID = length(idSet);
					for (unsigned i = 0; i < nID; i++) {
						BamAlignmentRecord bar = idSet[i];
						BamTagsDict tagDictid(bar.tags);
						bool ok = setTagValue(tagDictid, "NH", nID, 'i');
						//if(ok)std::cerr << "NH: OK\n";else std::cerr << "NH: NOTOK\n";
						ok = setTagValue(tagDictid, "IH", i, 'i');
						//if(ok)std::cerr << "iH: OK\n";else std::cerr << "iH: NOTOK\n";
						// CC reference name of the next hit; "=" for the same chromosome
						// TODO once Strings can be stored, store them as Strings.
						// i.e. currently SEQAN is not able to handle strings...!!!
						if (i < nID - 1) {
							if (idSet[i].rId == idSet[i + 1].rId) {
								setTagValue(tagDictid, "CC", '=');
							} else {
								setTagValue(tagDictid, "CC", idSet[i + 1].rId);
							}
						}
						if (options.inFormat == FORMAT_SAM) {
							writeFail = notZero(writeFail,
									write2(sFile, bar, context, tag));
						} else {
							writeFail = notZero(writeFail,
									write2(bFile, bar, context, tag));
						}

					}
					unsigned len = length(stats.nhHisto);
					resize(stats.nhHisto, std::max(len, nID + 1), 0);
					stats.nhHisto[nID] += 1;
				}
				clear(idSet);
				appendValue(idSet, record);
			}
		}
		oldRec = record;

	}
	// NH tag
	//finally write out last group
	int nID = length(idSet);
	for (int i = 0; i < nID; i++) {
		BamAlignmentRecord bar = idSet[i];
		BamTagsDict tagDictid(bar.tags);
		setTagValue(tagDictid, "NH", nID, 'i');
		if (options.inFormat == FORMAT_SAM) {
			writeFail = notZero(writeFail,
					write2(sFile, bar, context, tag));
		} else {
			writeFail = notZero(writeFail,
					write2(bFile, bar, context, tag));
		}

	}

	//rename (tmpnam, outFPname);
	if (options.inFormat == FORMAT_SAM) {
		fclose(sFile);
	} else {
		close(bFile);
	}
	if (options.outFileName != "" && writeFail == 0) {
		if (myrename(tmpFPnam, toCString(options.outFileName)) != 0) {
			perror("there was a problem renaming the output file...\n");
		}
	} else {
		if (writeFail != 0) {
			std::cerr
					<< "WARNINIG: there was a problem writing to tmp file ; removing tmp file\n";
		} else {
			std::cerr << "WARNINIG: no output file given; removing tmp file\n";
		}
		remove(tmpFPnam);
	}

	// Print stats
	printStatsIDsorted(stats, header, options, reads);
	return 0;
}

void clearROI(ROI &roi) {
	roi.min = INT_MAX;
	roi.max = 0;
	roi.rId = -2;
	resize(roi.count, 0, 0);
}

void writeOutROI(FILE * roiFile, ROI roi, BamHeader header,
		CharString const & dir) {
	std::stringstream roiIdStr; // declare the string
	roiIdStr << "region_" << roi.roiID;
	fprintf(roiFile, "%s\t%d\t%d\t%s\t%d\t%s\t%d",
			toCString(header.sequenceInfos[roi.rId].i1), roi.min + 1,
			roi.max + 1, toCString(dir), roi.len, roiIdStr.str().c_str(),
			roi.countMax);
	for (unsigned i = 0; i < roi.len; i++) {
		fprintf(roiFile, "\t%d", roi.count[i]);

	}
	fprintf(roiFile, "\n");
}

template<typename TStream, typename TFormat>
int writeOutRecordSet(TStream &fp, StringSet<BamAlignmentRecord> recSet,
		BamIOContext<StringSet<CharString> > & context, TFormat const & tag,
		ROI & roi) {

	unsigned maxC = 0;
	int writeFail = 0;
	for (unsigned i = 0; i < length(roi.count); i++) {
		if (maxC < roi.count[i])
			maxC = roi.count[i];
	}
	roi.countMax = maxC;
	for (unsigned i = 0; i < length(recSet); i++) {
		BamTagsDict tagDict(recSet[i].tags);
		setTagValue(tagDict, "rl", roi.len, 'i');
		setTagValue(tagDict, "ri", roi.min + 1, 'i');
		setTagValue(tagDict, "rx", roi.max + 1, 'i');
		setTagValue(tagDict, "ry", maxC, 'i');
		setTagValue(tagDict, "rd", roi.roiID, 'i');
		setTagValue(tagDict, "rn", length(recSet), 'i');
		writeFail = notZero(writeFail, write2(fp, recSet[i], context, tag));
		if (writeFail != 0) {
			std::cerr << "failed writing in writeOutRecordSet: "
					<< write2(std::cerr, recSet[i], context, tag) << std::endl;
		}
	}
	return maxC;
}
bool errorMessageSent = false;

template<typename TFormat>
void roiStats(ROI& roi, BamAlignmentRecord& record, int readStartPos,
		unsigned &roiN, const Options& options, FILE*& sFile,
		StringSet<BamAlignmentRecord>& recordSet,
		BamIOContext<StringSet<CharString> >& context, const TFormat& tag,
		Stream<Bgzf>& bFile, FILE*& roiFile, BamHeader& header, Stats& stats,
		int readEndPos) {

	if (roi.rId != record.rId || roi.max < readStartPos) {
		//rew ROI
		roiN++; // number of ROIs
		roi.roiID = roiCount++;
		if (options.verbosity >= 2) {
			std::cerr << "new roi: roi.roiID:" << roi.roiID << " roiN: " << roiN
					<< std::endl;
		}
		if (roi.max == 0) {
			//first roi
			if (options.verbosity >= 2) {
				std::cerr << "roi.max = " << roi.max << std::endl;
			}
		} else {
			// within writeOutRecordSet the roi.countMax is set!
			// This is writing out the previous one
			if (options.inFormat == FORMAT_SAM) {
				writeOutRecordSet(sFile, recordSet, context, tag, roi);
			} else {
				writeOutRecordSet(bFile, recordSet, context, tag, roi);
			}
			if (options.roiFileName != "") {
				CharString dir = "+";
				if (hasFlagRC(recordSet[0]))
					dir = "-";
				if (options.verbosity >= 2) {
					std::cerr << "directionality: " << dir << std::endl;
				}
				if (options.verbosity >= 2) {
					std::cerr << "writing out ROI" << std::endl;
				}

				writeOutROI(roiFile, roi, header, dir);
			}
			//resize(roiList, length(roiList) + 1);
			//roiList[length(roiList) - 1] = roi;
			unsigned len = length(stats.roiLenHisto);
			resize(stats.roiLenHisto, std::max(len, roi.len + 1), 0);
			stats.roiLenHisto[roi.len] += 1;
			len = length(stats.roicountMax);
			resize(stats.roicountMax, std::max(len, roi.countMax + 1), 0);
			stats.roicountMax[roi.countMax] += 1;
			clear(recordSet);
			clearROI(roi);
		}

		//initialize new ROI
		roi.min = readStartPos;
		roi.max = readEndPos;
		roi.len = readEndPos - readStartPos;
		roi.rId = record.rId;
		resize(roi.count, roi.len + 1, 1);
		appendValue(recordSet, record);
	} else {
		//TODO I guess if there are too many reads appended we run into a problem...
		if (options.verbosity >= 2) {
			std::cerr << "old roi: roi.roiID:" << roi.roiID << " roiN: " << roiN
					<< std::endl;
		}
		if (length(recordSet) <= options.maxCoverage) {
			appendValue(recordSet, record);
		} else {
			if (!errorMessageSent) {
				std::cerr << "WARNING: max Size of Record set is reached ("
						<< options.maxCoverage
						<< "). Closing Region and creating new one:"
						<< toCString(header.sequenceInfos[roi.rId].i1)
						<< " at: " << roi.min << " to " << roi.max << std::endl
						<< "all following warnings will not be displayed"
						<< std::endl;
				errorMessageSent = true;
			}
			appendValue(recordSet, record);
			if (options.inFormat == FORMAT_SAM) {
				writeOutRecordSet(sFile, recordSet, context, tag, roi);
			} else {
				writeOutRecordSet(bFile, recordSet, context, tag, roi);
			}
			clear(recordSet);
		}

		//read falls within current ROI
		if (roi.min > readStartPos) {
			std::cerr << "This file is not sorted by position" << std::endl;
			exit(-1);
		}
		if (readEndPos > roi.max) {
			roi.max = readEndPos;
			roi.len = readEndPos - roi.min;
			resize(roi.count, roi.len + 1, 0);
		}
		for (int i = readStartPos - roi.min; i < readEndPos - roi.min; i++) {
			roi.count[i]++;
		}
		if (options.verbosity >= 2) {
			for (unsigned i = 0; i < roi.len; i++) {
				std::cerr << roi.count[i] << " ";
			}
			std::cerr << std::endl;
		}

	}

}

/*
 * Analyze reads that are sorted by position
 */

template<typename TStreamOrReader, typename TSeqString, typename TSpec,
		typename TFormat>
int analyze_posSorted(TStreamOrReader & reader,
		StringSet<TSeqString, TSpec> & seqs, Options const & options,
		TFormat const & tag, BamHeader header,
		BamIOContext<StringSet<CharString> > & context) {

	// Create temporary file to ensure that the file is ordered by ID and not overwrite output file
	// if file is not properly sorted the header section will not be written correctly

	char tmpFPnam[2000];
	FILE *sFile;
	FILE *roiFile;
	Stream<Bgzf> bFile;
	createTmpFile(tmpFPnam, options, &sFile, bFile);
	int writeFail = 0;
	double programStartTime = sysTime();
	String<__uint64> qualSum;
	// Read alignments.
	if (options.verbosity >= 1)
		std::cerr << "Reading alignments" << std::endl;
	Align<Dna5String> align;
	__int64 reads = 0;

	Stats stats;
	resize(stats.flagsHisto, 65536, 0);

	// roi File
	if (options.roiFileName != "") {
		roiFile = fopen(toCString(options.roiFileName), "w");
		if (roiFile == NULL) {
			std::cerr << "Could not open " << options.roiFileName << std::endl;
			return 1;
		}
		fprintf(roiFile, "# 1.  column: reference name\n");
		fprintf(roiFile,
				"# 2.  column: starting position on reference genome\n");
		fprintf(roiFile, "# 3.  column: end position on reference genome\n");
		fprintf(roiFile, "# 4.  column: strand +/-(. if empty)\n");
		fprintf(roiFile, "# 5.  column: length of region\n");
		fprintf(roiFile, "# 6.  column: region_name\n");
		fprintf(roiFile, "# 7.  column: max_Count\n");
		fprintf(roiFile, "# 8+. columns: individual counts\n");
	}
	// If reference is given load it...
	StringSet<CharString> seqIds;
	StringSet<Dna5String> faSeqs;
	if (options.refFileName != "") {
		std::cerr << "Reading references from " << options.refFileName
				<< std::endl;
		SequenceStream seqIO(toCString(options.refFileName),
				SequenceStream::READ, SequenceStream::FASTA);
		int read2out = readAll(seqIds, faSeqs, seqIO);
		if (read2out != 0) {
			std::cerr << "Could not read reference from " << options.refFileName
					<< " because: " << read2out << std::endl;
			return 1;
		}
		if (options.verbosity >= 2) {
			std::cerr << "Finished reading references from "
					<< options.refFileName << std::endl;
			std::cerr << " this took: " << sysTime() - programStartTime << "\n";
		}
		std::cerr << " read FASTA sequence " << length(faSeqs) << "\n";
	}

	// check if already run
	// Since there is a bug in IGV that doesn't allow for multiple versions of a given program to be exectued we are removing all old records.
	BamHeader newHeader;
	newHeader.sequenceInfos = header.sequenceInfos;
	unsigned myIdx = 0;
	for (unsigned i = 0; i < length(header.records); i++) {
		String<char> porgID = "ID";
		if (findTagKey(myIdx, porgID, header.records[i])) {
			String<char> tagValue;
			getTagValue(tagValue, myIdx, header.records[i]);
			if (tagValue == PROGNAME_POS) {
				if (options.verbosity >= 2)
					std::cerr << "already run " << std::endl;
				if (options.forceRedo) {
					erase(header.records, myIdx);
				} else {
					std::cerr
							<< "Program has already been run, try using -f to force execution "
							<< "\n";
					exit(-1);
				}
			} else {
				appendValue(newHeader.records, header.records[i]);
			}
		} else {
			appendValue(newHeader.records, header.records[i]);
		}
	}

	//new bam header line
	typedef typename BamHeaderRecord::TTag TTag;
	BamHeaderRecord progRecord;
	progRecord.type = BAM_HEADER_PROGRAM;
	appendValue(progRecord.tags, TTag("ID", PROGNAME_POS));
	appendValue(progRecord.tags, TTag("VN", VERSION));
	appendValue(progRecord.tags, TTag("CL", toCString(options.cmd)));
	appendValue(newHeader.records, progRecord);
	header = newHeader;

	if (options.inFormat == FORMAT_SAM) {
		writeFail = notZero(writeFail, write2(sFile, newHeader, context, tag));
		//fflush(sFile);
	} else {
		writeFail = notZero(writeFail, write2(bFile, newHeader, context, tag));
		//fflush(bFile);
	}
	// Some sanity checks:
	if (options.refFileName != "") {
		if (length(faSeqs) < length(header.sequenceInfos)) {
			std::cerr
					<< "There are more sequences in the SAM/BAM file than there are in the fasta file"
					<< std::endl;
			exit(-1);
		}
		if (options.verbosity >= 2)
			std::cerr << "sequence lengths for header "
					<< length(header.sequenceInfos) << " and fasta "
					<< length(faSeqs) << std::endl;

		for (unsigned i = 0; i < length(header.sequenceInfos); i++) {
			if (unsigned(header.sequenceInfos[i].i2) != length(faSeqs[i])) {
				std::cerr << "sequence lengths for "
						<< header.sequenceInfos[i].i1 << " and " << seqIds[i]
						<< " don't match (" << header.sequenceInfos[i].i2
						<< " : " << length(faSeqs[i]) << ")" << std::endl;
				exit(-1);
			}
			if (options.verbosity >= 2)
				std::cerr << "sequence lengths for "
						<< header.sequenceInfos[i].i1 << " and " << seqIds[i]
						<< " do match (" << header.sequenceInfos[i].i2 << " : "
						<< length(faSeqs[i]) << ")" << std::endl;
		}
		if (options.verbosity >= 2)
			std::cerr << "fasta file and header seem to correspond"
					<< std::endl;
	}

	resize(stats.rIDs, length(header.sequenceInfos) + 1, 0);
	StringSet<BamAlignmentRecord> recordSet;
	StringSet<BamAlignmentRecord> recordSetRC;
	BamAlignmentRecord oldRec;
	BamAlignmentRecord record;

	// Read alignments.
	if (options.verbosity >= 2)
		std::cerr << "Reading alignments" << std::endl;
	unsigned maxLen = 1;
	unsigned readCounter = 0;
	ROI roi;
	ROI roiRC;

	//String<ROI> roiList;
	unsigned roiN = 0;
	unsigned roiNRC = 0;
	while (!atEnd(reader)) {
		readCounter++;
		if (!(readCounter % 100000)) {
			if (readCounter == 100000) {
				std::cerr
						<< "N(seq)\ttime[sec]\tRSS[kB]\tShared Memory[kB]\tPrivate Memory[kB]\n";
			}
			std::cerr << readCounter << "\t" << sysTime() - programStartTime
					<< "\t";
#if defined(_WIN32)
#else
			int tSize = 0, resident = 0, share = 0;
			std::ifstream buffer("/proc/self/statm");
			buffer >> tSize >> resident >> share;
			buffer.close();
			long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
			double rss = resident * page_size_kb;
			std::cerr << rss << "\t";
			double shared_mem = share * page_size_kb;
			std::cerr << shared_mem << "\t";
			std::cerr << rss - shared_mem;
#endif
			std::cerr << "\n";
		}

		// Read alignment record.
		if (readRecord(record, context, reader, tag) != 0) {
			std::cerr << "Could not read alignment!" << std::endl;
			return 1;
		}
		if (options.verbosity >= 2)
			write2(std::cerr, record, context, Sam());
		;

		stats.numRecords++;
		stats.alignedRecords += !hasFlagUnmapped(record);
		stats.flagsHisto[record.flag]++;

		BamTagsDict tagDict(record.tags);

		//Lengths histogram
		lengthHist(stats, record);

		// count reference ids
		{
			if (record.rId == BamAlignmentRecord::INVALID_REFID) {
				stats.rIDs[0]++;
			} else {
				stats.rIDs[record.rId + 1]++;
			}
		}

		//stats with fasta reference genome...
		fastaStats(options, record, maxLen, stats, faSeqs, align, qualSum,
				reads);

		// av average quality value
		averageQV(record, stats, tagDict);

		/*
		 * This part depends on the sorting...
		 */
		if (!hasFlagUnmapped(record)) {
			/*
			 * we basically can only handle one group at a time any other property has to be calculated for each record
			 */
			// MD,NM distance to reference (should already be calculated... We would need the reference genome for this and compare the it
			// therefore we assume that this has already been correctly calculated.
			// H0, H1, H2 this would be cool to have, but it is too expensive to calculate
			// ni number of introns
			intronCount(record, stats, tagDict);

			// il total intron length
			intronLength(record, stats, tagDict);

			//Only occurs the first time
			int readStartPos = record.pos;
			int readEndPos = 0;
			for (unsigned ti = 0; ti < length(tagDict); ++ti) {
				CharString tmp;
				tmp = getTagKey(tagDict, ti);
				if (tmp == "cl") {
					__int32 x = 0;
					bool res = extractTagValue(x, tagDict, ti);
					SEQAN_ASSERT_MSG(res, "Not a valid integer at pos %u!", ti);
					readEndPos = unsigned(x) + readStartPos;
				}
			}
			if (readEndPos == 0) {
				std::cerr
						<< "There was a problem getting the length of the read"
						<< std::endl;
				write2(std::cerr, record, context, Sam());
			}
			/*
			 * ROIs
			 */
			if (options.verbosity >= 2) {
				std::cerr << readStartPos << " <-start  " << readEndPos
						<< " <-end" << std::endl;
				std::cerr << roi.min << " " << roi.max << " " << roi.len << " "
						<< std::endl;
			}
			if (options.verbosity >= 2) {
				write2(std::cerr, record, context, Sam());
				if (hasFlagRC(record))
					std::cerr << "hasFlag" << std::endl;
				else
					std::cerr << "doesn't have Flag " << std::endl;
			}
			if ((!hasFlagRC(record)) && options.strandSpecific) {
				if (options.verbosity >= 2) {
					std::cerr << " not RC AND SS  " << std::endl;
				}
				//roiAnalyzer.pushRecord(record, ...);
				roiStats(roi, record, readStartPos, roiN, options, sFile,
						recordSet, context, tag, bFile, roiFile, header, stats,
						readEndPos);
			} else {
				if (options.verbosity >= 2) {
					std::cerr << " RC AND SS or ND  " << std::endl;
				}
				roiStats(roiRC, record, readStartPos, roiNRC, options, sFile,
						recordSetRC, context, tag, bFile, roiFile, header,
						stats, readEndPos);

			}
		} else { //unmapped
			if (roiRC.max != 0) {
				if (options.verbosity >= 1)
					std::cerr << "last group" << std::endl;
				//resize(roiList, length(roiList) + 1);
				if (options.inFormat == FORMAT_SAM) {
					writeOutRecordSet(sFile, recordSet, context, tag, roi);
					writeOutRecordSet(sFile, recordSetRC, context, tag, roi);
				} else {
					writeOutRecordSet(bFile, recordSet, context, tag, roi);
					writeOutRecordSet(bFile, recordSetRC, context, tag, roi);
				}

				//roiList[length(roiList) - 1] = roi;
				unsigned len = length(stats.roiLenHisto);
				resize(stats.roiLenHisto, std::max(len, roi.len + 1), 0);
				stats.roiLenHisto[roi.len] += 1;
				len = length(stats.roicountMax);
				resize(stats.roicountMax, std::max(len, roi.countMax + 1), 0);
				stats.roicountMax[roi.countMax] += 1;
				clear(recordSet);
				clear(recordSetRC);
				clearROI(roi);
				clearROI(roiRC);
			}
//			try {
			if (options.verbosity >= 1) {
				int tSize = 0, resident = 0, share = 0;
#if defined(_WIN32)
#else
				std::ifstream buffer("/proc/self/statm");
				buffer >> tSize >> resident >> share;
				buffer.close();
				long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
				double rss = resident * page_size_kb;
#endif
			}
			if (options.inFormat == FORMAT_SAM) {
				writeFail = notZero(writeFail,
						write2(sFile, record, context, tag));
			} else {
				writeFail = notZero(writeFail,
						write2(bFile, record, context, tag));
			}
//			} catch (int e) {
//				std::cerr << "ERROR: caught signal " << e << std::endl;
//				std::cerr << "optionsFormat" << options.inFormat << std::endl;
//				std::cerr << "record NR: " << stats.numRecords << std::endl;
//					std::cerr << "rss:" << rss << "\t";
//					double shared_mem = share * page_size_kb;
//					std::cerr << "shared mem: " << shared_mem << "\t";
//					std::cerr << "rss-shared: " << rss - shared_mem << "\n";
//				}
//				exit(-1);
//			}
		}

	}
	//finally write out last group
	if (roiRC.max != 0) {
		std::cerr << "last group (RC)" << std::endl;
		//resize(roiList, length(roiList) + 1);
		if (options.inFormat == FORMAT_SAM) {
			writeOutRecordSet(sFile, recordSetRC, context, tag, roiRC);
		} else {
			writeOutRecordSet(bFile, recordSetRC, context, tag, roiRC);
		}
		if (options.roiFileName != "") {
			if (options.verbosity >= 2) {
				std::cerr << "writing out ROI" << std::endl;
			}
			CharString dir = "+";
			if (hasFlagRC(recordSetRC[0]))
				dir = "-";
			writeOutROI(roiFile, roiRC, header, dir);
		}
		//roiList[length(roiList) - 1] = roi;
		unsigned len = length(stats.roiLenHisto);
		resize(stats.roiLenHisto, std::max(len, roi.len + 1), 0);
		stats.roiLenHisto[roi.len] += 1;
		len = length(stats.roicountMax);
		resize(stats.roicountMax, std::max(len, roi.countMax + 1), 0);
		stats.roicountMax[roi.countMax] += 1;
		clear(recordSetRC);
		clearROI(roiRC);
	}
	if (roi.max != 0) {
		std::cerr << "last group" << std::endl;
		//resize(roiList, length(roiList) + 1);
		if (options.inFormat == FORMAT_SAM) {
			writeOutRecordSet(sFile, recordSet, context, tag, roi);
		} else {
			writeOutRecordSet(bFile, recordSet, context, tag, roi);
		}
		if (options.roiFileName != "") {
			if (options.verbosity >= 2) {
				std::cerr << "writing out ROI" << std::endl;
			}
			CharString dir = "+";
			if (hasFlagRC(recordSet[0]))
				dir = "-";
			writeOutROI(roiFile, roi, header, dir);
		}
		//roiList[length(roiList) - 1] = roi;
		unsigned len = length(stats.roiLenHisto);
		resize(stats.roiLenHisto, std::max(len, roi.len + 1), 0);
		stats.roiLenHisto[roi.len] += 1;
		len = length(stats.roicountMax);
		resize(stats.roicountMax, std::max(len, roi.countMax + 1), 0);
		stats.roicountMax[roi.countMax] += 1;
		clear(recordSet);
		clearROI(roi);
	}

	//rename (tmpnam, outFPname);
	if (options.inFormat == FORMAT_SAM) {
		fclose(sFile);
	} else {
		close(bFile);
	}
	if (options.outFileName != "" && writeFail == 0) {
		if (myrename(tmpFPnam, toCString(options.outFileName)) != 0) {
			perror("there was a problem renaming the output file...\n");
		}
	} else {
		if (writeFail != 0) {
			std::cerr
					<< "WARNINIG: there was a problem writing the tmp file; removing tmp file\n";
		} else {
			std::cerr << "WARNINIG: no output file given; removing tmp file\n";
		}
		remove(tmpFPnam);
	}
	if (options.roiFileName != "") {
		fclose(roiFile);
	}

	printStatsPosSorted(stats, header, options, reads);

	//printROIs(roiList, options);

	return 0;
}

/*
 * Here we deal with the sorting and organizing different order of events...
 */
template<typename TStreamOrReader, typename TSeqString, typename TSpec,
		typename TFormat>
int doWork(TStreamOrReader & reader, StringSet<TSeqString, TSpec> & seqs,
		Options const & options, TFormat const & tag) {
	StringSet<CharString> refNames;
	NameStoreCache<StringSet<CharString> > refNamesCache(refNames);
	BamIOContext<StringSet<CharString> > context(refNames, refNamesCache);
	String<__uint64> qualSum;

// Read header.
	BamHeader header;
	if (options.verbosity >= 2)
		std::cerr << "Reading header" << std::endl;
	if (readRecord(header, context, reader, tag) != 0) {
		std::cerr << "Could not read header!" << std::endl;
		return 1;
	}
	if (length(header.sequenceInfos) == 0) {
		std::cerr
				<< "There was no information on sequences in the header or header is missing. Please correct\n"
				<< std::endl;
		return -1;

	}
	int ret = 0;
	switch (options.sortOrder) {
	case SORTORDER_ID:
		ret = analyze_idSorted(reader, seqs, options, tag, header, context);
		if (ret != 0) {
			//TODO remove temp file if not -vv
		}
		//TODO sort
		//sortAlignedReads(reader, SortBeginPos());
		//analyze_posSorted(reader, seqs, options, tag);
		break;
	case SORTORDER_POS:
		ret = analyze_posSorted(reader, seqs, options, tag, header, context);
		if (ret != 0) {
			//TODO remove temp file if not -vv
		} //TODO sort
		  //sortAlignedReads(reader, SortReadId());
		  //analyze_idSorted(reader, seqs, options, tag);
		break;
	case SORTORDER_NA:
		//TODO sort
		//sortAlignedReads(reader, SortReadId());
		//analyze_idSorted(reader, seqs, options, tag);
		//TODO sort
		//sortAlignedReads(reader, SortBeginPos());
		//analyze_posSorted(reader, seqs, options, tag);
		std::cerr
				<< "no sorting order has been given. Due to the current limitatios of seqan to not being able to sort SAM records we have to abord.";
		break;
	}

	return 0;
}

int mainWithOptions(Options & options) {
	typedef Iterator<String<CharString> >::Type TIterator;
	programStartTime = sysTime();
	std::cout << ">>Parameters\n";
	std::cout << "#Name\tvalue" << std::endl;
	std::cout << "Version\t" << VERSION << "\n";
	std::cout << "input file: \t\"" << options.inFileName << "\"" << std::endl;
	std::cout << "output file: \t\"" << options.outFileName << "\""
			<< std::endl;
	if (options.inFormat == FORMAT_SAM)
		std::cout << "file format: \tSAM" << std::endl;
	else
		std::cout << "file format: \tBAM" << std::endl;
	if (options.forceRedo)
		std::cout << "forceRedo\tyes" << std::endl;
	else
		std::cout << "forceRedo\tno" << std::endl;
	std::cout << "reference file\t\"" << options.refFileName << "\""
			<< std::endl;
	std::cout << "ROI filename \t\"" << options.roiFileName << "\""
			<< std::endl;
	std::cout << "maxCoverage\t" << options.maxCoverage << "" << std::endl;
	if (options.sortOrder == SORTORDER_ID)
		std::cout << "sortOrder\tID" << std::endl;
	else if (options.sortOrder == SORTORDER_POS)
		std::cout << "sortOrder\tPosition" << std::endl;
	else
		std::cout << "sortOrder\tnone" << std::endl;
	std::cout << "tmpDir\t" << options.tmpDir << "" << std::endl;
	if (options.strandSpecific)
		std::cout << "strand specific\tyes" << std::endl;
	else
		std::cout << "strand specific\tno" << std::endl;
	std::cout << ">>END_MODULE" << std::endl;

	StringSet<CharString> seqIds;
	StringSet<Dna5String> seqs;

    /*
    BamStream inStream(toCString(options.inFileName));
    if (!isGood(inStream))
    {
        std::cerr << "ERROR: Could not open " << options.inFileName << "\n";
        return 1;
    }
    doWork(inStream, seqs, options);
    */

	String<char, MMap<> > seqMMapString;

	if (options.inFormat == FORMAT_SAM) {
		if (options.verbosity >= 2)
			std::cerr << "Opening SAM file " << options.inFileName << std::endl;
		String<char, MMap<> > samMMapString;
		if (!open(samMMapString, toCString(options.inFileName), OPEN_RDONLY)) {
			std::cerr << "Could not open " << options.inFileName << std::endl;
			return 1;
		}
		RecordReader<String<char, MMap<> >, SinglePass<Mapped> > samReader(
				samMMapString);
		return doWork(samReader, seqs, options, Sam());
	} else // options.inFormat == FORMAT_BAM
	{
		if (options.verbosity >= 2)
			std::cerr << "Opening BAM file " << options.inFileName << std::endl;
		Stream<Bgzf> bamStream;
		if (!open(bamStream, toCString(options.inFileName), "r")) {
			std::cerr << "Could not open " << options.inFileName << std::endl;
			return 1;
		}
		return doWork(bamStream, seqs, options, Bam());
	}
	return 0;
}

/*
 //#else

 int main(int, char const **) {
 std::cerr << "bam_stats can only be compiled correctly with zlib."
 << std::endl;
 return 0;
 }

 //#endif  // #if SEQAN_HAS_ZLIB
 */
#endif  // #ifndef SANDBOX_JAGLA_APPS_SAMBAMSTAT_SAMBAMSTAT_H_
