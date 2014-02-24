// ==========================================================================
//                                   faUniq
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
// Author: Bernd Jagla <bernd.jagla@pasteur.fr>
// ==========================================================================

#ifndef SANDBOX_JAGLA_APPS_FAUNIQ_FAUNIQ_H_
#define SANDBOX_JAGLA_APPS_FAUNIQ_FAUNIQ_H_
#include <string>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <cstdio>
#include <fstream>
#if SEQAN_HAS_ZLIB
#include <zlib.h>
#endif  // #if SEQAN_HAS_ZLIB
#if SEQAN_HAS_BZIP2
#include <bzlib.h>
#endif  // #if SEQAN_HAS_BZIP2
#include <seqan/file.h>
#include <seqan/stream.h>
#include <seqan/misc/misc_cmdparser.h>

#include <seqan/misc/misc_cmdparser.h>

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================
const unsigned int BUFFER_SIZE = 1024 * 1024 * 4;

struct Options {
	bool showHelp;
	bool showVersion;
	CharString inputFileName;
	CharString outputFileName;
	String<CharString> texts;
	int verbosity;

	Options() {
		// Set defaults.
		showHelp = false;
		showVersion = false;
		verbosity = 0;
	}
};

// ============================================================================
// Metafunctions
// ============================================================================

///We use external Dna5 strings to construct the suffix array and lcp table in external memory.
///The concat-direct stringset stores all contigs in a single file.
typedef String<Dna5, External<> > TSequence;
typedef StringSet<CharString> TSequenceIds;
typedef StringSet<TSequence, Owner<ConcatDirect<> > > TSequences;

namespace seqan {

    // our stringset has not more than 4 billion chars
    // 32bit unsigned ints are sufficient to address a sequence
    template <>
    struct Size<TSequences>
    {
        typedef __uint32 Type;
    };

    // ... or an entry in the suffix array or lcp table
    template <>
    struct StringSetLimits<TSequences>
    {
        typedef String<__uint32> Type;
    };
}

// ============================================================================
// Functions
// ============================================================================

void setupCommandLineParser(CommandLineParser & parser,
		Options const & options) {
	addVersionLine(parser, "0.1");

	addTitleLine(parser, "**********************");
	addTitleLine(parser, "* faUniq *");
	addTitleLine(parser, "**********************");
	addTitleLine(parser, "");
	addTitleLine(parser, "(c) 2012 by Bernd Jagla <bernd.jagla@pasteur.fr>");

	addUsageLine(parser, "[OPTIONS] -if input");

	addSection(parser, "Main Options");
	addOption(parser,
			CommandLineOption("if", "input-file", "Input file (Fasta).",
					OptionType::String | OptionType::Label,
					options.inputFileName));
	addOption(parser,
			CommandLineOption("of", "output-file", "Output file (Wig).",
					OptionType::String | OptionType::Label,
					options.outputFileName));
	addOption(parser,
			CommandLineOption("v", "verbose", "Verbose output.",
					OptionType::Bool));
	addOption(parser,
			CommandLineOption("vv", "very-verbose", "Very verbose output.",
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

	getOptionValueLong(parser, "input-file", options.inputFileName);
	getOptionValueLong(parser, "output-file", options.outputFileName);

	// just changed values from 2/3 to 1/2
	// I guess it should be changed in the code not here...
	if (isSetLong(parser, "verbose"))
		options.verbosity = 1;
	if (isSetLong(parser, "very-verbose"))
		options.verbosity = 2;

	options.texts = getArgumentValues(parser);

	return 0;
}

int mainWithOptions(Options & options) {
	double programStartTime = sysTime();
	double startTime = 0;
	typedef Iterator<String<CharString> >::Type TIterator;
	if (options.verbosity >= 1) {
		std::cerr << "Option Arguments:" << std::endl;
		std::cerr << "  input file:  \"" << options.inputFileName << "\""
				<< std::endl;
		std::cerr << "  output file: \"" << options.outputFileName << "\""
				<< std::endl;
		std::cerr << "Non-option Arguments:" << std::endl;
		for (TIterator it = begin(options.texts); it != end(options.texts);
				++it) {
			std::cerr << "  " << *it << std::endl;
		}
	}
	startTime = sysTime();

	// -----------------------------------------------------------------------
	// Load FASTA
	// -----------------------------------------------------------------------
	double before = sysTime();

	typedef Iterator<TSequenceIds>::Type TSequenceIdsIter;
	typedef Iterator<TSequences>::Type TSequencesIter;

	MultiSeqFile multiSeqFile;
	if (!open(multiSeqFile.concat, toCString(options.inputFileName),
			OPEN_RDONLY)) {
		std::cerr << std::endl << "Could not open mmap file for reading."
				<< std::endl;
		return 1;
	}
	AutoSeqFormat format;
	guessFormat(multiSeqFile.concat, format);
	split(multiSeqFile, format);

	unsigned seqCount = length(multiSeqFile);
	TSequences seqs;
	TSequenceIds seqIDs;

	Dna5String seq;
	// CharString qual;
	CharString id;
	for (unsigned i = 0; i < seqCount; ++i) {
		assignSeq(seq, multiSeqFile[i], format); // read sequence
		assignSeqId(id, multiSeqFile[i], format); // read sequence id
		// we use reserve and append, as assign is not supported
		// by StringSet<..., Owner<ConcatDirect<> > >
		char * cid = toCString(id);
		strtok(cid, " \t");
		appendValue(seqs, seq, Generous());
		appendValue(seqIDs, cid, Generous());
	}

	double after = sysTime();
	//fprintf(stderr, "\t%f\n", after - before);
	for (unsigned i = 0; i < seqCount; i++) {
		std::cerr << seqIDs[i] << " " << length(seqs[i]) << "\n";
	}

//	StringSet<String<char> > mySet;
//	resize(mySet, 3);
//	mySet[0] = "SeqAn is a library for sequence analysis.";
//	mySet[1] = "The String class is the fundamental sequence type in SeqAn.";
//	mySet[2] = "Subsequences can be handled with SeqAn's Segment class.";

	///Then we create an @Class.Index@ of this @Class.StringSet@.
	typedef Index<TSequences> TMyIndex;
	TMyIndex myIndex(seqs);
	///Now we require that the index has a suffix array and a lcp table using the function @Function.indexRequire@.

// TODO needed to uncomment these line under windows....	
//	indexRequire(myIndex, EsaSA());
//	indexRequire(myIndex, EsaLcp());

	///We iterate over the lcp table and output the position, the value of the lcp table and the corresponding suffix using the functions @Function.saAt@ and @Function.lcpAt@.

	for (Size<TMyIndex>::Type i = 0; i < length(myIndex); ++i) {
		SAValue<TMyIndex>::Type p = saAt(i, myIndex);
		std::cout << seqIDs[p.i1] << "\t" << p.i2 + 1 << "\t"
				<< lcpAt(i, myIndex) << "\n";
	}
	return 0;
}

#endif  // #ifndef SANDBOX_JAGLA_APPS_FAUNIQ_FAUNIQ_H_
