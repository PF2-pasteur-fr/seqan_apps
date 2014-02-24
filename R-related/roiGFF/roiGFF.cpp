// ==========================================================================
//                                   roiPCA
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

#include <fstream>
#include <iostream>
#include <seqan/arg_parse.h>
#include <seqan/file.h>
#include <seqan/stream.h>

#include "roiGFF.h"

#define VERSION "0.3"

// Parse command line with ArgumentParser class.

seqan::ArgumentParser::ParseResult
parseCommandLine(Options & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("roiGFF");

    // Set short description, version, and date.
    setShortDescription(parser, "get ROIs based on GFF annotation");
    setVersion(parser, VERSION);
    setDate(parser, "March 2013");

    // Define usage line and long description.
    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] \\fB-if\\fP \\fIIN.roi\\fP"
				 "\\fB-ig\\fP \\fIIN.gff\fP"
                 "\\fB-of\\fP \\fIOUT.roi\\fP"
				 "[\\fB-ss\\fP] [\\fB-t\\fP \\fItype\\fP] "
				 "[\\fB-p\\fP \\fIparent attribute\\fP]]");
	addDescription(parser,
       "Takes regions from GFF file, calculates the overlap with the ROI file"
       "and creates a new ROI file based on the GFF regions.");
	
	// General Options

    addOption(parser, seqan::ArgParseOption("v", "verbose", "Verbose mode."));
    addOption(parser, seqan::ArgParseOption("vv", "vverbose", "Very verbose mode."));

    // Input / Output Options

    addSection(parser, "Input / Output Parameters");

	addOption(parser, seqan::ArgParseOption("if", "roi-file", "roi file", seqan::ArgParseOption::INPUTFILE));
    setValidValues(parser, "roi-file", "roi");
    setRequired(parser, "roi-file");

	addOption(parser, seqan::ArgParseOption("ig", "gff-file", "gff file", seqan::ArgParseOption::INPUTFILE));
    setValidValues(parser, "gff-file", "gff");
    setRequired(parser, "gff-file");

    addOption(parser, seqan::ArgParseOption("of", "output-file", "Output file", seqan::ArgParseOption::OUTPUTFILE));
    setValidValues(parser, "output-file", "roi");
    setRequired(parser, "output-file");

    // Optional Parameter
	addOption(parser,
		seqan::ArgParseOption("ss", "strandspecific",
			"calculate strandspecific stats (only position sorted)"));

	addOption(parser,
		seqan::ArgParseOption("t", "type",
		"Type to used (3. column in GFF file). If not set all types will be used."
		"There can only be one or none set.", seqan::ArgParseOption::STRING));
	setDefaultValue(parser, "type", "");

	addOption(parser,
		seqan::ArgParseOption("p", "parent",
			"tag to used to identify regions that should be concatednated."
			"They are sorted by start position and then concatenated. If empty/not set"
			"there will be no concatenation performed", 
			seqan::ArgParseOption::STRING));
	setDefaultValue(parser, "parent", "");


    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Extract option values.
	if (isSet(parser, "verbose"))
        options.verbosity = 1;
    if (isSet(parser, "vverbose"))
        options.verbosity = 2;

    getOptionValue(options.roiFileName, parser, "roi-file");
    getOptionValue(options.gffFileName, parser, "gff-file");
    getOptionValue(options.outputFileName, parser, "output-file");
    getOptionValue(options.strandSpecific, parser, "strandspecific");
    getOptionValue(options.type, parser, "type");
    getOptionValue(options.parent, parser, "parent");
	
    return seqan::ArgumentParser::PARSE_OK;

}

// Program entry point
int main(int argc, char const ** argv)
{
    // Parse the command line.
    Options options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.  Otherwise, exit with code 0 (e.g. help
    // was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    return mainWithOptions(options);
 }


