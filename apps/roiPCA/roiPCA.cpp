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

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "roiPCA.h"

using namespace seqan;

seqan::ArgumentParser::ParseResult
parseCommandLine(Options & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("roi_pca");

    // Set short description, version, and date.
    setShortDescription(parser, "Analysis of ROI data files.");
    setVersion(parser, "0.2");
    setDate(parser, "March 2013");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fB-if\\fP \\fIIN.roi\\fP \\fB-of\\fP \\fIOUT.roi\\fP");

    addDescription(parser,
                   "Read in the ROI file \\fIIN.roi\\fP, add data columns and write out as \\fIOUT.roi\\fP.");
    // General Options

    addSection(parser, "General Options");

    addOption(parser, seqan::ArgParseOption("v", "verbose", "Verbose mode."));
    addOption(parser, seqan::ArgParseOption("vv", "vverbose", "Very verbose mode."));

    // Input / Output Options

    addSection(parser, "Input / Output Parameters");

    addOption(parser, seqan::ArgParseOption("if", "input-file", "ROI file.", seqan::ArgParseOption::INPUTFILE));
    setValidValues(parser, "input-file", "roi");
    setRequired(parser, "input-file");

    addOption(parser, seqan::ArgParseOption("of", "output-file", "ROI output file", seqan::ArgParseOption::OUTPUTFILE));
    setValidValues(parser, "output-file", "roi");
    setRequired(parser, "output-file");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.

    options.verbosity = 1;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "vverbose"))
        options.verbosity = 3;

    getOptionValue(options.inputFileName, parser, "input-file");
    getOptionValue(options.outputFileName, parser, "output-file");

	return seqan::ArgumentParser::PARSE_OK;
}

// Program entry point.
int main(int argc, char const ** argv)
{
	// Parse the command line.
    Options options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.  Otherwise, exit with code 0 (e.g. help
    // was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    std::cerr << "__OPTIONS_____________________________________________________________________\n"
              << "\n"
              << "INPUT FILE \t" << options.inputFileName << "\n"
              << "OUTPUT FILE\t" << options.outputFileName << "\n"
              << "\n";

    ComputeRoiStatsApp app(options);
    return app.run();
}
