// ==========================================================================
//                                 Clean NGS
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#include <fstream>
#include <iostream>

#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/stream.h>

#include "clean_ngs.h"

// Parse command line with ArgumentParser class.

seqan::ArgumentParser::ParseResult
parseCommandLine(Options & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("clean_ngs");

    // Set short description, version, and date.
    setShortDescription(parser, "NGS Data Cleaning");
    setVersion(parser, VERSION);
    setDate(parser, "March 2013");

    // Define usage line and long description.
    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] \\fB-adf\\fP \\fIADAPTER.txt\\fP \\fB-if1\\fP \\fIIN.fq\\fP "
                 "\\fB-of1\\fP \\fIOUT.fq\\fP");
    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] \\fB-adf\\fP \\fIADAPTER.txt\\fP \\fB-if1\\fP \\fIIN_1.fq\\fP \\fB-if2\\fP \\fIIN_2.fq\\fP"
                 "\\fB -of1\\fP \\fIOUT_1.fq\\fP \\fB-of2\\fP \\fIOUT_2.fq\\fP");

    addDescription(parser,
                   "Use adapters definition file \\fIADAPTER.txt\\fP and remove the adapters from the input read"
                   "files.  The program writes out the cleansed reads and the rejected reads.");

    // General Options

    addOption(parser, seqan::ArgParseOption("v", "verbose", "Verbose mode."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Very verbose mode."));

    // Input / Output Options

    addSection(parser, "Input / Output Parameters");

    addOption(parser, seqan::ArgParseOption("adf", "adapter-file",
                                            "Adapter definition file.  See section Adapters File Format below for more"
                                            "information.", seqan::ArgParseOption::INPUTFILE));
    setValidValues(parser, "adapter-file", "txt dat");
    setRequired(parser, "adapter-file");

    char const * SEQ_FORMATS = "fastq fq fastq.gz fq.gz txt txt.gz dat";

    addOption(parser, seqan::ArgParseOption("if1", "input-file1", "Input file 1.", seqan::ArgParseOption::INPUTFILE));
    setValidValues(parser, "input-file1", SEQ_FORMATS);
    setRequired(parser, "input-file1");

    addOption(parser, seqan::ArgParseOption("if2", "input-file2", "Input file 2.", seqan::ArgParseOption::INPUTFILE));
    setValidValues(parser, "input-file2", SEQ_FORMATS);

    addOption(parser, seqan::ArgParseOption("of1", "output-file1", "Output file 1.", seqan::ArgParseOption::OUTPUTFILE));
    setValidValues(parser, "output-file1", SEQ_FORMATS);
    setRequired(parser, "output-file1");

    addOption(parser, seqan::ArgParseOption("of2", "output-file2", "Output file 2.", seqan::ArgParseOption::OUTPUTFILE));
    setValidValues(parser, "output-file2", SEQ_FORMATS);

    addOption(parser, seqan::ArgParseOption("rf1", "rejected-file1", "Rejected file 1.", seqan::ArgParseOption::OUTPUTFILE));
    setValidValues(parser, "rejected-file1", SEQ_FORMATS);
    setRequired(parser, "rejected-file1");

    addOption(parser, seqan::ArgParseOption("rf2", "rejected-file2", "Rejected file 2.", seqan::ArgParseOption::OUTPUTFILE));
    setValidValues(parser, "rejected-file2", SEQ_FORMATS);

    // Optional Parameter

    addSection(parser, "Optional Parameters");

    addOption(parser, seqan::ArgParseOption("l", "min-len",
                                            "Minimal length of sequencing.Set to 0 to turn off"
                                            "minimal length filter.", seqan::ArgParseOption::INTEGER));
    setMinValue(parser, "min-len", "0");
    setDefaultValue(parser, "min-len", "0");

    addOption(parser, seqan::ArgParseOption("L", "max-len",
                                            "Maximal length of sequencing.Set to 0 to turn off"
                                            "maximal length filter.", seqan::ArgParseOption::INTEGER));
    setMinValue(parser, "max-len", "0");
    setDefaultValue(parser, "max-len", "0");

    addOption(parser, seqan::ArgParseOption("q3", "qual-thresh3",
                                            "Quality threshold for 3-prime trimming.Set to"
                                            " 0 to turn off 3-prime trimming.", seqan::ArgParseOption::INTEGER));
    setMinValue(parser, "qual-thresh3", "0");
    setDefaultValue(parser, "qual-thresh3", "0");

    addOption(parser, seqan::ArgParseOption("q5", "qual-thresh5",
                                            "Quality threshold for 5-prime trimming.Set to"
                                            " 0 to turn off 5-prime trimming.", seqan::ArgParseOption::INTEGER));
    setMinValue(parser, "qual-thresh5", "0");
    setDefaultValue(parser, "qual-thresh5", "0");

    // Adapter Definition File Format Section

    addTextSection(parser, "Adapter Definition File Format");
    addText(parser, "The adapter sequences are defined in an adapter definition file with the extension \".txt\" (dat is only allowed for Galaxy)");
    addText(parser, "comment lines start with \"#\".");
    addText(parser, "Each non-comment(non-empty) line is tab delimited");
    addText(parser, "and has 7 fields (columns):");
	addText(parser, "1. name of the sequence");
    addText(parser, "2. sequence (N's are allowed and match any nucleotide, useful for multiplexing)");
    addText(parser, "3. fraction of identity between adapter sequence and target sequence (threshold)");
    addText(parser, "4. quality threshold");
    addText(parser, "5. minimum length of overlap between adapter and target sequence");
    addText(parser, "6. if set to >0 perform comparison with truncated adapter sequences (very time consuming)");
    addText(parser, "7. if set to 1 this is a leader sequence");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.

    if (isSet(parser, "verbose"))
        options.verbosity = 1;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 2;

    getOptionValue(options.adptFileName, parser, "adapter-file");
    getOptionValue(options.inputFileName1, parser, "input-file1");
    getOptionValue(options.inputFileName2, parser, "input-file2");
    getOptionValue(options.outputFileName1, parser, "output-file1");
    getOptionValue(options.outputFileName2, parser, "output-file2");
    getOptionValue(options.rejectedFileName1, parser, "rejected-file1");
    getOptionValue(options.rejectedFileName2, parser, "rejected-file2");

    // TODO(holtgrew): Check that if input-file2 is given that output-file2 is given?
    // TODO(holtgrew): Check that if paired input and output and rejected-file1 is given then rejected-file2 is given?

    getOptionValue(options.minLen, parser, "min-len");
    getOptionValue(options.maxLen, parser, "max-len");
    getOptionValue(options.qualityThreshold3, parser, "qual-thresh3");
    getOptionValue(options.qualityThreshold5, parser, "qual-thresh5");

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
