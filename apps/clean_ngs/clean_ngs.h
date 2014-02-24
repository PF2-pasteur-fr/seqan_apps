// ==========================================================================
//                                 Clean NGS
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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

#ifndef SANDBOX_MY_SANDBOX_APPS_ADAPTERREMOVAL_ADAPTERREMOVAL_H_
#define SANDBOX_MY_SANDBOX_APPS_ADAPTERREMOVAL_ADAPTERREMOVAL_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#define VERSION "0.9"

const unsigned int MAX_READ_LEN = 5000;

//Version 0.3: added Memory consumption to log
//Version 0.9 writing reads to rejected when maxLen constrain is met

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------

// Options for the Clean NGS Application.

struct Options
{
    // Verbosity of output: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // Paths to input read files.
    CharString inputFileName1;
    CharString inputFileName2;

    // Paths to output read files of accepted reads.
    CharString outputFileName1;
    CharString outputFileName2;

    // Paths to output read files of rejected reads.
    CharString rejectedFileName1;
    CharString rejectedFileName2;

    // Path to input file with adaptor information.
    CharString adptFileName;

    // Lengths for minimal/maximal length filter.
    int minLen;
    int maxLen;

    // Threshold for quality-based trimming on 5'/3' ends.
    int qualityThreshold5;
    int qualityThreshold3;

    // Whether or not to expect paired-end data, shortcut to !empty(inputFileName2).
    bool pairedEnd() const
    {
        return !empty(inputFileName2);
    }

    Options() :
        verbosity(0), minLen(0), maxLen(0), qualityThreshold5(0), qualityThreshold3(0)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

int mainWithOptions(Options const & options);

#endif  // #ifndef SANDBOX_MY_SANDBOX_APPS_ADAPTERREMOVAL_ADAPTERREMOVAL_H_
