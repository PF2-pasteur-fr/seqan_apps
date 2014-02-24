// ==========================================================================
//                              adaptor_record.h
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#include <seqan/sequence.h>
#include <seqan/stream.h>

#ifndef SANDBOX_JAGLA_APPS_ADAPTERREMOVAL_ADAPTOR_RECORD_H_
#define SANDBOX_JAGLA_APPS_ADAPTERREMOVAL_ADAPTOR_RECORD_H_

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct AdapterRecord
{
	seqan::CharString id;
	seqan::Dna5String seq;
	double similarityThreshold;
	int qualityThreshold;
	int minOverlap;
	bool trimAdapters;
	bool threePrime;
	unsigned int useCount;

	AdapterRecord() :
		similarityThreshold(0), qualityThreshold(0), minOverlap(0),
		trimAdapters(false), threePrime(false), useCount(0)
	{}
};

struct AdapterTsv_;
typedef seqan::Tag<AdapterTsv_> AdapterTsv;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function isHeader()
// ----------------------------------------------------------------------------

template <typename TStream>
bool isHeader(seqan::RecordReader<TStream, seqan::SinglePass<> > & reader,
              AdapterTsv const & /*tag*/)
{
    return value(reader) == '#';
}

// ----------------------------------------------------------------------------
// Function isRecord()
// ----------------------------------------------------------------------------

template <typename TStream>
bool isRecord(seqan::RecordReader<TStream, seqan::SinglePass<> > & reader,
              AdapterTsv const & /*tag*/)
{
    return value(reader) != '#';
}

// ----------------------------------------------------------------------------
// Function skipHeader()
// ----------------------------------------------------------------------------

template <typename TStream>
int skipHeader(seqan::RecordReader<TStream, seqan::SinglePass<> > & reader,
               AdapterTsv const & /*tag*/)
{
    using namespace seqan;
    
    int ret = skipLine(reader);
    if (ret != 0 && ret != EOF_BEFORE_SUCCESS)
    {
        std::cerr << "Error reading header\n";
        return 1;
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function readRecord()
// ----------------------------------------------------------------------------

template <typename TStream>
int readRecord(AdapterRecord & record,
               seqan::RecordReader<TStream, seqan::SinglePass<> > & reader,
               AdapterTsv const & /*tag*/)
{
    using namespace seqan;

    SEQAN_ASSERT(isRecord(reader, AdapterTsv()));

    CharString buffer;
    int res = 0;

    clear(record.id);
    clear(record.seq);

    // Read id.
    res = readUntilOneOf(record.id, reader, '\t', '\r', '\n');
    if (res != 0)
    {
        std::cerr << "Could not read ID field!\n";
        return 1;
    }
    if (value(reader) != '\t')
    {
        std::cerr << "No TAB after ID field!\n";
        return 1;
    }
    skipChar(reader, '\t');

    // Read sequence.
    res = readUntilOneOf(record.seq, reader, '\t', '\r', '\n');
    if (res != 0)
    {
        std::cerr << "Could not read SEQ field!\n";
        return 1;
    }
    if (value(reader) != '\t')
    {
        std::cerr << "No TAB after SEQUENCE field!\n";
        return 1;
    }
    skipChar(reader, '\t');

    // Read similarity threshold.
    clear(buffer);
    res = readUntilOneOf(buffer, reader, '\t', '\r', '\n');
    if (res != 0)
    {
        std::cerr << "Could not read SIMILARITY THRESHOLD field!\n";
        return 1;
    }
    if (value(reader) != '\t')
    {
        std::cerr << "No TAB after SIMILARITY THRESHOLD field!\n";
        return 1;
    }
    skipChar(reader, '\t');
    if (!lexicalCast2(record.similarityThreshold, buffer))
    {
        std::cerr << "Could not cast SIMILARITY THRESHOLD field to floating point number!\n";
        return 1;
    }

    // Read quality threshold.
    clear(buffer);
    res = readUntilOneOf(buffer, reader, '\t', '\r', '\n');
    if (res != 0)
    {
        std::cerr << "Could not read QUALITY THRESHOLD field!\n";
        return 1;
    }
    if (value(reader) != '\t')
    {
        std::cerr << "No TAB after QUALITY THRESHOLD field!\n";
        return 1;
    }
    skipChar(reader, '\t');
    if (!lexicalCast2(record.qualityThreshold, buffer))
    {
        std::cerr << "Could not cast QUALITY THRESHOLD field to integer!\n";
        return 1;
    }

    // Read min overlap.
    clear(buffer);
    res = readUntilOneOf(buffer, reader, '\t', '\r', '\n');
    if (res != 0)
    {
        std::cerr << "Could not read MIN OVEFRLAP field!\n";
        return 1;
    }
    if (value(reader) != '\t')
    {
        std::cerr << "No TAB after MIN OVERLAP field!\n";
        return 1;
    }
    skipChar(reader, '\t');
    if (!lexicalCast2(record.minOverlap, buffer))
    {
        std::cerr << "Could not cast MIN OVERLAP field to integer!\n";
        return 1;
    }

    // Read trim adaptors flag.
    char c = value(reader);
    if (c != '0' && c != '1' && c != 'Y' && c != 'N' && c != 'y' && c != 'n')
    {
        std::cerr << "TRIM ADAPTERS field has invalid value: " << c << "\n";
        return 1;
    }
    record.trimAdapters = (c == '1' || c == 'Y' || c == 'y');
    goNext(reader);
    res = skipChar(reader, '\t');
    if (res != 0)
    {
        std::cerr << "No TAB after TRIM ADAPTERS field!\n";
        return 1;
    }

    // Read direction left flag.
    // 0 = normal adapter
    // 1 = leader sequence
    c = value(reader);
    if (c != '0' && c != '1')
    {
        std::cerr << "THREE PRIME field has invalid value: " << c << "\n";
        return 1;
    }
    record.threePrime = (c == '0');
    goNext(reader);
    if (!atEnd(reader) && value(reader) != '\n' && value(reader) != '\r')
    {
        std::cerr << "Trailing characters!\n";
        return 1;
    }
    skipLine(reader);

    return 0;
}

// ----------------------------------------------------------------------------
// Function operator<<()
// ----------------------------------------------------------------------------

template <typename TStream>
TStream & operator<<(TStream & stream, AdapterRecord const & record)
{
    stream << record.id << '\t' << record.seq << '\t' << record.similarityThreshold << '\t'
           << record.qualityThreshold << '\t' << record.minOverlap << '\t'
           << record.trimAdapters << '\t' << record.threePrime;
    
    return stream;
}

#endif  // #ifndef SANDBOX_JAGLA_APPS_ADAPTERREMOVAL_ADAPTOR_RECORD_H_
