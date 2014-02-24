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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Implementation of the NGS Cleaning.
// ==========================================================================

#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>

#include "clean_ngs.h"

#include "adaptor_record.h"

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Record
// ----------------------------------------------------------------------------

// Stores information (sequence, qualities) for left and right read such that they can be treated at the same time more
// easily.  This is important for rejecting both reads if one fulfills the rejection criteria.

struct Record
{
    // Identifiers.
    CharString id1, id2;

    // Sequence.
    Dna5String seq1, seq2;

    // Qualities.
    CharString qual1, qual2;

    // Original sequence.
    Dna5String orgseq1, orgseq2;

    // Original qualities.
    CharString orgqual1, orgqual2;

    Record() {}
};

// ----------------------------------------------------------------------------
// Class Stats
// ----------------------------------------------------------------------------

struct Stats
{
    // TODO(holtgrew): Use a String/std::vector.
    // Histogram over read lengths.
    int length[MAX_READ_LEN];

    // Number of reads trimmed on 3'/5' end.
    int count3PrimeTrimmed;
    int count5PrimeTrimmed;

    // Number of reads rejected due to length constraints.
    int lengthConstrain;
    int maxLenConstrain;

    // Counters for good and bad reads.
    int goodReads;
    int badReads;

    Stats() :
        count3PrimeTrimmed(0), count5PrimeTrimmed(0), lengthConstrain(0), maxLenConstrain(0),
        goodReads(0), badReads(0)
    {
        for (int i = 0; i < (int) MAX_READ_LEN; i++)
            length[i] = 0;
    }

};

// ----------------------------------------------------------------------------
// Class CleanNgsApp
// ----------------------------------------------------------------------------

// Application class with global state.

class CleanNgsApp
{
public:
    CleanNgsApp(Options const & options) :
        options(options)
    {}

    // Run the cleaning process.
    int run();

private:
    Options const & options;

    // Statistics.
    Stats stats;

    // Adapter records read from options.adptFileName.
    String<AdapterRecord> adptRecords;
    // Adapter records for forward and reverse strand.
    String<AdapterRecord> ffRecords;
    String<AdapterRecord> rvRecords;

    // Streams for writing out good and rejected reads.
    seqan::SequenceStream outLeftGood, outLeftBad, outRightGood, outRightBad;

    // Print program header to out.
    void printHeader(std::ostream & out);

    // Open output files out{Left,Right}{Good,Bad}.
    int openOutputFiles();

    // Perform the actual trimming.
    int performTrimming(double programStartTime);

    // Print program footer to out.
    void printFooter(std::ostream & out, double programStartTime);

    // Read adapter records from path.
    int readAdapterRecords(char const * path);

    // Write reads to rejected files.
    // TODO(holtgrew): Make return void.
    bool removeRead(Record const & record);

    // Take a record with single-end or paired-end data and perform 5' quality trimming.
    bool qualityThreshold(Record & record);

    // TODO(holtgrew): Can we collapse qualityThreshold and qualityThresholdRev into one function using templates and Reversed Strings?

    // Take a record with single-end or paired-end data and perform 3' quality trimming.
    bool qualityThresholdRev(Record & record);

    // Perform actual adapter removal.
    bool adapterRemoval(Dna5String & seq1, CharString & qual1, Record record,
                        String<AdapterRecord> & ffRecords, bool rmRead);
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function mainWithOptions()
// ----------------------------------------------------------------------------

int mainWithOptions(Options const & options)
{
    CleanNgsApp app(options);
    return app.run();
}

// ----------------------------------------------------------------------------
// Member Function CleanNgsApp::removeRead()
// ----------------------------------------------------------------------------

bool CleanNgsApp::removeRead(Record const & record)
{
    //=============================
    // write bad
    //=============================
    if (!empty(options.rejectedFileName1))
        if (writeRecord(outLeftBad, record.id1, record.orgseq1, record.orgqual1) != 0)
        {
            std::cerr << "Could not write left bad read!\n" << " " << record.id1
                      << " " << record.orgseq1 << " " << record.orgqual1 << "\n";
            return 1;
        }
    if (!empty(options.rejectedFileName2))
        if (options.pairedEnd() && writeRecord(outRightBad, record.id2, record.orgseq2, record.orgqual2) != 0)
        {
            std::cerr << "Could not write right bad read!\n";
            return 1;
        }
    if (options.verbosity >= 2)
    {
        std::cout << "removed read\n\n";
    }
    stats.badReads++;
    return true;
}

// ----------------------------------------------------------------------------
// Member Function CleanNgsApp::qualityThreshold()
// ----------------------------------------------------------------------------

//we could take the pairedEnd part out.... but I think it is OK like that as it is still quite small

bool CleanNgsApp::qualityThreshold(Record & record)
{
    bool retVal = false;

    int startPos = 0, endPos = length(record.qual1);
    while (startPos < endPos && record.qual1[startPos] < options.qualityThreshold5)
        startPos++;

    if (startPos > 0)
    {
        stats.count5PrimeTrimmed++;
        record.seq1 = infix(record.seq1, startPos, endPos);
        record.qual1 = infix(record.qual1, startPos, endPos);
    }

    if (length(record.seq1) < (unsigned)(options.minLen))
    {
        stats.lengthConstrain++;
        return removeRead(record);
    }

    if (options.pairedEnd())
    {
        int startPos = 0, endPos = length(record.qual2);
        while (startPos < endPos && record.qual2[startPos] < options.qualityThreshold5)
            startPos++;

        if (startPos > 0)
        {
            stats.count5PrimeTrimmed++;
            record.seq2 = infix(record.seq2, startPos, endPos);
            record.qual2 = infix(record.qual2, startPos, endPos);
        }

        if (length(record.seq2) < (unsigned)(options.minLen))
        {
            stats.lengthConstrain++;
            return removeRead(record);
        }
    }

    return retVal;
}

// ----------------------------------------------------------------------------
// Member Function CleanNgsApp::qualityThresholdRev()
// ----------------------------------------------------------------------------

// TODO(berndj): cant we use reverse string????

bool CleanNgsApp::qualityThresholdRev(Record & record)
{
    bool retVal = false;

    int startPos = length(record.qual1) - 1;
    while (startPos >= 0 && record.qual1[startPos] < options.qualityThreshold3)
        startPos--;

    if ((unsigned) startPos < length(record.qual1) - (unsigned)1)
    {
        stats.count3PrimeTrimmed++;
        record.seq1 = infix(record.seq1, 0, startPos + 1);
        record.qual1 = infix(record.qual1, 0, startPos + 1);
    }

    if (length(record.seq1) < (unsigned)(options.minLen))
    {
        stats.lengthConstrain++;
        return removeRead(record);
    }

    if (options.pairedEnd())
    {
        int startPos = length(record.qual2) - 1;
        while (startPos >= 0 && record.qual2[startPos] < options.qualityThreshold3)
            startPos--;

        if ((unsigned) startPos < length(record.qual2) - 1)
        {
            stats.count3PrimeTrimmed++;
            record.seq2 = infix(record.seq2, 0, startPos);
            record.qual2 = infix(record.qual2, 0, startPos);
        }

        if (length(record.seq2) < (unsigned)(options.minLen))
        {
            stats.lengthConstrain++;
            return removeRead(record);
        }
    }

    return retVal;
}

// ----------------------------------------------------------------------------
// Member Function CleanNgsApp::adapterRemoval()
// ----------------------------------------------------------------------------

bool CleanNgsApp::adapterRemoval(Dna5String & seq1, CharString & qual1, Record record,
                                 String<AdapterRecord> & ffRecords, bool rmRead)
{
    int seqLen = length(seq1);
    int minSeqIdx = seqLen;
    bool usedAdpt = false;
    //loop over all adapters
    for (Iterator<String<AdapterRecord>, Rooted>::Type adpIter = begin(
             ffRecords); !atEnd(adpIter); goNext(adpIter))
    {
        seqLen = length(seq1);
        minSeqIdx = seqLen;

        // loop over all adapters
        int alen = length(adpIter->seq);
        Dna5String adapter = adpIter->seq;
        int partAdptLen = 0;
        if (adpIter->trimAdapters)
        {
            partAdptLen = alen - adpIter->minOverlap + 1;
        }
        // trim adapter (one iteration if no trimming
        for (int adpI = 0; adpI <= partAdptLen; adpI++)
        {
            // loop over all partial adapters
            Dna5String adapterCompStr = infix(adapter, adpI, length(adapter));
            alen = length(adapterCompStr);
            // Compare to adapter to sequence by sliding over read sequence
            // we add 1 because if minLen is set to 1 without it there would
            // be minimum of 2 bases being compared
            for (int seqIdx = 0; seqIdx < seqLen - adpIter->minOverlap + 1;
                 seqIdx++)
            {
                Dna5String seqPart = infix(seq1, seqIdx, length(seq1));
                CharString qual1A = infix(qual1, seqIdx, length(seq1));
                int pSim = 0;
                for (int seqPartIdx = 0;
                     (unsigned) (seqPartIdx) < length(seqPart);
                     seqPartIdx++)
                {
                    // adapter is shorter than sequence
                    if (seqPartIdx >= alen)
                    {
                        break;
                    }
                    // adapter can be compared with sequence
                    char qualC = qual1A[seqPartIdx];
                    if (seqPart[seqPartIdx] == adapterCompStr[seqPartIdx]
                       || qualC <= adpIter->qualityThreshold
                       || seqPart[seqPartIdx] == 'N')
                    {
                        pSim++;
                    }
                }

                int complen = alen;
                if (length(seqPart) < (unsigned) (alen))
                {
                    complen = length(seqPart);
                }
                //std::cout << adpIter->seq << " " << seq1 << " " << seqIdx << "\n";
                if (adpIter->similarityThreshold < pSim / (double) (complen))
                {
                    //std::cout << "Found" << "\n";
                    if (seqIdx < minSeqIdx)
                    {
                        minSeqIdx = seqIdx;
                        usedAdpt = true;
                        adpIter->useCount++;
                        if (minSeqIdx == 0)
                        {
                            rmRead = removeRead(record);
                            return rmRead;
                        }
                    }

                    break;
                }

            }

        }

        seq1 = infix(seq1, 0, minSeqIdx);
        qual1 = infix(qual1, 0, minSeqIdx);
    }
    if (length(seq1) < (unsigned) (options.minLen))
    {
        stats.lengthConstrain++;
        rmRead = removeRead(record);
    }
    return rmRead;
}

// ----------------------------------------------------------------------------
// Member Function CleanNgsApp::printHeader()
// ----------------------------------------------------------------------------

void CleanNgsApp::printHeader(std::ostream & out)
{
    if (options.verbosity >= 1)
    {
        out << "**********************\n"
            << "* adapterRemoval     *\n"
            << "**********************\n\n";
        out << "_____ARGUMENTS_____\n" << "\n";
        out << "adaptors file:  " << options.adptFileName << "\n";
        if (!empty(options.inputFileName2))
            out << "mode:           paired-end\n";
        else
            out << "mode:           single-end\n";
        out << "input left:     " << options.inputFileName1 << "\n";
        if (!empty(options.inputFileName2))
            out << "input right:    " << options.inputFileName2 << "\n";
        out << "output left:    " << options.outputFileName1 << "\n";
        if (!empty(options.outputFileName2))
            out << "output right:   " << options.outputFileName2 << "\n";
        if (!empty(options.rejectedFileName1))
            out << "rejected left:   " << options.rejectedFileName1
                << "\n";
        if (!empty(options.rejectedFileName2))
            out << "rejected right:   " << options.rejectedFileName2
                << "\n";
        out << "minLen:         " << options.minLen << "\n";
        out << "maxLen:         " << options.maxLen << "\n";
        out << "Quality threshold 3 prime:  " << options.qualityThreshold3
            << "\n";
        out << "Quality threshold 5 prime:  " << options.qualityThreshold5
            << "\n";
        out << "\n";
    }
}

// ----------------------------------------------------------------------------
// Member Function CleanNgsApp::openOutputFiles()
// ----------------------------------------------------------------------------

int CleanNgsApp::openOutputFiles()
{
    open(outLeftGood, toCString(options.outputFileName1), seqan::SequenceStream::WRITE, seqan::SequenceStream::FASTQ);
    if (!isGood(outLeftGood))
    {
        std::cerr << "Could not open file " << options.outputFileName1 << "\n";
        return 1;
    }

    if (!empty(options.inputFileName2))
    {
        open(outRightGood, toCString(options.outputFileName2), seqan::SequenceStream::WRITE, seqan::SequenceStream::FASTQ);
        if (!isGood(outRightGood))
        {
            std::cerr << "Could not open file " << options.outputFileName2 << "\n";
            return 1;
        }
    }

    if (!empty(options.rejectedFileName1))
    {
        open(outLeftBad, toCString(options.rejectedFileName1), seqan::SequenceStream::WRITE, seqan::SequenceStream::FASTQ);
        if (!isGood(outLeftBad))
        {
            std::cerr << "Could not open file " << options.rejectedFileName1 << "\n";
            return 1;
        }
    }

    if (!empty(options.rejectedFileName2))
    {
        open(outRightBad, toCString(options.rejectedFileName2), seqan::SequenceStream::WRITE, seqan::SequenceStream::FASTQ);
        if (!isGood(outRightBad))
        {
            std::cerr << "Could not open file " << options.rejectedFileName2 << "\n";
            return 1;
        }
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Member Function CleanNgsApp::readAdapterRecords()
// ----------------------------------------------------------------------------

int CleanNgsApp::readAdapterRecords(char const * path)
{
    std::ifstream adptIn(path, std::ios::binary | std::ios::in);
    if (!adptIn.good())
    {
        std::cerr << "Error: couldn't open Adapter file" << toCString(options.adptFileName) << "\n";
        return 1;
    }

    RecordReader<std::ifstream, SinglePass<> > adapterReader(adptIn);
    while (!atEnd(adapterReader))
    {
        int ret = 0;
        AdapterRecord record;

        // Skipping comments.
        if (isHeader(adapterReader, AdapterTsv()))
        {
            ret = skipHeader(adapterReader, AdapterTsv());
            if (ret != 0)
            {
                std::cerr << "ret!=0\n";
                return 1;
            }
        }
        else if (value(adapterReader) == '\n' || value(adapterReader) == '\r')
        {
            ret = skipLine(adapterReader);
            if (ret != 0)
            {
                std::cerr << "ERROR skipping empty line, should not happen!\n";
            }
            //return 1;
        }
        else
        {
            ret = readRecord(record, adapterReader, AdapterTsv());
            if (ret != 0)
            {
                std::cerr << "2ret!=0\n";
                return 1;
            }
            appendValue(adptRecords, record);
            if (options.verbosity >= 2)
                std::cerr << "RECORD: " << record << "\n";
        }
    }

    // Split adapters in 5' and 3' starting (normal adapter / leader sequence)

    for (Iterator<String<AdapterRecord>, Rooted>::Type it = begin(adptRecords); !atEnd(it); goNext(it))
    {
        if (it->threePrime)
        {
            appendValue(ffRecords, *it);
        }
        else
        {
            reverse(it->seq);
            appendValue(rvRecords, *it);
        }
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Member Function CleanNgsApp::performTrimming()
// ----------------------------------------------------------------------------

int CleanNgsApp::performTrimming(double programStartTime)
{
    seqan::SequenceStream inStream1(toCString(options.inputFileName1));
    if (!isGood(inStream1))
    {
        std::cerr << "Couldn't open file: " << options.inputFileName1 << "\n";
        return 1;
    }
    seqan::SequenceStream inStream2;
    if (options.pairedEnd())
    {
        open(inStream2, toCString(options.inputFileName2));
        if (!isGood(inStream2))
        {
            std::cerr << "Couldn't open file: " << options.inputFileName2 << "\n";
            return 1;
        }
    }

    // Read FASTA file and output "$id\t$seq".
    int readCounter = 0;
    Record record;
    while (!atEnd(inStream1))
    {
        readCounter++;
        if (!(readCounter % 100000))
        {
            if (readCounter == 100000)
            {
                std::cerr
                << "N(seq)\ttime[sec]\tRSS[kB]\tShared Memory[kB]\tPrivate Memory[kB]\n";
            }
            std::cerr << readCounter << "\t" << sysTime() - programStartTime << "\t";
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
        // Test modus
        /*
         if (readCounter > 5)
         break;
         */

        if (readRecord(record.id1, record.seq1, record.qual1, inStream1)
            != 0)
        {
            std::cerr << "ERROR reading FASTQ." << std::endl;
            return 1;
        }
        record.orgseq1 = record.seq1;
        record.orgseq2 = record.seq2;
        record.orgqual1 = record.qual1;
        record.orgqual2 = record.qual2;
        if (options.pairedEnd())
        {
            if (readRecord(record.id2, record.seq2, record.qual2, inStream2) != 0)
            {
                std::cerr << "ERROR reading 2nd FASTQ file." << std::endl;
                return 1;
            }
        }
        bool rmRead = false;

        if (options.verbosity >= 2)
        {
            std::cout << "Org:" << record.id1 << "\t" << record.seq1 << "\t"
                      << record.qual1 << "\n";
            if (options.pairedEnd())
                std::cout << "Org:" << "P" << record.id2 << "\t" << record.seq2
                          << "\t" << record.qual2 << "\n";
        }
        // Simply check whether read count balanced in case of paired-end data.
        // TODO ?? check ids??
        if (options.pairedEnd() && length(record.seq1) != length(record.seq2))
        {
            std::cerr << "ERROR: Unbalanced read files. "
                      << options.inputFileName1 << " has " << length(record.seq1)
                      << " alignments, " << options.inputFileName2 << " has "
                      << length(record.seq2) << " alignments.\n";
            return 1;
        }

        // =========================
        // 5 Prime quality trimming
        // =========================
        if (options.qualityThreshold5 > 0)
            if (qualityThreshold(record))
                continue;

        if (options.verbosity >= 2)
        {
            std::cout << "5pr:" << record.id1 << "\t" << record.seq1 << "\t"
                      << record.qual1 << "\n";
            if (options.pairedEnd())
                std::cout << "5pr:" << "P" << record.id2 << "\t" << record.seq2
                          << "\t" << record.qual2 << "\n";
        }
        // =========================
        // 3 Prime quality trimming
        // =========================
        if (options.qualityThreshold3 > 0)
            if (qualityThresholdRev(record))
                continue;

        if (options.verbosity >= 2)
        {
            std::cout << "3pr:" << record.id1 << "\t" << record.seq1 << "\t"
                      << record.qual1 << "\n";
            if (options.pairedEnd())
                std::cout << "3pr:" << "P" << record.id2 << "\t" << record.seq2
                          << "\t" << record.qual2 << "\n";
        }

        int seqLen = length(record.seq1);
        if (seqLen == 0)
            continue;  // shouldn't happen
        // =========================
        // Forward adapter (regular, non leader sequences)
        // seq1
        // =========================
        if (adapterRemoval(record.seq1, record.qual1, record, ffRecords,
                           rmRead))
        {
            continue;
        }

        // =========================
        // Forward adapter (regular, non leader sequences)
        // seq2
        // =========================
        if (options.pairedEnd())
        {
            int seqLen = length(record.seq2);
            if (seqLen == 0)
                continue;  // shouldn't happen
            if (adapterRemoval(record.seq2, record.qual2, record, ffRecords,
                               rmRead))
            {
                continue;
            }

        } //paired end
        if (options.verbosity >= 2)
        {
            std::cout << "AdF:" << record.id1 << "\t" << record.seq1 << "\t"
                      << record.qual1 << "\n";
            if (options.pairedEnd())
                std::cout << "AdF:" << "P" << record.id2 << record.seq2 << "\t"
                          << record.qual2 << "\n";
        }

        // =========================
        // Reverse adapter (leader sequences)
        // seq1
        // =========================
        reverse(record.seq1);
        reverse(record.qual1);
        if (adapterRemoval(record.seq1, record.qual1, record, rvRecords,
                           rmRead))
        {
            reverse(record.seq1);
            reverse(record.qual1);
            continue;
        }
        reverse(record.seq1);
        reverse(record.qual1);

        // =========================
        // Reverse adapter (leader sequences)
        // seq2
        // =========================
        if (options.pairedEnd())
        {
            int seqLen = length(record.seq2);
            if (seqLen == 0)
                continue;  // shouldn't happen
            reverse(record.seq2);
            reverse(record.qual2);
            if (adapterRemoval(record.seq2, record.qual2, record, rvRecords,
                               rmRead))
            {
                reverse(record.seq2);
                reverse(record.qual2);
                continue;
            }
            reverse(record.seq2);
            reverse(record.qual2);

        } //paired end
        if (options.verbosity >= 2)
        {
            std::cout << "AdR:" << record.id1 << "\t" << record.seq1 << "\t"
                      << record.qual1 << "\n";
            if (options.pairedEnd())
                std::cout << "AdR:" << "P" << record.id2 << "\t" << record.seq2
                          << "\t" << record.qual2 << "\n";
        }

        //=============================
        // last but not least check MaxLen if set
        //=============================
        if ((unsigned) options.maxLen > 0
           && length(record.seq1) > (unsigned) options.maxLen)
        {
            stats.maxLenConstrain++;
            removeRead(record);
            continue;
            if (options.pairedEnd() && length(record.seq2) > (unsigned) options.maxLen)
            {
                removeRead(record);
                stats.maxLenConstrain++;
                continue;
            }

        }

        //=============================
        // final result
        //=============================
        if (options.verbosity >= 2)
        {
            std::cout << "lst:" << record.id1 << "\t" << record.seq1 << "\t"
                      << record.qual1 << "\n";
            if (options.pairedEnd())
                std::cout << "Lst:" << "P" << record.id2 << "\t" << record.seq2
                          << "\t" << record.qual2 << "\n";
            std::cout << "\n";
        }
        //=============================
        // write good
        //=============================
        if (writeRecord(outLeftGood, record.id1, record.seq1, record.qual1) != 0)
        {
            std::cerr << "Could not write left good read!\n";
            return 1;
        }
        if (options.pairedEnd()
           && writeRecord(outRightGood, record.id2, record.seq2, record.qual2) != 0)
        {
            std::cerr << "Could not write right good read!\n";
            return 1;
        }
        stats.length[length(record.seq1)]++;
        stats.goodReads++;

    }

    return 0;
}

// ----------------------------------------------------------------------------
// Member Function CleanNgsApp::printFooter()
// ----------------------------------------------------------------------------

void CleanNgsApp::printFooter(std::ostream & out, double programStartTime)
{
    out << ">>Parameters\n";
    out << "#Name\tvalue" << "\n";
    out << "Version\t" << VERSION << "\n";
    out << "adaptors file:\t" << options.adptFileName << "\n";
    if (!empty(options.inputFileName2))
        out << "mode:\tpaired-end\n";
    else
        out << "mode:\tsingle-end\n";
    out << "input left:\t" << options.inputFileName1 << "\n";
    if (!empty(options.inputFileName2))
        out << "input right:\t" << options.inputFileName2 << "\n";
    out << "output left:\t" << options.outputFileName1 << "\n";
    if (!empty(options.outputFileName2))
        out << "output right:\t" << options.outputFileName2 << "\n";
    if (!empty(options.rejectedFileName1))
        out << "rejected left:\t" << options.rejectedFileName1 << "\n";
    if (!empty(options.rejectedFileName2))
        out << "rejected right:\t" << options.rejectedFileName2 << "\n";
    out << "minLen:\t" << options.minLen << "\n";
    out << "maxLen:\t" << options.maxLen << "\n";
    out << "Quality threshold 3 prime:\t" << options.qualityThreshold3
        << "\n";
    out << "Quality threshold 5 prime:\t" << options.qualityThreshold5
        << "\n";
    out << ">>END_MODULE" << "\n";

    out << ">>Statistics\n";
    out << "#Measure\tvalue" << "\n";
    out << "time\t" << sysTime() - programStartTime << "\n";
    out << "good reads\t" << stats.goodReads << "\n";
    out << "bad reads\t" << stats.badReads << "\n";
    out << "trimmed 3 prime\t" << stats.count3PrimeTrimmed << "\n";
    out << "trimmed 5 prime\t" << stats.count5PrimeTrimmed << "\n";
    out << "removed due to length\t" << stats.lengthConstrain << "\n";
    out << "removed due to max length\t" << stats.maxLenConstrain << "\n";
    out << ">>END_MODULE" << "\n";

    out << ">>sequence length distribution\n";
    out << "#seq length\tcount" << "\n";
    for (unsigned int i = 0; i < MAX_READ_LEN; i++)
        if (stats.length[i] > 0)
            out << i << "\t" << stats.length[i] << "\n";
    out << ">>END_MODULE" << "\n";
    out << ">>adapter sequences" << "\n";
    out << "#adapter ID\tcount" << "\n";
    for (Iterator<String<AdapterRecord>, Rooted>::Type it = begin(ffRecords); !atEnd(it); goNext(it))
        out << it->id << "\t" << it->useCount << "\n";
    out << ">>END_MODULE" << "\n";

    out << ">>leader sequences\n";
    out << "#adapter ID\tcount" << "\n";
    for (Iterator<String<AdapterRecord>, Rooted>::Type it = begin(rvRecords); !atEnd(it); goNext(it))
        out << it->id << "\t" << it->useCount << "\n";
    out << ">>END_MODULE" << "\n";
}

// ----------------------------------------------------------------------------
// Member Function CleanNgsApp::run()
// ----------------------------------------------------------------------------

int CleanNgsApp::run()
{
    double programStartTime = sysTime();

    // Print arguments in verbose mode.
    printHeader(std::cout);

    // Open output files.
    if (openOutputFiles() != 0)
        return 1;

    // Read Adapter Records.
    if (readAdapterRecords(toCString(options.adptFileName)) != 0)
        return 1;

    // Perform the actual trimming.
    if (performTrimming(programStartTime) != 0)
        return 1;

    // Stats
    printFooter(std::cout, programStartTime);

    return 0;
}
