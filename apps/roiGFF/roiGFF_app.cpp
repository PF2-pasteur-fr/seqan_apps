// ==========================================================================
//                                   class roiGFF
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================


#include "../bam2roi/roi.h"
#include "roiGFF.h"


// ----------------------------------------------------------------------------
// Class CleanNgsApp
// ----------------------------------------------------------------------------

// Application class with global state.

class roiGFFApp
{
public:

     roiGFFApp(Options const & options) : gffCount(0), options(options), globalMax(0), nRoi(0)
    {}

    // combine GFF and ROI.
    int run();

private:
	Options const & options;

	seqan::StringSet<ROI> roiSet;

	// FRAGMENT(read-record)
	unsigned int gffCount;
	// overall max count 
	unsigned int globalMax;
	//number of ROIs read
	unsigned int nRoi;

	// read ROI file
	void readROIFile(seqan::CharString inputFileName) ;

};



// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// read ROI file
void roiGFFApp::readROIFile(seqan::CharString inputFileName) 
{
	std::fstream stream(seqan::toCString(inputFileName),
			std::ios::binary | std::ios::in);
	if (!stream.good())
		std::cerr << "Error message";
	seqan::RecordReader<std::fstream, seqan::SinglePass<> > reader(stream);
	seqan::CharString result;
	seqan::StringSet<seqan::CharString> roiRow;
	seqan::StringSet<seqan::CharString> countStr;
	seqan::StringSet<seqan::CharString> countSplit;
	ROI currentROI;
	while (!atEnd(reader)) 
	{
		clear(result);
		int res = readLine(result, reader);
		if (res == 0) 
		{
			// TODO do something
		}
		if (length(result) > 0) 
		{
			// comment lines start with "#" and will be ignored
			if (result[0] != '#') 
			{
				clear(roiRow);
				clear(countStr);
				strSplit(roiRow, result, '\t');
				if (length(roiRow) < 6) 
				{
					std::cerr << "there is a problem with the following line: "
							<< std::endl << result << std::endl;
					exit(-1);
				}
				resize(roiSet, length(roiSet) + 1);
				currentROI = back(roiSet);
				back(roiSet).chr = suffix(roiRow[0], 0);
				back(roiSet).startPos = atoi(seqan::toCString(((seqan::CharString) roiRow[1])));
				back(roiSet).endPos = atoi(seqan::toCString(((seqan::CharString) roiRow[2])));
				back(roiSet).dir = roiRow[3];
				back(roiSet).len = (unsigned) atol(
						seqan::toCString(((seqan::CharString) roiRow[4])));
				back(roiSet).regionName=suffix(roiRow[5], 0);
				back(roiSet).countMax = (unsigned) atol(
						seqan::toCString(((seqan::CharString) roiRow[6])));
				
				strSplit(countSplit, roiRow[length(roiRow)-1], ',');

				for (int i = 0;	(unsigned) i < length(countSplit); i++) 
				{
					int val = atoi(seqan::toCString(((seqan::CharString)countSplit[i])));
					if (val > globalMax)
						globalMax = val;
					append(back(roiSet).count, val);
				}

				if (options.verbosity >= 2)
					std::cout << ((ROI) back(roiSet)).startPos << " "
							<< ((ROI) back(roiSet)).endPos << " "
							<< ((ROI) back(roiSet)).len << " "
							<< ((ROI) back(roiSet)).countMax << " "
							<< ((ROI) back(roiSet)).dir << " "
							<< ((ROI) back(roiSet)).chr << " NN "
							<< toCString(result) << std::endl;
			}
		}
	}
	nRoi = length(roiSet);
	std::cout << "found " << length(roiSet) << " regions of interest" << std::endl;
	std::cout << "global max: " << globalMax << std::endl;
	stream.close();
}

int roiGFFApp::run() {

	std::cout << "Option Arguments:" << std::endl;
	std::cout << "  input GFF file:  \"" << options.gffFileName << "\""
			<< std::endl;
	std::cout << "  input ROI file:  \"" << options.roiFileName << "\""
			<< std::endl;
	std::cout << "  output ROI file: \"" << options.outputFileName << "\""
			<< std::endl;
	std::cout << "Non-option Arguments:" << std::endl;

//read ROI
	readROIFile(options.roiFileName);

//read GFF
/*	std::fstream stream(toCString(options.inputGFF),
			std::ios::binary | std::ios::in);
	if (!stream.good())
		return 1;
	RecordReader<std::fstream, SinglePass<> > reader(stream);
	String<Gff2Record> gffRecords;
	int res = read2(gffRecords, reader, Gff2());
	if (res != 0)
		return res;
	/*  Print out GFF for testing

	 for (unsigned i = 0; i < length(gffRecords); ++i) {
	 std::cout << gffRecords[i].seqName << "\t" << gffRecords[i].source
	 << "\t" << gffRecords[i].feature << "\t" << gffRecords[i].start
	 << "\t" << gffRecords[i].end << "\t" << gffRecords[i].score
	 << "\t" << gffRecords[i].strand << "\t" << gffRecords[i].frame
	 << "\t" << gffRecords[i].attributes << "\t"
	 << "\t" << gffRecords[i].idName << "\t"
	 << gffRecords[i].comments << std::endl;

	 }
	 * /
	// search for each GFF entry
	StringSet<ROI> gffROI;
	ROI currentROI;
	for (unsigned i = 0; i < length(gffRecords); ++i) {
		resize(gffROI, length(gffROI) + 1);
		currentROI = back(gffROI);
		back(gffROI).chr = gffRecords[i].idName;
		back(gffROI).startPos = gffRecords[i].start;
		back(gffROI).endPos = gffRecords[i].end;
		back(gffROI).dir = gffRecords[i].strand;
		unsigned int rLen = gffRecords[i].end - gffRecords[i].start + 1;
		back(gffROI).len = rLen;
		back(gffROI).countMax = 0;
		for (unsigned int ci = 0; ci < rLen; ci++) {
			append(back(gffROI).count, 0);
		}
		for (unsigned r = 0; r < nRoi; r++) {
			if (roiSet[r].chr != gffRecords[i].seqName) {
				continue;
			}
			if (options.strandSpecific && roiSet[r].dir != gffRecords[i].strand)
				continue;
			if ((unsigned) (roiSet[r]).startPos < gffRecords[i].end
					&& (unsigned) roiSet[r].endPos > gffRecords[i].start) {
				int startOffset = roiSet[r].startPos - gffRecords[i].start;
				if (startOffset < 0)
					startOffset = 0;
				//std::cout << "working on roi Nr: " << r << std::endl;
				int minOffset = gffRecords[i].start - roiSet[r].startPos;
				if (minOffset < 0)
					minOffset = 0;
				int endOffset = roiSet[r].endPos - gffRecords[i].start;
				if (endOffset > rLen)
					endOffset = rLen;
				for (; startOffset < endOffset; startOffset++, minOffset++) {
					/*std::cout
					 << gffRecords[i].start << "\t"
					 << gffRecords[i].end << "\t"
					 << roiSet[r].startPos << "\t"
					 << roiSet[r].endPos << "\t"
					 << rLen << "\t"
					 << roiSet[r].len << "\t"
					 << startOffset << "\t"
					 << minOffset << "\t"
					 << endOffset << "\t"
					 << std::endl;* /
					back(gffROI).count[startOffset] += roiSet[r].count[minOffset];
					if((unsigned)(back(gffROI).countMax) < (unsigned)(back(gffROI).count[startOffset]))
					back(gffROI).countMax = (unsigned)(back(gffROI).count[startOffset]);
				}
			}
		}
	}
	std::ofstream myfile;
	myfile.open(toCString(options.outputROI));
	std::cout << "printing to output file " << toCString(options.outputROI) << std::endl;
	myfile << "# 1.  column: reference name" << std::endl;
	myfile << "# 2.  column: starting position on reference genome"
			<< std::endl;
	myfile << "# 3.  column: end position on reference genome" << std::endl;
	myfile << "# 4.  column: strand +/-(. if empty)" << std::endl;
	myfile << "# 5.  column: length of region" << std::endl;
	myfile << "# 6.  column: max Count" << std::endl;
	myfile << "# 7+. columns: individual counts" << std::endl;
	for (unsigned r = 0; r < length(gffRecords); r++) {
		myfile << ((ROI) (gffROI[r])).chr << "\t";
		myfile << ((ROI) (gffROI[r])).startPos << "\t";
		myfile << ((ROI) (gffROI[r])).endPos << "\t";
		myfile << ((ROI) (gffROI[r])).dir << "\t";
		myfile << ((ROI) (gffROI[r])).len << "\t";
		myfile << ((ROI) (gffROI[r])).countMax << "\t"; //TODO always 0 WHY!!! BUG!!!!
		for (unsigned li = 0; li < ((ROI) (gffROI[r])).len - 1; li++) {
			myfile << ((ROI) (gffROI[r])).count[li] << "\t";
		}
		myfile << ((ROI) (gffROI[r])).count[((ROI) (gffROI[r])).len - 1]
				<< std::endl;
	}
	myfile.close();
	*/
	return 0;
}


// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function mainWithOptions()
// ----------------------------------------------------------------------------

int mainWithOptions(Options const & options)
{
    roiGFFApp app(options);
    return app.run();
}
