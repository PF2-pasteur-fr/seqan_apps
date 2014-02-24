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

#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include "ap.h"
#include "interpolation.h"
#include "dataanalysis.h"
#include "linalg.h"
#include "alglibinternal.h"
#include <math.h>

#ifndef SANDBOX_JAGLA_APPS_ROIPCA_ROIPCA_H_
#define SANDBOX_JAGLA_APPS_ROIPCA_ROIPCA_H_

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

struct Options
{
    // Verbosity: 0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
	unsigned int verbosity;

    // Input and Output file names.
	CharString inputFileName;
	CharString outputFileName;

	Options() : verbosity(0)
    {}
};

// ----------------------------------------------------------------------------
// Class ROI
// ----------------------------------------------------------------------------

struct ROI
{
	int min;
	int max;
	CharString rId;
	CharString dir;
	CharString name;
	String<unsigned> count;
	unsigned len;
	unsigned countMax;

	ROI() : min(seqan::maxValue<int>()), max(0), rId(-2), dir("."), len(0), countMax(0)
    {}
};

// ----------------------------------------------------------------------------
// class RoiStats
// ----------------------------------------------------------------------------

// Collection of vectors that are used for the ROI computation.

class RoiStats
{
public:
	alglib::real_1d_array roiFirstnDimEigenValues, roiLengthNormalizedEigenValues;
	alglib::real_2d_array roiFirstnDimEigenVector, roiLengthNormalizedEigenVector;
	alglib::ae_int_t info;
	alglib::real_1d_array roiFirstnDimVec_prj_1;
	alglib::real_1d_array roiLengthNormalized_prj_1;
	alglib::real_1d_array roiFirstnDimVec_prj_2;
	alglib::real_1d_array roiLengthNormalized_prj_2;
	alglib::real_1d_array roiFirstnDimVector_1;
	alglib::real_1d_array roiLengthNormalizedVector_1;
	alglib::real_1d_array roiFirstnDimVector_2;
	alglib::real_1d_array roiLengthNormalizedVector_2;
};

// ----------------------------------------------------------------------------
// Class ComputeRoiStatsApp
// ----------------------------------------------------------------------------

class ComputeRoiStatsApp : public RoiStats
{
public:
    Options options;
    unsigned nDim;
    unsigned globalMax;
    unsigned nRoi;

    String<ROI> roi;

    ComputeRoiStatsApp(Options const & options) :
            options(options), nDim(100), globalMax(0), nRoi(0)
    {}

    int run()
    {
        int res = 0;

        double startLoad = sysTime();
        if (options.verbosity >= 1)
            std::cerr << "__PREPARATION_________________________________________________________________\n"
                      << "\n"
                      << "LOADING " << options.inputFileName << " ...";
        res = readRoiFile(options.inputFileName);
        if (res != 0)
            return res;
        if (options.verbosity >= 1)
            std::cerr << " OK\n"
                      << "\n"
                      << "  Loaded " << length(roi) << " regions of interest.\n"
                      << "  Loading took " << (sysTime() - startLoad) << " s\n"
                      << "\n";
        res = analyze();
        if (res != 0)
            return res;

        if (options.verbosity >= 1)
            std::cerr << "__WRITING OUTPUT______________________________________________________________\n"
                      << "\n";
        res = writeOutput();
        if (res != 0)
            return res;

        std::cerr << "DONE.\n";
        return 0;
    }

    // Read ROI file into this->roi.
    int readRoiFile(CharString const & inputFileName);

    // Print memory statistics on Linux.
    void printMemStats();

    // Perform analysis.
    int analyze();

    // Write output.
    int writeOutput();
};


// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function peakiness()
// ----------------------------------------------------------------------------

double peakiness(ROI const & roi)
{
	double sumDiff = 0.0;

	if (roi.len <= 1)
		return 0;

	for (unsigned i = 0; i < roi.len - 1; i++)
		sumDiff += fabs(
				(double) roi.count[i] / (double) roi.countMax
						- (double) roi.count[i + 1] / (double) roi.countMax);

    return sumDiff / (roi.len - 1);
}

// ----------------------------------------------------------------------------
// Function aoc()
// ----------------------------------------------------------------------------

double aoc(ROI const & roi)
{
	if (roi.countMax == 0 || roi.len == 0)
		return 0;
	double sumDiff = 0.0;
	for (unsigned i = 0; i < roi.len; i++)
		sumDiff += (double) roi.count[i] / (double) roi.countMax;

    return sumDiff / roi.len;
}

// ----------------------------------------------------------------------------
// Function changeSlope()
// ----------------------------------------------------------------------------

double changeSlope(ROI const & roi)
{
	unsigned changedSlope = 0;
	for (unsigned i = 0; i < roi.len - 2; i++)
    {
		if (roi.count[i] > roi.count[i + 1]
				&& roi.count[i + 1] < roi.count[i + 2])
			changedSlope++;
		else if (roi.count[i] < roi.count[i + 1]
				&& roi.count[i + 1] > roi.count[i + 2])
			changedSlope++;
	}

	return 1.0 * changedSlope / (roi.len - 2);
}

// ----------------------------------------------------------------------------
// Member Function ComputeRoiStatsApp::readRoiFile()
// ----------------------------------------------------------------------------

inline int ComputeRoiStatsApp::readRoiFile(CharString const & inputFileName)
{
    // Open file and create record reader.
	std::fstream stream(toCString(inputFileName), std::ios::binary | std::ios::in);
	if (!stream.good())
		std::cerr << "Could not open " << inputFileName << " for reading!\n";
	RecordReader<std::fstream, SinglePass<> > reader(stream);

    // Setup buffers.
	CharString buffer;
	StringSet<CharString> header;
	StringSet<CharString> countStr;
	StringSet<CharString> countSplit;
	ROI currentROI;

    // Read all records.
	while (!atEnd(reader))
    {
        clear(buffer);
        int res = 0;

        // Skip if this is a header/comment.
        if (value(reader) == '#')
        {
            skipLine(reader);
            continue;
        }

        // Read reference name.
        clear(currentROI.rId);
        if ((res = readUntilTabOrLineBreak(currentROI.rId, reader)) != 0)
            return res;

        // Skip TAB.
        if (skipChar(reader, '\t') != 0)
            return res;

        // Read and parse start position.
        clear(buffer);
        if ((res = readUntilTabOrLineBreak(buffer, reader)) != 0)
            return res;
        if (!lexicalCast2(currentROI.min, buffer))
            return 1;
        //currentROI.min -= 1;  // transform to 0-based

        // Skip TAB.
        if (skipChar(reader, '\t') != 0)
            return res;

        // Read and parse end position.
        clear(buffer);
        if ((res = readUntilTabOrLineBreak(buffer, reader)) != 0)
            return res;
        if (!lexicalCast2(currentROI.max, buffer))
            return 1;
        
        // Skip TAB.
        if (skipChar(reader, '\t') != 0)
            return res;

        // Read strand.
        clear(buffer);
        if ((res = readUntilTabOrLineBreak(buffer, reader)) != 0)
            return res;
        if (buffer[0] != '.' && buffer[0] != '+' && buffer[0] != '-')
            return 1;
        currentROI.dir = buffer;
        
        // Skip TAB.
        if (skipChar(reader, '\t') != 0)
            return res;

        // Read and parse region length.
        clear(buffer);
        if ((res = readUntilTabOrLineBreak(buffer, reader)) != 0)
            return res;
        if (!lexicalCast2(currentROI.len, buffer))
            return 1;

        // Skip TAB.
        if (skipChar(reader, '\t') != 0)
            return res;

        // Read and parse region name.
        clear(currentROI.name);
        if ((res = readUntilTabOrLineBreak(currentROI.name, reader)) != 0)
            return res;

        // Skip TAB.
        if (skipChar(reader, '\t') != 0)
            return res;

        // Read max count.
        clear(buffer);
        if ((res = readUntilTabOrLineBreak(buffer, reader)) != 0)
            return res;
        if (!lexicalCast2(currentROI.countMax, buffer))
            return 1;

        // Skip TAB.
        if (skipChar(reader, '\t') != 0)
            return res;

        // Individual counts.
        clear(currentROI.count);
        clear(buffer);
        for (; !atEnd(reader) && value(reader) != '\r' && value(reader) != '\n'; goNext(reader))
        {
            if (value(reader) != ',')
            {
                if (!isdigit(value(reader)))
                    return 1;  // Error parsing.
                appendValue(buffer, value(reader));
            }
            else
            {
                if (empty(buffer))
                    continue;
                unsigned count = 0;
                if (!lexicalCast2(count, buffer))
                    return 1;  // Error parsing.
                appendValue(currentROI.count, count);
                currentROI.countMax = std::max(currentROI.countMax, back(currentROI.count));
                clear(buffer);
            }
        }
        if (!empty(buffer))
        {
            unsigned count = 0;
            if (!lexicalCast2(count, buffer))
                return 1;  // Error parsing.
            appendValue(currentROI.count, count);
            currentROI.countMax = std::max(currentROI.countMax, back(currentROI.count));
        }

        // Skip over line.
        if (!atEnd(reader))
            skipLine(reader);
        appendValue(roi, currentROI);

        if (options.verbosity >= 2)
            std::cout << ((ROI) back(roi)).min << " "
                      << ((ROI) back(roi)).max << " "
                      << ((ROI) back(roi)).len << " "
                      << ((ROI) back(roi)).countMax << " "
                      << ((ROI) back(roi)).dir << " "
                      << ((ROI) back(roi)).name << " "
                      << ((ROI) back(roi)).rId << " NN "
                    /*<< result*/ << std::endl;
    }

	nRoi = length(roi);
    return 0;
}

// ----------------------------------------------------------------------------
// Member Function ComputeRoiStatsApp::printMemStats()
// ----------------------------------------------------------------------------

inline void ComputeRoiStatsApp::printMemStats()
{
#if defined(_WIN32)
#else
	int tSize = 0, resident = 0, share = 0;
	std::ifstream buffer("/proc/self/statm");
	buffer >> tSize >> resident >> share;
	buffer.close();
	long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
	double rss = resident * page_size_kb;
	std::cerr << "Memstats: " << rss << "\t";
	double shared_mem = share * page_size_kb;
	std::cerr << shared_mem << "\t";
	std::cerr << rss - shared_mem << "\n";
#endif
}

// ----------------------------------------------------------------------------
// Member Function ComputeRoiStatsApp::analyze()
// ----------------------------------------------------------------------------

inline int ComputeRoiStatsApp::analyze()
{
    if (options.verbosity >= 1)
        std::cerr << "__ANALYSIS____________________________________________________________________\n\n";

	alglib::real_2d_array roiFirstnDim, roiLengthNormalized;
	alglib::real_1d_array sumRoiFirstnDim;
	roiFirstnDim.setlength(nRoi, nDim);
	roiLengthNormalized.setlength(nRoi, nDim);
	sumRoiFirstnDim.setlength(nDim);
	for (unsigned li = 0; li < nDim; li++) {
		sumRoiFirstnDim(li) = 0;
	}
	std::cout << "PCA of first nDim " << std::endl;
	std::cout << "nDim: " << nDim << std::endl;
/*
	std::cout << "r2: " << r2.rows() << " " << r2.cols() << std::endl;
	std::cout << "r3: " << r3.rows() << " " << r3.cols() << std::endl;
	std::cout << "2r2: " << sr2.length() << std::endl;
	std::cout << "nDim: " << nDim << std::endl;
	std::cout << "nROI: " << nRoi << std::endl;
*/
	/* we are calculating a normalized vector of the first 100 nt of a given region (roiFirstnDim)
	 * NOT strandspecific
	 * sr2 will be used to recenter the data
	 */

    if (options.verbosity >= 1)
        std::cerr << "PCA of first " << nDim << " values ...";
    double startTime = sysTime();
	//PCA of first nDim values
	for (unsigned r = 0; r < nRoi; r++) {
		for (unsigned li = 0; li < nDim; li++) {
			//std::cout << r << " " << li << " " << roi[r].count[li] <<  " " << roi[r].countMax << " "  << std::endl;
			if (li < ((ROI) (roi[r])).len) {
				roiFirstnDim(r, li) = (double) (((ROI) (roi[r])).count[li])
						/ (double) ((ROI) (roi[r])).countMax;
				sumRoiFirstnDim(li) += roiFirstnDim(r, li);
			} else {
				roiFirstnDim(r, li) = 0.0;
			}
		}
	}
    if (options.verbosity >= 1)
        std::cerr << " OK\n"
                  << "  Took " << (sysTime() - startTime) << " s\n\n";
    if (options.verbosity >= 2)
        printMemStats();

	// PCA of resized values
    if (options.verbosity >= 1)
        std::cerr << "PCA of all values ...";
    startTime = sysTime();
	for (unsigned r = 0; r < nRoi; r++) {
		if (((ROI) (roi[r])).len < nDim) {
			int step = floor(
					(double) nDim / (double) ((ROI) (roi[r])).len + 0.5);
			unsigned idx = 0;
			for (unsigned li = 0; li < nDim; li++) {
				roiLengthNormalized(r, li) = (double) ((ROI) (roi[r])).count[idx]
						/ (double) ((ROI) (roi[r])).countMax;
				if (li > 0 && (li % step) == 0) {
					idx++;
					if (idx >= ((ROI) (roi[r])).len)
						idx = ((ROI) (roi[r])).len - 1;
				}
			}
		} else {
			unsigned li = 0;
			unsigned oldStep = 0;
			unsigned nStep = 0;
			double sum = 0.0;
			for (unsigned idx = 0; idx < ((ROI) (roi[r])).len; idx++) {
				int step = floor(
						(double) nDim / (double) ((ROI) (roi[r])).len
								* (double) idx);
				if ((int)oldStep != step) {
					roiLengthNormalized(r, li) = sum / (double) ((ROI) (roi[r])).countMax
							/ (double) nStep;
					sum = (double) ((ROI) (roi[r])).count[idx];
					nStep = 1;
					oldStep = step;
					li++;
				} else {
					nStep++;
					sum += (double) ((ROI) (roi[r])).count[idx];
				}
			}
			roiLengthNormalized(r, li) = sum / (double) ((ROI) (roi[r])).countMax
					/ (double) nStep;
		}
	}
    if (options.verbosity >= 1)
        std::cerr << " OK\n"
                  << "  Took " << (sysTime() - startTime) << " s\n\n";
    if (options.verbosity >= 2)
        printMemStats();

    if (options.verbosity >= 1)
        std::cerr << "recentering data ...";
    startTime = sysTime();
	//recenter data
	for (unsigned li = 0; li < nDim; li++) {
		sumRoiFirstnDim(li) = sumRoiFirstnDim(li) / nRoi;
	}
	for (unsigned r = 0; r < nRoi; r++) {
		for (unsigned li = 0; li < nDim; li++) {
			roiFirstnDim(r, li) = roiFirstnDim(r, li) - sumRoiFirstnDim(li);
		}
	}
    if (options.verbosity >= 1)
        std::cerr << " OK\n"
                  << "  Took " << (sysTime() - startTime) << " s\n\n";
    if (options.verbosity >= 2)
        printMemStats();
	//std::cout << "2.: roiFirstnDim " << roiFirstnDim.tostring(5).c_str() << std::endl;

    if (options.verbosity >= 1)
        std::cerr << "PCA ...";
    startTime = sysTime();

	roiFirstnDimVector_1.setlength(nDim);
	roiLengthNormalizedVector_1.setlength(nDim);
	roiFirstnDimVector_2.setlength(nDim);
	roiLengthNormalizedVector_2.setlength(nDim);
	roiFirstnDimVec_prj_1.setlength(nRoi);
	roiLengthNormalized_prj_1.setlength(nRoi);
	roiFirstnDimVec_prj_2.setlength(nRoi);
	roiLengthNormalized_prj_2.setlength(nRoi);

    if (options.verbosity >= 2)
        printMemStats();

	//PCA
	alglib::pcabuildbasis(roiFirstnDim, (alglib::ae_int_t) nRoi, (alglib::ae_int_t) nDim,
			info, roiFirstnDimEigenValues, roiFirstnDimEigenVector);
	alglib::pcabuildbasis(roiLengthNormalized, (alglib::ae_int_t) nRoi, (alglib::ae_int_t) nDim,
			info, roiLengthNormalizedEigenValues, roiLengthNormalizedEigenVector);
	//std::cout << "info: " << info << std::endl;

	for (unsigned i = 0; i < nDim; i++) {
		roiFirstnDimVector_1(i) = roiFirstnDimEigenVector(i, 0);
		roiLengthNormalizedVector_1(i) = roiLengthNormalizedEigenVector(i, 0);
		roiFirstnDimVector_2(i) = roiFirstnDimEigenVector(i, 1);
		roiLengthNormalizedVector_2(i) = roiLengthNormalizedEigenVector(i, 1);
	}
    if (options.verbosity >= 1)
        std::cerr << " OK\n\n"
                  << "  Took " << (sysTime() - startTime) << " s\n";

    if (options.verbosity >= 2)
        printMemStats();

    if (options.verbosity >= 1)
        std::cout << "Projection ...";
    startTime = sysTime();
	//Projection
	alglib::rmatrixmv(nRoi, nDim, roiFirstnDim, 0, 0, 0, roiFirstnDimVector_1, 0, roiFirstnDimVec_prj_1, 0);
	alglib::rmatrixmv(nRoi, nDim, roiLengthNormalized, 0, 0, 0, roiLengthNormalizedVector_1, 0, roiLengthNormalized_prj_1, 0);
	alglib::rmatrixmv(nRoi, nDim, roiFirstnDim, 0, 0, 0, roiFirstnDimVector_2, 0, roiFirstnDimVec_prj_2, 0);
	alglib::rmatrixmv(nRoi, nDim, roiLengthNormalized, 0, 0, 0, roiLengthNormalizedVector_2, 0, roiLengthNormalized_prj_2, 0);
    if (options.verbosity >= 1)
        std::cerr << " OK\n"
                  << "  Took " << (sysTime() - startTime) << " s\n\n";

    if (options.verbosity >= 2)
        printMemStats();

	return 0;
}

// ----------------------------------------------------------------------------
// Member Function ComputeRoiStatsApp::analyze()
// ----------------------------------------------------------------------------

inline int ComputeRoiStatsApp::writeOutput()
{
	std::cerr << "WRITING TO " << options.outputFileName << " ...";

	std::ofstream myfile;
	myfile.open(toCString(options.outputFileName));
	myfile << "# 1.  column: reference name" << std::endl;
	myfile << "# 2.  column: starting position on reference genome"
			<< std::endl;
	myfile << "# 3.  column: end position on reference genome" << std::endl;
	myfile << "# 4.  column: strand +/-(. if empty)" << std::endl;
	myfile << "# 5.  column: length of region" << std::endl;
	myfile << "# 6.  column: feature name" << std::endl;
	myfile << "# 7.  column: max_Count" << std::endl;
	myfile
			<< "# 8.  column: projection on 1 principle component of PCA of first "
			<< nDim << "roi values" << std::endl;
	myfile
			<< "# 9.  column: projection on 1 principle component of PCA of whole roi mapped to "
			<< nDim << " values" << std::endl;
	myfile
			<< "# 10.  column: projection on 1 principle component of PCA of first "
			<< nDim << "roi values" << std::endl;
	myfile
			<< "# 11.  column: projection on 1 principle component of PCA of whole roi mapped to "
			<< nDim << " values" << std::endl;
	myfile << "# 12.  column: peakiness: sum(abs(xi - xi+1))/N" << std::endl;
	myfile << "# 13. column: counts peaks (xi>xi+1<xi+2 | xi<xi+1>xi+2)"
			<< std::endl;
	myfile << "# 14. column: aoc : area under curve, sum(xi/xMax)/N"
			<< std::endl;
	myfile << "# 15+. columns: individual counts" << std::endl;
	myfile << "##chr\tstart\tend\tstrand\tlength\tfeature name\tmax count\tproj1(100)\tproj1(all)"
		      "\tproj2(100)\tproj2(all)\tpeakiness\tcountPeaks\taoc\tcounts" <<std::endl;
	for (unsigned r = 0; r < nRoi; r++) {
		myfile << ((ROI) (roi[r])).rId << "\t";
		myfile << ((ROI) (roi[r])).min << "\t";
		myfile << ((ROI) (roi[r])).max << "\t";
		myfile << ((ROI) (roi[r])).dir << "\t";
		myfile << ((ROI) (roi[r])).len << "\t";
		myfile << ((ROI) (roi[r])).name << "\t";
		myfile << ((ROI) (roi[r])).countMax << "\t";
		myfile << roiFirstnDimVec_prj_1(r) << "\t"; //projection of first nDim values
		myfile << roiLengthNormalized_prj_1(r) << "\t"; // projection of resized rois
		myfile << roiFirstnDimVec_prj_2(r) << "\t"; //projection of first nDim values
		myfile << roiLengthNormalized_prj_2(r) << "\t"; // projection of resized rois
		myfile << peakiness(((ROI) (roi[r]))) << "\t";
		myfile << changeSlope(((ROI) (roi[r]))) << "\t";
		myfile << aoc(((ROI) (roi[r]))) << "\t";
		for (unsigned li = 0; li < ((ROI) (roi[r])).len - 1; li++) {
			myfile << ((ROI) (roi[r])).count[li] << ",";
		}
		myfile << ((ROI) (roi[r])).count[((ROI) (roi[r])).len - 1] << std::endl;
	}
	myfile.close();
	/*std::cout << "sr2 " << sr2.tostring(5).c_str()<< std::endl;
	 std::cout << "3.: r2 " << r2.tostring(5).c_str()<< std::endl;
	 std::cout << "eigenV1 " << eigenV1.tostring(5).c_str()<< std::endl;
	 std::cout << "eigenMatrix " << eigenVectors2.tostring(5).c_str()<< std::endl;
	 std::cout << "1st PC: "
	 << eigenVectors2(0,0) << " "
	 << eigenVectors2(1,0) << " "
	 << eigenVectors2(2,0) << std::endl;
	 std::cout << "projection1" << projection1.tostring(5).c_str()<< std::endl;
	 std::cout << "eigenV1 " << eigenValues2[0]<< std::endl;
	 std::cout << "eigenV2 " << eigenValues2[1]<< std::endl;
	 std::cout << "eigenV3 " << eigenValues[2]<< std::endl;
	 */
    
    std::cerr << " OK\n";

    return 0;
}

#endif  // #ifndef SANDBOX_JAGLA_APPS_ROIPCA_ROIPCA_H_
