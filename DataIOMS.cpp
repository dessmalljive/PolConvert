/* DATAIO - FITS-IDI interface to PolConvert

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

*/

#include <sys/types.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dirent.h>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/ms/MeasurementSets/MSColumns.h>
// MSColumns.h automatically includes MSMainColumns.h, where many of the definitions live

#include <casacore/ms/MeasurementSets/MSPolarization.h>
#include <casacore/ms/MeasurementSets/MSDataDescription.h>
#include <casacore/ms/MeasurementSets/MSSpectralWindow.h>
#include <casacore/ms/MeasurementSets/MSSpWindowColumns.h>
#include <casacore/measures/Measures/MDirection.h>

#include <casacore/casa/Arrays/MatrixMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>

#include <casacore/scimath/Mathematics/SquareMatrix.h>
#include <casacore/tables/Tables/TableIter.h>
#include <casacore/ms/MeasurementSets/MSIter.h>
#include "./DataIOMS.h"



using namespace casacore;

void DataIOMS::log(char* message) {
    fprintf(logFile, "%s", message);
    std::cout << message;
    fflush(logFile);
}

/*

The DataIO.h header includes these many variables.

   FreqSetup *Freqs;
   AutoCorrelation *AutoCorrs;
   int *an1, *an2, *linAnts, *field; //, *sour;
   double *ParAng[2];
   double *UVDist;
   int NLinAnt, NIFs, Nband, status, Nants, Nfreqs, currFreq;
   long Nvis, currVis, *indexes;
   bool success, currConj;
   bool *conjugate, *is1, *is2, *is1orig, *is2orig;
   int Flux, *NAV;
   int nautos;

   // We have to care about:

   float ***averAutocorrs;
   long NLinVis;

   double day0;
   double *Freqvals[MAXIF], *doRange, *JDTimes;
   bool *conjugate, *is1, *is2;
   bool *is1orig, *is2orig;

   // Passed in and stored and not our problem otherwise
   ArrayGeometry *Geometry;

   // Not even used in DataIOSWIN so probably safe to ignore:

   int *Basels, *Freqids,
   bool *conjugate
   int *an1, *an2

*/

/*
  The following API calls occur in _PolConvert.cpp:

  DataIO *DifXData ;  // Polymorphism to SWIN or FITS-IDI.
  DifXData = new DataIOSWIN(nSWINFiles, SWINFiles, nALMA,
  DifXData = new DataIOFITS(outputfits, nALMA, almanums,
  DifXData->succeed())
  DifXData->getDay0());
  DifXData->getNfreqs();
  DifXData->getNchan(ii);
  DifXData->getMixedNvis());
  ((plRange[0]+DifXData->getDay0())-2400000.5)*86400.;
  DifXData->getAmpRatio(currAntIdx, ii, ij);
  DifXData->getAmpRatio(currAntIdx, im, ij);
  if(!DifXData->setCurrentIF(ii)){
  DifXData->getFrequencies(DifXFreqs);
  while(DifXData->getNextMixedVis(
  currFile = DifXData->getFileNumber();
  DifXData->succeed()
  DifXData->finish();
  DifXData->applyMatrix(Ktotal[currAntIdx],XYSWAP[currAntIdx],auxB2,
                        currAntIdx,plotFile[IFplot]);
  DifXData->zeroWeight();
  DifXData->setCurrentMixedVis();
  DifXData->finish();

*/

// const int secperday = 86400;

DataIOMS::DataIOMS(std::string msfn, int NLinAnt, int *LinAnts, double *Range, int nIF2Conv, int *IF2Conv,
                   bool Overwrite, bool doConj, bool doSolve,
                   ArrayGeometry *Geom, bool doPar, FILE *logF)
    :
    currRow(0),
    nextRow(0),
    currSPW(0),
    doneFirst(false),
    doWriteCirc(doSolve),
    logFile(logF),
    doConjugate(doConj),
    canPlot(false),
    doParang(doPar),
    isTwoLinear(false) {
    success = true;

    
    std::cerr << std::endl << "MS file name: " << msfn << std::endl;
    // Transcribe linAnts; somewhere they get one(1) added for AIPSness; we take it off again for CASAness
    std::cerr << "Linants: "; 
    for (int i=0; i < NLinAnt; i++) {
        linAnts.insert(LinAnts[i]-1);
        std::cerr << LinAnts[i]-1 << " ";
    }
    std::cerr << std::endl;
    // According to polconvert.py, if nDoIF is 0 then we should do _all_ IFs.

    for (int i=0; i < nIF2Conv; i++) {
        IFsToDo.push_back(IF2Conv[i]);
        // We get the IFs labelled from 1, FITS style
    }
    Geometry = Geom;
    //currSPW(0) // Should we really not set this?!
    try {
        ms = MeasurementSet(msfn);
        std::cerr << "Opened the MS at least" << std::endl;
    } catch (casacore::TableNoFile t) {
        std::cout << "Measurement set " << msfn << " not found" << std::endl;
        success = false;
    }

    if (!success) {
        return;
    }

    MSColumns msc(ms);
    // FIXME: is this really the first day?
    day0 = floor(msc.time()(0)/86400.0) + 2400000.5;

    Double nRange0 = ((Range[0] + day0) - 2400000.5)*86400.;
    Double nRange1 = ((Range[1] + day0) - 2400000.5)*86400.;
    doRange.push_back(nRange0);
    doRange.push_back(nRange1);

    std::cerr << "Time range originally " << Range[0] << " -- " << Range[1] << std::endl;
    std::cerr << "Time range now " << nRange0 << " -- " << nRange1 << std::endl;
    std::cerr << "day0 " << day0 << std::endl;

    // Set Nfreqs, which in the FITS case is the number of spectral
    // windows times the number of "frequency setups"
    MSSpectralWindow& swtable = ms.spectralWindow();
    // NFreqs is the DataIO name for the number of spectral windows
    Nfreqs = swtable.nrow();
    Freqs = new FreqSetup[Nfreqs];

    const MSSpWindowColumns& spwc(swtable);
    for (int i=0; i < Nfreqs; i++) {
        Freqs[i].Nchan = spwc.numChan()(i);
        Freqs[i].BW = spwc.totalBandwidth()(i);
        Freqs[i].Nu0 = spwc.refFrequency()(i); // In what units?
        Freqs[i].SB = spwc.netSideband()(i); // Encoded how?
        Freqvals[i] = new double[Freqs[i].Nchan];
        Array<Double> freqs(spwc.chanFreq()(i));
        for (Int j=0; j < Freqs[i].Nchan; j++) {
            Freqvals[i][j] = freqs(IPosition(1, j));
        }
    }

    for (int i=0; i < Nfreqs; i++) {
        std::cerr << "Freqs[" << i << "]" << std::endl;
        std::cerr << "   Nchan " << Freqs[i].Nchan << std::endl;
        std::cerr << "   BW " << Freqs[i].BW << std::endl;
        std::cerr << "   Nu0 " << Freqs[i].Nu0 << std::endl;
        std::cerr << "   SB " << Freqs[i].SB << std::endl;
    }

    std::cerr << "IFsToDo = [";
    for (Int i: IFsToDo ) {
        std::cerr << i << ' ';
    }
    std::cerr << "]" << std::endl;

    // We have:
    // int NLinAnt, int *LinAnts,

    NLinVis = 12; // FIXME; if this is used we would actually have to count I guess.

    // An array of NLinAnt * NFreqs * Nchans
    // Prepare memory for average autocorrs.
    // DataIO.h defines averAutcorrs as float*** so we have to stick to that.

    NAV = new int[NLinAnt];
    averAutocorrs = new float**[NLinAnt];

    for (int i=0; i<NLinAnt; i++){
        NAV[i] = 0;
        averAutocorrs[i] = new float*[Nfreqs];
        for (int j=0; j<Nfreqs; j++){
            averAutocorrs[i][j] = new float[Freqs[j].Nchan];
            for (int k=0; k < Freqs[j].Nchan; k++){
                averAutocorrs[i][j][k] = 1.0;
            }
        }
    }
    std::cerr << "Saving circulars " << std::endl;
    saveCirculars();

}

DataIOMS::~DataIOMS() {
  delete[] NAV;
  delete[] Freqs;
}

bool DataIOMS::setCurrentIF(int i)
{
    std::cerr << "DataIOMS::setCurrentIF(" << i <<")" << std::endl;
    currSPW = i;
    currRow = 0;
    nextRow = 0;
    currAnt = 0;
    std::cerr << "DataIOMS::setCurrentIF currRow=" << currRow << std::endl;
    std::cerr << "Linants currently: ";
    for (auto a : linAnts) std::cerr << a << " ";
    std::cerr << std::endl;
    return true;  // Risky I guess, but what else?
}

/*
int DataIO::getNfreqs() {return Nfreqs;};  // NUMBER OF IFs
int DataIO::getNant() {return Nants;};  // TOTAL NUMBER OF ANTENNAS
int DataIO::getNchan(int freqid) {return Freqs[freqid].Nchan;};  // # OF CHANNEL IN freqid IF
long DataIO::getMixedNvis() {
    std::cerr << "Someone asked DataIO::getMixedNvis()" << std::endl;
    return NLinVis;
}  // NUMBER OF MIXED-POLARIZATION VISIBILITIES FOUND

So we need to set: Nfreqs, Nants and NLinVis
 */


void DataIOMS::getFrequencies(double* output) {

    MSSpectralWindow& swtable = ms.spectralWindow();
    const MSSpWindowColumns& spwc(swtable);
    Nfreqs = swtable.nrow();
    const Array<Double>& Freqvals = spwc.chanFreq()(currSPW);

    size_t i = 0;
    Array<Double>::const_iterator endIter = Freqvals.end();
    for (Array<Double>::const_iterator iter=Freqvals.begin(); iter != endIter; ++iter) {
        output[i] = *iter;
        i++;
    }
    std::cerr << "DataIOMS::getFrequencies: i= " << i << std::endl;

    // memcpy(output, Freqvals[currFreq], Freqs[currFreq].Nchan*sizeof(double));
}


// NB: everything is a reference here; these are all the *return* values!
bool DataIOMS::getNextMixedVis(double &JDTime, int &antenna, int &otherAnt, bool &conj, int &calField)
{

    MSColumns msc(ms);
    canPlot = false;
    currRow = nextRow;
    while (1) {
        if (!(currRow < ms.nrow())) {
            std::cerr << "Ran out of rows after " << currRow << std::endl;
            return false;
        }
        while (1) {
            Int dd = msc.dataDescId()(currRow);
            Int ispw = msc.dataDescription().spectralWindowId()(dd);
            if (ispw != currSPW) {
                currRow++;
                if ((currRow % 1000001) == 0) {
                    std::cerr << "Looking for SPW " << currSPW << " got " << ispw << " row " << currRow << std::endl;

                }
                if (!(currRow < ms.nrow())) {
                    std::cerr << "Ran out of rows after " << currRow << std::endl;
                    std::cerr << "Found " << nmixed << " mixed visibilities" << std::endl;
                    return false;
                }
            } else { break; }
        }
        // Now we know we're on the right SPW. We didn't check if we
        // already were; in that case we just skip zero rows!
        JDTime = msc.time()(currRow);
        calField = msc.fieldId()(currRow); // Check!
        bool lins[2];
        const Int ant1 = msc.antenna1()(currRow);
        const Int ant2 = msc.antenna2()(currRow);

        Bool lin1 = linAnts.find(ant1)!=linAnts.end();
        Bool lin2 = linAnts.find(ant2)!=linAnts.end();
        
        if ((ant1 != ant2) && lin1 && lin2) {
            isTwoLinear = true;
            canPlot = false;
        } else {
            canPlot = true;
        }
        if (currAnt==0) {
            if (lin1) {
                // Return values must be AIPSed too!
                antenna = ant1+1;
                otherAnt = ant2+1;
                conj = true;
                // For next time round
                nextRow = currRow;
                currAnt=1;
                // std::cerr << "Got ant1 = " << ant1 <<std::endl;
                currConj = conj;
                nmixed++;
                return true;
            }
            else {
                // Whether we found a match for ant1, we look at ant2 next.
                currAnt=1;
                // Don't increment row or return!
            }
        }
        // currAnt == 1
        if (lin2) {
            antenna = ant2+1;
            otherAnt = ant1+1;
            conj = false;
            // And for when we come back:
            // FIXME: I think this is wrong.
            nextRow = currRow+1;
            currAnt = 0;
            currConj = conj;
            nmixed++;
            return true;
        }
        else {
            // We've checked both antennas; now we advance to the next row and
            // look for the next SPW match.
            currAnt = 0;
            currRow++;
        }
    }
    // Mollify a compiler warning, although I think this is unreachable:
    currConj = conj;
    return false;
}

bool DataIOMS::setCurrentMixedVis()
{
    MSColumns msc(ms);
    IPosition s(bufferVis.shape());
    // I would like to get the time, but digging the actual time out of the CASA format is a challenge. Can I think of another way?
    const Array<Complex>& oldRow(msc.data()(currRow));
    msc.correctedData()(currRow) = bufferVis;
    // std::cerr << "DataIOMS::setCurrentMixedVis: set row " << currRow << std::endl;
    // std::cerr << "Done setting row " << currRow << std::endl;
    bool isChanged = !allNear(oldRow, bufferVis, 1e-9);
    return true;
}

int DataIOMS::getFileNumber() { return 0; }
// ::zeroWeight is empty in DataIOFITS.cpp, so I guess I can leave it empty here too?
void DataIOMS::zeroWeight() {}

// Modify the visibilities read by "getNextMixedVis" by the calibration matrix supplied
void DataIOMS::applyMatrix(std::complex<float> *Mvec[2][2], bool swap, bool print,
                           int thisAnt, FILE *plotFile)
{
    MSColumns msc(ms);
    Matrix<Complex> currentVis = msc.data()(currRow);
    IPosition visShape(currentVis.shape());
    // std::cerr << "Shape is " << visShape << std::endl;
    bufferVis.resize(visShape);
    Matrix<Complex> zeros(visShape);
    zeros = 0;
    
    // We didn't store the parallactic angles so we have to get them again.
    double parAng[2];
    double AuxPA1, AuxPA2;

    // It turns out we need to store antennas, uvw and time to file
    // below anyway.
    Int a1 = msc.antenna1()(currRow);
    Int a2 = msc.antenna2()(currRow);

    Bool lin1 = linAnts.find(a1)!=linAnts.end();
    Bool lin2 = linAnts.find(a2)!=linAnts.end();
        
    Array<Double> uvw(msc.uvw()(currRow));
    double duvw[3];
    for (int i=0; i<3; i++) {
        duvw[i] = uvw(IPosition(1, i));
    }
    double uvdist = (uvw(IPosition(1, 0))*uvw(IPosition(1, 0)) +
                     uvw(IPosition(1, 1))*uvw(IPosition(1, 1)));
    Int souidx = msc.fieldId()(currRow);
    Double t = msc.time()(currRow);
    getParAng(souidx, a1, a2, duvw, t, AuxPA1, AuxPA2);

    // If I actually knew C++ I could probably do better than this.
    parAng[0] = AuxPA1;
    parAng[1] = AuxPA2;

    if (print && canPlot) {
        int zero = 0;
        int pos = ftell(plotFile);
        std::cerr << "DataIOMS a1=" << a1 << " a2=" << a2
                  << " lin1=" << lin1 << " lin2=" << lin2 << endl;
        std::cerr << "DataIOMS writing antennas " << a1+1 << " " << a2+1
                  << " for time " << t << " at pos " << pos << " with currConj " << currConj
                  << " currentRow " << currRow << " currSPW " << currSPW
                  << std::endl;
        fwrite(&zero, sizeof(int), 1, plotFile);
        fwrite(&t, sizeof(double), 1, plotFile);
        int aips1 = a1+1;
        int aips2 = a2+1;
        fwrite(&aips1, sizeof(int), 1, plotFile);
        fwrite(&aips2, sizeof(int), 1, plotFile);
        fwrite(&parAng[0], sizeof(double), 1, plotFile);
        fwrite(&parAng[1], sizeof(double), 1, plotFile);
        fwrite(&uvdist, sizeof(double), 1, plotFile);

        if (a1 == 1 && a2 == 12) {
            std::cerr << "Angles for baseline 2-13 ";
        }
        for (int k=0; k<Freqs[currSPW].Nchan; k++) {
            /* The default ordering of correlation pairs in FITS-IDI is (probably): RR,LL,RL,LR;
               for Measurement Set the convention is RR, LL, RL, LR. */
            // Visibility matrix:
            Matrix<Complex> e(IPosition(2,2,2));
            Matrix<Complex> e0(IPosition(2,2,2));

            // Calibration Jones matrix
            Matrix<Complex> M(IPosition(2,2,2));
            // Matrix to store back in buffer
            Matrix<Complex> b(IPosition(2,2,2));
            // Parallactic angle Jones matrix
            Matrix<Complex> pa(IPosition(2,2,2));

            e = 0;
            e0 = 0;
            M = 0;
            b = 0;
            pa = 0; // I have trust issues

            e(IPosition(2, 0, 0)) = currentVis(IPosition(2, 0, k));
            e(IPosition(2, 1, 1)) = currentVis(IPosition(2, 1, k));
            e(IPosition(2, 0, 1)) = currentVis(IPosition(2, 2, k));
            e(IPosition(2, 1, 0)) = currentVis(IPosition(2, 3, k));


            e0(IPosition(2, 0, 0)) = currentVis(IPosition(2, 0, k));
            e0(IPosition(2, 1, 1)) = currentVis(IPosition(2, 1, k));
            e0(IPosition(2, 0, 1)) = currentVis(IPosition(2, 2, k));
            e0(IPosition(2, 1, 0)) = currentVis(IPosition(2, 3, k));

            for (int i=0; i<2; i++) {
                for (int j=0; j<2; j++) {
                    M(IPosition(2, i, j)) = Mvec[i][j][k];
                }
            }
            if (currConj) {
                b = product(M, e);
            } else {
                M = adjoint(M);
                b = product(e, M);
            }
            // Parallactic angle:
            pa = 0;
            int paind0, paind1;
            if (currConj) {
                paind0 = 0;
                paind1 = 1;
            } else {
                paind0 = 1;
                paind1 = 0;
            }
            pa(0, 0) = std::polar(1.0, parAng[paind0]);
            pa(1, 1) = std::polar(1.0, -parAng[paind1]);
            if (currConj) {
                if (doParang && parAng[0]>-1.e8) {
                    b = product(b, pa);
                }
            } else {
                if (doParang && parAng[0]>-1.e8) {
                    // Are we leaking memory when we do this?
                    pa = adjoint(pa);
                    b = product(pa, b);
                }
            }
            // Now transcribe b back into bufferVis
            bufferVis(IPosition(2, 0, k)) = b(0, 0);
            bufferVis(IPosition(2, 1, k)) = b(1, 1);
            bufferVis(IPosition(2, 2, k)) = b(0, 1);
            bufferVis(IPosition(2, 3, k)) = b(1, 0);
            if (!currConj) {
                e = adjoint(e);
                b = adjoint(b);
                // We have already taken the adjoint of M above.
            }
            if (a1 == 1 && a2 == 12) {
                
                std::complex<float> c = currentVis(IPosition(2, 0, k));
                std::cerr << arg(c) << " ";
            }
            for (int i=0; i<2; i++) {
                for (int j=0; j<2; j++) {
                    fwrite(&e(IPosition(2, i, j)), sizeof(std::complex<float>), 1, plotFile);
                }
            }
            for (int i=0; i<2; i++) {
                for (int j=0; j<2; j++) {
                    fwrite(&b(IPosition(2, i, j)), sizeof(std::complex<float>), 1, plotFile);
                }
            }
            for (int i=0; i<2; i++) {
                for (int j=0; j<2; j++) {
                    fwrite(&M(IPosition(2, i, j)), sizeof(std::complex<float>), 1, plotFile);
                }
            }
        } // End of k loop
        if (a1 == 1 && a2 == 12) {
            std::cerr << std::endl;
        }
        if (isTwoLinear) {
            // FIXME?
        }
    }
    // if (allNear(bufferVis, zeros, 1e-9)) {
    //     std::cerr << "bufferVis zeroed!" << std::endl;
    //     std::cerr << "parAng: " << parAng[0] << ", " << parAng[1] << std::endl;
    // }
}

// Close all files.
void DataIOMS::finish()
    // Measurement Set is closed by its destructor, so we don't do anything here.
{ }


void DataIOMS::saveCirculars() {
    std::cerr << "DataIOMS::saveCirculars()" << std::endl;
    std::cerr << "DataIOMS::saveCirculars(): IFsToDo again = [";
    for (Int i: IFsToDo) {
        std::cerr << i << ' ';
    }
    std::cerr << "]" << std::endl;
    
    char buffer[256];
    MSColumns msc(ms);
    Int saveSource = 1; // FIXME!
    if (!doWriteCirc) {
        std::cerr << "Not writing circulars this time?" << std::endl;
        return;
    }
    // AUXILIARY BINARY FILES TO STORE CIRCULAR VISIBILITIES:
    ofstream circFile[40];
    // = ofstream[IFsToDo.size()];

    int iif=0;
    std::cerr << "About to open some h*cking files!" << std::endl;

    for (Int spw : IFsToDo) {
        std::cerr << "About to open file for spw " << iif << std::endl;
        sprintf(buffer, "POLCONVERT.FRINGE/OTHERS.FRINGE_IF%d", iif+1);
        circFile[iif].open(buffer, ios::out | ios::binary);
        if (!circFile[iif].is_open()) {
            std::cerr << "Well, that went well!" << std::endl;
            exit(1);
        }
        circFile[iif].write((char*) &Freqs[spw].Nchan, sizeof(int));
        for (int j=0; j < Freqs[spw].Nchan; j++) {
            circFile[iif].write((char*) &Freqvals[spw][j], sizeof(double));
        }
        iif++;
    }
    std::cerr << "Times here: " << doRange[0] << ", " << doRange[1] << std::endl;
    int low = 0;
    int high = 0;
    int good = 0;
    for (casacore::rownr_t currRow=0; currRow<ms.nrow(); currRow++) {
        Double t = msc.time()(currRow);
        /* if ((currRow % 10000) == 0) {
            std::cerr << "Low " << low
                      << " High " << high
                      << " Good " << good
                      << " Now " << t
                      << " (Now - begin) " << t - doRange[0]
                      << " (end - Now) " << t - doRange[1]
                      << std::endl;
        } */
        if (t < doRange[0]) { low++; continue; }
        if (t > doRange[1]) { high++; continue; }
        good++;
        Int a1 = msc.antenna1()(0);
        Int a2 = msc.antenna2()(0);
        Bool lin1 = linAnts.find(a1)!=linAnts.end();
        Bool lin2 = linAnts.find(a2)!=linAnts.end();
        if (lin1 || lin2) {
            std::cerr << "Bailing because linants" << std::endl;
            continue;
        }
        Int dd = msc.dataDescId()(currRow);
        Int ispw = msc.dataDescription().spectralWindowId()(dd);
        bool found = false;
        int ifidx = 0;
        for (auto iif : IFsToDo) {
            if (iif==ispw) {
                found = true;
                break;
            }
            ifidx++;
        }
        if (!found) {
            std::cerr << "Bailing because spw" << std::endl;
            continue;
        }
        Int souidx = msc.fieldId()(currRow);
        /*
        if (souidx != saveSource && saveSource >= 0) {
            std::cerr << "Bailing because source" << std::endl;
            continue;
            } */
        ofstream& outfile(circFile[ifidx]);
        const Array<Double>& uvw = msc.uvw()(currRow);
        double duvw[3];
        for (int i=0; i<3; i++) {
            duvw[i] = uvw(IPosition(1, i));
        }
        std::cerr << "Actually writing" << std::endl;
        double AuxPA1, AuxPA2;
        getParAng(souidx, a1, a2, duvw, t, AuxPA1, AuxPA2);
        outfile.write(reinterpret_cast<char*>(&t), sizeof(double));
        int aips1 = a1+1;
        int aips2 = a2+1;
        outfile.write(reinterpret_cast<char*>(&aips1), sizeof(int));
        outfile.write(reinterpret_cast<char*>(&aips2), sizeof(int));
        outfile.write(reinterpret_cast<char*>(&AuxPA1), sizeof(double));
        outfile.write(reinterpret_cast<char*>(&AuxPA2), sizeof(double));
            
        Matrix<Complex> bufferVis = msc.data()(currRow);

        // Oh, what are the dimensions of Freqs?
        for (int j=0; j<Freqs[ispw].Nchan; j++){
            size_t sz = sizeof(std::complex<float>);
            IPosition a11(3, 0, 0, j);
            IPosition a22(3, 1, 0, j);
            IPosition a12(3, 2, 0, j);
            IPosition a21(3, 3, 0, j);
            
            outfile.write(reinterpret_cast<char*>(&bufferVis(a11)), sz);
            outfile.write(reinterpret_cast<char*>(&bufferVis(a12)), sz);
            outfile.write(reinterpret_cast<char*>(&bufferVis(a21)), sz);
            outfile.write(reinterpret_cast<char*>(&bufferVis(a22)), sz);
        }
    }
    // CLOSE AUXILIARY BINARY FILES:
    iif = 0;
    for (Int spw : IFsToDo){
        circFile[iif].close();
        iif++;
    }
}
