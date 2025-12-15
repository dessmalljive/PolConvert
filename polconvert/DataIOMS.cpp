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

   // Passed in and stored, and now very much our problem!
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
    doConjugate(doConj), // Not used in this class
    canPlot(false),
    doParang(doPar), // Referred to but never set in this class
    isTwoLinear(false) // Not used in this class
{

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
        // I'm just going to assume we want it writeable all the time.
        ms.reopenRW();
        std::cerr << "Opened the MS at least" << std::endl;
    } catch (casacore::TableNoFile t) {
        std::cout << "Measurement set " << msfn << " not found" << std::endl;
        success = false;
    }

    if (!success) {
        return;
    }

    MSColumns msc(ms);
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
  
    std::cerr << "Time for destruction time " << std::endl;
}

bool DataIOMS::setCurrentIF(int i)
{
    std::cerr << "DataIOMS::setCurrentIF(" << i <<")" << std::endl;
    currSPW = i;
    currRow = 0;
    nextRow = 0;
    currAntInd = 0;
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
    
    std::cerr << "DataIOMS::getFrequencies: currSPW= " << currSPW << std::endl;
    std::cerr << "DataIOMS::getFrequencies: i= " << i << std::endl;
}

// NB: everything is a reference here; these are all the *return* values!
bool DataIOMS::getNextMixedVis(double &JDTime, int &antenna, int &otherAnt, bool &conj, int &calField)
{
    MSColumns msc(ms);
    canPlot = true;
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
        const Int ant1 = msc.antenna1()(currRow);
        const Int ant2 = msc.antenna2()(currRow);

        Bool lin1 = linAnts.find(ant1)!=linAnts.end();
        Bool lin2 = linAnts.find(ant2)!=linAnts.end();
        
        if ((ant1 != ant2) && lin1 && lin2 && currAntInd==0) {
            cerr << "[getNextMixedVis] returning once from " << currRow << " for lin-lin " << endl;
            canPlot = False;
            isTwoLinear = True;
            antenna = ant1+1;
            otherAnt = ant2+1;
            conj = true;
            // For next time round
            // OK, here goes!
            // We now take the view that we have to return from this row twice
            // and let applyMatrix apply twice and generally have an interesting time
            nextRow = currRow;
            currAntInd=1;
            currConj = conj;
            nmixed++;
            return true;

        }
        // We handle the case that both are the same linear antenna in this branch
        if (currAntInd==0) {
            if (lin1) {
                isTwoLinear = false; // Because we handled that case above and exited with currAntInd==1
                cerr << "[getNextMixedVis] returning linear ant1 from "<< currRow << endl;
                antenna = ant1+1;
                otherAnt = ant2+1;
                conj = true;
                // For next time round
                nextRow = currRow;
                currAntInd=1;
                currConj = conj;
                nmixed++;
                return true;
            } else {
                // Whether we found a match for ant1, we look at ant2 next.
                currAntInd=1;
                // Don't increment row or return!
            }
        }
        // currAntInd==1 from here on
        if (lin2 && ant2 != ant1) {
            cerr << "[getNextMixedVis] returning linear ant2 from "<< currRow << endl;
            antenna = ant2+1;
            otherAnt = ant1+1;
            conj = false;
            nextRow = currRow+1;
            currAntInd = 0;
            currConj = conj;
            nmixed++;
            return true;
        } else {
            // We've checked both antennas; now we advance to the next row and
            // look for the next SPW match.
            currAntInd = 0;
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

    if (isTwoLinear && currConj) {
        std::cerr << "First time through we just don't" << std::endl;
    }
    
    // I would like to get the time, but digging the actual time out of the CASA format is a challenge.
    // Can I think of another way?
    const Array<Complex>& oldRow(msc.data()(currRow));
    msc.correctedData().put(currRow, bufferVis);
    if (1) {
        const Int dd = msc.dataDescId()(currRow);
        const Int ant1 = msc.antenna1()(currRow);
        const Int ant2 = msc.antenna2()(currRow);
        const Int ispw = msc.dataDescription().spectralWindowId()(dd);
        if (ispw==0 && ant1==10 && ant2==12) {
            std::cerr << "[setCurrentMixedVis] Row " << currRow
                      << " Before:"  << oldRow(IPosition(3, 10, 0)) << std::endl
                      << " after: "  << bufferVis(IPosition(3, 10, 0)) << std::endl;
        }
    }
    return true;
}

int DataIOMS::getFileNumber() { return 0; }
// ::zeroWeight is empty in DataIOFITS.cpp, so I guess I can leave it empty here too?
void DataIOMS::zeroWeight() {}

// Modify the visibilities read by "getNextMixedVis" by the calibration matrix supplied
// in Mvec
void DataIOMS::applyMatrix(std::complex<float> *Mvec[2][2], bool swap, bool print,
                           int thisAnt, FILE *plotFile)
{
    MSColumns msc(ms);
    const Matrix<Complex>& currentVis (msc.data()(currRow));
    IPosition visShape(currentVis.shape());
    // std::cerr << "bufferVis before " << bufferVis << std::endl;
    bufferVis.resize(visShape);
    // std::cerr << "bufferVis after " << bufferVis << std::endl; 
    Matrix<Complex> zeros(visShape);
    zeros = 0;
    
    Int a1 = msc.antenna1()(currRow);
    Int a2 = msc.antenna2()(currRow);

    const Double C = 299792458.0;
    Array<Double> uvw(msc.uvw()(currRow));
    double duvw[3];
    for (int i=0; i<3; i++) {
        duvw[i] = uvw(IPosition(1, i))/C; // Convert m to lightseconds
    }
    double uvdist = (duvw[0]*duvw[0] + duvw[1]*duvw[1]);
    Int souidx = msc.fieldId()(currRow);
    Double t = msc.time()(currRow);

    // We didn't store the parallactic angles so we have to get them again.
    double parAng[2];
    double AuxPA1, AuxPA2;
    getParAng(souidx, a1, a2, duvw, t, AuxPA1, AuxPA2);
    int aips1, aips2;

    // Sometimes we have to flip antenna AIPS orderings
    if (!currConj || isTwoLinear) {
        aips1 = a2+1;
        aips2 = a1+1;
        parAng[0] = AuxPA2;
        parAng[1] = AuxPA1;
    } else {
        aips1 = a1+1;
        aips2 = a2+1;
        parAng[0] = AuxPA1;
        parAng[1] = AuxPA2;
    }
    if (isTwoLinear) {
        std::cerr << "[ApplyMatrix] Both lin " << aips1 << " " << aips2
                  << " swap " << swap << " print " << print
                  << " currConj " << currConj 
                  << " canPlot " << canPlot << std::endl;
    }
    // OK there are going to be some things in here that lack elegance
    // But if we have a lin-lin combination then currConj is true the
    // first time and false the second, while neither of swap or
    // canPlot is actually set to anything useful in the MS version of
    // the code. That's probably entirely my fault but I am going to
    // start from where we are instead of going back to figure out
    // better values for those.
    bool firstOfBoth = (isTwoLinear && currConj);
    bool secondOfBoth = (isTwoLinear && !currConj);
    bool writeIt = (print && canPlot && (!firstOfBoth));
    if (writeIt) {
        int zero = 0;
        fwrite(&zero, sizeof(int), 1, plotFile);
        fwrite(&t, sizeof(double), 1, plotFile);
        fwrite(&aips1, sizeof(int), 1, plotFile);
        fwrite(&aips2, sizeof(int), 1, plotFile);
        fwrite(&parAng[0], sizeof(double), 1, plotFile);
        fwrite(&parAng[1], sizeof(double), 1, plotFile);
        fwrite(&uvdist, sizeof(double), 1, plotFile);
    }
    for (int k=0; k<Freqs[currSPW].Nchan; k++) {
        /* The default ordering of correlation pairs in FITS-IDI is (probably): RR,LL,RL,LR;
           for Measurement Set the convention is RR, RL, LR, LL. */
        // Visibility matrix:
        Matrix<Complex> e(IPosition(2,2,2));
        // Calibration Jones matrix
        Matrix<Complex> M(IPosition(2,2,2));
        // Matrix to store back in buffer
        Matrix<Complex> b(IPosition(2,2,2));
        e = 0;
        M = 0;
        b = 0;
        Matrix<Complex> e0(e);
        if (!secondOfBoth) {
            e0(IPosition(2, 0, 0)) = currentVis(IPosition(2, 0, k));
            e0(IPosition(2, 0, 1)) = currentVis(IPosition(2, 1, k));
            e0(IPosition(2, 1, 0)) = currentVis(IPosition(2, 2, k));
            e0(IPosition(2, 1, 1)) = currentVis(IPosition(2, 3, k));
            // The conjugation below also needs explaining! Maybe it's an AIPS/CASA sign
            // convention thing?
            e = conj(e0); 
        } else {
            // Assume we don't conjugate here
            e(IPosition(2, 0, 0)) = bufferVis(IPosition(2, 0, k));
            e(IPosition(2, 0, 1)) = bufferVis(IPosition(2, 1, k));
            e(IPosition(2, 1, 0)) = bufferVis(IPosition(2, 2, k));
            e(IPosition(2, 1, 1)) = bufferVis(IPosition(2, 3, k));
        }
           
        for (int i=0; i<2; i++) {
            for (int j=0; j<2; j++) {
                M(IPosition(2, i, j)) = Mvec[i][j][k];
            }
        }
        Matrix<Complex> e_out, M_out, b_out;
        // Parallactic angle Jones matrix
        Matrix<Complex> P0(IPosition(2,2,2)), P1(IPosition(2,2,2));
        P0 = 0;
        P1 = 0;
        std::complex<double> ph0 = std::polar(1.0, AuxPA1);
        P0(0, 0) = ph0;
        P0(1, 1) = 1.0/ph0;
        std::complex<double> ph1 = std::polar(1.0, AuxPA2);
        P1(0, 0) = ph1;
        P1(1, 1) = 1.0/ph1;
        // Why would parallactic angles be -1e8 anyway?!
        bool doParang0 = doParang && AuxPA1 > -1.e8;
        bool doParang1 = doParang && AuxPA2 > -1.e8;
        if (currConj) {
            b = product(M, e);
            if (doParang0) {
                b = product(P0, b);
            }
            b_out = b;
            e_out = e;
            M_out = M;
        } else { // !currConj
            b = product(e, adjoint(M));
            if (doParang1) {
                b = product(b, adjoint(P1));
            }
            b_out = adjoint(b);
            e_out = adjoint(e);
            M_out = adjoint(M);
        }
                

        // small@jive.eu: bufferVis is used for output in 
        // setCurrentMixedVis. Comparing with AIPS shows we need the
        // adjoint of b for this. Not currently sure why.
        Matrix<Complex> b_adjoint = adjoint(b);

        // This is a knob we can twiddle
        if (!firstOfBoth) {
            bufferVis(IPosition(2, 0, k)) = b_adjoint(0, 0);
            bufferVis(IPosition(2, 1, k)) = b_adjoint(0, 1);
            bufferVis(IPosition(2, 2, k)) = b_adjoint(1, 0);
            bufferVis(IPosition(2, 3, k)) = b_adjoint(1, 1);
        } else {
            bufferVis(IPosition(2, 0, k)) = b(0, 0);
            bufferVis(IPosition(2, 1, k)) = b(0, 1);
            bufferVis(IPosition(2, 2, k)) = b(1, 0);
            bufferVis(IPosition(2, 3, k)) = b(1, 1);
        }
            
        if (k == 10) {
            const Int dd = msc.dataDescId()(currRow);
            const Int ant1 = msc.antenna1()(currRow);
            const Int ant2 = msc.antenna2()(currRow);
            const Int ispw = msc.dataDescription().spectralWindowId()(dd);
            
            if (ispw==0 && ((ant1==10 && ant2==12) ||
                            (ant1==12 && ant2==10))) {
                std::cerr << "[applyMatrix] linAnts ";
                for (Int const& a : linAnts) { std::cerr << a << ' '; }
                std::cerr << std::endl;
                std::cerr <<" isTwoLinear " << isTwoLinear << std::endl;
                std::cerr << "[applyMatrix] Row " << currRow << " spw " << ispw << " currConj " << currConj 
                          << " baseline " << ant1 << "-" << ant2 << " M " << M << std::endl;
                std::cerr << " b " << b << " e " << e << " doParang1 " << doParang1 << std::endl;
                std::cerr << "b_adjoint " << b_adjoint << endl;
                std::cerr << "bufferVis() ";
                for (int i=0; i<4; i++) {
                    std::cerr << bufferVis(IPosition(2, i, k)) << " ";
                }
                std::cerr << std::endl;
                    
            }
        }

        
        if (writeIt) {
            for (int i=0; i<2; i++) {
                for (int j=0; j<2; j++) {
                    fwrite(&e_out(IPosition(2, i, j)), sizeof(std::complex<float>), 1, plotFile);
                }
            }
            for (int i=0; i<2; i++) {
                for (int j=0; j<2; j++) {
                    fwrite(&b_out(IPosition(2, i, j)), sizeof(std::complex<float>), 1, plotFile);
                }
            }
            for (int i=0; i<2; i++) {
                for (int j=0; j<2; j++) {
                    fwrite(&M_out(IPosition(2, i, j)), sizeof(std::complex<float>), 1, plotFile);
                }
            }
        }
    } // End of k loop
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
    // Int saveSource = 1; // FIXME!
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
        if (t < doRange[0]) { low++; continue; }
        if (t > doRange[1]) { high++; continue; }
        Int a1 = msc.antenna1()(currRow);
        Int a2 = msc.antenna2()(currRow);
        Bool lin1 = linAnts.find(a1)!=linAnts.end();
        Bool lin2 = linAnts.find(a2)!=linAnts.end();
        if (lin1 || lin2) {
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
        //     std::cerr << "Bailing because spw" << std::endl;
            continue;
        }
        Int souidx = msc.fieldId()(currRow);
        // FIXME: Do we really want to check source id?
        // if (souidx != saveSource && saveSource >= 0) {
        //    //std::cerr << "Bailing because source" << std::endl;
        //     continue;
        // } 
        good++;
        ofstream& outfile(circFile[ifidx]);
        const double C = 3e8;
        const Array<Double>& uvw = msc.uvw()(currRow);
        double duvw[3];
        for (int i=0; i<3; i++) {
            duvw[i] = uvw(IPosition(1, i))/C;
        }
        double AuxPA1, AuxPA2;
        getParAng(souidx, a1, a2, duvw, t, AuxPA1, AuxPA2);
        outfile.write(reinterpret_cast<char*>(&t), sizeof(double));
        int aips1 = a1+1;
        int aips2 = a2+1;
        outfile.write(reinterpret_cast<char*>(&aips1), sizeof(int));
        outfile.write(reinterpret_cast<char*>(&aips2), sizeof(int));
        outfile.write(reinterpret_cast<char*>(&AuxPA1), sizeof(double));
        outfile.write(reinterpret_cast<char*>(&AuxPA2), sizeof(double));
            
        const Matrix<Complex>& currentVis = msc.data()(currRow);

        for (int j=0; j<Freqs[ispw].Nchan; j++){
            size_t sz = sizeof(std::complex<float>);
            IPosition a11(2, 0, j);
            IPosition a12(2, 1, j);
            IPosition a21(2, 2, j);
            IPosition a22(2, 3, j);

            // Do conjugates here!
            Complex c11 = std::conj(currentVis(a11));
            Complex c12 = std::conj(currentVis(a12));
            Complex c21 = std::conj(currentVis(a21));
            Complex c22 = std::conj(currentVis(a22));
                        
            outfile.write(reinterpret_cast<char*>(&c11), sz);
            outfile.write(reinterpret_cast<char*>(&c12), sz);
            outfile.write(reinterpret_cast<char*>(&c21), sz);
            outfile.write(reinterpret_cast<char*>(&c22), sz);
        }
    }
    // CLOSE AUXILIARY BINARY FILES:
    iif = 0;
    for (Int spw : IFsToDo){
        circFile[iif].close();
        iif++;
    }
    std::cerr << "saveCircs good: " << good << std::endl;
}
