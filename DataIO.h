/* DATAIO - FITS-IDI interface to PolConvert

             Copyright (C) 2013-2022  Ivan Marti-Vidal
             Nordic Node of EU ALMA Regional Center (Onsala, Sweden)
             Max-Planck-Institut fuer Radioastronomie (Bonn, Germany)
             University of Valencia (Spain)  

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
#include <math.h>
#include <complex>
#include <vector>
#include "fitsio.h"
#include "XTHelpers.h"

#ifndef __DATAIO_H__
#define __DATAIO_H__


typedef std::complex<float> cplx32f;

typedef struct {
double Nu0;
int Nchan;
double BW;
bool SB;
} FreqSetup;


typedef struct {
xt::xtensor<double,1> BaseLine[3];
xt::xtensor<double,1> SinDec;
xt::xtensor<double,1> CosDec;
xt::xtensor<double,1> RA;
xt::xtensor<double,1> AntLon;
xt::xtensor<int,1> Mount;
xt::xtensor<double,1> Lat;
int NtotAnt;
int NtotSou;
xt::xtensor<int,2> BasNum;
} ArrayGeometry;


typedef struct {
    int AntIdx;
    int IF;
    int Pol;
    double JD;
    xt::xtensor<double,1> AC;
    
} AutoCorrelation;


//   int NtotAnt, NtotSou;
//  double *BaseLine[3], *SinDec, *CosDec, *AntLon, *TanLat, **BasNum;
//   double *ParAngs[2];


/* Class to read FITS-IDI or SWIN files, setup the data streams,
and iterate over all the baselines with mixed polarization. */

class DataIO {
  public:

   ~DataIO();

   DataIO(); 

   int getNfreqs();
   int getNant();

   int getNchan(int freqid);
   long getMixedNvis();
   bool succeed();

   double getDay0();

   void getFrequencies(double* Freqarray);

// Compute parallactic angle.
  void getParAng(int sidx, int Ant1, int Ant2, double*UVW, double &MJD, double &P1, double &P2);

// Estimate amplitude ratios from the autocorrelations:
  std::complex<float> getAmpRatio(int ant, int spw, int chan);
  
  
// Set the current IF (and reset the mixed-vis. counter):
   virtual bool setCurrentIF(int i) = 0;

// Get the file number of the current visibility (only useful for SWIN files; always returns 0 for FITS-IDI):
   virtual int getFileNumber() = 0;


/* Very important function. Finds the next combination of the 4 correlation
   products. Returns the visibilities as an array of pointers;
   It also returns the time of the visibility and the ids of the antennas in the baseline. 
   If more than one lin-pol antenna are present, this function will still work, and will 
   return the pure-linear visibilities twice, one iteration for each antenna (the same is 
   true for the auto-correlations of each lin-pol antenna). 
   Returns false if no more mixed-pol visibilities are found for the corresponding IF. */
   virtual bool getNextMixedVis(double &JDTime, int &antenna, int &otherAnt, bool &conj, int &calField) = 0;


/* Complementary to the function above. Once the visibilities have been 
   corrected, this function puts them back into the right place. */
   virtual bool setCurrentMixedVis() = 0;


 // Modify the visibilities read by "getNextMixedVis" by the calibration matrix supplied
 // Saves the result in the "bufferVis" pointer  
  virtual void applyMatrix(std::complex<float> *M[2][2], bool swap, bool print, int thisAnt, FILE *plotFile) = 0;


 // Flag bad (unconvertable) data:
  virtual void zeroWeight() = 0;

// Close all files.
   virtual void finish() = 0;



   static const int MAXIF = 256;
   std::vector<FreqSetup> Freqs;
   ArrayGeometry *Geometry;
   std::vector<AutoCorrelation> AutoCorrs;
// averAutocorrs[ant][freq] is a length-Nchan(freq) tensor (ragged across freqs).
   std::vector<std::vector<xt::xtensor<float,1>>> averAutocorrs;
// Freqvals[freq] is a length-Nchan(freq) tensor (ragged across freqs).
   std::vector<xt::xtensor<double,1>> Freqvals;
   double *doRange; // borrowed (not owned by DataIO/subclasses)
   xt::xtensor<double,1> JDTimes;
   xt::xtensor<int,1> Basels, an1, an2, linAnts, field;
   xt::xtensor<double,1> ParAng[2];
   xt::xtensor<double,1> UVDist;
   int NLinAnt, NIFs, Nband, status, Nants, Nfreqs, currFreq;
   long NLinVis, Nvis, currVis;
   xt::xtensor<long,1> indexes;
   bool success, currConj;
   xt::xtensor<bool,1> is1, is2, is1orig, is2orig;
   int Flux;
   xt::xtensor<int,1> NAV;
   int nautos;
   double day0;


};




#endif 

