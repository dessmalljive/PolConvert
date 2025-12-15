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
#include <stdlib.h>  
#include <string.h>
#include <math.h>
#include <dirent.h>
#include "./DataIO.h"
#include "fitsio.h"





DataIO::~DataIO() {

};

DataIO::DataIO() :
    currConj(false)
{ printf("\nCreating VLBI data structure"); nautos=0; };


// SELF-EXPLANATORY FUNCTIONS:
int DataIO::getNfreqs() {return Nfreqs;};  // NUMBER OF IFs
int DataIO::getNant() {return Nants;};  // TOTAL NUMBER OF ANTENNAS
int DataIO::getNchan(int freqid) {
    int n = Freqs[freqid].Nchan; 
    std::cerr << "Someone asked DataIO::getMixedNchan() :" <<  n  << std::endl;
    return n;

};  // # OF CHANNEL IN freqid IF


long DataIO::getMixedNvis() {
    std::cerr << "Someone asked DataIO::getMixedNvis()" << std::endl;
    return NLinVis;
}  // NUMBER OF MIXED-POLARIZATION VISIBILITIES FOUND


// GET FREQUENCIES OF CURRENT IF:
void DataIO::getFrequencies(double* output){
    std::cerr << "Mysteriously here " << std::endl;
    memcpy(output, Freqvals[currFreq], Freqs[currFreq].Nchan*sizeof(double));
};

// CHECK WHETHER SOMETHING FAILED:
bool DataIO::succeed(){return success;};

// GET FIRST DAY OF OBSERVATION (MODIFIED JULIAN DATE)
double DataIO::getDay0(){return day0;};

cplx32f DataIO::getAmpRatio(int ant, int spw, int chan){
    return cplx32f(averAutocorrs[ant][spw][chan],0.);
};

void DataIO::getParAng(int sidx, int Ant1, int Ant2, 
                       double*UVW, double &MJD, double &P1, double &P2) {

    double GMST;
    double CT1, CT2, HAng, H1, H2, Elev1, Elev2;

    Elev1 = 0.0; Elev2 = 0.0;

    const double pi = 3.1415926535;
    double days = MJD/86400.;
    double t = (days-51544.0)/36525.;
    double Hh = days - floor(days);
    // Where the h*ck does this come from?
    double GMsec = 24110.54841 + 8640184.812866*t + 0.093104*t*t - 0.0000062*t*t*t;
    GMST = (GMsec/86400. + Hh)*2.*pi;
    
    if (sidx < Geometry->NtotSou && Ant1 < Geometry->NtotAnt && Ant2 < Geometry->NtotAnt){
        CT1 = Geometry->CosDec[sidx]*tan(Geometry->Lat[Ant1]);
        CT2 = Geometry->CosDec[sidx]*tan(Geometry->Lat[Ant2]);

        HAng = GMST - Geometry->RA[sidx];
        H1 = HAng + Geometry->AntLon[Ant1];
        H2 = HAng + Geometry->AntLon[Ant2];

        if (Geometry->Mount[Ant1] > 3){
            Elev1 = asin(sin(Geometry->Lat[Ant1])*Geometry->SinDec[sidx] +
                         cos(Geometry->Lat[Ant1])*Geometry->CosDec[sidx]*cos(H1));};

        if (Geometry->Mount[Ant2] > 3){
            Elev2 = asin(sin(Geometry->Lat[Ant2])*Geometry->SinDec[sidx] +
                         cos(Geometry->Lat[Ant2])*Geometry->CosDec[sidx]*cos(H2));};
 

        switch (Geometry->Mount[Ant1]){
        case 0:  P1 = atan2(sin(H1), CT1 - Geometry->SinDec[sidx]*cos(H1)); break; // ALT-AZ
        case 1: P1 = 0.; break; // EQ
        case 2: P1 = 0.; break; // ORBITAL (NO WAY!)
        case 3: P1 = atan2(cos(H1), Geometry->SinDec[sidx]*sin(H1)); break; // X-Y (E-W?)
        case 4: P1 = atan2(sin(H1), CT1 - Geometry->SinDec[sidx]*cos(H1)) + Elev1; break; // NA-R
        case 5: P1 = atan2(sin(H1), CT1 - Geometry->SinDec[sidx]*cos(H1)) - Elev1; break; // NA-L
        default: P1 = 0.0;
        };

        switch (Geometry->Mount[Ant2]){
        case 0:  P2 = atan2(sin(H2), CT2 - Geometry->SinDec[sidx]*cos(H2));  break; // ALT-AZ
        case 1: P2 = 0.;  break; // EQ
        case 2: P2 = 0.;  break; // ORBITAL (NO WAY!)
        case 3: P2 = atan2(cos(H2), Geometry->SinDec[sidx]*sin(H2));  break; // X-Y (E-W?)
        case 4: P2 = atan2(sin(H2), CT2 - Geometry->SinDec[sidx]*cos(H2)) + Elev2;  break; // NA-R
        case 5: P2 = atan2(sin(H2), CT2 - Geometry->SinDec[sidx]*cos(H2)) - Elev2;  break; // NA-L
        default: P2 = 0.0;
        };
        // std::cerr << "Ant 1= " << Ant1 << " Ant2= " << Ant2 << " sidx= " << sidx << "Mount[Ant1]= " << Geometry->Mount[Ant1] << " Mount[Ant2]= " << Geometry->Mount[Ant2] << std::endl;
        // std::cerr << "P1= " << P1 << " P2= " << P2 << std::endl;
        // std::cerr << "Time is now: " << t 
        //           << " Hh " << Hh 
        //           << " GMsec " << GMsec
        //           << std::endl;

        // std::cerr << "Geometry->Lat[Ant1] " << Geometry->Lat[Ant1] << std::endl 
        //           << "Geometry->AntLon[Ant1] " << Geometry->AntLon[Ant1] << std::endl
        //           << "Geometry->Lat[Ant2] " << Geometry->Lat[Ant2] << std::endl 
        //           << "Geometry->AntLon[Ant2] " << Geometry->AntLon[Ant2] << std::endl
        //           << "Geometry->RA[sidx] " << Geometry->RA[sidx] << std::endl; 
        // for (int s=0; s < Geometry->NtotSou; s++) {
        //     std::cerr << "Source " << s << " RA " << Geometry->RA[s] << std::endl;
        // }

    } else {
        P1 = 0.0;
        P2 = 0.0;
    }
}


