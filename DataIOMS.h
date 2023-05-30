/* DATAIO - Measurement Set interface to PolConvert

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
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/measures/Measures/MeasTable.h>
#include <set>
#include "DataIO.h"


/* Class to read Measurement Sets, setup the data streams,
   and iterate over all the baselines with mixed polarization. */

class DataIOMS: public DataIO {
public:

     // Override:
     virtual void getFrequencies(double* output);

     ~DataIOMS();

     DataIOMS(std::string msfn, int Nant, int *Ants, double *doRange, int nIF2Conv,
              int *IF2Conv, bool Overwrite, bool doConj, bool doSolve, 
              ArrayGeometry *Geom, bool doPar, FILE *logF);

     bool setCurrentIF(int i);

     /* Very important function. Finds the next combination of the 4 correlation
        products. Returns the visibilities as an array of pointers;
        It also returns the time of the visibility and the ids of the antennas in the baseline. 
        If more than one lin-pol antenna are present, this function will still work, and will 
        return the pure-linear visibilities twice, one iteration for each antenna (the same is 
        true for the auto-correlations of each lin-pol antenna). 
        Returns false if no more mixed-pol visibilities are found for the corresponding IF. */
     bool getNextMixedVis(double &JDTime, int &antenna, int &otherAnt, bool &conj, int &calField);

     /* Complementary to the function above. Once the visibilities have been 
        corrected, this function puts them back into the right place. */
     bool setCurrentMixedVis();

     // Dummy function:
     int getFileNumber();

     // Modify the visibilities read by "getNextMixedVis" by the calibration matrix supplied
     // Saves the result in the "bufferVis" pointer  
     void applyMatrix(std::complex<float> *M[2][2], bool swap, bool print, int thisAnt, FILE *plotFile);

     // Flag bad (unconvertible) data:
     void zeroWeight();

     // Close all files.
     void finish();

     void log(char*);


     
private:
     casacore::Array<casacore::Complex> bufferVis;
     casacore::rownr_t currRow, nextRow;
     casacore::Int currSPW;
     casacore::Int currAnt;
     std::vector<casacore::Int> IFsToDo;
     bool doneFirst;
     bool doWriteCirc;
     
     std::set<casacore::Int> linAnts;
     std::vector<casacore::Double> doRange;
     FILE *logFile;
     long *Vis2Save;
     long NVis2Save;
     int nmixed;
     
     bool doConjugate, canPlot, doParang, isTwoLinear;
     casacore::MeasurementSet ms;

     
     void saveCirculars(void);   
     
};
