/* getAntInfo -
             Copyright (C) 2017-2022  Ivan Marti-Vidal
             Centro Astronomico de Yebes (Spain)
             Universitat de Valencia (Spain)
  
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
#include <Python.h>

#if PY_MAJOR_VERSION >= 3
#define NPY_NO_DEPRECATED_API 0x0
#endif
#include <numpy/npy_common.h>
#include <numpy/arrayobject.h>

#include "fitsio.h"
#include "./SWIG.h"
#include "./_getAntInfo.h"

int Nants = 1;
double *Coords = new double[3];
int *Mounts = new int[1];

PyObject *getAntInfo(PyObject *self, PyObject *args){
  PyObject *IDI;
  int status, ist, imt, iAux;
  long lAux, jj;

  status = 0;

  if (!PyArg_ParseTuple(args, "O", &IDI)){
     printf("Failed initialization of getAntInfo! Check inputs!\n"); 
     PyObject *ret = Py_BuildValue("i", -1);
     return ret;
  }

  std::string fname;
  fname = PyString_AsString(IDI);

  fitsfile *ifile;
  char ARRAY[] = "ARRAY_GEOMETRY";
  char STABXYZ[] = "STABXYZ";
  char MNTSTA[] = "MNTSTA";

  fits_open_file(&ifile, fname.c_str(), READONLY, &status);

  if (status){
     printf("\n\nPROBLEM OPENING FILE!  ERR: %i\n\n",status);
     return Py_BuildValue("i",1);
  }
  fits_movnam_hdu(ifile, BINARY_TBL, ARRAY, 1, &status);

  if (status){
     printf("\n\nPROBLEM ACCESSING ARRAY INFO!  ERR: %i\n\n",status); 
     return Py_BuildValue("i",1);
  };

  fits_get_num_rows(ifile, &lAux, &status);
  if(status){
    printf("\n\nPROBLEM ACCESSING ARRAY INFO!  ERR: %i\n\n",status);
    return Py_BuildValue("i",1);
  }

  Nants = (int) lAux;

// Allocate memory:
  delete Coords;
  delete Mounts;

  Coords = new double[3*Nants];
  Mounts = new int[Nants];
 
// Read coordinates:
  fits_get_colnum(ifile, CASEINSEN, STABXYZ, &ist, &status);
  if(status){
    printf("\n\nPROBLEM ACCESSING ARRAY INFO!  ERR: %i\n\n",status);
    return Py_BuildValue("i",1);
  }
  fits_get_colnum(ifile, CASEINSEN, MNTSTA, &imt, &status);
  if(status){
    printf("\n\nPROBLEM ACCESSING ARRAY INFO!  ERR: %i\n\n",status);
    return Py_BuildValue("i",1);
  }

  for (jj=0;jj<Nants;jj++){
    fits_read_col(ifile, TDOUBLE, ist, jj+1, 1, 3, NULL, &Coords[3*jj], &iAux, &status);
    if(status){
      printf("\n\nPROBLEM ACCESSING ARRAY COORDINATES DATA!  ERR: %i\n\n",status);
      return Py_BuildValue("i",2);
    }
    fits_read_col(ifile, TINT, imt, jj+1, 1, 1, NULL, &Mounts[jj], &iAux, &status);
    if(status){
      printf("\n\nPROBLEM ACCESSING ARRAY MOUNTS DATA!  ERR: %i\n\n",status);
      return Py_BuildValue("i",2);
    }
    printf("ANTENNA %2li: X=%-9.1f Y=%-9.1f Z=%-9.1f | MOUNT: %i\n",
       jj,Coords[3*jj],Coords[3*jj+1],Coords[3*jj+2],Mounts[jj]);
  }

  fits_close_file(ifile, &status);

  if(status){
    printf("\n\nPROBLEM CLOSING FITS-IDI!  ERR: %i\n\n",status);
    return Py_BuildValue("i",3);
  }
  return Py_BuildValue("i",0);
}

PyObject *getCoords(PyObject *self, PyObject *args){
  // Build numpy array:
  PyObject *CoordArr;
  int nd = 2;
  npy_intp* CD = new npy_intp[nd];
  CD[0] = Nants; CD[1] = 3;
  CoordArr = PyArray_SimpleNewFromData(nd, CD, NPY_DOUBLE, (void*) Coords);
  return CoordArr;
}

PyObject *getMounts(PyObject *self, PyObject *args){
  // Build numpy array:
  PyObject *MountArr;
  int nd = 1;
  npy_intp* MD = new npy_intp[nd];
  MD[0] = Nants;
  MountArr = PyArray_SimpleNewFromData(nd, MD, NPY_INT, (void*) Mounts);
  return MountArr;
}


// eof
