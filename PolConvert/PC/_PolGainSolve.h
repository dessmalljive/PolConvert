/* Available functions */

// This may be for internal use only
bool solveSystem(int Neq, double *HessianOrig, double *ResidualsOrig, double *Solution, double *Errors);

// Methods for the Python object we masquerade as:
PyObject *PolGainSolve(PyObject *self, PyObject *args);
PyObject *ReadData(PyObject *self, PyObject *args);
PyObject *GetChi2(PyObject *self, PyObject *args);
PyObject *GetIFs(PyObject *self, PyObject *args);
PyObject *GetNchan(PyObject *self, PyObject *args);
PyObject *DoGFF(PyObject *self, PyObject *args);
PyObject *SetFringeRates(PyObject *self, PyObject *args);
PyObject *GetNScan(PyObject *self, PyObject *args);
PyObject *FreeData(PyObject *self, PyObject *args);
PyObject *SetFit(PyObject *self, PyObject *args);

