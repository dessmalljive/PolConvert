#include <Python.h>
#include "./_PolConvert.h"
#include "./_PolGainSolve.h"
#include "./_getAntInfo.h"
#include "./_XPCal.h"
#include "./_XPCalMF.h"

// _PolConvert docstrings 
static char module_docstring[] =
    "This module provides an interface for PolConvert.";

static char PolConvert_docstring[] =
    "Converts mixed-polarization visibilities into pure circular-polarization basis";
static char setPCMode_docstring[] =
    "Sets the running mode to either ALMA (true; the default) or non-ALMA (false).";

// _PolGainCal docstrings 
static char PolGainSolve_docstring[] =
    "Solves for cross-polarization gains from mixed-polarization visibilities";
static char ReadData_docstring[] =
    "Reads data (one IF) for PolGainSolve";
static char FreeData_docstring[] =
    "Releases the data pointers of PolGainSolve";
static char GetChi2_docstring[] =
    "Computes the Chi2 for a given set of cross-pol gains";
static char GetIFs_docstring[] =
    "Returns the array of frequencies for a given IF";
static char DoGFF_docstring[] =
    "Performs a simplified GFF (delays and rates). The reference antenna is set by not adding it to the list of fittable antennas";
static char SetFringeRates_docstring[] =
    "Forces the antenna fringe-rates to the list give, before the GCPFF is performed";
static char GetNScan_docstring[] =
    "Returns the number of scans for the given IF.";
static char SetFit_docstring[] =
    "Allocates memory for the GCPFF.";
static char GetNchan_docstring[] =
    "Returns the number of channels for the given IF.";

// _XPCal docstrings 
static char XPCal_docstring[] =
    "Reads and interpolates phasecal-tone phase values, given in DiFX format.\n The first argument is the name of the phasecal file.\n The second argument is an integer: \n 0 = return the scan-averaged cross-polarization phases. \n 1 = return the scan-averaged phases in the firs pol. channel, referred to the tone of lowest frequency. \n 3 = return the phasecals in the first pol. channel for each tone and integration time.\n\n The output is a list of numpy arrays: frequency (MHz), phase (deg), rate (Hz), reference time (MJD), time offsets (s), phasecal matrix (i.e., phases for each tone and time offset).\n\nThe third argument is a group delay (in ns) to be substracted from the phasecals before the phase connection among tones.";

static char XPConvert_docstring[] =
    "Overwrites the pcal DiFX file with a new version, where X and Y are substitued by R and L. In this version, both the R and L entries are exact copies of X.";

// _XPCalMF docstrings
static char XPCalMF_docstring[] =
    "Reads and interpolates phasecal-tone phase values, given in DiFX format. \n Supports multi-file. It can either return the time average of the cross-polarization phase \n or overwrite the pcal file with the polconverted solutions.\n";

// _getAntInfo docstrings
static char getAntInfo_docstring[] =
    "Returns antenna information";
static char getCoords_docstring[] =
    "Returns antenna coordinates";
static char getMounts_docstring[] =
    "Returns the mount types";


/* Module specification */
static PyMethodDef module_methods[] = {
    // _PolConvert
    {"PolConvert", PolConvert, METH_VARARGS, PolConvert_docstring},
    {"setPCMode", setPCMode, METH_VARARGS, setPCMode_docstring},
    // _PolGainSolve
    {"PolGainSolve", PolGainSolve, METH_VARARGS, PolGainSolve_docstring},
    {"ReadData", ReadData, METH_VARARGS, ReadData_docstring},
    {"GetChi2", GetChi2, METH_VARARGS, GetChi2_docstring},
    {"GetIFs", GetIFs, METH_VARARGS, GetIFs_docstring},
    {"DoGFF", DoGFF, METH_VARARGS, DoGFF_docstring},
    {"SetFringeRates", SetFringeRates, METH_VARARGS, SetFringeRates_docstring},
    {"GetNScan",GetNScan, METH_VARARGS, GetNScan_docstring},
    {"GetNchan",GetNchan, METH_VARARGS, GetNchan_docstring},
    {"FreeData", FreeData, METH_VARARGS, FreeData_docstring},
    {"SetFit", SetFit, METH_VARARGS, SetFit_docstring},
    // _XPCal 
    {"XPCal", XPCal, METH_VARARGS, XPCal_docstring},
    {"XPConvert", XPConvert, METH_VARARGS, XPConvert_docstring},
    // _XPCalMF
    {"XPCalMF", XPCalMF, METH_VARARGS, XPCalMF_docstring},
    // _getAntInfo
    {"getAntInfo", getAntInfo, METH_VARARGS, getAntInfo_docstring},
    {"getCoords", getCoords, METH_VARARGS, getCoords_docstring},
    {"getMounts", getMounts, METH_VARARGS, getMounts_docstring}, 
    {NULL, NULL, 0, NULL}   /* terminated by list of NULLs, apparently */
    
};

/* Initialize the module */

static struct PyModuleDef polconvert_module_def = {
    PyModuleDef_HEAD_INIT,
    "polconvert",          /* m_name */
    module_docstring,       /* m_doc */
    -1,                     /* m_size */
    module_methods,         /* m_methods */
    NULL,NULL,NULL,NULL     /* m_reload, m_traverse, m_clear, m_free */
};


PyMODINIT_FUNC PyInit_polconvert(void)
{
    PyObject *m = PyModule_Create(&polconvert_module_def);
    PCMode = true;
    return(m);
}
