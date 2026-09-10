from distutils.core import setup, Extension
import numpy as np
import glob

printM  = f"""
#######################################################################
# Compiling with numpy version {np.__version__}
#                              {(np.__file__)}
#######################################################################
"""
print(printM)

## CHANGE IF NEEDED:
cfitsio='/usr/include/cfitsio'


sourcefiles0 = ['polconvert.cpp',
               'CalTable.cpp', 'DataIO.cpp', 'DataIOFITS.cpp',
               'DataIOSWIN.cpp', 'DataIOMS.cpp', 'Weighter.cpp', 'util.cpp',
               '_PolConvert.cpp', '_PolGainSolve.cpp',
               '_getAntInfo.cpp', # FIXME: this is only for FITS and it has a shared definition of message
               '_XPCal.cpp',
               '_XPCalMF.cpp']

PCPath =  'PolConvert/PC/'
sourcefiles = [PCPath + fn for fn in sourcefiles0]


casacore_includes = [
  '/usr/include/'
  ]

casacore_libs = ['casa_casa', 'casa_ms',
                 'casa_tables', 'casa_measures',
                 'casa_scimath', 'casa_scimath_f']


compile_args = ["-Wno-deprecated","-O0","-std=c++11", "-g"]
link_args = ["-Xlinker", "-export-dynamic"]


include_dirs = [PCPath] + [np.get_include()] + casacore_includes

print(f"{include_dirs=}")
print(f"{sourcefiles=}")
setup(
  ext_modules=[
      Extension(
          name="PolConvert.polconvert",
          sources=sourcefiles,
          libraries=['cfitsio', 'fftw3'] + casacore_libs,
          include_dirs=include_dirs,
          extra_compile_args=compile_args,
          extra_link_args=link_args
      ),
  ]
)

