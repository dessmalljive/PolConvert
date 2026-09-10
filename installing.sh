conda config --add channels conda-forge
conda create --name PCONV python=3.9 scipy=1.6.3 numpy=1.20.2 matplotlib=3.4.1
conda activate PCONV

conda create --name PCONV2 python=3.10

# I am having a miserable time, as usual, with this stuff
conda create --name PCONV3 python=3.9 scipy=1.6.3 numpy=1.20.2 matplotlib=3.4.1
conda activate PCONV3
conda install pkgw-forge::casa-tools
conda install astropy

# The above won't compile the PolConvert code I wrote, which isn't great.


# In file included from DataIOMS.cpp:42:
# ./DataIOMS.h:78:16: error: ‘rownr_t’ in namespace ‘casacore’ does not name a type
#    78 |      casacore::rownr_t currRow, nextRow;
#       |                ^~~~~~~

conda create --name PCONV4 python=3.10 scipy 'numpy<2' matplotlib
conda activate PCONV4
conda install tomli
conda install python-casacore

# There is already a complete set of data in 
/home/small/miniconda310/lib/casa/data/


# So of course it doesn't work, because it can't find it.
# Specifically it can't find the TAI_UTC table, which is installed in MANY places, including but very much not limited to:

# /home/small/miniconda310/lib/casa/data/geodetic/TAI_UTC
# /home/small/miniconda310/lib/python3.10/site-packages/casadata/__data__/geodetic/TAI_UTC

# The data needs to be in 
# /home/small/miniconda310/envs/PCONV4/lib/casa/data/
#
# So...

cp -r /home/small/miniconda310/lib/casa /home/small/miniconda310/envs/PCONV4/lib/

# Now it runs until it core dumps (still in the first pass I think). Good* times!


