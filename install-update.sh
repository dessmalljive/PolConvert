#! /usr/bin/env sh
conda activate py3.9
rm -rf *.so build
python3 setup.py build_ext --inplace
