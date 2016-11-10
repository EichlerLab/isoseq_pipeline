#! /bin/bash
# Copies the files and scripts needed by the barcode trimmer to the working directory
#
# Copes from ~mldough/my_scripts/barcode/
# Files include:
# 	classify_barcode.py
# 	custom_barcode_primers_forward.fa
# 	custom_barcode_primers_reverse.fa
# 	fluid_barcode_identification.py_log
# 	data/PBMATRIX.txt
cp ~mldough/my_scripts/barcode/classify_barcode.py ~mldough/my_scripts/barcode/custom_barcode_primers_forward.fa ~mldough/my_scripts/barcode/custom_barcode_primers_reverse.fa ~mldough/my_scripts/barcode/fluid_barcode_identification.py_log `pwd`
cp -r ~mldough/my_scripts/barcode/data `pwd`
