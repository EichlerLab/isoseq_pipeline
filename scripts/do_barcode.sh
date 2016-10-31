#! /bin/bash
# 
# 1) Copies the files and scripts needed by the barcode trimmer to the working directory
# 2) Runs the barcode identifier
#
source ~mldough/my_scripts/barcode/drop_barcode_files.sh && python fluid_barcode_identification.py_log
