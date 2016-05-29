#!/bin/bash
set -e

mnimaster="/home/chris/MRI/Manchester/data/CommonBrains/MNI_EPI_funcRes.nii"

for INFILE in "$@"
do
  OUTFILE=$(basename ${INFILE%.*})
  3dUndump -master "$mnimaster" -xyz -datum float -prefix $OUTFILE $INFILE 2>tmp
done
