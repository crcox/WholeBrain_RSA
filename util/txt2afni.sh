#!/bin/bash -e

INFILE=$1
tmp=${INFILE%.*}
OUTFILE=${tmp##*/}
mnimaster="/home/chris/Manchester/CommonBrains/MNI_EPI_funcRes.nii"

3dUndump -master "$mnimaster" -xyz -datum float -prefix $OUTFILE $INFILE 2>tmp
