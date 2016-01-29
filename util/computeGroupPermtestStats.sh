#!/bin/bash
set -e

TMPDIR=$(mktemp -dt "$(basename $0).XXXXXXXXXX")

# Takes a list of directories as arguments.
# Directories contain permutation results for a group of subjects whole all had
# the similarity matrix permuted in the same way (i.e., based on the same
# random seed).

TOP=$(pwd)
for d in "$@"
do
  echo $d
  cd $d
  ${HOME}/src/WholeBrain_RSA/util/txt2afni.sh *_nodestrength.mni
  3dmerge -prefix group.mean.nodestrength.b4 -1blur_fwhm 4 *_nodestrength+tlrc.HEAD
  cd $TOP
  if [ -f group.mean.nodestrength.b4+tlrc.HEAD ]; then
    3dbucket -fbuc -glueto group.mean.nodestrength.b4+tlrc ${d}/group.mean.nodestrength.b4+tlrc
  else
    3dbucket -fbuc -prefix group.mean.nodestrength.b4+tlrc ${d}/group.mean.nodestrength.b4+tlrc
  fi
done

3dTstat -mean -sos -l2norm -prefix perm.group.stats.nodestrength.b4 group.mean.nodestrength.b4+tlrc
