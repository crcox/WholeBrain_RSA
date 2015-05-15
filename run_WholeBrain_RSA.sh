#!/bin/sh
# script for execution of deployed applications
#
# Sets up the MCR environment for the current $ARCH and executes 
# the specified command.
#

## Download the runtime environment from SQUID
wget -q "http://proxy.chtc.wisc.edu/SQUID/r2013b.tar.gz"
tar xzf r2013b.tar.gz
rm r2013b.tar.gz

## Download all large data files listed in URLS from SQUID
mkdir data/
cd data/
while read url; do
  wget -q "http://proxy.chtc.wisc.edu/SQUID/${url}"
done < URLS
while read url; do
  wget -q "http://proxy.chtc.wisc.edu/SQUID/${url}"
done < URLS_SHARED
cd ../

echo "BEFORE RUNNING"
ls -l

exe_name=$0
exe_dir=`dirname "$0"`
echo "------------------------------------------"
echo Setting up environment variables
MCRROOT="v82"
echo ---
LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
XAPPLRESDIR=${MCRROOT}/X11/app-defaults ;
export XAPPLRESDIR;
export LD_LIBRARY_PATH;
echo LD_LIBRARY_PATH is ${LD_LIBRARY_PATH};
eval "${exe_dir}/WholeBrain_RSA"

echo "AFTER RUNNING"
ls -l

exit
