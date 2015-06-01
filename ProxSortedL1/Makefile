MCC=/usr/local/MATLAB/R2013b/bin/mcc
MEX=/usr/local/MATLAB/R2013b/bin/mex
MFLAGS=-m -R -singleCompThread -R -nodisplay -R -nojvm
SRCDIR=.
IDIRS=
.PHONEY: clean all

all: proxSortedL1Mex.mexa64

proxSortedL1Mex.mexa64: proxSortedL1.c proxSortedL1Mex.c
	$(MEX) $<

clean:
	-rm *.dmr
	-rm _condor_std???
	-rm readme.txt
	-rm mccExcludedFiles.log
