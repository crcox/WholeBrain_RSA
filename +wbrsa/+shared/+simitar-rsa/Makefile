MATLABDIR ?= /usr/local/MATLAB/R2013b
MCC:=/usr/local/MATLAB/R2013b/bin/mcc
MEX:=/usr/local/MATLAB/R2013b/bin/mex
MFLAGS:=-m -R -singleCompThread -R -nodisplay -R -nojvm
SRCDIR= ./
INCL= -I$(MATLABDIR)/extern/include
PROGS= findNeighbours simitar fastscoring
.PHONEY: all clean

all: $(PROGS)

$(PROGS): %: %.c
	$(MEX) $(INCL) $<

clean:
	-rm *.dmr
	-rm _condor_std???
	-rm readme.txt
	-rm mccExcludedFiles.log
