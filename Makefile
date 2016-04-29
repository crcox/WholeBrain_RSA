export MATLABDIR = /usr/local/MATLAB/R2013b
export MCC=$(MATLABDIR)/bin/mcc
export MEX=$(MATLABDIR)/bin/mex
export MFLAGS=-m -R -singleCompThread -R -nodisplay -R -nojvm -N
TOP := $(shell pwd)
SRCTAR=source_code.tar.gz
SRC=src
DEP=dependencies
JSON=$(DEP)/jsonlab
SIMITAR=$(DEP)/simitar-rsa
INCL= -I $(SRC) -I $(JSON) -I $(SIMITAR) 
.PHONEY: all clean-all clean-simitar clean-postbuild simiart sdist

all: setup simitar WholeBrain_RSA clean-postbuild

setup: $(SRC) $(DEP)

$(SRC) $(DEP):
	tar xzvf $(SRCTAR)

simitar:
	$(MAKE) -C $(SIMITAR)

WholeBrain_RSA: $(SRC)/WholeBrain_RSA.m $(SIMITAR)/simitar.mexa64
	$(MCC) -v $(MFLAGS) $(INCL) -o $@ $<

clean-postbuild:
	rm *.dmr
	rm mccExcludedFiles.log
	rm readme.txt
	rm run_WholeBrain_RSA.sh

sdist:
	tar czhf $(SRCTAR) src dependencies

clean-all:
	-rm WholeBrain_RSA
	-rm $(SIMITAR)/simitar.mexa64
	-rm $(LIBSVM)/findneighbours.mexa64
	-rm $(LIBSVM)/fastscoring.mexa64
