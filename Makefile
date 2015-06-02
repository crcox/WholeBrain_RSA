MCC=/usr/local/MATLAB/R2013b/bin/mcc
MEX=/usr/local/MATLAB/R2013b/bin/mex
MFLAGS=-m -R -singleCompThread -R -nodisplay -R -nojvm
SRCDIR=.
IDIRS=
.PHONEY: clean clean-all all source_code.tar.gz extract

all: WholeBrain_RSA binaries.tar.gz

sdist:
	tar czhf source_code.tar.gz src dependencies ProxSortedL1

extract:
	-mkdir source_code/
	tar xzf source_code.tar.gz -C ./source_code/

WholeBrain_RSA: $(SRCDIR)/WholeBrain_RSA.m
	$(MCC) $(MFLAGS) $(IDIRS) -o $@ WholeBrain_RSA.m 

binaries.tar.gz: WholeBrain_RSA run_WholeBrain_RSA.sh
	mkdir bin/
	mv WholeBrain_RSA WholeBrain_RSA.sh bin/
	tar czvf $@ bin/

clean:
	-rm *.dmr
	-rm _condor_std???
	-rm readme.txt
	-rm mccExcludedFiles.log

clean-all: binaries.tar.gz
	rm -rf src/ bin/
