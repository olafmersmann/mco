.PHONEY: clean test check build install pkg data

install: clean
	R CMD INSTALL --no-multiarch pkg

#test: install
#	Rscript pkg/inst/unittests/runner.r

check: clean
	R CMD check pkg && rm -fR pkg.Rcheck

clean:
	rm -fR pkg/src/*.o pkg/src/*.so pkg.Rcheck .RData .Rhistory

pkg: clean
	git log --no-merges -M --date=iso --format=medium pkg/ > pkg/ChangeLog
	R CMD build pkg
	R CMD build --binary pkg
	rm -f pkg/ChangeLog
