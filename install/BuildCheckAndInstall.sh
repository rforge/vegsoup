R CMD build /Users/roli/Documents/vegsoup/pkg
R CMD check vegsoup_0.1-3.tar.gz
R CMD INSTALL -l /Users/roli/Library/R/2.13/library vegsoup_0.1-3.tar.gz
