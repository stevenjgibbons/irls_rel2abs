#
# Makefile for irls_rel2abs.
# Steven J Gibbons 2026/01/11 (Oslo)
# Before trying to compile, please ensure that the following
# lines have the correct locations of LAPACK, and BLAS
# BINDIR must be set to the directory in which you want the executable to reside
#
LAPACK= /usr/lib/x86_64-linux-gnu/liblapack.a
BLAS= /usr/lib/x86_64-linux-gnu/libblas.a
TOPDIR=   /home/stg/SRC
# BINDIR=  $(TOPDIR)/BIN
BINDIR=  .
#
# PLEASE CHECK ALL THE ABOVE LINES FOR YOUR SYSTEM ----
#
PROGNAME= irls_rel2abs
ALLSOURCECODE=  \
   $(PROGNAME).f 
#
SOURCES= \
        $(ALLSOURCECODE)
#
OPTIM=	  -O3
EXEFILE= $(BINDIR)/$(PROGNAME)
FORTRAN= gfortran
#
LIBS=    $(LAPACK) $(BLAS)
#
backup:
	cp -ip $(ALLSOURCECODE) ./BACKUP ; \
	cd ./BACKUP ; \
	\rm -f *.gz ; \
	gzip $(ALLSOURCECODE) ; \
	cd ../
#
$(PROGNAME):	$(ALLSOURCECODE) $(LIBS)
	$(FORTRAN) -o $(EXEFILE) $(ALLSOURCECODE) $(LIBS) $(OPTIM)
