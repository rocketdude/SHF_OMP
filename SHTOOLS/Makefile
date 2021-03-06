###################################################################################
#
#	make all
#		Compile program in the current directory. Optionally, one
#		can specify the parameters F95="my compiler" and 
#		F95FLAGS="my compiler flags". The default is to use "f95".
#
#	make all2
#		Variant of make all: LAPACK subroutine names have
#		an underscore appended to them in the source files in order to 
#		use FFTW and LAPACK	libraries with conflicting underscore 
#		conventions.
#
#	make all3
#		Variant of make all: FFTW subroutine names have
#		an underscore appended to them in the source files in order to 
#		use FFTW and LAPACK	libraries with conflicting underscore 
#		conventions.
#
#	make install
#		To be run after successfully making "all". Move the entire 
#		directory to usr/local/SHTOOLS2.x. Generally requires root priviledges. 
#
#	make clean
#		Remove the compile lib, module, and object files.
#
#	make examples
#		Compile example programs. Optionally, one
#		can specify the parameters F95="my compiler" and 
#		F95FLAGS="my compiler flags", which should be identical to
#		those used to make "all".
#
#	make remove-examples
#		Delete compiled example programs.
#
#	make run-examples
#		Run all examples programs.
#
#	make doc
#		Create the man and html-man pages from input POD files.
#		These are PRE-MADE in the distribution, so it shouldn't
#		be necessary to recreate these unless there is some kind
#		of problem.
#
#	make remove-doc
#		Remove the man and html-man pages.
#
#
#	Written by Mark Wieczorek (July 2012).
#
#####################################################################################

VERSION = 2.7

F95 = f95

SHELL=/bin/tcsh
MAKE = make
DOCDIR = src/doc
SRCDIR = src
EXDIR = examples


.PHONY: all all2 all3 install doc remove-doc clean getflags examples remove-examples run-examples
	
all: getflags
	$(MAKE) -C $(SRCDIR) -f Makefile all VERSION=$(VERSION) F95=$(F95) F95FLAGS="$(F95FLAGS)"
	@echo
	@echo MAKE SUCCESSFUL!
	@echo	
	@echo ---------------------------------------------------------------------------------------------------
	@echo Compile your code with the following flags:
	@echo
	@echo $(F95) $(MODFLAG) $(F95FLAGS) -Llibpath -lSHTOOLS$(VERSION) -lfftw3 -lm -O3 -s -m64
	@echo
	@echo where modpath and libpath are replaced with their respective paths.
	@echo ---------------------------------------------------------------------------------------------------
	@echo

all2: getflags
	$(MAKE) -C $(SRCDIR) -f Makefile all2 VERSION=$(VERSION) F95=$(F95) F95FLAGS="$(F95FLAGS)"
	@echo
	@echo MAKE SUCCESSFUL!
	@echo
	@echo ---------------------------------------------------------------------------------------------------
	@echo Compile your code with the following flags:
	@echo
	@echo $(F95) $(MODFLAG) $(F95FLAGS) -L<libpath> -lSHTOOLS$(VERSION) -lfftw3 -lm -O3 -s -m64
	@echo
	@echo where modpath and libpath are replaced with their respective paths.
	@echo ---------------------------------------------------------------------------------------------------
	@echo
	
all3: getflags
	$(MAKE) -C $(SRCDIR) -f Makefile all3 VERSION=$(VERSION) F95=$(F95) F95FLAGS="$(F95FLAGS)"
	@echo
	@echo MAKE SUCCESSFUL!
	@echo	
	@echo ---------------------------------------------------------------------------------------------------
	@echo Compile your code with the following flags:
	@echo
	@echo $(F95) $(MODFLAG) $(F95FLAGS) -L<libpath> -lSHTOOLS$(VERSION) -lfftw3 -lm -O3 -s -m64
	@echo
	@echo where modpath and libpath are replaced with their respective paths.
	@echo ---------------------------------------------------------------------------------------------------
	@echo
	
getflags:
ifeq ($(F95),f95)
# Default Absoft Pro Fortran flags
F95FLAGS ?= -m64 -O3 -YEXT_NAMES=LCS -YEXT_SFX=_
MODFLAG = -p modpath
endif

ifeq ($(F95),g95)
# Default g95 flags.
F95FLAGS ?= -m64 -O3 -fno-second-underscore 
MODFLAG = -Imodpath
endif

ifeq ($(F95),gfortran)
# Default gfortran flags
F95FLAGS ?= -m64 -O3
MODFLAG = -Imodpath
endif

ifeq ($(F95),ifort)
# Default intel fortran flags
F95FLAGS ?= -m64 -free -O3 -Tf
MODFLAG = -Imodpath
endif

ifeq ($(origin F95FLAGS), undefined)
F95FLAGS = -m64 -O3
MODFLAG = -Imodpath
endif


install:
	@cp -R ../SHTOOLS /usr/local/SHTOOLS$(VERSION)/
	@echo
	@echo ENTIRE DIRECTORY MOVED TO /usr/local/SHTOOLS$(VERSION)
	

doc: 
	$(MAKE) -C $(DOCDIR) -f Makefile VERSION=$(VERSION)
	@echo DOCUMENTATION SUCCESSFULLY CREATED


remove-doc:
	@rm -f man/man1/*.1
	@rm -f www/man/*.html
	@echo
	@echo REMOVED MAN AND HTML-MAN FILES
	

clean:
	$(MAKE) -C $(SRCDIR) -f Makefile VERSION=$(VERSION) clean
	@echo
	@echo REMOVED LIB, MODULE, AND OBJECT FILES

	
examples: getflags
	$(MAKE) -C $(EXDIR) -f Makefile all VERSION=$(VERSION) F95=$(F95) F95FLAGS="$(F95FLAGS)"
	@echo
	@echo MAKE OF EXAMPLES SUCCESSFUL

remove-examples:
	$(MAKE) -C $(EXDIR) -f Makefile clean
	@echo
	@echo REMOVED EXAMPLE EXECUTABLES

run-examples:
	$(MAKE) -C $(EXDIR) -f Makefile run-examples
	@echo
	@echo RAN ALL EXAMPLE PROGRAMS

	