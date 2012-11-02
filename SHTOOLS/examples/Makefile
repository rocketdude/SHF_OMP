######################################################################
#
#	make all
#		Compile example programs in their directories.
#
# 	make remove-examples
#		Delete the compile executibles.
#
#	make run-examples
#
#	Written by Mark Wieczorek (June 2012).
#
#######################################################################

F95 = f95

SHELL=/bin/tcsh

MAKE = make

.PHONY: all getflags clean


all: getflags

	@echo
	@echo -------------------------------------------------------------------------------------------------------
	@echo Compiling example code in SHTOOLS/examples with the following flags:
	@echo
	@echo $(F95) $(MODFLAG) $(F95FLAGS) -L../lib -lSHTOOLS$(VERSION) -lfftw3 -lm -O3 -s -m64
	@echo -------------------------------------------------------------------------------------------------------
	@echo
	
	$(F95) $(MODFLAG) -L../lib -lSHTOOLS$(VERSION) -lfftw3 -lm -O3 -s $(F95FLAGS) fftw_wisdom/create_wisdom.f95 -o fftw_wisdom/create_wisdom
	
	$(F95) $(MODFLAG) -L../lib -lSHTOOLS$(VERSION) -lfftw3 -lm -O3 -s $(F95FLAGS) MarsCrustalThickness/MarsCrustalThickness.f95 -o MarsCrustalThickness/MarsCrustalThickness
	
	$(F95) $(MODFLAG) -L../lib -lSHTOOLS$(VERSION) -lfftw3 -lm -O3 -s $(F95FLAGS) SHCilmPlus/TestCilmPlus.f95 -o SHCilmPlus/TestCilmPlus
	
	$(F95) $(MODFLAG) -L../lib -lSHTOOLS$(VERSION) -lfftw3 -lm -O3 -s $(F95FLAGS) SHExpandDH/TestExpandDH.f95 -o SHExpandDH/TestExpandDH
	
	$(F95) $(MODFLAG) -L../lib -lSHTOOLS$(VERSION) -llapack -lblas -lfftw3 -lm -O3 -s $(F95FLAGS) SHExpandLSQ/TestSHExpandLSQ.f95 -o SHExpandLSQ/TestSHExpandLSQ
	
	$(F95) $(MODFLAG) -L../lib -lSHTOOLS$(VERSION) -llapack -lblas -lfftw3 -lm -O3 -s $(F95FLAGS) SHLocalizedAdmitCorr/LocalizedAdmitCorr.f95 -o SHLocalizedAdmitCorr/LocalizedAdmitCorr
	
	$(F95) $(MODFLAG) -L../lib -lSHTOOLS$(VERSION) -lfftw3 -lm -O3 -s $(F95FLAGS) SHMag/SHMag.f95 -o SHMag/SHMag
	
	$(F95) $(MODFLAG) -L../lib -lSHTOOLS$(VERSION) -lfftw3 -lm -O3 -s $(F95FLAGS) SHRotate/TestSHRotate.f95 -o SHRotate/TestSHRotate
	
	$(F95) $(MODFLAG) -L../lib -lSHTOOLS$(VERSION) -lfftw3 -lm -O3 -s $(F95FLAGS) TimingAccuracy/TimingAccuracyDH.f95 -o TimingAccuracy/TimingAccuracyDH
	
	$(F95) $(MODFLAG) -L../lib -lSHTOOLS$(VERSION) -lfftw3 -lm -O3 -s $(F95FLAGS) TimingAccuracy/TimingAccuracyDHC.f95 -o TimingAccuracy/TimingAccuracyDHC
	
	$(F95) $(MODFLAG) -L../lib -lSHTOOLS$(VERSION) -lfftw3 -lm -O3 -s $(F95FLAGS) TimingAccuracy/TimingAccuracyGLQ.f95 -o TimingAccuracy/TimingAccuracyGLQ
	
	$(F95) $(MODFLAG) -L../lib -lSHTOOLS$(VERSION) -lfftw3 -lm -O3 -s $(F95FLAGS) TimingAccuracy/TimingAccuracyGLQC.f95 -o TimingAccuracy/TimingAccuracyGLQC
	
make run-examples:
	 @echo CREATE_WISDOM
	 cd fftw_wisdom ; create_wisdom < input.txt 
	 @echo
	 @echo MARS CRUSTAL THICKNESS
	 cd MarsCrustalThickness ; MarsCrustalThickness < input.txt
	 @echo
	 @echo TEST CILMPLUS
	 cd SHCilmPlus ; TestCilmPlus
	 @echo
	 @echo TEST EXPANDDH
	 cd SHExpandDH ; TestExpandDH
	 @echo
	 @echo TEST SHEXPANDLSQ
	 cd SHExpandLSQ ; TestSHExpandLSQ < input.txt
	 @echo
	 @echo LOCALIZEDADMITCORR
	 cd SHLocalizedAdmitCorr ; LocalizedAdmitCorr < input.txt
	 @echo
	 @echo SHMAG
	 cd SHMag ; SHMag < input.txt
	 @echo
	 @echo TEST SHROTATE
	 cd SHRotate ; TestSHRotate < input.txt
	 @echo
	 @echo *** The following timing and accuracy tests will take almost a day to complete ***
	 @echo
	 @echo TIMING ACCURACY GLQ
	 cd TimingAccuracy ;	TimingAccuracyGLQ < input1 > output1
	 @echo
	 @echo TIMING ACCURACY GLQ
	 cd TimingAccuracy ; TimingAccuracyGLQ < input2 > output2
	 @echo
	 @echo TIMING ACCURACY DH
	 cd TimingAccuracy ; TimingAccuracyDH < input3 > output3
	 @echo
	 @echo TIMING ACCURACY DH
	 cd TimingAccuracy ; TimingAccuracyDH < input4 > output4
	 @echo
	 @echo TIMING ACCURACY DHC
	 cd TimingAccuracy ; TimingAccuracyDHC < input5 > output5
	 @echo
	 @echo TIMING ACCURACY DHC
	 cd TimingAccuracy ; TimingAccuracyDHC < input6 > output6
	 @echo
	 @echo TIMING ACCURACY GLQC
	 cd TimingAccuracy ; TimingAccuracyGLQC < input7 > output7
	 @echo
	 @echo TIMING ACCURACY GLQC
	 cd TimingAccuracy ; TimingAccuracyGLQC < input8 > output8

getflags:

ifeq ($(F95),f95)
F95FLAGS ?= -m64 -O3 -YEXT_NAMES=LCS -YEXT_SFX=_
MODFLAG = -p ../modules/
endif

ifeq ($(F95),g95)
F95FLAGS ?= -m64 -O3 -fno-second-underscore 
MODFLAG = -I../modules/
endif

ifeq ($(F95),gfortran)
F95FLAGS ?= -m64 -O3
MODFLAG = -I../modules/
endif

ifeq ($(F95),ifort)
F95FLAGS ?= -m64 -free -O3 -Tf
MODFLAG = -I../modules/
endif

ifeq ($(origin F95FLAGS), undefined)
F95FLAGS = -m64 -O3
MODFLAG = -I../modules/
endif


clean:
	@rm -f SHCilmPlus/TestCilmPlus
	@rm -f fftw_wisdom/create_wisdom
	@rm -f MarsCrustalThickness/MarsCrustalThickness
	@rm -f SHExpandDH/TestExpandDH
	@rm -f SHExpandLSQ/TestSHExpandLSQ
	@rm -f SHLocalizedAdmitCorr/LocalizedAdmitCorr
	@rm -f SHMag/SHMag
	@rm -f SHRotate/TestSHRotate
	@rm -f TimingAccuracy/TimingAccuracyDH
	@rm -f TimingAccuracy/TimingAccuracyDHC
	@rm -f TimingAccuracy/TimingAccuracyGLQ
	@rm -f TimingAccuracy/TimingAccuracyGLQC
	


