SHTOOLS COMPILATION INSTRUCTIONS

In order to compile the SHTOOLS archive into a static library, in most
cases, it is necessary to type only the command

	make

in the directory SHTOOLS, or

	make F95="MyCompiler"

where "MyCompiler" is replaced by the compiler name. If the compiler
name is not specified, the default compiler will be assumed to be "f95".
Two free Fortran 90/95 compilers are gfortran (http://gcc.gnu.org/) and
g95 (http://www.g95.org/). To modify the default compiler flags set in
the makefile, use

	make F95="MyCompiler" F95FLAGS="MyCompilerFlags"

Successful compilation will create the library file libSHTOOLS2.x.a in
the directory lib and will place a few compiled module files in the
directory modules. Following compilation, the entire directory can be
moved to /usr/local (if desired) by typing

	make install

(this will require administrator privileges). To compile and then run
the example files, use

	make examples F95="MyCompiler" F95FLAGS="MyCompilerFlags"

and

	make run-examples

Note that the timing tests could take almost a day to complete. To
delete the compiled archive, module files, object files, and example
executables, use

	make clean

and

	make remove-examples

If you need to recompile SHTOOLS a second time using a different set of
compiler flags, it will be necessary to first remove all the previously
compiled object files by typing "make clean".

The man and html web pages are included pre-built in the SHTOOLS
directory (if it is necessary to rebuild these, use "make doc").
However, in order to access the unix man pages, it will be necessary to
add "shtoolslocation"/man to your man path, where "shtoolslocation" is
the path of the root directory of SHTOOLS. The link to the local html
web page is "shtoolslocation"/SHTOOLS.html.


COMPILER FLAGS

Default compiler options are specified in the Makefile for a few common
compilers (Absoft f95, gfortran, g95, and ifort). If it is necessary to
change these, consider the following guidelines:

One should always use some form of optimization when compiling SHTOOLS,
such as by specifying the option

	-O3

The biggest difficulty in compiling SHTOOLS is setting the compiler
flags so that the external subroutine names are in a format that is
compatible with the FFTW and LAPACK libraries. In general, it is
necessary to ensure that the SHTOOLS subroutine names are in lower case
and have the right number of underscores appended to them.

For Absoft ProFortran, this is achieved by setting

	-YEXT_NAMES=LCS -YEXT_SFX=_

For g95, it will be necessary to use either

	-fno-second-underscore</tt> (most likely)

or

	-fno-underscoring

For gfortran, it is generally not necessary to use any special flags,
though it could arise that either

	-fno-underscoring

or

	-fsecond-underscore

might be necessary.

For the Intel Fortran compiler ifort, it will be necessary to use

	-free -Tf

in order that the compiler recognizes files with the extension .f95 as
fortran 95 files. In this case, the f95 file should come directly after
the option -Tf.

Setting the right compiler flags is more complicated when the FFTW and
LAPACK libraries have different naming and underscoring conventions. In
order to accommodate this case, underscores have been explicitly added
to LAPACK and FFTW subroutine names in an alternative set of source
files. In order to compile SHTOOLS with underscores explicitly appended
to LAPACK routine names, use

	make all2

In order to compile SHTOOLS with underscores explicitly appended to FFTW
routine names, use

	make all3

For both cases, compiler flags should be set so that underscores are not
appended to routine names. See the FAQ for further information.

To generate 64 bit code, use the compiler option

	-m64

For this case, it will be necessary to use 64-bit compiled FFTW and
LAPACK libraries.
