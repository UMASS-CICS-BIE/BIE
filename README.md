========================================
=      Bayesian Inference Engine       =
========================================

Files and directories:

README		This file

TODO		Wish list of features to add, things to fix, big and small

doc		DOC++ documentation templates (cd doc; make html). Directory.

include		Include file for all classes.  Directory.

libVector	Utility library 1 (vector class, Gaussian integration
		developed by MDW).  Directory.

libsrnd		Utility library 2 (random number classes, GNU string class,
		mostly part of the old libg++, no longer maintained by gcc
		folks).  Directory.

libcutils	(old directory)

src		Source for simulation classes.  Directory of directories.

testing		Test routines. Directory.
		(e.g. cd testing; make multilevelN)

cli		Command line interface.  Directory.

Galphat         GAlaxy PHotometric ATtributes: probabilistically rigorous 
                galaxy image analysis

PopModel	Realistic galaxy stellar population model (currently under 
                development)

==================================================
The main routines are in the directory testing.
To compile everything you will need to:

First: make sure that you have a MPIHOME environment
variable set.  E.g. in csh

	set MPIHOME /usr/local/openmpi

or in sh

	MPIHOME=/usr/local/openmpi

Then,

	aclocal -I m4
	autoheader
	libtoolize --force
	automake --add-missing
	autoconf
	./configure

To make the entire source tree type

	make

And if desired

        make install

The build may be run in parallel (e.g. make -j4 <other options>) but
the installiation should be done serially.  Alternatively, to compute
only one of the test routines the Gaussian splatting model for
example, you could:

	cd testing
	make multilevelS

Or,

	./autoconf.sh

to do the build in one step.  You may follow this with a "make
install" if desired.


You will need the CommonC++ libraries in order to compile the
command line interface (which uses the wrappers to POSIX threads).
If these are not installed in the standard system locations, you
can specify the location of this library with the

	--with-commoncpp=<top level ComonC++ directory>

flag to either "autogen.sh" or "configure".

You will also need the boost library version >=1.35 (system,
serialization, mpi and file_system). You can specify the location of
this library with the

    --with-boost=<top level Boost library>

or

    --with-boost-libdir=<.so file location>


In order to compile the documentation in the doc directory, you will
need the doxygen utility, available from

	http://www.doxygen.org

or most of Linux distribution repositories.

==================================================
The default test input files are now for
the Gaussian splat model (multilevelS).

Good choices for parameters for first attempts might be:

./multilevelS -T 1 -n 1000 -t 100000 -N 10 -I 7 -f`pwd`

where

-T 1		use kd tessellation with bounds
-n 1000		take up to 1000 steps per level
-t 100000	maximum temperature increase for tempered annealing
-I 7		use 7 integration knots in each dimension per tile
-f <datadir>    defines the directory containing the data file (e.g. `pwd`)

for -T 2 (MappedGrid)
-N 6		use 6x6 tiling

for -T 3 (QuadGrid)
-N 6		6-level tree

NB: not a good idea to use -T 0 (kd tessellation without bounds) 
    because perimeter bins will be unbounded and cause numerical
    integrals to be grossly inaccurate

==================================================

A good first attempt for simple galaxy model:

./multilevelN -T 3 -n 1000 -t 100000 -N 2 -I 7 -f`pwd`

where

-T 3		use quad grid
-n 1000		take up to 1000 steps per level
-t 100000	maximum temperature increase for tempered annealing
-I 7		use 7 integration knots in each dimension per tile
-N 2            use 2-level tree to start
-f <datadir>    defines the directory containing the data file (e.g. `pwd`)

for -T 2 (MappedGrid)
-N 6		use 6x6 tiling

for -T 3 (QuadGrid)
-N 6		6-level tree

NB: not a good idea to use -T 0 (kd tessellation without bounds) 
    because perimeter bins will be unbounded and cause numerical
    integrals to be grossly inaccurate


==================================================

Notes on the command line interface (CLI):

The CLI should have all the functionality of the routines in
testing but allow the user set parameters, instantiate classes
and run a simulation interactively.  Plans are to multithread the
control and running simulation but this had not been implemented.
See cli/simuclass_script for an example script.  The cli also as
interactive help.
