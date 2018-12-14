# SPLAT!

A Terrestrial RF Path and Terrain Analysis Tool for Unix/Linux

## About

This version is a work in progress. It does exactly the same calculations as the basic
Splat! system, but uses CPU threading to speed things up.

Future version will either use OpenCL or Vulkan to hand computation off to a graphics
card in the hopes of even more speed improvements. In preparation for this, itwom3.0 was
made fully C99-compliant, as all the current implementations of OpenCL drivers require
that. (Later versions of OpenCL allow C++, but none of the common GPU drivers support that).

## Getting Started

Build instructions are in the file README.

For this version, you must have either clang or gcc installed, and it must be a version that supports at
least C++11 .

You also need libzip and libbzip. You can often install these from packages. For instance, on Debian:

apt-get install zlib1g-dev libbz2-dev


## Changes

* Build system

  * The build system has partially been converted to Gnu Make. More work needs to be done here to support
    building the subdirectories.

* ITWOM 3.0

  * itwom3.0.cpp was renamed to itwom3.0.c and is now compiled with the C compiler rather than the C++ one.
    In addition, every effort was made to make it fully C99-compliant rather than the peculiar amalgam of C
    and C++ syntax that it was using.
  
    itwom3.0 no longer uses the C++ complex number templates (obtained by doing #include <complex>). Instead
    a small complex-number library is introduced. This was done as part of the effort to make it fully
    C99 compliant.

  * All the static variables were removed from itwom3.0.c in an effort to make the code fully reentrant. In
    some cases this means that we recalculate variables multiple times. In other cases, new "state" structs
    have been introduced that act as local contexts for repeated calls to the same function. While this does
    slow things down a bit over using statics, the difference is most noticeable on older systems (on an Atom
    D510 it increases run time from 4 1/2 minutes to 5 minutes). On newer systems it seems to make little to
    no difference. At any rate, making that code reentrant means that the code can be multithreaded, and the
    gains from that far offset any loss.

  * A number of minor fixes were made to itwom3.0 that were disguised by using the C++ compiler. For instance,
    in a number of places abs() was called when fabs() was meant.

* splat.cpp

  * The PlotLOSMap() and PlotLRMap() functions have been converted to run multithreaded if a "-mt" flag is
    passed on the command line.

## To Do

* Make the memory allocation for the arrays in splat.cpp be dynamic so there's no need for two different
  versions of splat.

* Add the OpenCL parallelizing.

