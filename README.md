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

You also need several utility libraries:
* libbzip2
* zlib
* libpng
* libjpeg

You can generally get these via system packages. For instance:
    
Centos 7:

`yum install bzip2-devel zlib-devel libpng-devel libjpeg-turbo-devel`

Debian (Buster):

`apt-get install libbz2-dev zlib1g-dev libjpeg-dev libpng-dev`

OSX (High Sierra):

`brew install jpeg libpng`

## Changes

* Build system

  The build system has been converted to Gnu Make (gmake). There is no need to run "configure" any more.
  The next release may see a CMake build system.
  
* SDF file names

  In order to accommodate Windows filesystems, the basic format of the SDF file naming has changed.
  The colon has been changed to an underline. For instance, instead of having a file named `46:47:122:123.sdf`
  it should now be `46_47_122_123.sdf`. The various unix versions (including OSX) will handle either version,
  but by default srtm2sdf will now create names with the underline. If you don't want this, a "-c" flag has
  been added to srtm2sdf to preserve your colon.

* splat.cpp

  * Incorporate John's antenna height changes from SPLAT 1.4.3 (unreleased).
  
  * Revert to using the ITM model by default, in accordance with SPLAT 1.4.3.
  
  * The PlotLOSMap() and PlotLRMap() functions have been converted to run multithreaded ~~if a "-mt" flag is
    ~~passed on the command line~~. If you want to run single-threaded, use "-st" on the command line.

  * WritePPM(), WritePPMSS(), etc were converted to WriteImage(), WriteImageSS(), etc, and functionality
    was added to allow them to emit png or jpg images instead of pixmaps. png's are now the default. Add "-ppm"
    or "-jpg" to the command line if you want to generate the others. The generated jpg's are smaller but the text
    can be hard to read in some instances. The png's are, of course, lossless, and nice and crisp, but they are
    larger and take slightly longer to generate. Both are an order of magnitude smaller than the pixmaps though.
    
  * All the DEMs and Paths are allocated on the heap rather than using pre-sized arrays and stack allocations.
    This means the same code can be used for 1x1 through 8x8 grids, and/or 3-deg (Normal) or 1-deg (High Definition)
    modes without recompiling. It also means that the code can adjust itself (or error out nicely) depending on the
    available memory of the host machine.
    
  * I did some minor fixup to the antenna azimuth (.az) and elevation (.el) reading so that you can now put a value
    of 0 in those files to indicate your antenna blocks that area perfectly.
    
* ITWOM 3.0

  * itwom3.0.cpp was renamed to itwom3.0.c and is now compiled with the C compiler rather than the C++ one.
    In addition, every effort was made to make it fully C99-compliant rather than the peculiar amalgam of C
    and C++ syntax that it was using.
  
    itwom3.0 no longer uses the C++ complex number templates (obtained by doing #include <complex>). Instead
    a small complex-number library is introduced. This was done as part of the effort to make it fully
    C99 compliant.

  * All the static variables were removed from function scopes within itwom3.0.c in an effort to make the code
    fully reentrant. This had resulted in a number of changes:
        - There are a few cases where these were actually consts; those have been defined as such
        and moved to the global scope.
        - In a few cases this means that we recalculate variables multiple times.
        - In other cases, new "state" structs have been introduced that act as local contexts for repeated
        calls to the same function.
    All these changes are geared towards making the code multithread-safe.
    
  * A few variables that are not modified by the various functions have been declared as const. This includes
    pfl[] arrays.

  * A number of minor fixes were made to itwom3.0 that were disguised by using the C++ compiler. For instance,
    in a number of places abs() was called when fabs() was meant.
    
  * Much much code documenting was done. There remains a lot to do though.

## To Do

* Since we have to link to zlib for the pngs, we might as well create kmz files if asked.
* More code cleanup.
* Split this into multiple files and c++ify things.
* Reformat tabs to four spaces.
