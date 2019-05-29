# SPLAT! for Windows (Command line)

### NOTE: SDF files under Windows *must* be named so that they have '_'s instead of ':'s in them.
  For instance, a file that was formerly "44:45:120:121.sdf" must now be "44_45_120_121.sdf".
   

## Visual Studio version
   SPLAT! requires the C++11 standard, which means Visual Studio 2015 or better.
   The Community Edition will work fine.

## Compiling SPLAT! for Windows, the TL;DR version:
  Checkout files from git. Open the sln. Choose "Build->Solution". Run Splat.exe with
  the same options you'd use with unix.

## Details
* This creates a Win32 Console Application. No fancy windows here!

* Settings
    * Unicode is not supported, nor multibyte character set. General->Character set should be left as "unset".
    * The preprocessor for both Debug and Release needs to have _CRT_SECURE_NO_WARNINGS
    * The optimization should be /Od for Debug and /O2 for Release.
    * Precompiled Headers should be set to "Not Using Precompiled Headers".
    
* Auxiliary Libraries
    You need bzip2, zlib, and png libraries at the minimum. jpeg is optional.

    The project files includes the necessary libraries under <Splat>/windows/x64_lib, but only for the
    x64 platform. If you want to compile the libraries yourself you can use the project file as
    a guideline. Here're some tips on where to get them and how to build them:

    * libbzip2: You can get the original code from https://sourceforge.net/projects/bzip2. At the time
        of writing, the release is 1.0.6. To build...well, it's sort of ugly. You'll need to create your
        own sln and vcxproj files and include blocksort.c, bzlib.c, compress.c, crctable.c, decompress.c,
        huffman.c, and randtable.c. Build them into a static library, no unicode. Make sure you use /MT
        and /MTd. There is a minimal version that does this at https://github.com/hoche/bzip2. Once built,
        copy the following:
        - bzip2-1.0.6/bzip.h to splat/bzip.h
        - bzip2-1.0.6/x64/Release/bzip.lib to splat/bzip.lib

            
    * zlib: You don't actually need to build this, but you can if you like. libpng will build it for you.
        Either way you will need to have the source files. Get them from http://www.zlib.net/. At the time
        of writing, the release is 1.2.11. Once you've unpacked it, copy the following:
        - zlib-1.2.11\zlib.h to splat\zlib.h
        - zlib-1.2.11\zconf.h to splat\zconf.h
        
        You don't need to, but if you want to build it independently of libpng, go into the contrib/vstudio directory
        and choose the subdirectory corresponding to your rev of Visual Studio. For instance VS2015 is vc14. Make sure
        the /MT setting is used (not /MD). You only need to build the static library (zlibstat). I had to use the
        "ReleaseWithoutAsm" configuration. When it's done, copy the following:
        - zlib-1.2.11\contrib\vstudio\vc14\x86\ZlibStatReleaseWithoutAsm\zlibstat.lib to splat\zlib.lib

            
    * libpng: Get http://www.libpng.org/pub/png/libpng.html. You'll need a version later than 1.5.x. At the
        time of writing, the release version is 1.6.36. To build, go into projects/vstudio. Edit zlib.props
        so that the <ZLibSrcDir> value points to the zlib you just downloaded - for instance:
            `<ZLibSrcDir>..\..\..\..\zlib-1.2.11</ZLibSrcDir>`
  
        Then open the .sln file. It is for VS2010; if you open it with a later version like VS2015 it will ask
        if you want to convert the project. Go ahead and allow it.
        
        * By default, the project will compile for 32-bit Windows. You probably don't want this; it must match
        the platform Splat is compiled with, which is x64. To set this:
            - In the Visual Studio menubar, click the "Solution Platforms" pulldown. If you don't have an "x64"
              option, click "Configuration Manager".
            - Under "Active Solution Platform", choose "<New>"
            - Choose "x64" as the platform type, and use the "Copy settings from Win32" option. Make sure the
              "Create new platform projects" checkbox is selected. Click "Ok".
        Set your "Solution Platform" to the newly-created "x64" library.
        
        * Next, change the configuration to the "Release Library", which builds the library for statically
        linking (not the dll). At minimum, build the following two projects:
          - libpng
          - pnglibconf
          - zlib (which should happen automatically when you build libpng)
        
        Note: You may get an error "unreferenced formal parameter" when building libpng. If so, try going to
        Configuration properties->C/C++->Advanced and adding 4100 to the "Disable Specific Warnings" line:
            `$(DisableSpecificWarnings);4100`
            
        * When done, copy the following:
          - lpng1636\png.h to splat\png.h
          - lpngl635\pngconf.h to splat\pngconf.h
          - lpng1636\pnglibconf.h to splat\png.h         (if you don't have this you didn't build pnglibconf properly)
          - lpng1636\projects\vstudio\Release Library\libpng16.lib to splat\libpng.lib
          - lpng1636\projects\vstudio\Release Library\zlib.lib to splat\zlib.lib
          
* Auxiliary Programs
    You need gnuplot to get the terrain and elevation profile plots. It can be installed either by downloading from
    the gnuplot site:
        `https://sourceforge.net/projects/gnuplot/`
        
    or by using chocolatey (the Windows package manager available from https://chocolatey.org/):
        `choco install gnuplot`
        
## Troubleshooting
    * Splat can't find png_write_row(), png_set_compression_level, etc:
        Make sure that when you compiled libpng, the Solution Platform matches Splat's (x64).
  
  
Last updated on 14 Jan 2019. hoche@grok.com
     
