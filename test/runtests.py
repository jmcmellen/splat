#!/usr/bin/env python3

import sys
import shutil
import unittest

err = False

# we requires python3 as a minimum
if sys.version_info[0] < 3:
    print('The testsuite requires Python 3.\n')
    exit(-1)
	
needed_modules=""

# check for numpy.
try:
    import numpy
except ImportError:
	needed_modules+="numpy "
	err = True

# check for scipy. This implies having numpy
try:
    import scipy
except ImportError:
	needed_modules+="scipy "
	err = True

# check for pillow, which is imported as PIL
try:
    import PIL
except ImportError:
	needed_modules+="pillow "
	err = True

try:
    import ssim
except ImportError:
	needed_modules+="pyssim "
	err = True

if len(needed_modules) > 0:
	print('The testsuite requires certain python modules. Try \"pip3 install %s\".\n' % needed_modules)
	
if shutil.which('gnuplot') is None:
	print('gnuplot not found. gnuplot is needed for some of the tests.\n')
	err = True
	
if err:
	exit(-1)


# now search the subdirectories under the "tests" directory
suite = unittest.TestLoader().discover('tests')
unittest.TextTestRunner().run(suite)
