#!/usr/bin/env python3

import sys
import unittest

errors=0

# we requires python3 as a minimum
if sys.version_info[0] < 3:
    print('The testsuite requires Python 3.\n')
    exit(-1)

# check for numpy.
try:
    import numpy
except ImportError:
    print('The testsuite requires the numpy module. Try \"pip3 install numpy\".\n')
    errors+=1

# check for scipy. This implies having numpy
try:
    import scipy
except ImportError:
    print('The testsuite requires the scipy module. Try \"pip3 install scipy\".\n')
    errors+=1

# check for pillow, which is imported as PIL
try:
    import PIL
except ImportError:
    print('The testsuite requires the pillow module. Try \"pip3 install pillow\".\n')
    errors+=1

try:
    import ssim 
except ImportError:
    print('The testsuite requires pyssim. Try \"pip3 install pyssim\".\n')
    errors+=1

if errors > 0:
    exit(-1)


# now search the subdirectories under the "tests" directory
suite = unittest.TestLoader().discover('tests')
unittest.TextTestRunner().run(suite)
