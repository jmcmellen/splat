#!/usr/bin/env python

# compare_ano.py <file1> <file2>
# 
# compares two .ano files emitted from splat (via the -ano flag). The
# contents don't have to be in the same order, but they should be from
# runs with:
#  a) the same transmitter site
#  b) the same receiver site
#  c) the same settings for the splat compilation (8x8, 2x2, etc)
#  d) the same terrain elevation files
#
# This is primarily to test for differences between linear processing,
# multithreaded-CPU processing, and GPU processing.
#
# 
# Bounds:
#
# In regular Splat, ipp=1200 and that gives a bounds of
#   (max_west - min_west) / 1200, or (for typical use) 2/1200 = 0.00167 decimal degrees
# as the maximum positional accuracy for latitude. Longitude is similar, but of course
# it compresses as you get farther from the equator.
#
# In HD Splat, ipp=3600, and the accuracy goes up to 0.00055.
#
# For Reference: 
# decimal degrees           meters (at equator)
#     0.0001                     11
#     0.00001                     1.1      
#     0.000001                     .11     
#
# https://en.wikipedia.org/wiki/Decimal_degrees
#

import sys
import math
from collections import defaultdict

ano1 = defaultdict(dict)

ppd = 1200
mpi = (1200 - 1)

max_west = -1
min_west = 360
max_north = 90
min_north = -90

def lonDiff(lon1, lon2):
    """cribbed from splat.cpp"""
    diff = lon1-lon2
    if diff <= -180:
        diff += 360

    if diff >= 180:
        diff -= 360

    return diff

def getElemCoords(lat, lon):
    x = round(ppd * (lat - min_north))
    y = mpi - round(ppd*(lonDiff(max_west, lon)))
    return (x,y)

def loadAnoIntoDict(anofilepath, anodict):
    """Read an anofile and load it into a dictionary of dictionaries, indexed by the same 
       coordinates that dem uses in splat.cpp. Note that this only works for regular Splat.
       Splat HD requires a different division factor.
    """
    with open(anofilepath) as fp:  
        line = fp.readline()
        cnt = 1

        while line:
            data = line.strip().split()
            if ";" in line:
                if (cnt == 1):
                    max_west = float(data[0].rstrip(","))
                    min_west = float(data[1])
                elif (cnt == 2):
                    max_north = float(data[0].rstrip(","))
                    min_north = float(data[1])
            else:
                lat = float(data[0].rstrip(","))
                lng = float(data[1].rstrip(","))
                azi = float(data[2].rstrip(","))
                elev = float(data[3].rstrip(","))
                loss = float(data[4])
                if len(data) == 6:
                    blocked = True
                else:
                    blocked = False

                #print("%0.7f,%0.7f,%0.7f,%0.7f,%0.7f,%s" % (lat, lng, azi, elev, loss, blocked))
                (xcoord, ycoord) = getElemCoords(lat, lng)
                anodict[xcoord][ycoord] = (lat, lng, xcoord, ycoord, azi, elev, loss, blocked) 

            line = fp.readline()
            cnt += 1

            if (cnt % 20000) == 0:
                sys.stdout.write('.')
                sys.stdout.flush()

        print("")


def findInAnoDict(findlat, findlng, anodict):
    """Search in the dictionary of lists of locations for a location within maxdist range"""
    """Returns none or the location data"""

    (x, y) = getElemCoords(findlat, findlng)
    entries = anodict.get(x)
    if entries == None:
        return None

    return entries.get(y)


def compareAnos(anofilepath, anodict):
    """Read an anofile and compare it to the data in anodict"""
    with open(anofilepath) as fp:  
        line = fp.readline()
        cnt = 1
        needHeader = True

        miss = 0

        while line:
            if ";" not in line:
                data = line.strip().split()
                lat = float(data[0].rstrip(","))
                lng = float(data[1].rstrip(","))
                azi = float(data[2].rstrip(","))
                elev = float(data[3].rstrip(","))
                loss = float(data[4])
                if len(data) == 6:
                    blocked = True
                else:
                    blocked = False

                data = findInAnoDict(lat, lng, ano1)
                (x, y) = getElemCoords(lat, lng)
                if (data == None):
                    print("not found")
                else:
                    if (abs(data[6] - loss) > 0.1):
                        pct = data[6]/loss
                        if (pct > 5):
                            if needHeader:
                                print(" latitude , longitude  ( demX , demY ): dbloss      latitude , longitude  ( demX , demY ): dbloss  %change")
                                needHeader = False
                            print("%0.7f,%0.7f (%d,%d): %0.4f vs %0.7f,%0.7f (%d,%d): %0.4f  %0.1f%%" % 
                                    ( data[0], data[1], data[2], data[3], data[6],
                                      lat, lng, x, y, loss, pct))
                            miss += 1

            line = fp.readline()
            cnt += 1

        print("misses: %d" % (miss))


def main():
    if len(sys.argv) < 3:
        print("compare_ano.py <file1> <file2>\n")
        return

    print("loading %s" % (sys.argv[1]) )
    loadAnoIntoDict(sys.argv[1], ano1)
    compareAnos(sys.argv[2], ano1)

if __name__ == "__main__":
    main()
