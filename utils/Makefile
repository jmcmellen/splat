ARCH := $(shell getconf LONG_BIT)
CLANG := $(shell command -v clang 2> /dev/null)
GXX := $(shell command -v g++ 2> /dev/null)
OS:=$(shell uname)

# prefer gcc/g++, if available
# there's really no good reason for using this over clang and this test should probably
# be reversed.
ifdef GXX
  CC=gcc
  CXX=g++
  CPPFLAGS_32:= 
  GCC_CFLAGS:=$(CPPFLAGS_$(ARCH))
else
  CC=clang
  CXX=clang++
  CLANG_CFLAGS:=
endif

#CPPFLAGS= -g -Wall -ffast-math $(CLANG_CFLAGS) $(GCC_CFLAGS)
CPPFLAGS= -O3 -Wall -ffast-math $(CLANG_CFLAGS) $(GCC_CFLAGS)

LDFLAGS = -lm -lbz2 -lz

all: citydecoder usgs2sdf srtm2sdf srtm2sdf-hd fontdata bearing

citydecoder: citydecoder.o
	$(CXX) $(CPPFLAGS) -o $@ $< $(LDFLAGS)

usgs2sdf: usgs2sdf.o
	$(CXX) $(CPPFLAGS) -o $@ $< $(LDFLAGS)

srtm2sdf: srtm2sdf.o
	$(CXX) $(CPPFLAGS) -o $@ $< $(LDFLAGS)

srtm2sdf-hd: srtm2sdf
	ln -sf srtm2sdf srtm2sdf-hd

fontdata: fontdata.o
	$(CXX) $(CPPFLAGS) -o $@ $< $(LDFLAGS)

bearing: bearing.o
	$(CXX) $(CPPFLAGS) -o $@ $< $(LDFLAGS)

.PHONY: clean
clean:
	@rm -f *.o citydecoder usgs2sdf srtm2sdf fontdata bearing

.SUFFIXES: .c .cpp .o
.c.o:
	$(CC) $(CPPFLAGS) -std=c99 -pedantic $(INCLUDES) -c $<

.cpp.o:
	$(CXX) $(CPPFLAGS) -std=c++11 $(INCLUDES) -c $<

