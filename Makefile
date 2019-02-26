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

SRCS = splat.cpp
OBJS = $(SRCS:.cpp=.o)
C_SRCS = itwom3.0.c
C_OBJS = $(C_SRCS:.c=.o)

LDFLAGS = -lm -lpthread -lz -lbz2 -ljpeg -lpng

all: splat utils splat-hd

splat: $(SRCS) $(C_OBJS)
	$(CXX) $(CPPFLAGS) -std=c++11 $(INCLUDES) -c $(SRCS)
	$(CXX) $(CPPFLAGS) -o $@ $(OBJS) $(C_OBJS) $(LDFLAGS)

splat-hd: splat
	ln -sf splat splat-hd

utils:
	cd utils && $(MAKE)

install:
	./install all

.PHONY: clean utils install
clean:
	cd utils && $(MAKE) clean
	@rm -f *.o splat splat-hd

.SUFFIXES: .c .cpp .o
.c.o:
	$(CC) $(CPPFLAGS) -std=c99 -pedantic $(INCLUDES) -c $<


