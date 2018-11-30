ARCH := $(shell getconf LONG_BIT)
CLANG := $(shell command -v clang 2> /dev/null)
OS:=$(shell uname)

ifdef CLANG
  CC=clang
  CXX=clang++
  CLANG_CFLAGS:=
else
  CC=gcc
  CXX=g++
  CPPFLAGS_32:= 
  ifneq "$(OS)" "Darwin"
    CPPFLAGS_64:=-march=x86-64 -mcmodel=medium
  endif

  #GCC_CFLAGS:=-fomit-frame-pointer -Wno-stringop-truncation -Wno-format-truncation -Wno-format-overflow $(CPPFLAGS_$(ARCH))
  GCC_CFLAGS:=-fomit-frame-pointer $(CPPFLAGS_$(ARCH))
endif

#CPPFLAGS= -g -Wall -ffast-math -pipe $(CLANG_CFLAGS) $(GCC_CFLAGS)
CPPFLAGS= -O2 -Wall -ffast-math -pipe $(CLANG_CFLAGS) $(GCC_CFLAGS)
#CPPFLAGS= -O2 -Wall -ffast-math -pipe $(CLANG_CFLAGS) $(GCC_CFLAGS) -DHT_MODE

SRCS = splat.cpp
C_SRCS = itwom3.0.c

OBJS = $(SRCS:.cpp=.o) $(C_SRCS:.c=.o)

#LDFLAGS = -lm -lbz2 -lpthread -Wl,-stack_size -Wl,0x1000000
LDFLAGS = -lm -lbz2 -lpthread

all: splat

splat: $(OBJS)
	$(CXX) $(CPPFLAGS) -o $@ $^ $(LDFLAGS)

.PHONY: clean
clean:
	@rm -f *.o splat splat-hd splat-mt

.SUFFIXES: .c .cpp .o
.c.o:
	$(CC) $(CPPFLAGS) -std=c99 -pedantic $(INCLUDES) -c $<

.cpp.o:
	$(CXX) $(CPPFLAGS) -std=c++11 $(INCLUDES) -c $<

