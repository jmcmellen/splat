CC=g++

ARCH := $(shell getconf LONG_BIT)
CPPFLAGS_32 := 
CPPFLAGS_64 := -march=x86-64 -mcmodel=medium

CPPFLAGS=-O2 -s -fomit-frame-pointer -ffast-math -pipe $(CPPFLAGS_$(ARCH))

SRCS = itwom3.0.cpp \
	  splat.cpp
OBJS = $(SRCS:.cpp=.o)

LDFLAGS = -lm -lbz2

all: splat splat-hd

splat: $(SRCS)
	$(CC) $(CPPFLAGS) -o $@ $^ $(LDFLAGS)

splat-hd: $(SRCS)
	$(CC) $(CPPFLAGS) -DHD_MODE=1 -o $@ $^ $(LDFLAGS)

.PHONY: clean
clean:
	@rm -f splat splat-hd

