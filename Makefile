CC=gcc
CXX=g++

ARCH := $(shell getconf LONG_BIT)
CPPFLAGS_32 := 
CPPFLAGS_64 := -march=x86-64 -mcmodel=medium

#CPPFLAGS=-g -std=c++11 -fomit-frame-pointer -ffast-math -pipe $(CPPFLAGS_$(ARCH)) -DUSE_MT_WORKQUEUE
CPPFLAGS=-O2 -fomit-frame-pointer -ffast-math -pipe $(CPPFLAGS_$(ARCH))

SRCS = splat.cpp
C_SRCS = itwom3.0.c

OBJS = $(SRCS:.cpp=.o) $(C_SRCS:.c=.o)

LDFLAGS = -lm -lbz2

all: splat splat-hd splat-mt

splat: $(OBJS)
	$(CXX) $(CPPFLAGS) -o $@ $^ $(LDFLAGS)

splat-hd: $(OBJS)
	$(CXX) $(CPPFLAGS) -DHD_MODE -o $@ $^ $(LDFLAGS)

splat-mt: $(OBJS)
	$(CXX) $(CPPFLAGS) -DUSE_MT_WORKQUEUE -o $@ $^ $(LDFLAGS) -lpthread

.PHONY: clean
clean:
	@rm -f *.o splat splat-hd splat-mt

.SUFFIXES: .c .cpp .o
.c.o:
	$(CXX) $(CPPFLAGS) -std=c99 -pedantic $(INCLUDES) -c $<

.cpp.o:
	$(CXX) $(CPPFLAGS) -std=c++11 $(INCLUDES) -c $<

