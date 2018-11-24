CC=gcc -O2 -Wall
CXX=g++ -O2 -Wall

ARCH := $(shell getconf LONG_BIT)
CPPFLAGS_32 := 
CPPFLAGS_64 := -march=x86-64 -mcmodel=medium

CPPFLAGS= -fomit-frame-pointer -ffast-math -pipe $(CPPFLAGS_$(ARCH))
#CPPFLAGS= -fomit-frame-pointer -ffast-math -pipe $(CPPFLAGS_$(ARCH)) -DUSE_MT_WORKQUEUE
#CPPFLAGS= -fomit-frame-pointer -ffast-math -pipe $(CPPFLAGS_$(ARCH)) -DUSE_MT_WORKQUEUE -DHT_MODE

SRCS = splat.cpp
C_SRCS = itwom3.0.c

OBJS = $(SRCS:.cpp=.o) $(C_SRCS:.c=.o)

LDFLAGS = -lm -lbz2
#LDFLAGS = -lm -lbz2 -lpthread

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

