CLANG := $(shell command -v clang 2> /dev/null)

ifdef CLANG
CC=clang
CXX=clang++
else
ARCH := $(shell getconf LONG_BIT)
CPPFLAGS_32 := 
CPPFLAGS_64 := -march=x86-64 -mcmodel=medium

GCC_CFLAGS:=-Wno-stringop-truncation -Wno-format-truncation -Wno-format-overflow $(CPPFLAGS_$(ARCH))

CC=gcc $(GCC_CFLAGS)
CXX=g++ $(GCC_CFLAGS)
endif

CPPFLAGS= -O2 -Wall -fomit-frame-pointer -ffast-math -pipe
#CPPFLAGS= -O2 -Wall -fomit-frame-pointer -ffast-math -pipe -DHT_MODE

SRCS = splat.cpp
C_SRCS = itwom3.0.c

OBJS = $(SRCS:.cpp=.o) $(C_SRCS:.c=.o)

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

