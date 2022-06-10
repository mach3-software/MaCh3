SRC := $(wildcard *.cpp)
TRG := $(SRC:.cpp= )
INC := $(wildcard *.h)

$(info Making: $(SRC) to $(TRG) with $(INC))

# Compiler
CC=g++ 
CXXFLAGS = -g -O3 -Wall -Werror -Wextra
ifdef MULTITHREAD
CXXFLAGS += -fopenmp
CXXFLAGS += -DMULTITHREAD
endif

# ROOT CFLAGS
ROOT_CFLAGS=$(shell root-config --cflags)
ROOT_LIBS=$(shell root-config --libs) -lEG -lGeom -lXMLIO
ROOT_INCLUDES=-I$(shell root-config --incdir)

# MaCh3 includes and libs
LIBDIR=../lib

LIBS=-L$(LIBDIR) ${ROOT_LIBS}
INCLUDES=$(ROOT_INCLUDES) $(INC)


all:	$(TRG)

%:	%.cpp
	$(CC) $(CXXFLAGS) $(ROOT_CFLAGS) $(INCLUDES) -o $@ $< $(LIBS)

clean:
	rm -rfv $(TRG)
