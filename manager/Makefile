CXX=g++
AR=ar
LD_SHARED=g++

ROOTLIBS := $(shell root-config --glibs --libs) 
ROOTINCLUDES := $(shell root-config --prefix)/include 
ROOTCFLAGS := $(shell root-config --cflags) -I../

INCDIR = -I. 
LIBDIR = ../lib

UNAME := $(shell uname)

ifneq ($(UNAME), Darwin)
EXPORT_DYNAMIC="-Wl,-export-dynamic"
endif

ifdef MULTITHREAD
OMP_DEFINES=-DMULTITHREAD
CXXFLAGS=-fopenmp
endif

ifdef CUDAPATH
GPU_DEFINES=-DCUDA
endif

CXXFLAGS   = -Wall -fPIC $(ROOTCFLAGS) $(INCDIR) $(OMP_DEFINES) $(GPU_DEFINES) -Werror -std=c++11
ifndef DEBUG
CXXFLAGS+=-O3 -g
endif

SYSLIB     = -lm -lg2c
LINK_ARGS_BIN = $(SYSLIB) $(ROOTLIBS)
SOFLAGS= -shared $(EXPORT_DYNAMIC) $(ROOTCFLAGS)
LDFLAGS=
ARFLAGS= rcsv

%.o:%.cpp
	g++ -c $(LINK_ARGS_BIN) $(CXXFLAGS) $<

libManager: libManager.a $(LIBDIR)/libManager.so

libManager.a: manager.o
	$(AR) $(ARFLAGS) $@ $^

$(LIBDIR)/libManager.so: libManager.a manager.o
	$(LD_SHARED) $(SOFLAGS) $^ -o $@ $(LIB_OBJECTS) $(ROOTLIBS) -L../lib -lconfig++

clean:
	rm -f ../lib/libManager.so	
	rm -f *.o lib*

distclean:
	rm -fv *~ *.o *.so *.a lib*
