CXX=g++
AR=ar
LD_SHARED=g++

ROOTLIBS := $(shell root-config --glibs --libs)  -lMinuit
ROOTINCLUDES := $(shell root-config --prefix)/include 
ROOTCFLAGS := $(shell root-config --cflags) -I../

INCDIR = -I. 
LIBDIR = ../lib

UNAME := $(shell uname)

ifneq ($(UNAME), Darwin)
EXPORT_DYNAMIC="-Wl,-export-dynamic"
endif

CXXFLAGS   = -Wall -O3 -g -fPIC $(ROOTCFLAGS) $(INCDIR) -Werror -std=c++11
SYSLIB     = -lm -lg2c
LINK_ARGS_BIN = $(SYSLIB) $(LIBDIR)
SOFLAGS= -shared $(EXPORT_DYNAMIC) $(ROOTCFLAGS)
LDFLAGS=
ARFLAGS= rcsv

ifdef COMBINED_KIN
CXXFLAGS += -DCOMBINED_KIN
endif

ifdef MULTITHREAD
CXXFLAGS += -fopenmp
OMP_DEFINES = -DMULTITHREAD
endif

OBJECTS := covarianceBase.o covarianceOsc.o covarianceXsec.o

%.o:%.cpp
	g++ -c $(CXXFLAGS) $< $(OMP_DEFINES)

libCovariance: libCovariance.a $(LIBDIR)/libCovariance.so

libCovariance.a: $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $^

$(LIBDIR)/libCovariance.so: libCovariance.a $(OBJECTS)
	$(LD_SHARED) $(SOFLAGS) $^ -o $@ $(LIB_OBJECTS) $(ROOTLIBS) -L../lib -lThrowParms

clean:
	rm -f ../lib/libCovariance.so	
	rm -f *.o lib*

distclean:
	rm -fv *~ *.o *.so *.a lib*
