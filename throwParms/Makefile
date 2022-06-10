CXX=g++
AR=ar
LD_SHARED=g++

ROOTLIBS   = -Wl,-rpath,$(shell root-config --libdir) $(shell root-config --libs) -lThread 
ROOTCFLAGS = $(shell root-config --cflags)
INCDIR = -I.
LIBDIR = ../lib

UNAME := $(shell uname)

ifneq ($(UNAME), Darwin)
EXPORT_DYNAMIC="-Wl,-export-dynamic"
endif

#CLHEPINC = $(CLHEP_BASE_DIR)/include
#CLHEPLIB = -L$(CLHEP_BASE_DIR)/lib -l$(CLHEP_LIB)



CXXFLAGS   = -Wall -O3 -g -fPIC $(ROOTCFLAGS) $(INCDIR) -Werror -std=c++11
SYSLIB     = -lm -lg2c
LINK_ARGS_BIN = $(SYSLIB) $(ROOTLIBS) $(LIBDIR)
SOFLAGS= -shared $(EXPORT_DYNAMIC) $(ROOTCFLAGS) $(ROOTLIBS)
LDFLAGS=
ARFLAGS= rcsv


%.o:%.cpp
	g++ -c $(CXXFLAGS) $<

libThrowParms: libThrowParms.a $(LIBDIR)/libThrowParms.so

libThrowParms.a: ThrowParms.o
	$(AR) $(ARFLAGS) $@ $^

$(LIBDIR)/libThrowParms.so: libThrowParms.a ThrowParms.o
	$(LD_SHARED) $(SOFLAGS) $^ -o $@ $(LIB_OBJECTS)

clean:
	rm -f *.o lib*

distclean:
	rm -fv *~ *.o *.so *.a lib*
