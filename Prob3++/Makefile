ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)

CXXFLAGS = -I. -Wall -fPIC -O3 -g -std=c++11

%.o : %.c
	$(RM) $@
	$(CC) -c -fPIC $(CFLAGS) -O2 -o $@ $<
	
%.o : %.cc
	$(RM) $@
	$(CXX) -c $(CXXFLAGS) -o $@ $*.cc


OBJS    = EarthDensity.o BargerPropagator.o mosc.o mosc3.o 


LIBBASE   = ThreeProb
VER       = 1.00
TAG       = 
LIBALIAS  = $(LIBBASE)$(TAG)
LIBNAME   = $(LIBALIAS)_$(VER)

lib3p     = lib$(LIBNAME).a


targets = $(lib3p) probRoot probLinear


$(lib3p) : $(OBJS)
	$(RM) $@
	ar clq $@ $(OBJS)
	ranlib $@
	

probRoot: probRoot.o $(lib3p) 
	$(RM) $@
	$(CXX) -o $@ $(CXXFLAGS) -L. $^ $(ROOTLIBS)


.PHONY: probRoot.o
probRoot.o: 
	$(CXX) -o probRoot.o $(ROOTCFLAGS) $(CXXFLAGS) -c probRoot.cc



probLinear: probLinear.o $(lib3p) 
	$(RM) $@
	$(CXX) -o $@ $(CXXFLAGS) -L. $^ $(ROOTLIBS)


.PHONY: probLinear.o
probLinear.o: 
	$(CXX) -o probLinear.o $(ROOTCFLAGS) $(CXXFLAGS) -c probLinear.cc

muOsc: muOsc.o $(lib3p) 
	$(RM) $@
	$(CXX) -o $@ $(CXXFLAGS) -L. $^ $(ROOTLIBS)


.PHONY: muOsc.o
muOsc.o: 
	$(CXX) -o muOsc.o $(ROOTCFLAGS) $(CXXFLAGS) -c muOsc.cc

	
.PHONY: all
all: $(targets)
	
	
.PHONY: clean
clean:
	$(RM) $(targets) *.o
	
emptyrule:: $(lib3p)

shared: $(OBJS)
	$(CXX) -shared  -o lib$(LIBNAME).so $(ROOTCFLAGS) $(CXXFLAGS) $(OBJS)
	mv lib$(LIBNAME).so ../lib/.
