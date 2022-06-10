CXX=g++
AR=ar
LD_SHARED=g++

ROOTLIBS := $(shell root-config --glibs --libs)
ROOTINCLUDES := $(shell root-config --prefix)/include 
ROOTCFLAGS := $(shell root-config --cflags)

INCDIR = -I. -I..
LIBDIR = ../lib

UNAME := $(shell uname)

# Darwin Mac OS X
ifneq ($(UNAME), Darwin)
EXPORT_DYNAMIC="-Wl,-export-dynamic"

else

# Flags to detect 32-bit or 64-bit OS platform
OS_SIZE = $(shell uname -m | sed -e "s/i.86/32/" -e "s/x86_64/64/")
OS_ARCH = $(shell uname -m | sed -e "s/i386/i686/")

# These flags will override any settings
ifeq ($(i386),1)
        OS_SIZE = 32
        OS_ARCH = i686
endif

ifeq ($(x86_64),1)
        OS_SIZE = 64
        OS_ARCH = x86_64
endif

SIZE_FLAG=-m$(OS_SIZE)

endif

# Usual CXX flags
CXXFLAGS   = -Wall -g -O3 -fPIC $(ROOTCFLAGS) $(INCDIR) -Werror -std=c++11
# If we want multi-threading
ifdef MULTITHREAD
CXXFLAGS += -fopenmp
OMP_DEFINES=-DMULTITHREAD
endif

SYSLIB     = -lm -lg2c
LINK_ARGS_BIN = $(SYSLIB) $(ROOTLIBS)
SOFLAGS= -shared $(EXPORT_DYNAMIC) $(ROOTCFLAGS)
LDFLAGS=
ARFLAGS= rcsv

# Define all the CUDA stuff here                                       
ifdef CUDAPATH
# Include the CUDA libraries and helper functions
CUDA_LIBS := -L$(CUDAPATH)/lib64 -I$(CUDAPATH)/include -I$(CUDAPATH)/common/inc -lcudart -I$(CUDAPATH)/samples/common/inc

# The CUDA objects we want to make
CUDA_OBJECTS := gpuSplineUtils.o SplineMonolith.o

# CUDA PTX options
CUDA_PTX="-dlcm=ca -allow-expensive-optimizations=true -fmad=true -O3 -warn-lmem-usage -warn-spills"

# CUDA optimisations
CUDA_OPTI=-prec-sqrt=false -use_fast_math -O3 -Werror cross-execution-space-call

# cross-compilar options
CUDA_XCOMP="-fpic -O3 -Wall -Wextra -Werror" -c

# Enable error checking on CUDA code
# Slows down GPU code slightly but on microsecond scale
CUDA_FLAGS=-DCUDA_ERROR_CHECK

CUDAVER=$(shell nvcc --version | grep -o 'V[0-9].*' |  cut -d. -f1 |sed 's:V::')


# Additional flags for CUDA debug (-g to debug g++, -G to debug CUDA)
ifdef DEBUG_CUDA_ND280
DEBUG_DEFINES+=-DDEBUG
CUDA_DEFINES+=-DDEBUG
CUDA_FLAGS+=-g -G
# Remove optimisations when debugging (having both is not supported by many CUDA versions)
CUDA_PTX="-dlcm=ca -warn-lmem-usage -warn-spills -v"
CUDA_OPTI=
endif

ifdef DEBUG_ND280_DUMP
DEBUG_DEFINES += -DDEBUG_DUMP
endif

#-Xptxas "-dlcm=ca -v --allow-expensive-optimizations=true -fmad=true"
#-Xcompiler "-fpic" -c
#-prec-sqrt=false -use_fast_math

#CWRET comment: could change arch here. also, prec-div is only used in (seg+khig)/2; could replace by *0.5
# -prec-sqrt is not needed (no sqrts in program code!), -use_fast_math forces intrinsic
#  tried this on GTX 660 and no speed increase, buu and P6 6c
#  don't use fastmath!

#CWret: -g and -G for debug flags to use cuda-gdb; slows stuff A LOT
#-pxtas-options=-v, -maxregcount=N

#-Xptxas -dlcm=ca increases L1 cache from 16kB to 48kB for Fermi
#-Xptxas -dlcm=cg caches global memory into shared memory and turns off the L1 cache entirely
#

%.o:%.cu

ifeq ($(shell expr $(CUDAVER) \>= 11), 1)
	nvcc $(SIZE_FLAG) \
		-arch=sm_52 \
		-gencode arch=compute_52,code=sm_52 \
		-gencode arch=compute_60,code=sm_60 \
		-gencode arch=compute_61,code=sm_61 \
		-gencode arch=compute_70,code=sm_70 \
		-gencode arch=compute_75,code=sm_75 \
		-gencode arch=compute_80,code=sm_80 \
		-gencode arch=compute_86,code=sm_86 \
		-gencode arch=compute_86,code=compute_86 \
                $(CUDA_OPTI) \
                -Xptxas $(CUDA_PTX) \
                -Xcompiler $(CUDA_XCOMP) \
                $(CUDA_FLAGS) $(CUDA_LIBS) $< $(CUDA_DEFINES)

else ifeq (($(shell expr $(CUDAVER) \>= 9)&&($(shell expr $(CUDAVER) \<= 10))), 1)
	nvcc $(SIZE_FLAG) \
		-gencode arch=compute_30,code=sm_30 \
		-gencode arch=compute_32,code=sm_32 \
		-gencode arch=compute_35,code=sm_35 \
		-gencode arch=compute_35,code=compute_35 \
		$(CUDA_OPTI) \
		-Xptxas $(CUDA_PTX) \
		-Xcompiler $(CUDA_XCOMP) \
		$(CUDA_FLAGS) $(CUDA_LIBS) $< $(CUDA_DEFINES) 
else
	nvcc $(SIZE_FLAG) \
		-gencode arch=compute_30,code=sm_30 \
		-gencode arch=compute_32,code=sm_32 \
		-gencode arch=compute_35,code=sm_35 \
		-gencode arch=compute_35,code=compute_35 \
		$(CUDA_OPTI) \
		-Xptxas $(CUDA_PTX) \
		-Xcompiler $(CUDA_XCOMP) \
		$(CUDA_FLAGS) $(CUDA_LIBS) $< $(CUDA_DEFINES) 
endif
endif
#primary access global memory and uses less that 16kB of shared memory, that's us!
#can explicitly cache global memory into shared memory and turn off L1 cachewith -Xptxas -dlcm=cg

%.o:%.cpp
	g++ -c $(CXXFLAGS) $< $(OMP_DEFINES) $(DEBUG_DEFINES)

libSplines: libSplines.a $(LIBDIR)/libSplines.so

# APPEND HERE
libSplines.a: splineBase.o splineFDBase.o $(CUDA_OBJECTS)
	$(AR) $(ARFLAGS) $@ $^

# AND HERE
$(LIBDIR)/libSplines.so: splineBase.o splineFDBase.o $(CUDA_OBJECTS)
	$(LD_SHARED) $(SOFLAGS) $^ -o $@ $(LIB_OBJECTS) $(ROOTLIBS)


# Cleaners
clean:
	rm -f ../lib/libSplines.so	
	rm -f *.o lib*

distclean:
	rm -fv *~ *.o *.so *.a lib*

