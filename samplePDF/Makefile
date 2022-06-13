CXX=g++
AR=ar
LD_SHARED=g++

UNAME := $(shell uname)
OS_SIZE = $(shell uname -m | sed -e "s/i.86/32/" -e "s/x86_64/64/")

ifneq ($(UNAME), Darwin)
EXPORT_DYNAMIC="-Wl,-export-dynamic"
else # Darwin Mac OS X

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
endif # End Mac OS X BS

# Debugs event rates, checks GPU to CPU weights
# DO NOT RUN FOR PRODUCTION
ifdef DEBUG_ND280
DEBUG_DEFINES = -DDEBUG=${DEBUG_ND280}
endif

# Dumps ROOT file of CPU and GPU weights
# DO NOT RUN FOR PRODUCTION
ifdef DEBUG_ND280_DUMP
DEBUG_DEFINES += -DDEBUG_DUMP
endif


ROOTLIBS := $(shell root-config --glibs --libs) -lMinuit
ROOTINCLUDES := $(shell root-config --prefix)/include 
ROOTCFLAGS := $(shell root-config --cflags) -I../ -I../utils/agf/include

INCDIR = -I. -I..
LIBDIR = ../lib

CXXFLAGS   = -Wall -O3 -g -fPIC $(ROOTCFLAGS) $(INCDIR) -Werror -std=c++11
SYSLIB     = -lm
LINK_ARGS_BIN = $(SYSLIB) $(ROOTLIBS)
SOFLAGS= -shared $(EXPORT_DYNAMIC) $(ROOTCFLAGS)
LDFLAGS=
ARFLAGS= rcsv

# openMP defines
ifdef MULTITHREAD
CXXFLAGS += -fopenmp
OMP_DEFINES=-DMULTITHREAD
ifeq ($(DEBUG_ND280),2)
ifneq ($(OMP_NUM_THREADS),1)
$(error I WAS TOLD TO DO DEBUG IN MP MODE, NEED OMP_NUM_THREADS=1)
endif # End the OMP NUM THREADS
endif # End the DEBUG check
endif # End the multithreads

LIB_OBJECTS += $(ROOTLIBS)

OBJECTS := samplePDFBase.o samplePDFFDBase.o samplePDFFDBase.o interfacePDFEbE.o Structs.o

ifdef COMBINED_KIN
CXXFLAGS += -DCOMBINED_KIN
endif

# Define all the CUDA stuff here
ifdef CUDAPATH
CUDA_LIBS := -L$(CUDAPATH)/lib64 -I$(CUDAPATH)/include -lcudart # -lcuda
LIB_OBJECTS += $(CUDA_LIBS)

CUDAVER=$(shell nvcc --version | grep -o 'V[0-9].*' |  cut -d. -f1 |sed 's:V::')

#CUDA_OBJECTS := samplePDFSKBaseGPU.o probGpu.o samplePDFND2019GPU.o kdeGpu.o
CUDA_OBJECTS := kdeGpu.o


# Note, -gencode arch, code works by
# 		arch = architecture for the virtual GPU; what CC does your code require?
# 		code = architecture for the actual GPU; what CC does your card run?
# 		Note: if you only want to make one single gencode, 	make arch small (first stage virtual)
# 									make code for whatever CC card you have
#
#			putting a virtual architecture (e.g. compute_30) in the code field leads to just-in-time compilation (see 5.7.2.1 of CUDA toolkit 7.5, nvcc section)
#			omitting gencode compiles for compute_20 and sm_20
#		-gencode arch=compute_20,code=sm_20 

%.o:%.cu
ifeq ($(shell expr $(CUDAVER) \< 9), 1)
	nvcc $(SIZE_FLAG) \
                -gencode arch=compute_30,code=sm_30 \
                -gencode arch=compute_32,code=sm_32 \
                -gencode arch=compute_35,code=sm_35 \
                -gencode arch=compute_35,code=compute_35 \
                -prec-sqrt=false -use_fast_math -O3 -Werror cross-execution-space-call \
                -Xptxas "-allow-expensive-optimizations=true -fmad=true -O3 -warn-lmem-usage -warn-spills" \
                -Xcompiler "-fpic -O3 -Wall -Wextra -Werror -Wno-error=unused-parameter" -c $<
else
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
                -prec-sqrt=false -use_fast_math -O3 -Werror cross-execution-space-call \
                -Xptxas "-allow-expensive-optimizations=true -fmad=true -O3 -warn-lmem-usage -warn-spills" \
                -Xcompiler "-fpic -O3 -Wall -Wextra -Werror -Wno-error=unused-parameter" -c $<
else
	nvcc $(SIZE_FLAG) \
		-gencode arch=compute_30,code=sm_30 \
		-gencode arch=compute_32,code=sm_32 \
		-gencode arch=compute_35,code=sm_35 \
		-gencode arch=compute_35,code=compute_35 \
                -prec-sqrt=false -use_fast_math -O3 -Werror cross-execution-space-call \
                -Xptxas "-allow-expensive-optimizations=true -fmad=true -O3 -warn-lmem-usage -warn-spills" \
                -Xcompiler "-fpic -O3 -Wall -Wextra -Werror -Wno-error=unused-parameter" -c $<
endif
endif
endif

%.o:%.cpp
	g++ -c $(LINK_ARGS_BIN) $(CXXFLAGS) $< $(OMP_DEFINES) $(DEBUG_DEFINES)

libSamplePDF: libSamplePDF.a $(LIBDIR)/libSamplePDF.so

libSamplePDF.a: $(OBJECTS) $(CUDA_OBJECTS) 
	$(AR) $(ARFLAGS) $@ $^

$(LIBDIR)/libSamplePDF.so: libSamplePDF.a $(OBJECTS) $(CUDA_OBJECTS)
	$(LD_SHARED) $(SOFLAGS) $^ -o $@ $(LIB_OBJECTS) $(ROOTLIBS) -L$(LIBDIR) -lThreeProb_1.00 -lSplines -lCovariance

clean:
	rm -f ../lib/libSamplePDF.so       
	rm -f *.o lib*

distclean:
	rm -fv *~ *.o *.so *.a lib*

