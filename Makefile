UNAME := $(shell uname)

# Because Apple computers are special
ifeq ($(UNAME), Darwin)
LIB_EXT="dylib"
else
LIB_EXT="so"
endif

# Then we also need the MACH3 environement variable to exist
ifndef MACH3
$(error MACH3 is not set! Have you sourced setup.sh?)
endif

##############################
# ----------------------------
# The default make
# ----------------------------
all:

# Perform a check of libconfig so we don't have to make config && make
ifeq "$(wildcard libagf-g.a)" ""
# Perform a second check to see if we need to rebuild the external entrirely
# Or just copy the already made libraries to our lib folder
ifeq "$(wildcard libconfig/lib/.libs/*)" ""
	make config
else
	make cpconfig
endif
endif

# Check that bin directory is made
ifeq "$(wildcard bin)" ""
	make links
endif

ifeq "$(wildcard lib)" ""
	make links
endif
	./setup.sh

	make -j4 -C Prob3++ shared
	make -j4 -C throwParms
	make -j4 -C manager
	make -j4 -C covariance
	make -j4 -C splines
	make -j4 -C samplePDF
	make -j4 -C mcmc
	make -j4 -C Diagnostics

# The make for libconfig
config:
	cd libconfig && ./configure --disable-examples && make 
	cp libconfig/lib/.libs/lib*.$(LIB_EXT)* lib/.

# The simpler case wher we have already made the externals but need to copy them
cpconfig:
	cp libconfig/lib/.libs/lib*.$(LIB_EXT)* lib/.

# The make for the links
links:
	@ [ -d lib ]     || mkdir lib
	@ [ -d bin ]     || mkdir bin

# The usual clean
clean:
	make clean -C Prob3++
	make clean -C throwParms
	make clean -C manager
	make clean -C covariance
	make clean -C splines
	make clean -C samplePDF
	make clean -C mcmc
	make clean -C Diagnostics

	rm -rf lib/lib*

# Removes more directories ontop of the usual clean
vclean:
	make clean
	make clean -C libconfig
