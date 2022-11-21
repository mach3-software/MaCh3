#!/bin/bash
#DO NOT DELETE THIS COMMENT OR THINGS WILL BREAK - SetMeUp performs a search and replace and will remove the following line otherwise

# e.g. this:
#export PATH=/home/cwret/procmail-3.22/new:${PATH}
#source /home/cwret/CMT/setup.sh
#source /home/cwret/root/bin/thisroot.sh
#module load 
#export CUDAPATH=${CUDA_HOME}
#export MACH3_DATA=/home/cwret/P6Data
#export MACH3_MC=/home/cwret/P6MC

#source /vols/t2k/users/ea2817/Software/ND280/highland2_install/ROOT/v5r34p34n00/Linux-x86_64/bin/thisroot.sh
#module load 
#export CUDAPATH=${CUDA_HOME}
export MACH3=${PWD}
#export MACH3_DATA=/vols/t2k/users/ea2817/MaCh3_OA2019_inputs/asg2019oa/Data_nd280Psyche_v3r47 
#export MACH3_MC=/vols/t2k/users/ea2817/MaCh3_OA2021_inputs/ND280_inputs


#################################################################

if [[ $HOSTNAME == *gpu.rl.ac.uk ]]; then
    # Add lockfile to path (required on emerald!) (commented out for those not on EMERALD)
    export PATH=$PATH:/home/oxford/eisox159/procmail/bin

    # setup root on EMERALD (commented out for those not on EMERALD)
    source /home/stfc/eisext13/root/bin/thisroot.sh
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/stfc/eisext13/root/lib
    export PATH=$PATH:/home/stfc/eisext13/root/bin
    export ROOTLIBDIR=/home/stfc/eisext13/root/lib
    export ROOTINCLUDEDIR=/home/stfc/eisext13/root/include
    module load cuda/8.0.44
    export CUDAPATH=${CUDA_HOME}
fi

# Cuda directory
#   Unset if you don't have a CUDA GPU
#   If you want GPU point CUDAPATH to your CUDA install and makefile will pick this up
#if [[ $HOSTNAME == *hep.ph.ic.ac.uk ]]; then
#  export CUDAPATH=/usr/local/cuda-8.0
#fi

# Multithreading? 
# MULTITHREAD=1 means MP is on, if environment variable doesn't exist it's off
unset MULTITHREAD
export MULTITHREAD=1

# Automatically set CUDA for Emerald
if [[ $HOSTNAME == *gpu.rl.ac.uk ]]; then
  module load cuda/8.0.44
  export CUDAPATH=${CUDA_HOME}
# Automatically set up for ComputeCanada
elif [[ $HOSTNAME == lg-1r[47]* ]]; then
  module load CUDA/7.5.18
  export CUDAPATH=${CUDA_HOME}

  # Also need GSL for some external dependencies
  module load GSL
elif [[ $HOSTNAME == *hep.ph.ic.ac.uk ]]; then
  
  #export CUDAPATH=/usr/local/cuda-10.2
  echo $CUDAPATH
    if [ -z $CUDAPATH ]; then
	source /vols/software/cuda/setup.sh #Can give this an argument to get a particular CUDA version, e.g. 10.2.2
	export CUDAPATH=$CUDA_PATH
    fi
#export CUDAPATH=$CUDA_PATH
fi

# Set the MaCh3 directory
if [ -z $MACH3 ]; then
  export MACH3=$(pwd)
fi

# Add the CUDA libraries and executables if we have CUDA enabled
if [ -d "$CUDAPATH" ]; then
  # Add the CUDA binaries (e.g. nvcc!)
  export PATH=${CUDAPATH}/bin:${PATH}
  # Add the CUDA libraries (e.g. cudarand)
  export LD_LIBRARY_PATH=${CUDAPATH}/lib64:/usr/lib64/nvidia:${LD_LIBRARY_PATH}
  # Report on what CUDA settings and GPUs we're running
  echo "******************************************"
  echo Running GPU on: $(nvidia-smi --query-gpu=name --format=csv,noheader)
  echo with CUDAPATH = $CUDAPATH
  echo "******************************************"
else
  # Report that we aren't running with GPU
  echo "******************************************"
  echo DID NOT SET CUDAPATH!
  echo This is OK if you do not want to run with a GPU
  echo "******************************************"
fi

mkdir -pv ${MACH3}/lib

echo ========================================================
echo "Are you resetting to change GPU/CPU options?"
echo "If so, don't forget to fully remake!"
echo "'make clean && make'"
echo ========================================================
echo "Finished setup"
