# Known Clusters {#known-clusters}

## Centrum Informatyczne Świerk (CIS), Warsaw Cluster
![image](https://github.com/mach3-software/MaCh3/assets/45295406/9c8b2900-d3c2-4a04-a6a5-99ff2c7d10b0)

```
source /mnt/opt/spack/0.20/share/spack/setup-env.sh
spack load gcc@8.5.0%gcc@=8.5.0
spack load cmake@3.26.3%gcc@=8.5.0
spack load python@3.8.12%gcc@=8.5.0
spack load binutils@2.30%gcc@=8.5.0
spack load vdt@0.4.3
spack load git@2.40.0%gcc@=8.5.0
spack load imake
spack load makedepend
spack load xz@5.4.1%gcc@=8.5.0 #This holds LZMA
spack load libxpm@3.5.12%gcc@=8.5.0
spack load ncurses@6.4%gcc@8.5.0
spack load davix@0.8.1%gcc@8.5.0
spack load openssl@1.1.1t%gcc@=8.5.0
source /mnt/home/share/t2k/kskwarczynski/CERN_ROOT/root-6.24.08/build_root/bin/thisroot.sh
export LD_LIBRARY_PATH=/mnt/opt/spack/0.20/opt/spack/linux-centos7-ivybridge/gcc-8.5.0/vdt-0.4.3-s6mhw6a6pl7zdfbnaekxl5ad3x7jayvn/lib:$LD_LIBRARY_PATH
```
On GPU nodes
```
export CUDAPATH=/usr/local/cuda-11.4
```
```
System: CentOS Linux 7
```

## LXPLUS, CERN Cluster
![image](https://github.com/mach3-software/MaCh3/assets/45295406/1fb468ab-c923-4869-8225-d5ae2473aeb3)

Works with default configuration:
```
System: Red Hat Enterprise Linux 9.3 
GCC:    11.4.1
CMake:  3.20.2
ROOT:   6.30/04
```

## Cedar, Digital Research Alliance of Canada
![image](https://github.com/mach3-software/MaCh3/assets/45295406/bfd11868-31b2-4f71-82d5-c4503c919434)
Only need
```
module load root
```
Then configuration looks like:
```
System: CentOS Linux 7
GCC:    12.3.1
CMake:  3.27.7
ROOT:   6.28/06 
```

## Your Cluster?