# FAQ {#faq}

## How to start?
We recomend starting MaCh3 by using [MaCh3 Tutorial](https://github.com/mach3-software/MaCh3Tutorial)

## Error: Need MACH3 environment variable
Some functionalities rely on setting `Env{MACH3}` which should point to path experiment specific MaCh3. This way MaCh3 can easily find `Env{MACH3}/inputs/SomeInput.root` for example. Safest approach is to inlcude this in experiemnt specyfic cmake

This is example of T2K **setup.MaCh3T2K.sh.in**
```
export MACH3T2K_ROOT=@CMAKE_INSTALL_PREFIX@
export MACH3=@CMAKE_INSTALL_PREFIX@
export MaCh3T2K_VERSION=@MaCh3T2K_VERSION@

add_to_PATH "${MACH3T2K_ROOT}/bin"
add_to_LD_LIBRARY_PATH ${MACH3T2K_ROOT}/lib
```

## Do I need GPU to use MaCh3
No, MaCh3 works with a CPU-only setup. In fact, MaCh3 is predominantly a CPU library. Only most performance intensive operations have GPU support at accelerate. We make rigorous checks to make sure the CPU and GPU produce the same results up to numeric precision.

## Version
MaCh3 is version controlled. You can find all releases [here](https://github.com/mach3-software/MaCh3/releases) you can also find easily what changed in each version [here](https://github.com/mach3-software/MaCh3/wiki/0.1.-History)

## Code Structure
<img width="879" alt="Structure" src="https://github.com/user-attachments/assets/47a4ec31-4a60-4603-ac84-359281c715b3">

Last updated: 12.12.2025