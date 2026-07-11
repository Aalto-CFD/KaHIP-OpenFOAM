![OpenFOAM v14](https://img.shields.io/badge/OpenFOAM-14-brightgreen)
# KaHIP-OpenFOAM
The [KaHIP](https://github.com/KaHIP/KaHIP) distributor wrapper for [OpenFOAM](https://github.com/OpenFOAM/OpenFOAM-dev).

## Overview
A load balancing distributor that repartitions the mesh with the parallel
[ParHIP](https://github.com/KaHIP/KaHIP) graph partitioner.
It is a drop-in `distributor` model, so any solver that redistributes the mesh
at run time can use it without code changes.

Partitioning is done through the C interface of the ParHIP shared library,
which is provided by an external KaHIP installation.

## Usage
Select the `kahip` distributor in `system/decomposeParDict` and load the
library there:
```cpp
distributor     kahip;
libs            ( "libkahipDecomp_aalto.so" );

kahip
{
    // Optional. ParHIP preconfiguration (default: fast)
    // One of: ultrafast, fast, eco, ultrafastSocial, fastSocial, ecoSocial
    method      fast;

    // Optional. Random seed (default: 0)
    seed        0;

    // Optional. Print the ParHIP output (default: no)
    verbose     no;
}
```
As a distributor it redistributes onto the running processors only, so
`numberOfSubdomains` must equal the number of processors the solver runs on.

### Limitations
- ParHIP supports a single weight per vertex, so multiple cell weights are
  combined by summation.
- Edge weights are not passed to ParHIP.
- `processorWeights` has no ParHIP equivalent and raises a fatal error.

## Compilation
Clone the repository
```sh
git clone git@github.com:Aalto-CFD/KaHIP-OpenFOAM.git
cd KaHIP-OpenFOAM
```
with OpenFOAM sourced, then point the build at a KaHIP installation, e.g.:
```sh
export KAHIP_INCLUDE_DIR=$(spack location -i kahip)/include
export KAHIP_LIB_DIR=$(spack location -i kahip)/lib
```
and compile with one of the two options below.

#### Option A: Compile into your OpenFOAM `$WM_PROJECT_SITE` (default)
```sh
wmake -all $PWD
```

#### Option B: Compile into your local user directory
```sh
FOAM_SITE_LIBBIN=$FOAM_USER_LIBBIN wmake -all $PWD
```
This places `libkahipDecomp_aalto.so` and the symlinked KaHIP libraries under
your `$FOAM_USER_LIBBIN` directory.
