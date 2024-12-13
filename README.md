[![Build Status](https://github.com/TEOS-10/GibbsSeaWater.jl/workflows/CI/badge.svg)](https://github.com/TEOS-10/GibbsSeaWater.jl/actions)
[![Build status Windows](https://ci.appveyor.com/api/projects/status/77kj4lug424x20y9/branch/master?svg=true)](https://ci.appveyor.com/project/Alexander-Barth/gibbsseawater-jl-ojx2d/branch/master)
[![codecov](https://codecov.io/github/TEOS-10/GibbsSeaWater.jl/graph/badge.svg?token=gQtGTvWWLK)](https://codecov.io/github/TEOS-10/GibbsSeaWater.jl)

# GibbsSeaWater.jl

GibbsSeaWater.jl is a Julia wrapper for [GSW-C](https://github.com/TEOS-10/GSW-C/), which is the C implementation of the Thermodynamic Equation of Seawater 2010 (TEOS-10).

## Installation

Start Julia and issue the following commands:

```julia
using Pkg
Pkg.add("GibbsSeaWater")
```

## Example

For arrays, one should use the [vectorized "dot" operator](https://docs.julialang.org/en/v1/manual/mathematical-operations/#man-dot-operators-1):

```julia
C = [45.8;34.7]
T = [28.9;22.8]
P = [10.0;50.0]
SP = gsw_sp_from_c.(C,T,P)
```

## About TEOS-10

Please check the [official site](http://www.teos-10.org) and [official repository](https://github.com/TEOS-10), which provide the official implementations (C/Fortran/Matlab/PHP) and the wrappers.
