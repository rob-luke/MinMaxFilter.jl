# MinMaxFilter

[![Build Status](https://travis-ci.org/codles/MinMaxFilter.jl.svg?branch=master)](https://travis-ci.org/codles/MinMaxFilter.jl)
[![Coverage Status](https://coveralls.io/repos/codles/MinMaxFilter.jl/badge.png?branch=master)](https://coveralls.io/r/codles/MinMaxFilter.jl?branch=master)


## Overview

This is a port of the Lemire min max filter as implemented by Bruno Luong.  
http://arxiv.org/abs/cs.DS/0610046  
http://lemire.me/  
http://www.mathworks.com/matlabcentral/fileexchange/24705-min-max-filter

## Install

```julia
Pkg.clone("https://github.com/codles/MinMaxFilter.jl.git")

```


## Example usage

### Three dimensional peak finding

Read in three dimensional data, plot and find peaks.  
Uses Winston and EEG packages for plotting.

```julia
using Winston
using EEG

filen = matopen(joinpath(dirname(@__FILE__), "data", "3d_array.mat"))
A = read(filen, "image")
A = convert(Array{FloatingPoint}, A)
close(filen)

p = plot_dat(A)

minval, maxval = minmax_filter(A, 6, verbose=false)

matching = A[2:size(maxval)[1]+1, 2:size(maxval)[2]+1, 2:size(maxval)[3]+1]
matching = matching .== maxval

peaks = maxval[matching]
peaks = peaks[peaks .>= 0.1 * maximum(peaks)]

for l = 1:length(peaks)
    x, y, z = ind2sub(size(A), find(A .== peaks[l]))
    println("x=$x, y=$y, z=$z, val=$(A[x,y,z])")
    oplot(p, x, y, z)
end
```

![Dipoles](/test/images/3D-dipole.png)
