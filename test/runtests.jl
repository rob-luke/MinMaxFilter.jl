using MinMaxFilter
using Base.Test
using MAT


# #######
#
# 1 dimension
#
# Compare to matlab slowminmaxfilt_algo.m
#
# t = [1:1024] ./1024; d = sin(2*pi*4*t);
# [minval, maxval] = slowminmaxfilt_algo(d, 100)
# dlmwrite('minmax_output.txt',[minval; maxval], 'delimiter', '\t', 'precision', '%.12f')
#
# ######

matlab = readdlm(joinpath(dirname(@__FILE__), "data", "minmax_output.txt"),'\t')

t = (1:1024)./1024
d = sin(2*pi*4*t)
minval, maxval = minmax_filter(d, 100, verbose=true)

@test_approx_eq minval matlab[1,:]
@test_approx_eq maxval matlab[2,:]


# #######
#
# 2 dimension
#
# Compare to matlab minmaxfilt.m
#
# [X,Y,Z] = peaks(100);
# smax = minmaxfilt(Z, 11, 'max')
#
# ######

filen = matopen("data/2d_array.mat")
A = read(filen, "Z")
A = convert(Array{FloatingPoint}, A)
close(filen)

filen = matopen("data/2d_array_max11.mat")
max_matlab = read(filen, "smax")
max_matlab = convert(Array{FloatingPoint}, max_matlab)
close(filen)

maxval = minmax_filter(A, 11, verbose=false)

@test_approx_eq maxval max_matlab


# #######
#
# 3 dimension
#
# Compare to matlab minmaxfilt.m
#
# amax=minmaxfilt(image,5,'max','valid');
#
# ######


filen = matopen("data/3d_array.mat")
A = read(filen, "image")
A = convert(Array{FloatingPoint}, A)
close(filen)

filen = matopen("data/3d_array_max5.mat")
max_matlab = read(filen, "amax")
max_matlab = convert(Array{FloatingPoint}, max_matlab)
close(filen)

maxval = minmax_filter(A, 5, verbose=false)

@test_approx_eq maxval max_matlab

