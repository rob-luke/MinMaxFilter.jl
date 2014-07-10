using MinMaxFilter
using Base.Test


# #######
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
