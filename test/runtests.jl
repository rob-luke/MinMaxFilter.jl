using MinMaxFilter
using Base.Test

#
# small example
#

A = zeros(5,5)
A[2,2] = 0.8
A[4,4] = 0.6
minval, maxval = extrema_filter(A, 2)
matching = [vec(A[2:end]) .!= vec(maxval); false]
@test A[reshape(matching, size(A))] == [0.8, 0.6]

A = zeros(5,5,5)
A[2,2,2] = 0.7
A[4,4,3] = 0.5
minval, maxval = extrema_filter(A, 2)
matching = [vec(A[2:end]) .!= vec(maxval); false]
@test A[reshape(matching, size(A))] == [0.7, 0.5]

A = zeros(5,5,5,5)
A[2,2,2,2] = 0.7
A[4,4,3,1] = 0.4
A[3,4,3,2] = 0.5
minval, maxval = extrema_filter(A, 2)
matching = [vec(A[2:end]) .!= vec(maxval); false]
@test A[reshape(matching, size(A))] == [0.4,0.7,0.5]
x, y, z, t = ind2sub(size(A), find(A .== 0.4))
@test x[1] == 4
@test y[1] == 4
@test z[1] == 3
@test t[1] == 1
