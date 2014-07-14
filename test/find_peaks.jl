using MAT
using MinMaxFilter


# ######
#
# 2D
#
# ######

filen = matopen(joinpath(dirname(@__FILE__), "data", "2d_array.mat"))
A = read(filen, "Z")
A = convert(Array{FloatingPoint}, A)
close(filen)

minval, maxval = minmax_filter(A, 11, verbose=false)

matching = A[2:size(maxval)[1]+1, 2:size(maxval)[2]+1]
matching = matching .== maxval

peaks = maxval[matching]

for l = 1:length(peaks)
    x, y = ind2sub(size(A), find(A .== peaks[l]))
    println("x=$x, y=$y, val=$(A[x,y])")
end


# ######
#
# 3D
#
# #####

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

file(p, "images/3D-dipole.png", width=1000, height=1000)
