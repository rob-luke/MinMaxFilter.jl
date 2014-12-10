module MinMaxFilter

# This is a port of the Lemire min max filter as implemented by Bruno Luong
# http://arxiv.org/abs/cs.DS/0610046
# http://lemire.me/
# http://www.mathworks.com/matlabcentral/fileexchange/24705-min-max-filter

export minmax_filter

using Base.Cartesian

type Wedge
    buffer::AbstractArray
    size::Int
    n::Int
    first::Int
    last::Int
    mxn::Int
end


for N = 2:4
    @eval begin
    function minmax_filter{T <: Number}(A::Array{T, $N}, window::Array{Int, 1}; verbose::Bool=true)

        if verbose; println("Min max filter on $(length(size(A))) dimensions"); end

        maxval_temp = copy(A); minval_temp = copy(A)

        for dim = 1:$N

            # For all but the last dimension
            @nloops $(N-1) i maxval_temp begin

                # Create index for full array (fa) length
                @nexprs $(N)   j->(fa_{j} = 1:size(maxval_temp)[j])
                @nexprs $(N-1) j->(fa_{j} = i_j)

                # Create index for short array (sa) length
                @nexprs $(N)   j->(sa_{j} = 1:size(maxval_temp)[j] - window[dim] + 1)
                @nexprs $(N-1) j->(sa_{j} = i_j)

                # Filter the last dimension
                (@nref $N minval_temp sa) = min_filter(vec( @nref $N minval_temp fa), window[dim])
                (@nref $N maxval_temp sa) = max_filter(vec( @nref $N maxval_temp fa), window[dim])

            end

            # Circular shift the dimensions
            maxval_temp = permutedims(maxval_temp, mod([1:$N], $N)+1)
            minval_temp = permutedims(minval_temp, mod([1:$N], $N)+1)

        end

        # The dimensions to extract
        @nexprs $N j->(a_{j} = 1:size(A, j)-window[j]+1)

        # Extract set dimensions
        maxval_out = @nref $N maxval_temp a
        minval_out = @nref $N minval_temp a

        return minval_out, maxval_out
    end
    end
end


function minmax_filter{T <: Number}(a::AbstractArray{T}, window::Int; verbose::Bool=false)

    if verbose; println("Running min max filter on array of length $(length(a)) with window length $window"); end

    n = length(a)

    # Initialise the output variables
    # This is the running minimum and maximum over the specified window length
    minval = zeros(T, 1, n-window+1)
    maxval = zeros(T, 1, n-window+1)

    # Initialise the internal wedges
    # U[1], L[1] are the location of the global maximum and minimum
    # U[2], L[2] are the maximum and minimum over (U1, inf)
    L = Wedge(zeros(1,window+1), window+1, 0, 1, 0, 0)          # Min
    U = Wedge(zeros(1,window+1), window+1, 0, 1, 0, 0)

    for i = 2:n
        if i > window
            if ~wedgeisempty(U)
                maxval[i-window] = a[getfirst(U)]
            else
                maxval[i-window] = a[i-1]
            end
            if ~wedgeisempty(L)
                minval[i-window] = a[getfirst(L)]
            else
                minval[i-window] = a[i-1]
            end
        end # window

        if a[i] > a[i-1]
            L = pushback(L, i-1)
            if i==window+getfirst(L); L=popfront(L); end
            while ~wedgeisempty(U)
                if a[i] <= a[getlast(U)]
                    if i == window+getfirst(U); U = popfront(U); end
                    break
                end
                U = popback(U)
            end

        else

            U = pushback(U, i-1)
            if i==window+getfirst(U); U=popfront(U); end

            while ~wedgeisempty(L)
                if a[i] >= a[getlast(L)]
                    if i == window+getfirst(L); L = popfront(L); end
                    break
                end
                L = popback(L)
            end

        end  # a>a-1

        #=if verbose=#
            #=println("---- i=$i")=#
            #=println(L.buffer[mod([L.first+(-1:L.n-2)], L.size)+1])=#
            #=println(U.buffer[mod([U.first+(-1:U.n-2)], U.size)+1])=#
        #=end=#

    end # for i

    i = n+1
    if ~wedgeisempty(U)
        maxval[i-window] = a[getfirst(U)]
    else
        maxval[i-window] = a[i-1]
    end

    if ~wedgeisempty(L)
        minval[i-window] = a[getfirst(L)]
    else
        minval[i-window] = a[i-1]
    end

    return minval, maxval
end


function min_filter(a::AbstractArray, window::Int; verbose::Bool=false)

    minval, maxval = minmax_filter(a, window, verbose=verbose)

    return minval
end


function max_filter(a::AbstractArray, window::Int; verbose::Bool=false)

    minval, maxval = minmax_filter(a, window, verbose=verbose)

    return maxval
end


function wedgeisempty(X::Wedge)
    X.n <= 0
end

function pushback(X::Wedge, v)
    X.last = mod(X.last, X.size) + 1
    X.buffer[X.last] = v
    X.n = X.n+1
    X.mxn = maximum([X.mxn, X.n])
    return X
end

function getfirst(X::Wedge)
    X.buffer[X.first]
end

function getlast(X::Wedge)
    X.buffer[X.last]
end

function popfront(X::Wedge)
    X.n = X.n-1
    X.first = mod(X.first, X.size) + 1
    return X
end

function popback(X::Wedge)
    X.n = X.n-1
    X.last = mod(X.last-2, X.size) + 1
    return X
end



end # module




