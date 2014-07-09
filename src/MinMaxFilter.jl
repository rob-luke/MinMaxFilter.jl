module MinMaxFilter

# This is a port of the Lemire min max filter as implemented by Bruno Luong
# http://arxiv.org/abs/cs.DS/0610046
# http://lemire.me/
# http://www.mathworks.com/matlabcentral/fileexchange/24705-min-max-filter

export minmax_filter

type Wedge
    buffer::AbstractArray
    size::Int
    n::Int
    first::Int
    last::Int
    mxn::Int
end

function minmax_filter(a::AbstractArray, window::Int; verbose::Bool=false)

    if verbose; println("Running min max filter on array of length $(length(a))"); end

    n = length(a)

    # Initialise the output variables
    # This is the running minimum and maximum over the specified window length
    minval = zeros(FloatingPoint, 1, n-window+1)  # TODO: can you initialise as empty or nans?
    maxval = zeros(FloatingPoint, 1, n-window+1)

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

        if verbose
            println("---- i=$i")
            println(L.buffer[mod([L.first+(-1:L.n-2)], L.size)+1])
            println(U.buffer[mod([U.first+(-1:U.n-2)], U.size)+1])
        end

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



