module Ergos

export bins!, intervaldynamic

function bins!(x,X)
    n, N = size(X)
    for i in 1:n
        for j in 1:N
            if (j-1)/N ≤ x[i] < (j)/N
                X[i,j] += 1
            end
        end
    end
    X
end

function intervaldynamic(f,N,n,iters)
    avgVisTime = zeros(Int,n,N)
    x0 = rand(n)
    bins!(x0, avgVisTime)
    for j in 1:iters
        x0 = f.(x0)
        bins!(x0, avgVisTime)
    end
    avgVisTime / iters
end

end # module