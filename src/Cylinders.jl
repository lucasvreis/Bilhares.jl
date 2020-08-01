module Cylinders

using OffsetArrays, StaticArrays
import Base: show

export Cylinder, intersect, univ, set, σ, σi

Cylinder{Usize,N} = OffsetArray{Bool,1,<:StaticArray{Tuple{N,Usize},Bool}}

Cylinder{u}(i, m::Vararg{Set{Int},N}) where {u,N} = Cylinder{u,N}(i, OffsetVector(SVector{N,Set{Int}}(m...), i-1))
Cylinder{u,N}(i, m::Vararg{Set{Int},N}) where {u,N} = Cylinder{u,N}(i, OffsetVector(SVector{N,Set{Int}}(m...), i-1))

size(c::Cylinder) = size(c.data)
getindex(c::Cylinder, i...) = getindex(c.data, i...)
firstindex(c::Cylinder) = firstindex(c.data)
lastindex(c::Cylinder) = lastindex(c.data)	
iterate(c::Cylinder) = iterate(c.data)
iterate(c::Cylinder,t) = iterate(c.data,t)

Base.show(io::IO, c::Cylinder) = print(io, "[$(c.i),$(join("{".*join.(c.data,",").*"}",","))]")


@inline σ(c::Cylinder{u,n}) where {u,n} =
                c.i == 0 ? 
                    Cylinder{u,n-1}(0, c[1:end]...) : 
                    Cylinder{u,n}(c.i-1, c...)

@inline σi(c::Cylinder{u,n}) where {u,n} = Cylinder{u,n}(c.i+1,c...)

@inline inv(::typeof(σ)) = σi
@inline (::typeof(σ))^(n::Int) = n == 0 ? c -> c : c -> (σ^(n-1))(σ(c))
@inline function ^(::typeof(σi), n::Int)
    function ξ(c::Cylinder{u,k}) where {u,k} 
        Cylinder{u,k}(c.i+n,c...) 
    end
end

set(x...) = Set(x)
univ(u) = Set(1:u)

# p(::Cylinder) = ∑()

function intersect(A::Cylinder{u,n}, B::Cylinder{u,m}) where {u,n,m}
    ai = A.i; af = ai + n - 1
    bi = B.i; bf = bi + m - 1

    #1º bloco - os cilindros não se intersectam
    if af < bi
        return Cylinder{u}(ai, A..., repeat([univ(u)],bi-af)..., B...)
    elseif bf < ai
        return Cylinder{u}(bi, B..., repeat([univ(u)],ai-bf)..., A...)
    end

    #2º bloco = os cilindros se intersectam
    # há várias maneiras disso acontecer. mas todos os casos
    # podem ser vistos como H I J

    hfim = max(ai,bi) - 1 # menor ou igual a af
    ii = max(ai,bi) # menor ou igual a af
    ifim = min(bf,af)
    ji = min(bf,af) + 1
    @show hfim, ii, ifim, ji
    return Cylinder{u}(min(ai,bi), A[begin:hfim]..., B[begin:hfim]...,
                           (A[ii:ifim] .∩ B[ii:ifim])..., 
                           A[ji:end]..., B[ji:end]...)
end

end