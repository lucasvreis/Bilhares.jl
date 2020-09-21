export Measure, HyperMeasure, HyperMeasure2, integrate, ∫, correlation, 𝟙

using Strided

# sampler : Int ⟶ Particle[]
struct Measure{B<:Billiard,G<:Function,F<:Function}
    sampler::G
    bd::B
    density::F
end

const 𝟙 = x -> 1

Measure(sampler, bd) = Measure(sampler, bd, 𝟙)

# f : Particle[] ⟶ Real
function integrate(f::Function, μ::Measure, N::Int)
    ps = μ.sampler(N)
    d  = μ.density.(ps)
    return @strided sum(f.(ps) .* d) / sum(d)
end

function integrate(f::Function, μ::Measure, ps::AbstractVector)
    d  = μ.density.(ps)
    return @strided sum(f.(ps) .* d) / sum(d)
end

function integrate(f::Function, μ::Measure{B,G,typeof(𝟙)}, N) where {B,G}
    ps = μ.sampler(N)
    return @strided mean(f.(ps))
end

∫(f,μ,N) = integrate(f,μ,N)

@inline function bla(bd, p)
    v = normalvec(bd[p.obs], p.pos) |> perp
    sqrt(v'*g_klein(p.pos)*v)
end

HyperMeasure(bd::Billiard{T}) where T = 
Measure(bd, p -> bla(bd, p)) do N
    random_on_border(HyperBKParticle, bd, N)
end

# Versão antiga - sampler com densidade da medida

# HyperMeasure(bd::Billiard{T}) where T = 
# Measure(bd) do N
#     l  = htotallength(bd)
#     rφ = 2rand(N) .- 1.
#     rξ = l*rand(N)
#     ps = Vector{HyperBKParticle{T}}(undef,N)
#     intervals = harcintervals(bd)
#     Threads.@threads for i in 1:N
#         p, v, _ = hfrom_bcoords(T(rξ[i]),T(rφ[i]),bd,intervals)
#         ps[i] = HyperBKParticle(p,v)
#     end
#     return ps
# end

function correlation(φ,ψ, μ::Measure, N::Int, n::Int)
    ps = μ.sampler(N)
    d  = μ.density.(ps)
    sd = sum(d)
    rψ = ψ.(ps) .* d
    ∫ψ = sum(rψ) / sd

    Rφ = MVector{n+1,Float64} |> zero
    Threads.@threads for i in 1:N
        p = ps[i]
        rφ = φ(p)
        Rφ[1] += rφ*rψ[i] - rφ*∫ψ*d[i]
        for j in 2:n+1
            bounce!(p, μ.bd)
            Rφ[j] += φ(p)*rψ[i] - rφ*∫ψ*d[i]
        end
    end
    
    Rφ ./ sd
end