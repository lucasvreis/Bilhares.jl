__precompile__()
module BilharesMarkov

export prob_estacionária, mat_estocástica, map_markov, média_temporal, HyperbolicParticle
export hbilliard_sq

using DynamicalBilliards
using StaticArrays, LinearAlgebra, Arpack, SparseArrays
using ProgressMeter, PyPlot
const SV = SVector{2}
const SM = SMatrix{2,2}

"""
    prob_estacionária(g::SparseMatrixCSC)
Probabilidade estacionária associada à matriz estocástica esparsa `g`,
calculada como um autovetor associado ao autovalor `1` de `g'`, cujas
entradas somam `1`. Caso a matriz não seja irredutível, a probabilidade
pode não ser única. 
"""
prob_estacionária(g::SparseMatrixCSC) = normalize(eigs(g';nev=1,which=:LR)[2][:,1],1) .|> abs


"""
    mat_estocástica(maps::AbstractMatrix)

Retorna uma cópia da matriz `A` com as linhas normalizadas de forma
que as somas de cada linha vale `1`, ou `0` se a linha inteira for nula.
"""
function mat_estocástica(A::AbstractMatrix) 
    P = copy(A)
    for m in eachrow(P)
        if sum(m) != 0
            m ./= sum(m)
        else
            m .*= 0
        end
    end
    P
end


"""
       map_markov(bd::Billiard{T}, Nbounces::Int,
                  np::Int, Nξ, Nφ=Nξ;
                  intervals=arcintervals(bd)) where {T}

Estima um processo de Markov gerado pela partição do bilhar `bd` em
`Nξ` comprimentos de arco iguais e `Nφ` intervalos de ângulos iguais,
resultando em `Nξ*Nφ` estados. As probabilidades de transição são estimadas
simulando `np` partículas cada uma por `Nbounces` reflexões no bordo.
"""
function map_markov(bd::Billiard{T}, Nbounces::Int,
                           np::Int, Nξ, Nφ=Nξ;
                           intervals=arcintervals(bd)) where {T}
    δξ = totallength(bd) / Nξ # Tamanho dos compartimentos-arco
    δφ = 2 / Nφ               # Tamanho dos compartimentos-ângulo

    # As partículas aleatórias
    ps = [randominside(bd) for _ in 1:np] 

    # Vetor para armazenar os dados de cada partícula em paralelo.
    pP = [spzeros(Int, Nξ * Nφ, Nξ * Nφ) for i in 1:np]

    # Transforma (arco,ângulo) em índice.
    vectoind(ξ,sφ) = floor(Int, (sφ + 1) / δφ) * Nξ +
                     floor(Int, ξ / δξ) + 1

    prog = Progress(np) # Barra de progresso

    # Itera pelas partículas em paralelo
    Threads.@threads for ip in 1:np
        # Partícula atual
        p = ps[ip]
        # Índice da reflexão anterior
        lastind = nothing
        for N in 1:Nbounces
        i, = bounce!(p,bd) # Reflete a partícula na próxima borda

        ξ, sφ = to_bcoords(p, bd[i]) # Coordenadas
        ξ += intervals[i] # Transforma em coordenadas globais

        ind = vectoind(ξ,sφ) # Índice
        
        if N > 1
            # Registra que a partícula `ip` saiu da região `lastind`
            # e foi para a região `ind` mais uma vez.
            pP[ip][lastind,ind] += 1
        end
        lastind = ind
        end
        next!(prog) # atualiza a barra de progresso
  end

  mat_estocástica(sum(pP)./1) # transforma em matriz estocástica
end

function média_temporal(bd::Billiard{T}, Nbounces::Int,
                       np::Int, Nξ, Nφ=Nξ;
                       intervals=arcintervals(bd)) where {T}
    δξ = totallength(bd) / Nξ # Tamanho dos compartimentos-arco
    δφ = 2 / Nφ               # Tamanho dos compartimentos-ângulo

    # As partículas aleatórias
    ps = [randominside(bd) for _ in 1:np] 

    # Vetor para armazenar os dados de cada partícula em paralelo.
    pP = [spzeros(Int, Nξ * Nφ) for _ in 1:np]
    firsti = zeros(Int, np)

    # Transforma (arco,ângulo) em índice.
    vectoind(ξ,sφ) = floor(Int, (sφ + 1) / δφ) * Nξ +
                     floor(Int, ξ / δξ) + 1

    prog = Progress(np) # Barra de progresso

    # Itera pelas partículas em paralelo
    Threads.@threads for ip in 1:np
        # Partícula atual
        p = ps[ip]
        for N in 1:Nbounces
            i, = bounce!(p,bd) # Reflete a partícula na próxima borda

            ξ, sφ = to_bcoords(p, bd[i]) # Coordenadas
            ξ += intervals[i] # Transforma em coordenadas globais

            ind = vectoind(ξ,sφ) # Índice
            
            pP[ip][ind] += 1

            if N == 1
                firsti[ip] = ind
            end
        end
        next!(prog) # atualiza a barra de progresso
    end
    
    tret = zeros(Nξ * Nφ)
    mret = zeros(Int, Nξ * Nφ)
    for (p,fi) in zip(pP,firsti)
        tret[fi] += Nbounces/p[fi]
        mret[fi] += 1
    end

    vcat(pP'...)./Nbounces, [mret[i] ≠ 0 ? tret[i]/mret[i] : 0 for i in 1:(Nξ * Nφ)] # transforma em matriz estocástica
end

normsq(v) = (x = zero(eltype(v)); for e in v x += e^2 end; x)

function hyperreflect(x,n,i)
    T = eltype(x)
    xsq = normsq(x)
    x⊥ = x[1]^2 - x[2]^2
    denom = (xsq + 1)^2
    J = SM{T}(2(1-x⊥),  -4x[1]x[2],
              -4x[1]x[2], 2(1+x⊥)) ./ denom
    J * (i - 2*n⋅i*n)
end

import DynamicalBilliards: normalvec, specular!, collision, extrapolate, cossin, nocollision,
                           resolvecollision!, propagate!, ispinned, accuracy, plot,
                           obcolor, obls

mutable struct HyperbolicParticle{T<:AbstractFloat} <: AbstractParticle{T}
    pos::SVector{2,T}
    vel::SVector{2,T}
    current_cell::SVector{2,Int}
    function HyperbolicParticle(
        pos::SVector{2,T}, vel::SVector{2,T}) where{T<:AbstractFloat}
        new{T}(pos, normalize(vel), SV{Int}(0,0))
    end
end
Base.copy(p::HyperbolicParticle) = HyperbolicParticle(p.pos, p.vel)

function HyperbolicParticle(ic::AbstractVector{S}) where {S<:Real}
    T = S<:Integer ? Float64 : S
    φ0 = ic[3]
    pos = SVector{2,T}(ic[1:2]); vel = SVector{2,T}(cossin(φ0)...)
    return HyperbolicParticle(pos, vel)
end
HyperbolicParticle(x::Real, y::Real, φ::Real) = HyperbolicParticle(collect(promote(x,y,φ)))
# HyperbolicParticle() = HyperbolicParticle(rand(), rand(), rand()*2π)
function HyperbolicParticle(pos::SV{T}, vel::SV{T}) where {T}
    S = T<:Integer ? Float64 : T
    return HyperbolicParticle(pos, vel)
end
show(io::IO, p::HyperbolicParticle{T}) where {T} =
print(io, "HyperbolicParticle{$T}\n",
"position: $(p.pos)\nvelocity: $(p.vel)")

function collision(p::HyperbolicParticle{T}, w::Wall{T}) where {T}
    n = normalvec(w, p.pos)
    denom = dot(p.vel, n)
    if denom ≥ 0.0
        return nocollision(T)
    else
        t = dot(w.sp - p.pos, n)/denom #??
        return t, p.pos + t*p.vel
    end
end
ispinned(p::HyperbolicParticle,bd) = false
@inline function specular!(p::HyperbolicParticle{T}, o::Obstacle{T})::Nothing where {T}
    n = normalvec(o, p.pos)
    p.vel = hyperreflect(p.pos,n,-p.vel)
    return nothing
end
@inline resolvecollision!(p::HyperbolicParticle, o::Obstacle) = specular!(p, o)
@inline function propagate!(p::HyperbolicParticle{T}, t::Real) where {T}
    p.pos += SV{T}(p.vel[1]*t, p.vel[2]*t)
end
@inline propagate!(p::HyperbolicParticle, newpos::SV, t::Real) = (p.pos = newpos)
@inline resolvecollision!(p::HyperbolicParticle, o::Obstacle) = specular!(p, o)

extrapolate(p::HyperbolicParticle, a, b, c, d, ω) = extrapolate(p::HyperbolicParticle, a, b, c, d)
function extrapolate(p::HyperbolicParticle{T}, prevpos::SV{T}, prevvel::SV{T}, ct::T,
                     dt::T) where {T <: AbstractFloat}

    tvec = collect(0:dt:ct)
    x = Vector{T}(undef, length(tvec))
    y = Vector{T}(undef, length(tvec))
    vx = Vector{T}(undef, length(tvec))
    vy = Vector{T}(undef, length(tvec))

    @inbounds for (i,t) ∈ enumerate(tvec)
        px = prevpos[1] + t*prevvel[1]
        py = prevpos[2] + t*prevvel[2]
        x′, y′ = hyper_transf(px,py)
        x[i] = x′
        y[i] = y′
        vx[i], vy[i] = prevvel
    end

    # finish with ct
    @inbounds if tvec[end] != ct
        push!(tvec, ct)
        x′, y′ = hyper_transf(p.pos...)
        push!(x, x′)
        push!(y, y′)
        push!(vx, p.vel[1]); push!(vy, p.vel[2])
    end

    return x, y, vx, vy, tvec
end

function hyper_transf(x,y)
    s² = x*x + y*y
    k = 1/(1+sqrt(1 - s²))
    return k*x, k*y
    # return x,y
end

function plot(w::Wall; kwargs...)
    if typeof(w) <: FiniteWall &&  w.isdoor
       PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]];
       color="black", linestyle = "-", lw = 2.0, kwargs...)
       PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]];
       color=(0, 0.9, 0.9), linestyle = "--", lw = 2.0, kwargs...)
    else
        x = range(w.sp[1],w.ep[1],length=100)
        y = range(w.sp[2],w.ep[2],length=100)
        p = (hyper_transf(a,b) for (a,b) in zip(x,y))
        x′, y′ = collect.(zip(p...))
        # circle1 = PyPlot.plt.Circle([0,0], 1;
        # edgecolor = "black", facecolor = "none",
        # linestyle = "-", lw = 2.0, kwargs...)
        # PyPlot.gca().add_artist(circle1)
        PyPlot.plot(x′,y′;
        color=obcolor(w),
        linestyle = obls(w), lw = 2.0, kwargs...)
    end
end

function hbilliard_sq(s)
    s = convert(AbstractFloat, s)
    o = typeof(s)(0.0)
    sp = [-s,s]; ep = [-s, -s]; n = [s,o]
    leftw = InfiniteWall(sp, ep, n, "Left wall")
    sp = [s,-s]; ep = [s, s]; n = [-s,o]
    rightw = InfiniteWall(sp, ep, n, "Right wall")
    sp = [s,s]; ep = [-s, s]; n = [o,-s]
    topw = InfiniteWall(sp, ep, n, "Top wall")
    sp = [-s,-s]; ep = [s, -s]; n = [o,s]
    botw = InfiniteWall(sp, ep, n, "Bottom wall")
    return Billiard(botw, rightw, topw, leftw)
end

function collision(p::HyperbolicParticle{T}, e::Ellipse{T}) where {T}
    # First check if particle is "looking at" eclipse if it is outside
    if e.pflag
        # These lines may be "not accurate enough" but so far all is good
        dotp = dot(p.vel, normalvec(e, p.pos))
        dotp ≥ 0.0 && return nocollision(T)
    end

    # http://www.ambrsoft.com/TrigoCalc/Circles2/Ellipse/EllipseLine.htm
    a = e.a; b = e.b
    # Translate particle with ellipse center (so that ellipse lies on [0, 0])
    pc = p.pos - e.c
    # Find μ, ψ for line equation y = μx + ψ describing particle
    μ = p.vel[2]/p.vel[1]
    ψ = pc[2] - μ*pc[1]

    # Determinant and intersection points follow from the link
    denomin = a*a*μ*μ + b*b
    Δ² = denomin - ψ*ψ
    Δ² ≤ 0 && return nocollision(T)
    Δ = sqrt(Δ²); f1 = -a*a*μ*ψ; f2 = b*b*ψ # just factors
    I1 = SV(f1 + a*b*Δ, f2 + a*b*μ*Δ)/denomin
    I2 = SV(f1 - a*b*Δ, f2 - a*b*μ*Δ)/denomin

    d1 = norm(pc - I1); d2 = norm(pc - I2)
    if e.pflag
        return d1 < d2 ? (d1, I1 + e.c) : (d2, I2 + e.c)
    else # inside the ellipse: one collision is _always_ valid
        if d1 < d2
            dmin, Imin = d1, I1
            dmax, Imax = d2, I2
        else
            dmin, Imin = d2, I2
            dmax, Imax = d1, I1
        end

        if dmin < accuracy(T) # special case for being very close to ellipse
            dotp = dot(p.vel, normalvec(e, Imin))
            dotp ≥ 0 && return (dmax, Imax + e.c)
        end
         # check which of the two points is ahead or behind the obstacle
        z1 = dot(pc - Imax, p.vel)
        return z1 < 0 ? (dmax, Imax + e.c) : (dmin, Imin + e.c)
    end
end


end
