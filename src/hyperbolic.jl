import DynamicalBilliards: normalvec, specular!, collision, extrapolate, cossin,
                           nocollision, ispinned, accuracy, plot, obcolor, obls

export HyperBKParticle, hbilliard_sq

"""
        HyperbolicParticle <: AbstractParticle

Partícula no modelo hiperbólico de Beltrami-Klein.
"""
mutable struct HyperBKParticle{T<:AbstractFloat} <: AbstractParticle{T}
    pos::SVector{2,T}
    vel::SVector{2,T}
    function HyperBKParticle(
        pos::SVector{2,T}, vel::SVector{2,T}) where{T<:AbstractFloat}
        new{T}(pos, normalize(vel))
    end
end

# HyperBKParticles propagate linearly
@traitimpl PropagatesLinearly{HyperBKParticle}

# How to copy HyperBKParticles
Base.copy(p::HyperBKParticle) = HyperBKParticle(p.pos, p.vel)

function HyperBKParticle(ic::AbstractVector{S}) where {S<:Real}
    T = S<:Integer ? Float64 : S
    φ0 = ic[3]
    pos = SVector{2,T}(ic[1:2]); vel = SVector{2,T}(cossin(φ0)...)
    return HyperBKParticle(pos, vel)
end

HyperBKParticle(x::Real, y::Real, φ::Real) = HyperBKParticle(collect(promote(x,y,φ)))

function HyperBKParticle(pos::SV{T}, vel::SV{T}) where {T}
    S = T<:Integer ? Float64 : T
    return HyperBKParticle(pos, vel)
end
Base.show(io::IO, p::HyperBKParticle{T}) where {T} = print(io, 
        "HyperbolicParticle{$T}\n", 
        "position: $(p.pos)\nvelocity: $(p.vel)"
)


# Reflection -- the elephant trouble

@inline function hyperreflect(x,n,i)
    sq1x  = 1 + sqrt(1 - normsq(x))
    PTK = SM(
        sq1x-x[1]^2,  -x[1]x[2],
        -x[1]x[2],  sq1x-x[2]^2
    )
    KTP = inv(PTK)
    t′ = KTP * SV(-n[2], n[1])
    i′ = KTP * i
    o′ =  2 * t′⋅i′ * t′ / normsq(t′) - i′ 
    normalize(PTK * o′)
end


@inline function specular!(p::HyperBKParticle{T}, o::Obstacle{T})::Nothing where {T}
    n = normalvec(o, p.pos)
    p.vel = hyperreflect(p.pos, n, p.vel)
    return nothing
end

extrapolate(p::HyperBKParticle, a, b, c, d, ω) = extrapolate(p::HyperBKParticle, a, b, c, d)
function extrapolate(p::HyperBKParticle{T}, prevpos::SV{T}, prevvel::SV{T}, ct::T,
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
    x = big(x); y = big(y)
    s² = x*x + y*y
    k = typeof(x)(1/(1+sqrt(1 - s²)))
    return k*x, k*y
    # return x,y
end

# NOT WORKING!!
function plot(w::Wall; kwargs...)
    if typeof(w) <: FiniteWall &&  w.isdoor
       PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]];
       color="black", linestyle = "-", lw = 2.0, kwargs...)
       PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]];
       color=(0, 0.9, 0.9), linestyle = "--", lw = 2.0, kwargs...)
    else
        x = range(w.sp[1],w.ep[1],length=50)
        y = range(w.sp[2],w.ep[2],length=50)
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



