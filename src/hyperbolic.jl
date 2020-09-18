
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

@inline hyper_rtransf(s²) = 1 / (1+sqrt(1 - s²))

@inline function hyper_transf(x,y)
    s² = x*x + y*y
    k = hyper_rtransf(s²)
    return k*x, k*y
end

function PyPlot.plot(w::Wall, ::Type{HyperBKParticle}; kwargs...)
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
        color = DynamicalBilliards.obcolor(w),
        linestyle = DynamicalBilliards.obls(w), lw = 2.0, kwargs...)
    end
end

function PyPlot.plot(d::Semicircle, ::Type{HyperBKParticle}; kwargs...)
    theta1 = atan(d.facedir[2], d.facedir[1])*180/π + 90
    theta2 = theta1 + 180
    edgecolor = DynamicalBilliards.obcolor(d)
    s1 = PyPlot.matplotlib.patches.Arc(d.c, 2d.r*hyper_rtransf(d.r^2), 2d.r*hyper_rtransf(d.r^2),
    theta1 = theta1, theta2 = theta2, edgecolor = edgecolor, lw = 2.0, kwargs...)
    PyPlot.gca().add_artist(s1)
end

PyPlot.plot(d::Circular, ::Type{HyperBKParticle}; kwargs...) = 
    _oc_plot(hyper_rtransf, d; kwargs...)


PyPlot.plot(ps::AbstractVector{<:AbstractParticle}, ::Type{HyperBKParticle}, colrs; ax=PyPlot.gca()) =
    _p_plot(hyper_transf, ps, colrs, ax)

DynamicalBilliards.timeseries!(p::HyperBKParticle{T}, bd::Billiard{T}, f) where T =
    timeseries!(p, bd, f; transf = hyper_transf)
