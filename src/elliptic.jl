
export EllipticParticle, ebilliard_sq

"""
        HyperbolicParticle <: AbstractParticle

Partícula no modelo hiperbólico de Beltrami-Klein.
"""
mutable struct EllipticParticle{T<:AbstractFloat} <: AbstractParticle{T}
    pos::SVector{2,T}
    vel::SVector{2,T}
    function EllipticParticle(
        pos::SVector{2,T}, vel::SVector{2,T}) where{T<:AbstractFloat}
        new{T}(pos, normalize(vel))
    end
end

# EllipticParticles propagate linearly
@traitimpl PropagatesLinearly{EllipticParticle}

# How to copy EllipticParticles
Base.copy(p::EllipticParticle) = EllipticParticle(p.pos, p.vel)

function EllipticParticle(ic::AbstractVector{S}) where {S<:AbstractFloat}
    T = S<:Integer ? Float64 : S
    φ0 = ic[3]
    pos = SVector{2,T}(ic[1:2]); vel = SVector{2,T}(cossin(φ0)...)
    return EllipticParticle(pos, vel)
end

EllipticParticle(x::AbstractFloat, y::AbstractFloat, φ::AbstractFloat) = EllipticParticle(collect(promote(x,y,φ)))


Base.show(io::IO, p::EllipticParticle{T}) where {T} = print(io, 
        "HyperbolicParticle{$T}\n", 
        "position: $(p.pos)\nvelocity: $(p.vel)"
)


# Reflection -- the elephant trouble

@inline function ellreflect(x,n,i)
    sqx2um = sqrt(normsq(x)+1)
    denom = sqx2um*(sqx2um+1)^2
    mx = x[1]x[2]
    GTS = 2 * SM(
        sqx2um+x[2]^2+1,  -mx,
        -mx,  sqx2um+x[1]^2+1
    ) ./ denom
    STG = inv(GTS)
    t′ = GTS * SV(-n[2], n[1])
    i′ = GTS * i
    o′ =  2 * t′⋅i′ * t′ ./ normsq(t′) - i′ 
    normalize(STG * o′)
end

@inline function specular!(p::EllipticParticle{T}, o::Obstacle{T})::Nothing where {T}
    n = normalvec(o, p.pos)
    p.vel = ellreflect(p.pos, n, p.vel)
    return nothing
end

@inline ell_rtransf(s²) = 2 / (1 + sqrt(1 + s²))

@inline function ell_transf(x,y)
    s² = x*x + y*y
    k = ell_rtransf(s²)
    return k*x, k*y
end

function PyPlot.plot(w::Wall, ::Type{EllipticParticle}; kwargs...)
    if typeof(w) <: FiniteWall &&  w.isdoor
       PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]];
       color="black", linestyle = "-", lw = 2.0, kwargs...)
       PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]];
       color=(0, 0.9, 0.9), linestyle = "--", lw = 2.0, kwargs...)
    else
        x = range(w.sp[1],w.ep[1],length=50)
        y = range(w.sp[2],w.ep[2],length=50)
        p = (ell_transf(a,b) for (a,b) in zip(x,y))
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

function PyPlot.plot(d::Semicircle, ::Type{EllipticParticle}; kwargs...)
    theta1 = atan(d.facedir[2], d.facedir[1])*180/π + 90
    theta2 = theta1 + 180
    edgecolor = DynamicalBilliards.obcolor(d)
    s1 = PyPlot.matplotlib.patches.Arc(d.c, 2d.r*ell_rtransf(d.r^2), 2d.r*ell_rtransf(d.r^2),
    theta1 = theta1, theta2 = theta2, edgecolor = edgecolor, lw = 2.0, kwargs...)
    PyPlot.gca().add_artist(s1)
end

PyPlot.plot(ps::AbstractVector{<:AbstractParticle}, ::Type{EllipticParticle}, colrs; ax=PyPlot.gca()) =
    _p_plot(ell_transf, ps, colrs, ax)

DynamicalBilliards.timeseries!(p::EllipticParticle{T}, bd::Billiard{T}, f) where T =
    timeseries!(p, bd, f; transf = ell_transf)
