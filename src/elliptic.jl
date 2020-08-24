import DynamicalBilliards: normalvec, specular!, collision, extrapolate, cossin,
                           nocollision, ispinned, accuracy, plot, obcolor, obls

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

function EllipticParticle(ic::AbstractVector{S}) where {S<:Real}
    T = S<:Integer ? Float64 : S
    φ0 = ic[3]
    pos = SVector{2,T}(ic[1:2]); vel = SVector{2,T}(cossin(φ0)...)
    return EllipticParticle(pos, vel)
end

EllipticParticle(x::Real, y::Real, φ::Real) = EllipticParticle(collect(promote(x,y,φ)))

function EllipticParticle(pos::SV{T}, vel::SV{T}) where {T}
    S = T<:Integer ? Float64 : T
    return EllipticParticle(pos, vel)
end
Base.show(io::IO, p::EllipticParticle{T}) where {T} = print(io, 
        "HyperbolicParticle{$T}\n", 
        "position: $(p.pos)\nvelocity: $(p.vel)"
)


# Reflection -- the elephant trouble

@inline function hyperreflect(x,n,i)
    mx = x[1]x[2]
    dx = x[1]^2 - x[2]^2
    denom = (normsq(x) - 4)^2
    STG = 4 * SM(
        4 + dx,  2mx,
        2mx,  4 - dx
    ) ./ denom
    GTS = inv(STG)
    t′ = GTS * SV(-n[2], n[1])
    i′ = GTS * i
    o′ =  2 * t′⋅i′ * t′ / normsq(t′) - i′ 
    normalize(STG * o′)
end


@inline function specular!(p::EllipticParticle{T}, o::Obstacle{T})::Nothing where {T}
    n = normalvec(o, p.pos)
    p.vel = hyperreflect(p.pos, n, p.vel)
    return nothing
end

extrapolate(p::EllipticParticle, a, b, c, d, ω) = extrapolate(p::EllipticParticle, a, b, c, d)
function extrapolate(p::EllipticParticle{T}, prevpos::SV{T}, prevvel::SV{T}, ct::T,
                     dt::T) where {T <: AbstractFloat}

    tvec = collect(0:dt:ct)
    x = Vector{T}(undef, length(tvec))
    y = Vector{T}(undef, length(tvec))
    vx = Vector{T}(undef, length(tvec))
    vy = Vector{T}(undef, length(tvec))

    @inbounds for (i,t) ∈ enumerate(tvec)
        px = prevpos[1] + t*prevvel[1]
        py = prevpos[2] + t*prevvel[2]
        x′, y′ = ellip_transf(px,py)
        x[i] = x′
        y[i] = y′
        vx[i], vy[i] = prevvel
    end

    # finish with ct
    @inbounds if tvec[end] != ct
        push!(tvec, ct)
        x′, y′ = ellip_transf(p.pos...)
        push!(x, x′)
        push!(y, y′)
        push!(vx, p.vel[1]); push!(vy, p.vel[2])
    end

    return x, y, vx, vy, tvec
end

function ellip_transf(x,y)
    x = big(x); y = big(y)
    s² = x*x + y*y
    k = typeof(x)(2 / (1+sqrt(1 + s²)))
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
        p = (ellip_transf(a,b) for (a,b) in zip(x,y))
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

function ebilliard_sq(s)
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


function DynamicalBilliards.timeseries!(p::EllipticParticle{T}, bd::Billiard{T}, f, raysplitters = nothing;
                     dt = typeof(p) <: Particle ? T(Inf) : T(0.01),
                     warning::Bool = true) where {T}

    ts = [zero(T)]
    x  = [p.pos[1]]; y  = [p.pos[2]]
    vx = [p.vel[1]]; vy = [p.vel[2]]
    ismagnetic = p isa MagneticParticle
    prevω = ismagnetic ? p.ω : T(0)
    isray = !isa(raysplitters, Nothing)
    isray && acceptable_raysplitter(raysplitters, bd)
    raysidx = raysplit_indices(bd, raysplitters)
    n, i, t, flight = 0, 0, zero(T), zero(T)
    prevpos = p.pos
    prevvel = p.vel

    @inbounds while check_condition(f, n, t, i, p) # count < t
        i, ct = bounce!(p, bd, raysidx, raysplitters)
        flight += ct
        if flight ≤ dt # push collision point only
            push!(ts, flight + t)
            push!(x, p.pos[1])
            push!(y, p.pos[2])
            push!(vx, p.vel[1]); push!(vy, p.vel[2])
        else
            nx, ny, nvx, nvy, nts = extrapolate(p, prevpos, prevvel, flight, dt, prevω)
            append!(ts, nts[2:end] .+ t)
            append!(x, nx[2:end]); append!(vx, nvx[2:end])
            append!(y, ny[2:end]); append!(vy, nvy[2:end])
        end
        prevpos = p.pos
        prevvel = p.vel
        t += flight; n += 1
        flight = zero(T)
        ismagnetic && isray && (prevω = p.ω)
    end
    return x, y, vx, vy, ts
end

