export billiard_csq, billiard_cmushroom, billiard_crectangle, billiard_csinai,
       random_on_border

using DynamicalBilliards: cellsize

include("plotting.jl")
include("billiards.jl")
include("gui.jl")

@inline normsq(v) = v[1]v[1] + v[2]v[2]

function DynamicalBilliards.randominside(T::Type{<:AbstractParticle}, bd::Billiard, N::Int) 
    [T(randominside_xyφ(bd)...) for _ in 1:N]
end

function random_on_border(P::Type{<:AbstractParticle}, bd::Billiard{T}, N::Int) where T
    ξ = totallength(bd)
    rξ = rand(N)*ξ
    rϕ = 2*rand(N).-1.
    [   begin
        p, v, _ = from_bcoords(T(rξ[i]),T(rϕ[i]),bd)
        P(p,v)
    end for i in 1:N]
end

function DynamicalBilliards.extrapolate(p::AbstractParticle{T}, prevpos::SV{T}, 
        prevvel::SV{T}, ct::T, dt::T, transf::Function) where {T <: AbstractFloat}

    tvec = collect(0:dt:ct)
    x = Vector{T}(undef, length(tvec))
    y = Vector{T}(undef, length(tvec))
    vx = Vector{T}(undef, length(tvec))
    vy = Vector{T}(undef, length(tvec))

    @inbounds for (i,t) ∈ enumerate(tvec)
        px = prevpos[1] + t*prevvel[1]
        py = prevpos[2] + t*prevvel[2]
        x′, y′ = transf(px,py)
        x[i] = x′
        y[i] = y′
        vx[i], vy[i] = prevvel
    end

    # finish with ct
    @inbounds if tvec[end] != ct
        push!(tvec, ct)
        x′, y′ = transf(p.pos...)
        push!(x, x′)
        push!(y, y′)
        push!(vx, p.vel[1]); push!(vy, p.vel[2])
    end

    return x, y, vx, vy, tvec
end

function DynamicalBilliards.timeseries!(p::AbstractParticle{T}, bd::Billiard{T}, f;
                     transf = (x...) -> x,
                     dt = typeof(p) <: Particle ? T(Inf) : T(0.01)) where T
    ts = [zero(T)]
    p0 = transf(p.pos...)
    x  = [p0[1]]; y  = [p0[2]]
    vx = [p.vel[1]]; vy = [p.vel[2]]
    n, i, t = 0, 0, zero(T)
    prevpos = p.pos
    prevvel = p.vel

    @inbounds while check_condition(f, n, t, i, p) # count < t
        i, ct = bounce!(p, bd)
        if ct ≤ dt # push collision point only
            push!(ts, ct + t)
            pp′ = transf(p.pos...)
            push!(x, pp′[1])
            push!(y, pp′[2])
            push!(vx, p.vel[1]); push!(vy, p.vel[2])
        else
            nx, ny, nvx, nvy, nts = extrapolate(p, prevpos, prevvel, ct, dt, transf)
            append!(ts, nts[2:end] .+ t)
            append!(x, nx[2:end]); append!(vx, nvx[2:end])
            append!(y, ny[2:end]); append!(vy, nvy[2:end])
        end
        prevpos = p.pos
        prevvel = p.vel
        t += ct; n += 1
    end
    return x, y, vx, vy, ts
end


