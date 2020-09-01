export billiard_csq, billiard_cmushroom, billiard_rpolygon, random_on_border

using DynamicalBilliards: cellsize

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


function extrapolate(p::AbstractParticle{T}, prevpos::SV{T}, prevvel::SV{T}, ct::T,
                     dt::T, transf::Function) where {T <: AbstractFloat}

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

function PyPlot.plot(bd::Billiard{T}, P::Type{<:AbstractParticle};
    ax = (PyPlot.figure(); PyPlot.gca())) where {T}
    PyPlot.sca(ax)
    for obst in bd; plot(obst, P); end
    # xmin, ymin, xmax, ymax = cellsize(bd)
    # dx = xmax - xmin; dy = ymax - ymin
    # ax.set_aspect("equal")
    # if !isinf(xmin) && !isinf(xmax)
    #     PyPlot.xlim(xmin - 0.1dx, xmax + 0.1dx)
    # end
    # if !isinf(ymin) && !isinf(ymax)
    #     PyPlot.ylim(ymin - 0.1dy, ymax + 0.1dy)
    # end
    return nothing
end

PyPlot.plot(o::Obstacle, ::Type{Particle}; kwargs...) = plot(o; kwargs...)

function _p_plot(transf, ps, colrs, ax)
    x = []; y = []
    for p in ps
        x′, y′ = transf(p.pos...)
        push!(x,x′)
        push!(y,y′)
    end
    ax.scatter(x, y; color=colrs, s=20.0, zorder=99)
end


PyPlot.plot(ps::AbstractVector{<:AbstractParticle}, ::Type{Particle}, colrs; kwargs...) = 
    plot(ps, colrs; kwargs...)
PyPlot.plot(ps::AbstractVector{<:AbstractParticle}, colrs; ax=PyPlot.gca()) =
    _p_plot((x...) -> x, ps, colrs, ax)

function billiard_csq(s)
    o = typeof(s)(0)
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

function billiard_cmushroom(stem_length = 1.0, stem_width=0.2, cap_radius=1.0,
    stem_location = 0.0; sl = stem_length, sw = stem_width, cr = cap_radius,
    sloc = stem_location, door = false, scale=1.0)

    abs(sloc) + sw/2 > cr && error("Stem is outside the mushroom cap!")

    sl, sw, cr, sloc = scale .* promote(sl, sw, cr, sloc)
    T = eltype(sl)
    leftcorn = SV(-sw/2 + sloc, -sl)
    rightcorn = SV(sw/2 + sloc, -sl)
    upleftcorn = SV(-sw/2 + sloc, 0)
    uprightcorn = SV(sw/2 + sloc, 0)

    stembot = FiniteWall(leftcorn, rightcorn, SV(0, sw), door, "Stem bottom")
    stemleft = FiniteWall(upleftcorn, leftcorn, SV(sw, 0), false, "Stem left")
    stemright = FiniteWall(rightcorn, uprightcorn, SV(-sw, 0), false, "Stem right")

    farleft = SV(-cr, 0)
    farright = SV(cr, 0)

    capbotleft = FiniteWall(
    farleft, upleftcorn, SV(0, sw), false, "Cap bottom left")
    capbotright = FiniteWall(
    uprightcorn, farright, SV(0, sw), false, "Cap bottom right")

    cap = Semicircle([0, 0], cr, [0, -T(1.0)], "Mushroom cap")
    return Billiard(stembot, stemright, capbotright, cap, capbotleft, stemleft)
end