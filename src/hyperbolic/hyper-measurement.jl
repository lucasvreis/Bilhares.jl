
using DynamicalBilliards: cross2D

export hfrom_bcoords, htotallength

# @inline function inter_line_circle(a,b)
#     d² = normsq(a-b)
#     a² = normsq(a)
#     m  = a[1]b[1] + a[2]b[2]
#     Δ  = √(d² - (cross2D(a,b))^2)
#     t₁ = (a² - m - Δ)/d²
#     t₂ = (a² - m + Δ)/d²
#     (1-t₁)*a + t₁*b, (1-t₂)*a + t₂*b
# end
# @inline function hdistance(u,v)
#     a, b = inter_line_circle(u,v)
#     1/2 * (log(norm(v-a)*norm(b-u)) - log(norm(u-a)*norm(b-v)))
# end

o_hdistance(r) = 1 / 2 * (log(abs(1 + r)) - log(abs(1 - r)))
hdistance(u, v) = ManifoldsBase.distance(ℍ, BKtoHyperboloid(u), BKtoHyperboloid(v))


htotallength(o::Wall) = hdistance(o.ep, o.sp)
htotallength(o::Semicircle) = π * sinh(o_hdistance(o.r))
htotallength(o::Circular) = 2π * sinh(o_hdistance(o.r))

htotallength(bd::Billiard) = sum(htotallength(x) for x in bd.obstacles)

function harcintervals(bd::Billiard{T,D}) where {T,D}
    intervals = SVector{D + 1,T}(0, map(x -> htotallength(x), bd.obstacles)...)
    return cumsum(intervals)
end
#=
to_bcoords(p::AbstractParticle, o::Obstacle) = to_bcoords(p.pos, p.vel, o)

function to_bcoords(pos::SV, vel::SV, o::Obstacle)
    n = normalvec(o, pos)
    sinφ = cross2D(vel, n)
    ξ = _ξ(pos, o)
    return ξ, sinφ
end

_ξ(pos::SV, o::Wall) = norm(pos - o.sp)

function _ξ(pos::SV{T}, o::Semicircle{T}) where {T <: AbstractFloat}
    # project pos on open face
    chrd = SV{T}(-o.facedir[2], o.facedir[1])
    d = (pos - o.c) / o.r
    x = dot(d, chrd)
    r =  acos(clamp(x, -1, 1)) * o.r
    return r
end

function _ξ(pos::SV{T}, o::Circular{T}) where {T <: AbstractFloat}
    d = (pos - o.c) / o.r
    if d[2] > 0
        r = acos(clamp(d[1], -1, 1))
    else
        r = π + acos(-clamp(d[1], -1, 1))
    end
    return r * o.r
end
=#

function hfrom_bcoords(ξ::T, sφ::T, o::Obstacle{T}) where {T}
    z = zero
    pos = hreal_pos(ξ, o)
    hpos = BKtoHyperboloid(pos)
    cφ = sqrt((1 - sφ) * (1 + sφ)) # = cos(asin(sφ))
    cφ′, sφ′ = D_h_BK(hpos) * vector_transport_to(
        ℍ, SVector(0.,0.,1.), SVector(cφ, sφ, 0.), hpos) |> normalize
    n = normalvec(o, pos)
    vel = SV{T}(n[1] * cφ′ + n[2] * sφ′, -n[1] * sφ′ + n[2] * cφ′)
    return pos, vel
end

function hreal_pos(ξ, o::Wall) 
    spH = BKtoHyperboloid(o.sp)
    epH = BKtoHyperboloid(o.ep)
    l   = log(ℍ, spH, epH)
    HyperboloidtoBK(exp(ℍ, spH, l * ξ / ManifoldsBase.distance(ℍ, spH, epH)))
end

function hreal_pos(ξ, o::Semicircle{T}) where T
    sξ, cξ = sincos(ξ / sinh(o_hdistance(o.r)))
    chrd   = SV{T}(-o.facedir[2], o.facedir[1])
    return -o.r * (sξ * o.facedir - cξ * chrd)
end

hreal_pos(ξ, o::Circular{T}) where T = o.r * SV{T}(cossin(ξ / sinh(o_hdistance(o.r))))

function hfrom_bcoords(ξ, sφ, bd::Billiard{T}, intervals=harcintervals(bd)) where T

    for (i, obst) ∈ enumerate(bd)
        if ξ <= intervals[i + 1]
            pos, vel  = hfrom_bcoords(ξ - intervals[i], sφ, obst)
            return pos, vel, i
        end
    end
    throw(DomainError(ξ, "ξ is too large for this billiard!"))
end
