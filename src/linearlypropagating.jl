using SimpleTraits

import DynamicalBilliards: collision, resolvecollision!, propagate!, Obstacle,
                            distance

export PropagatesLinearly

"""
        PropagatesLinearly{P}

Trait for linearly propagating particles. Must define `specular!`, [TODO]...
"""
@traitdef PropagatesLinearly{P<:AbstractParticle}

@traitfn function collision(p::P, w::Wall{T}) where {T,P<:AbstractParticle{T}; 
                                                        PropagatesLinearly{P}}
    n = normalvec(w, p.pos)
    denom = dot(p.vel, n)
    if denom ≥ 0.0
        return nocollision(T)
    else
        t = dot(w.sp - p.pos, n)/denom #??
        return t, p.pos + t*p.vel
    end
end

@traitfn function collision(p::P, w::FiniteWall{T}) where {T,P<:AbstractParticle{T}; 
                                                              PropagatesLinearly{P}}
    n = normalvec(w, p.pos)
    denom = dot(p.vel, n)
    # case of velocity pointing away of wall:
    denom ≥ 0.0 && return nocollision(T)
    posdot = dot(w.sp-p.pos, n)
    # Case of particle starting behind finite wall:
    posdot ≥ 0.0 && return nocollision(T)
    colt = posdot/denom
    i = p.pos + colt * p.vel
    dfc = norm(i - w.center)
    if dfc > w.width/2
        return nocollision(T)
    else
        return colt, i
    end
end

@traitfn function collision(p::P, d::Circular{T}) where {T,P<:AbstractParticle{T};
                                                            PropagatesLinearly{P}}

    dotp = dot(p.vel, normalvec(d, p.pos))
    dotp ≥ 0.0 && return nocollision(T)

    dc = p.pos - d.c
    B = dot(p.vel, dc)           #pointing towards circle center: B < 0
    C = dot(dc, dc) - d.r*d.r    #being outside of circle: C > 0
    Δ = B*B - C

    Δ ≤ 0.0 && return nocollision(T)
    sqrtD = sqrt(Δ)

    # Closest point:
    t = -B - sqrtD
    return t, p.pos + t*p.vel
end

@traitfn function collision(p::P, d::Semicircle{T}) where {T,P<:AbstractParticle{T}; 
                                                              PropagatesLinearly{P}}
    dc = p.pos - d.c
    B = dot(p.vel, dc)         #velocity towards circle center: B > 0
    C = dot(dc, dc) - d.r*d.r    #being outside of circle: C > 0
    Δ = B^2 - C

    Δ ≤ 0 && return nocollision(T)
    sqrtD = sqrt(Δ)

    nn = dot(dc, d.facedir)
    if nn ≥ 0 # I am NOT inside semicircle
        # Return most positive time
        t = -B + sqrtD
    else # I am inside semicircle:
        t = -B - sqrtD
        # these lines make sure that the code works for ANY starting position:
        if t ≤ 0 || distance(p, d) ≤ accuracy(T)
            t = -B + sqrtD
        end
    end
    # This check is necessary to not collide with the non-existing side
    newpos = p.pos + p.vel * t
    if dot(newpos - d.c, d.facedir) ≥ 0 # collision point on BAD HALF;
        return nocollision(T)
    end
    # If collision time is negative, return Inf:
    t ≤ 0.0 ? nocollision(T) : (t, p.pos + t*p.vel)
end

# TODO AFTER THIS LINE
@inline @traitfn function resolvecollision!(p::::PropagatesLinearly, o::Obstacle)
    specular!(p, o)
end

@inline @traitfn function propagate!(p::P, t::Real) where {T,P;PropagatesLinearly{P}}
    p.pos += SV{T}(p.vel[1]*t, p.vel[2]*t)
end

@inline @traitfn propagate!(p::::PropagatesLinearly, newpos::SV, t::Real) = (p.pos = newpos)

# add collisions for other obstacles