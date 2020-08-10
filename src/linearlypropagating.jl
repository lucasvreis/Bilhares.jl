using SimpleTraits

import DynamicalBilliards: collision, resolvecollision!, propagate!, Obstacle

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

# TODO AFTER THIS LINE
@inline @traitfn function resolvecollision!(p::::PropagatesLinearly, o::Obstacle)
    specular!(p, o)
end

@inline @traitfn function propagate!(p::P, t::Real) where {T,P;PropagatesLinearly{P}}
    p.pos += SV{T}(p.vel[1]*t, p.vel[2]*t)
end

@inline @traitfn propagate!(p::::PropagatesLinearly, newpos::SV, t::Real) = (p.pos = newpos)

# add collisions for other obstacles