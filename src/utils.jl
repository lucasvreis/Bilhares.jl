
@inline normsq(v) = v[1]v[1] + v[2]v[2]

## Does it belong here?
function DynamicalBilliards.randominside(T::Type{<:AbstractParticle}, bd::Billiard, N::Int) 
    [T(randominside_xyφ(bd)...) for _ in 1:N]
end
