using Manifolds
using ManifoldsBase

const ℍ = Hyperbolic(2)

BKtoHyperboloid(x) = SV3(x[1], x[2], 1) / sqrt((1 - norm(x)) * (1 + norm(x)))
HyperboloidtoBK(x) = SV(x[1], x[2]) / x[3]

g_klein(x) = I / (1 - normsq(x)) + Diagonal(x.^2) / (1 - normsq(x))^2

dg_klein(x) = (
    d = (1 - normsq(x));
    (1 + x[1]^2 / d) * (1 + x[2]^2 / d) / d^2
)

# sqrt_dg_klein(x) = (
#     sqrt(normsq(x) + 1) / ((1 - norm(x))*(1 + norm(x)))^2
# )

D_h_BK(x) = @SMatrix[
    1 0 -x[1] / x[3];
    0 1 -x[2] / x[3]
] / x[3]

D_BK_h(x) = @SMatrix[
    1-x[2]^2 x[1]x[2];
    x[1]x[2] 1-x[1]^2;
    x[1]     x[2]
] / ((1-norm(x))*(1+norm(x)))^(3/2)

norm3m(x) = sqrt(x[3]^2 - x[1]^2 - x[2]^2)

sqrt_dg_klein(x) = (
    (d = D_BK_h(x); norm3m(cross(d[:,1],d[:,2])) |> sqrt)
)

@inline function hyperreflect(x, n, i)
    sq1x  = 1 + sqrt(1 - normsq(x))
    PTK = @SMatrix[
        sq1x - x[1]^2  -x[1]x[2];
        -x[1]x[2]  sq1x - x[2]^2
    ]
    KTP = inv(PTK)
    t′ = KTP * SV(-n[2], n[1])
    i′ = KTP * i
    o′ =  2 * t′ ⋅ i′ * t′ / normsq(t′) - i′ 
    normalize(PTK * o′)
end

@inline hyper_rtransf(s²) = 1 / (1 + sqrt(1 - s²))

@inline function hyper_transf(x, y)
    s² = x * x + y * y
    k = hyper_rtransf(s²)
    return k * x, k * y
end