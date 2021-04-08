include("constants.jl")

struct Precipitator
    r₁::AbstractFloat
    r₂::AbstractFloat
    l::AbstractFloat
    T::AbstractFloat
    p::AbstractFloat
    q::AbstractFloat
end

relerr(n::Int, eps::AbstractFloat, f::List) =
    abs(takenth(f, n)[1] / takenth(f, n + 1)[1] - 1) < eps ? 
    takenth(f, n + 1)[1] :
    relerr(n + 1, eps, f)

λ(Λ) = 6.6e-8 * (101315.0 / Λ.p) * (Λ.T / 293.15)
cc(Λ, d) = 1.0 + λ(Λ) / d * (2.34 + 1.05 * exp(-0.39 * d / λ(Λ)))
η(Λ) = 1.83245e-5 * exp(1.5 * log(Λ.T / 296.1)) * (406.55) / (Λ.T + 110.4)
vtoz(Λ, v) = Λ.q ./ (2.0π * Λ.l * v) * log(Λ.r₂ / Λ.r₁)
ztov(Λ, z) = Λ.q ./ (2.0π * Λ.l * z) * log(Λ.r₂ / Λ.r₁)
dtoz(Λ, d) = ec * cc(Λ, d) / (3.0π * η(Λ) * d)
ztod(Λ, z, d) = ec * cc(Λ, d) / (3.0π * η(Λ) * z)
ztod(Λ, z) = relerr(1, 1e-11, iterated(d -> ztod(Λ, z, d), 10e-9))
vtod(Λ, v) = @>> v vtoz(Λ) ztod(Λ)
dtov(Λ, d) = @>> d dtoz(Λ) ztov(Λ)

Er(Λ, v, r) = v / (r * log(Λ.r₂ / Λ.r₁))
function Er(Λ, v, r, x) 
    out = (x > 0.012) ? v / (r * log(Λ.r₂ / Λ.r₁)) : 0
end

function ux(Λ, r)
    r₁ = Λ.r₁
    r₂ = Λ.r₂
    A1 =
        -2.0 *
        Λ.q *
        log(r₂ / r₁) *
        (r^2 + (r₁^2.0 * log.(r / r₂) - r₂^2.0 * log(r / r₁)) / log(r₂ / r₁))
    A2 =
        pi * (
            2.0 * log(r₁ / r₂) * r₁^4.0 + log(r₂ / r₁) * r₁^4.0 + log(r₂ / r₁) * r₂^4.0 -
            r₁^4.0 + 2.0 * r₁^2 * r₂^2 - r₂^4.0
        )

   return A1 / A2
end

function get_weights(Λ)
    rrs = range(Λ.r₁, stop = Λ.r₂, length = 1000)
    rmid = (rrs[2:end] - rrs[1:end-1])./2 .+ rrs[1:end-1]
    us = map(x->ux(Λ, x), rmid)

    frac = map(1:length(rrs)-1) do i
        pi * us[i] * (rrs[i+1]^2 -rrs[i]^2)./1.666e-5 / (rrs[i+1]-rrs[i])
    end
    frac = frac ./ maximum(frac)

    @> interpolate((rmid,), frac, Gridded(Linear())) extrapolate(0)
end