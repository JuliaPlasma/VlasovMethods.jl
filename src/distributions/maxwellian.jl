@doc raw"""
    MaxwellianDistribution(v)

The normalised Maxwellian of unit temperature and zero mean in `d = length(v)` velocity
dimensions,

```math
f_M(v) = (2\pi)^{-d/2} \, \exp\!\left(-\tfrac{1}{2} |v|^2\right) ,
\qquad \int_{\mathbb{R}^d} f_M \, dv = 1 .
```

!!! note "The normalisation is dimension-dependent"
    It was fixed at ``1/(2\pi)``, which is the ``d = 2`` value. That was correct where the
    function was used — `project_Maxwellian` existed only for a two-dimensional velocity
    space — but it is wrong for any other `d`, and silently so: in one dimension it is short by
    a factor ``\sqrt{2\pi}``, so a "normalised" Maxwellian integrated to ``0.399`` rather than
    to ``1``. The Lenard-Bernstein manuscript's ``f_M(v) = e^{-v^2/2}/\sqrt{2\pi}`` is the
    ``d = 1`` case of the expression above.
"""
function MaxwellianDistribution(v)
    d = length(v)
    return (2π)^(-d / 2) * exp(-dot(v, v) / 2)
end
