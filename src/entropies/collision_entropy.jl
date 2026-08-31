struct CollisionEntropy{XD, VD, DT <: DistributionFunction{XD, VD}} <: Entropy
    dist::DT
    # entropy::ET 

    function CollisionEntropy(dist::DistributionFunction{XD, VD}) where {XD, VD}
        new{XD, VD, typeof(dist)}(dist)
    end
end

function (ent::CollisionEntropy)(nquad = 5)
    params = (sdist = ent.dist,)
    gauss_quad_2d((v, params) -> ent.dist.spline(v) * log(ent.dist.spline(v)), ent.dist.basis, nquad, params)
end

# ## TODO: add functions for computing the entropy given a distribution
# function compute_entropy!(entropy, dist <: DistributionFunction{XD, VD}) where {XD, VD}

# end
