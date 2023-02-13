struct CollisionEntropy{XD, VD, DT <: DistributionFunction{XD, VD}} <: Entropy

    dist::DT
    cache::SplineDistributionCache{DT}
    # entropy::ET 
    
    function CollisionEntropy(dist::DistributionFunction{XD,VD}) where {XD,VD}
        new{XD, VD, typeof(dist)}(dist, SplineDistributionCache(dist))
    end
end

# ## TODO: add functions for computing the entropy given a distribution
# function compute_entropy!(entropy, dist <: DistributionFunction{XD, VD}) where {XD, VD}
    
# end