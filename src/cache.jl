
abstract type Cache{T} end

_cachehash(T) = hash(Threads.threadid(), hash(T))

struct CacheDict{parentType}
    parent::parentType
    caches::Dict{UInt64, Cache}

    function CacheDict(p)
        caches = Dict{UInt64, Cache}()
        caches[_cachehash(eltype(p))] = p
        new{typeof(p)}(p, caches)
    end
end

Base.parent(cd::CacheDict) = cd.parent

@inline function Base.getindex(c::CacheDict, T::DataType)
    key = _cachehash(T)
    if haskey(c.caches, key)
        c.caches[key]
    else
        c.caches[key] = Cache(T, parent(c))
    end::CacheType(T, parent(c))
end

# Every collision cache is a particle/spline pair with operator-specific buffers added by its
# own inner constructor. The two methods that re-type a cache for a new element type are
# identical for every cache, so they are generated once rather than repeated per cache, and so
# is the element type, which is that of the cache's spline: `CacheDict(p)` stores its parent `p`
# under `eltype(p)`, so that a lookup of the parent's type finds the parent. It is defined on the
# type, which `eltype` of an instance reads, so a cache may also define it on an instance.
macro collision_cache(CacheT)
    esc(quote
        Base.eltype(::Type{<:$CacheT{T}}) where {T} = T
        function Cache(AT, c::$CacheT)
            $CacheT{AT}(c.pdist, similar(AT, c.sdist))
        end
        function CacheType(AT, c::$CacheT{<:Any, PT}) where {PT}
            $CacheT{AT, PT, similar_type(AT, c.sdist)}
        end
    end)
end
