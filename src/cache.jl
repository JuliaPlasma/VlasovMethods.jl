
abstract type Cache{T} end

# The caches of one model, keyed by element type. Like each cache, the dictionary belongs to one
# task: a model evaluated on several threads at once shares its buffers and the dictionary.
struct CacheDict{parentType}
    parent::parentType
    caches::Dict{DataType, Cache}

    function CacheDict(p)
        caches = Dict{DataType, Cache}()
        caches[eltype(p)] = p
        new{typeof(p)}(p, caches)
    end
end

Base.parent(cd::CacheDict) = cd.parent

@inline function Base.getindex(c::CacheDict, T::DataType)
    if haskey(c.caches, T)
        c.caches[T]
    else
        c.caches[T] = Cache(T, parent(c))
    end::CacheType(T, parent(c))
end

# Every collision cache is a particle/spline pair with operator-specific buffers added by its
# own inner constructor. The two methods that re-type a cache for a new element type are
# identical for every cache, so they are generated once rather than repeated per cache, and so
# is the element type, which is that of the cache's spline: `CacheDict(p)` stores its parent `p`
# under `eltype(p)`, so that a lookup of the parent's type finds the parent. It is defined on the
# type, which `eltype` of an instance reads, so a cache may also define it on an instance.
#
# A cache with a fourth type parameter, the type of its `values` buffers, passes the word `buffers`.
macro collision_cache(CacheT, buffers = nothing)
    celltype = if buffers === :buffers
        :($CacheT{AT, PT, similar_type(AT, c.sdist), _buffer_type(AT, c.sdist)})
    else
        :($CacheT{AT, PT, similar_type(AT, c.sdist)})
    end
    esc(quote
        Base.eltype(::Type{<:$CacheT{T}}) where {T} = T
        function Cache(AT, c::$CacheT)
            $CacheT{AT}(c.pdist, similar(AT, c.sdist))
        end
        function CacheType(AT, c::$CacheT{<:Any, PT}) where {PT}
            $celltype
        end
    end)
end
