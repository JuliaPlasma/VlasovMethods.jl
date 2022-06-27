
using Plots, LaTeXStrings

function remap_to_domain(x::AbstractVector, xmin::Number, xmax::Number)
    X = copy(x)
    L = xmax - xmin
    
    for i in eachindex(X)
        while X[i] < xmin
            X[i] += L
        end
        while X[i] > xmax
            X[i] -= L
        end
    end

    return X
end


function plot_distribution(F::AbstractMatrix{T}, xGrid::AbstractVector{T}, vGrid::AbstractVector{T}) where {T}
    heatmap(xGrid, vGrid, F',
            xlims = (minimum(xGrid), maximum(xGrid)),
            ylims = (minimum(vGrid), maximum(vGrid)),
            title = "Vlasov-Poisson",
            xlabel = L"$x$",
            ylabel = L"$v$",
            legend = :none,
            grid   = :none)
end

function plot_distribution(f::Function, xGrid::AbstractVector{T}, vGrid::AbstractVector{T}) where {T}
    F = zeros(length(xGrid), length(vGrid))

    for i in eachindex(xGrid)
        for j in eachindex(vGrid)
            F[i,j] = f(xGrid[i], vGrid[j])
        end
    end

    plot_distribution(F, xGrid, vGrid)
end

function plot_distribution(f::Function, xmin::T, xmax::T, vmin::T, vmax::T; nx = 100, nv = 100) where {T}
    xGrid = collect(LinRange(xmin, xmax, nx))
    vGrid = collect(LinRange(vmin, vmax, nv))
    plot_distribution(f, xGrid, vGrid)
end


function plot_distribution(X::AbstractVector, V::AbstractVector, W::AbstractVector, xmin, xmax, vmin, vmax; nx = 100, nv = 100)
    xGrid = collect(LinRange(xmin, xmax, nx))
    vGrid = collect(LinRange(vmin, vmax, nv))

    X = remap_to_domain(X, xmin, xmax)
    F = zeros(nx, nv)

    for k in eachindex(X,V)
        if xmin < X[k] < xmax && vmin < V[k] < vmax
            i = Int(floor((X[k]-xmin)/(xmax - xmin)*nx))
            j = Int(floor((V[k]-vmin)/(vmax - vmin)*nv))
            F[i+1,j+1] += W[k]
        end
    end

    plot_distribution(F, xGrid, vGrid)
end


function plot_particles(X::AbstractVector, V::AbstractVector, W::AbstractVector, xmin, xmax, vmin, vmax)
    scatter(remap_to_domain(X, xmin, xmax), V,
            xlims = (xmin, xmax),
            ylims = (vmin, vmax),
            title = "Vlasov-Poisson",
            xlabel = L"$x$",
            ylabel = L"$v$",
            legend = :none,
            grid   = :none)
end

function plot_particles(X::AbstractVector, V::AbstractVector, W::AbstractVector, xGrid::AbstractVector, vGrid::AbstractVector)
    plot_particles(X, V, W, minimum(xGrid), maximum(xGrid), minimum(vGrid), maximum(vGrid))
end
