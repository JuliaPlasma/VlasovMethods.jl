
function IM_rule!(g, x, y, f, Δt, params)
    g .= x .- y .- Δt .* f((y .+ x) ./ 2, params)
end
