
"""
save spline solver parameters
"""
function h5save(h5::H5DataStore, P::PoissonSolverPBSplines; path::AbstractString = "/")
    group = _create_group(h5, path)
    attributes(group)["p"] = P.p
    attributes(group)["L"] = P.L
end

function PoissonSolvers.PoissonSolverPBSplines(h5::H5DataStore, path::AbstractString = "/")
    group = h5[path]
    p = read(attributes(group)["p"])
    L = read(attributes(group)["L"])
    # κ = read(group["parameters/κ"])

    group = group["integrator"]
    n = read(attributes(group)["nh"])

    PoissonSolverPBSplines(p, n, L)
end

function PoissonSolvers.PoissonSolverPBSplines(fpath::AbstractString, path::AbstractString = "/")
    h5open(fpath, "r") do file
        PoissonSolverPBSplines(file, path)
    end
end
