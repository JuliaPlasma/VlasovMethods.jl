using ParallelTestRunner
using VlasovMethods

const TESTS = [
    "core" => [
        "Spline Distribution" => "distributions/spline_distribution.jl",
        "Grid Distribution" => "distributions/grid_distribution.jl",
        "Tensor Spline Distribution" => "distributions/tensor_spline_distribution.jl",
        "Low-Rank Distribution" => "lowrank/distribution.jl",
        "Low-Rank Fields" => "lowrank/fields.jl",
        "Low-Rank BUG" => "lowrank/bug.jl",
        "Projector Splitting" => "lowrank/projector_splitting.jl",
        "Particle Distribution" => "distributions/particle_distribution.jl",
        "Projections" => "projections/projections.jl",
        "Reduced Tensors" => "gridbased/reduced_tensors.jl",
        "Geometric Integrator" => "methods/geometric_integrator.jl",
        "Splitting Method" => "methods/splitting.jl",
        "Parameter Order" => "integration/parameter_order.jl",
        "Cache Dict" => "cache.jl",
        "Vlasov–Poisson" => "integration/vlasov_poisson.jl",
        "Vlasov with External Fields" => "models/vlasov_with_external_fields.jl",
        "ElectromagneticFields Extension" => "integration/electromagnetic_fields_ext.jl",
        "SODE Ensemble" => "integration/sodeensemble.jl",
        "Names" => "integration/names.jl",
        "Moments" => "gridbased/moments.jl",
        "Vlasov–Poisson Advection" => "models/vlasov_poisson.jl",
        "Grid 1D1V Vlasov–Poisson" => "gridbased/vlasov_poisson_1d1v.jl",
        "Spline 1D1V Vlasov–Poisson" => "splinebased/vlasov_poisson_1d1v.jl",
        "Snapshots" => "particles/snapshots.jl",
        "Electric Field" => "electric_field.jl",
        "Reduced Electric Field" => "particles/electric_field.jl",
        "Reduced March" => "particles/time_marching.jl",
        "Diagnostics" => "particles/diagnostics.jl",
        "Reduced Equations" => "particles/reduced_equations.jl",
        "Bump-on-tail Scripts" => "integration/bump_on_tail_scripts.jl"
    ],
    "quality" => [
        "Aqua" => "quality/aqua.jl",
        "JET" => "quality/jet.jl",
        "ExplicitImports" => "quality/explicit_imports.jl",
        "Plotting" => "plotting.jl",
        "Plotting extension" => "integration/makie_extension.jl"
    ],
    "slow" => [
        "Bump-on-tail Runs" => "integration/bump_on_tail_runs.jl",
        "HNN" => "hnn.jl",
        "Multi-modal Vlasov Frontends" => "rom/multimodal_vlasov_frontends.jl",
        "GeometricIntegrators Scaling" => "integration/geometric_integrators_scaling.jl",
        "Spline Galerkin Energy and Rates" => "verification/spline_galerkin.jl"
    ],
    "collisions" => [
        "Landau Solver" => "integration/landau_solver.jl",
        "Metriplectic Solve" => "integration/metriplectic_solve.jl",
        "Collision Operators" => "integration/collision_operators.jl",
        "Lenard–Bernstein Buffers" => "models/lenard_bernstein_buffers.jl",
        "Metriplectic Lenard–Bernstein Brackets" => "models/lenard_bernstein_metriplectic.jl",
        "Landau Buffers" => "models/landau.jl",
        "Pairwise Landau" => "models/landau_pairwise.jl",
        "Discrete Gradient" => "methods/discrete_gradient.jl",
        "Reference Cases" => "examples/reference_cases.jl",
        "Collisions" => "gridbased/collisions.jl"
    ],
    "verification" => [
        "Reduced March at T = 1000" => "verification/time_marching.jl",
        "Reduced Equations at T = 100" => "verification/reduced_equations.jl",
        "Bounded Velocity Closure Runs" => "verification/bounded_closure.jl",
        "Low-Rank BUG Landau Rates" => "verification/lowrank_bug_landau_rates.jl",
        "Low-Rank BUG Two-Stream Rates" => "verification/lowrank_bug_two_stream_rates.jl",
        "Low-Rank BUG Two-Stream Rates, Float32" => "verification/lowrank_bug_two_stream_rates_float32.jl",
        "Projector Splitting Rates" => "verification/projector_splitting_rates.jl"
    ]
]

const SELECTED = filter(!startswith("--"), ARGS)
const GROUPS = isempty(SELECTED) ?
               ["core", "quality", "slow", "collisions", "verification"] : SELECTED
const UNKNOWN = setdiff(SELECTED, first.(TESTS))
isempty(UNKNOWN) ||
    error("no test group ", join(UNKNOWN, ", "), " in runtests.jl; its groups are ", join(first.(TESTS), ", "))

runtests(VlasovMethods, filter(startswith("--"), ARGS);
    testsuite = Dict(label => :(include($(joinpath(@__DIR__, path))))
    for (group, files) in TESTS if group in GROUPS for (label, path) in files),
    serial = [label for (group, files) in TESTS if group in ("metal", "cuda")
              for (label, _) in files],
    serial_position = :after,
    retries = 0,
    exeflags = ["--threads=$(Threads.nthreads())"],
    max_default_jobs = max(1, Sys.CPU_THREADS ÷ Threads.nthreads()))
