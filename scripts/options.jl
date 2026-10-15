# The command line and the output directories of the drivers in this directory that write files.
# `electric_field_timing.jl` writes none and does not include this: its one argument is the name
# of the path its cold subprocess times.
#
# Code is in `scripts/`; a driver writes its numbers to `runs/` and its figures and reports to
# `results/`, both at the repository root and ignored by git. `--runs-dir` and `--results-dir`
# redirect them, which is why `RUNS_DIR` and `RESULTS_DIR` are `Ref`s and not constants.
#
# The parser is shared so that every driver rejects an unknown flag the same way. A driver that
# ignored one would run at its defaults and exit 0, and a mistyped regeneration command would then
# look like it worked.

const RUNS_DIR = Ref(joinpath(@__DIR__, "..", "runs"))
const RESULTS_DIR = Ref(joinpath(@__DIR__, "..", "results"))

flag(name) = "--" * replace(String(name), "_" => "-")

parse_value(::T, value) where {T <: Number} = parse(T, value)
parse_value(::AbstractString, value) = String(value)
parse_value(::Vector{T}, value) where {T} = parse.(T, split(value, ","))

"""
    parse_arguments(args; defaults...) -> NamedTuple

The options of a driver: one `--flag value` pair for each keyword of `defaults`, the flag spelled
with `-` for `_` (`final_time` is `--final-time`), plus `--runs-dir` and `--results-dir`. Each
value is parsed to the type of its default, a `Vector` default from a comma-separated list. An
unknown flag, a flag given twice, a missing value and a value that begins with `--` each throw an
`ArgumentError`.

Sets `RUNS_DIR` and `RESULTS_DIR` and creates both directories, so a driver calls this before it
writes anything.
"""
function parse_arguments(args; defaults...)
    defaults = (; runs_dir = RUNS_DIR[], results_dir = RESULTS_DIR[], defaults...)
    known = Dict(flag(name) => name for name in keys(defaults))
    given = Dict{Symbol, Any}()

    for pair in Iterators.partition(args, 2)
        argument = first(pair)
        haskey(known, argument) || throw(ArgumentError(
            "unknown argument `$(argument)`; this script takes " *
            join(sort!(collect(keys(known))), ", ")))
        name = known[argument]
        haskey(given, name) && throw(ArgumentError("`$(argument)` is given twice"))
        length(pair) == 2 && !startswith(pair[2], "--") ||
            throw(ArgumentError("`$(argument)` needs a value"))
        given[name] = parse_value(defaults[name], pair[2])
    end

    options = merge(defaults, NamedTuple(given))
    RUNS_DIR[] = mkpath(options.runs_dir)
    RESULTS_DIR[] = mkpath(options.results_dir)

    return options
end

"""
    choose(choices, argument, value)

The entry of the dictionary `choices` under `value`, the value of `argument`; an `ArgumentError`
that names the accepted keys when there is none.
"""
function choose(choices, argument, value)
    return get(choices, value) do
        throw(ArgumentError("unknown $(argument) `$(value)`; this script takes " *
                            join(sort!(collect(keys(choices))), ", ")))
    end
end

runs_path(name) = joinpath(RUNS_DIR[], name)
results_path(name) = joinpath(RESULTS_DIR[], name)
