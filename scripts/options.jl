# The shared layer under every driver in this directory: where output goes and how the command line
# is parsed. No work happens at top level, so a driver can `include` this and pay only the load.
#
# Code in `scripts/`, data in `runs/`, figures in `results/` — the tree-wide rule in
# `Packages/CLAUDE.md`. `RUNS_DIR` and `RESULTS_DIR` are `Ref`s rather than constants so that
# `--runs-dir` and `--results-dir` can redirect them; a constant derived from `@__DIR__` is what
# forces a caller to *copy* a script in order to write somewhere else. `parse_arguments` creates
# both, so a run writes into a directory that did not exist before it.
#
# The two `Ref`s do execute, so including this twice in one session redefines them and resets a
# `--runs-dir` already set. Each driver includes it once.
#
# The parser is shared, and the rejection of an unknown flag is the reason. A driver that parses
# its own arguments disagrees with the next one about what an unrecognised flag means, and the
# dangerous shape is the one that renders nothing, prints "done" and exits 0: a regeneration
# recipe with a mistyped flag then looks like it worked.

const RUNS_DIR = Ref(joinpath(@__DIR__, "..", "runs"))
const RESULTS_DIR = Ref(joinpath(@__DIR__, "..", "results"))

const COMMON_OPTIONS = ("--runs-dir", "--results-dir")

"""
    parse_arguments(args, extra = ()) -> (positional, options)

Split `args` into positional names and `--flag value` pairs, accepting `--runs-dir` and
`--results-dir` everywhere plus whatever `extra` a driver adds, and **rejecting anything else**
with an `ArgumentError`.

Sets `RUNS_DIR` and `RESULTS_DIR` as a side effect and creates both directories, so a driver calls
this before it touches either.
"""
function parse_arguments(args, extra = ())
    known = (COMMON_OPTIONS..., extra...)
    options = Dict{String, String}()
    positional = String[]

    i = firstindex(args)
    while i ≤ lastindex(args)
        argument = args[i]
        if startswith(argument, "--")
            argument in known || throw(ArgumentError(
                "unknown option `$(argument)`; this script takes " * join(known, ", ")))
            i < lastindex(args) || throw(ArgumentError("`$(argument)` needs a value"))
            options[argument] = args[i + 1]
            i += 2
        else
            push!(positional, argument)
            i += 1
        end
    end

    haskey(options, "--runs-dir") && (RUNS_DIR[] = options["--runs-dir"])
    haskey(options, "--results-dir") && (RESULTS_DIR[] = options["--results-dir"])

    mkpath(RUNS_DIR[])
    mkpath(RESULTS_DIR[])

    return positional, options
end

"""
    no_positional(positional, script)

Throw an `ArgumentError` for a script that takes no positional argument, naming the ones it was
given.
"""
function no_positional(positional, script)
    isempty(positional) || throw(ArgumentError(
        "`$(script)` takes no positional argument, got `" * join(positional, "`, `") * "`"))

    return nothing
end

"""
    option_string(options, name, default) -> String

The `--<name>` value, or `default`.
"""
option_string(options, name, default) = get(options, name, default)

"""
    option_float(options, name, default) -> Float64

The `--<name>` value, or `default`.
"""
function option_float(options, name, default)
    haskey(options, name) || return default

    return parse(Float64, options[name])
end

"""
    option_int(options, name, default) -> Int

The `--<name>` value, or `default`.
"""
function option_int(options, name, default)
    haskey(options, name) || return default

    return parse(Int, options[name])
end

"""
    option_choice(options, name, default, allowed) -> String

The `--<name>` value, or `default`, rejecting a value outside `allowed`.
"""
function option_choice(options, name, default, allowed)
    value = get(options, name, default)
    value in allowed || throw(ArgumentError(
        "unknown $(name) `$(value)`; this script takes " * join(allowed, ", ")))

    return value
end

"""
    option_ints(options, name, default) -> Vector{Int}

The `--<name>` value as a comma-separated list of integers, or `default`.
"""
function option_ints(options, name, default)
    haskey(options, name) || return default

    return [parse(Int, field) for field in split(options[name], ",")]
end

"""
    runs_path(name) -> String

`name` under `RUNS_DIR`.
"""
runs_path(name) = joinpath(RUNS_DIR[], name)

"""
    results_path(name) -> String

`name` under `RESULTS_DIR`.
"""
results_path(name) = joinpath(RESULTS_DIR[], name)
