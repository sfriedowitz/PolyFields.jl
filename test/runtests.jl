using Test
using PolyFields

#==============================================================================#

add_jl(x) = endswith(x, ".jl") ? x : x*".jl"

if length(ARGS) > 0
    tests = map(add_jl, ARGS)
else
    tests = [
        "test_cell.jl",
        "test_species.jl",
        "test_interactions.jl",
        "test_system.jl",
        "test_updaters.jl"
    ]
end

println("Testing...")

for test in tests
    include(test)
end