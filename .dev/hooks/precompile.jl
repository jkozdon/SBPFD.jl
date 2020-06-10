using JuliaFormatter

include(joinpath(@__DIR__, "..", "formatter_options.jl"))

format(@__FILE__; formatter_options...)
