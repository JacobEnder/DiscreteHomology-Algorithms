#===================================================================#
# precompute_templates.jl
#
# Offline script to build the full Q^3 quotient template catalog and
# an automatically minimized universal generating subcatalog.
#
# Usage:
#   julia precompute_templates.jl
#
# Outputs (relative to this directory):
#   data/templates_full.json
#   data/templates_min.json
#===================================================================#

include("src/GraphTypes.jl")
include("src/Cube.jl")
include("src/Isomorphism.jl")
include("src/Cycles.jl")
include("src/Embeddings.jl")
include("src/Templates.jl")
include("src/TemplateMinimize.jl")

out_full = joinpath(@__DIR__, "data", "templates_full.json")
out_min  = joinpath(@__DIR__, "data", "templates_min.json")

println("Building full template catalog from Q^3 partitions...")
full = build_template_catalog(; drop_degenerate=true, drop_zero_boundary=true)
mkpath(dirname(out_full))
write_template_catalog(out_full, full)
println("  wrote $out_full  (|T|=$(length(full)))")

println("Minimizing by universal redundancy elimination...")
min = minimize_template_catalog(full)
write_template_catalog(out_min, min)
println("  wrote $out_min  (|T_min|=$(length(min)))")
