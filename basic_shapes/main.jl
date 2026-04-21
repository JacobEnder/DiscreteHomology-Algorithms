#===================================================================#
# main.jl
# Entry point for the v2 discrete H_2 computation pipeline.
#
# Usage:
#   julia main.jl path/to/graphs.json
#
# If no json file is provided, we run a small built-in smoke test.
#===================================================================#

include("src/GraphTypes.jl")
include("src/Cube.jl")
include("src/Isomorphism.jl")
include("src/Cycles.jl")
include("src/Embeddings.jl")
include("src/Templates.jl")
include("src/TemplateMinimize.jl")
include("src/ChainComplex.jl")
include("src/ComputeH2.jl")

const DEFAULT_FULL_TEMPLATES = joinpath(@__DIR__, "data", "templates_full.json")
const DEFAULT_MIN_TEMPLATES  = joinpath(@__DIR__, "data", "templates_min.json")

"""Load a template catalog from disk, or build it (and write it) if missing."""
function load_or_build_templates(; minimized::Bool=true, verbose::Bool=true)
    want = minimized ? DEFAULT_MIN_TEMPLATES : DEFAULT_FULL_TEMPLATES
    if isfile(want)
        println("in here")
        verbose && println("Loading templates from $want ...")
        return read_template_catalog(want)
    end

    verbose && println("No template file found; building template catalog from Q^3 partitions...")
    full = build_template_catalog(; drop_degenerate=true, drop_zero_boundary=true)
    mkpath(dirname(DEFAULT_FULL_TEMPLATES))
    write_template_catalog(DEFAULT_FULL_TEMPLATES, full)
    verbose && println("Wrote full catalog to $DEFAULT_FULL_TEMPLATES (|T|=$(length(full)))")

    if minimized
        verbose && println("Minimizing catalog by universal redundancy elimination...")
        min = minimize_template_catalog(full)
        write_template_catalog(DEFAULT_MIN_TEMPLATES, min)
        verbose && println("Wrote minimized catalog to $DEFAULT_MIN_TEMPLATES (|T_min|=$(length(min)))")
        return min
    else
        return full
    end
end

"""Run the pipeline on a vector of graphs."""
function run_pipeline(graphs::Vector{SimpleGraph}; templates::Vector{Template})
    results = []
    for (i, G) in enumerate(graphs)
        println("="^70)
        println("Graph $i: $(G.name)")
        println("  |V|=$(length(G.vertices))  |E|=$(length(non_loop_edges(G)))")
        println("="^70)

        # --- Timed homology computation ---
        local cc, res
        elapsed = @elapsed begin
            cc = build_chain_complex(G; templates=templates, verbose=true)
            res = compute_h2(cc)
        end
        # ----------------------------------

        println("\nResult:")
        println("  dim(H_2) = $(res.h2)")
        println("  |C_1|=$(res.n_1cells), |C_2|=$(res.n_2cells), |C_3|=$(res.n_3cells)")
        println("  rank(∂_2)=$(res.rank_d2), rank(∂_3)=$(res.rank_d3)")
        println("  nullity(∂_2)=$(res.nullity_d2)")
        println("  Homology computation time: $(round(elapsed; digits=6)) s")

        push!(results, (G.name, res, elapsed))
    end

    println("\n" * "="^70)
    println("Summary")
    println("="^70)
    for (name, r, elapsed) in results
        println("  $name: dim(H_2)=$(r.h2)  [$(round(elapsed; digits=6)) s]")
    end
end

# ---------------------------------------------------------------
# Main
# ---------------------------------------------------------------

templates = load_or_build_templates(; minimized=true, verbose=true)

if length(ARGS) >= 1
    filename = ARGS[1]
    graphs = parse_graph_json(filename)
    run_pipeline(graphs; templates=templates)
else
    println("No input json provided; running a smoke test...")

    # Triangle graph K3 (H_2 should be 0 in this chain model)
    V = [1,2,3]
    rawE = Tuple{Int,Int}[(1,2),(2,3),(1,3)]
    E = make_symmetric_reflexive(rawE, V)
    Gtri = SimpleGraph("K3", V, E, build_adjacency(V, E))

    # Square C4
    V4 = [1,2,3,4]
    rawE4 = Tuple{Int,Int}[(1,2),(2,3),(3,4),(4,1)]
    E4 = make_symmetric_reflexive(rawE4, V4)
    Gsq = SimpleGraph("C4", V4, E4, build_adjacency(V4, E4))

    run_pipeline([Gtri, Gsq]; templates=templates)
end