using LinearAlgebra

"""
    energy(lap_matrix)

Computes the Laplacian energy of the graph which laplacian matrix is the given matrix.
"""
function energy(lap_matrix::Matrix)
    return sum(eigvals(lap_matrix) .^ 2)
end

function _laplacian_vertex_centrality(
    lap_matrix::SparseMatrixCSC, vertex_index::Integer, full_energy::Number, normalized::Bool
)
    vertex_removed_lap_matrix = lap_matrix[1:end .!= vertex_index, 1:end .!= vertex_index]
    new_diag = diag(lap_matrix) .- abs.(lap_matrix[:, vertex_index])
    vertex_removed_lap_matrix[diagind(vertex_removed_lap_matrix)] = new_diag[1:end .!= vertex_index]

    new_energy = energy(Matrix(vertex_removed_lap_matrix))
    lapl_cent = full_energy - new_energy

    if normalized
        lapl_cent = lapl_cent / full_energy
    end

    return lapl_cent
end

"""
    laplacian_centrality(g, normalized=true, vs=nothing)

Calculate the laplacian centrality of a graph `g` across all vertices or
a specified subset of vertices `vs`. Return a vector representing the
centrality calculated for each node in `g`.

The Laplacian Centrality of a node ``i`` is measured by the drop in the
Laplacian Energy after deleting node ``i`` from the graph. The Laplacian Energy
is the sum of the squared eigenvalues of a graph's Laplacian matrix.

### Optional Arguments
- `vs=nothing`: subset of vertices to compute centrality for.
- `normalize=true`: If true, normalize the centrality values by the
full energy of the original graph, making centralities sum equal to 1.
If  false, the centrality score for each node is the drop in Laplacian
energy when that node is removed.

### References
- Qi, X., Fuller, E., Wu, Q., Wu, Y., and Zhang, C.-Q. (2012). 
Laplacian centrality: A new centrality measure for weighted networks.

# Examples
```jldoctest
julia> using Graphs

julia> laplacian_centrality(star_graph(3))
3-element Vector{Float64}:
 1.0
 0.6
 0.6

julia> laplacian_centrality(path_graph(4))
4-element Vector{Float64}:
 0.3749999999999997       
 0.7499999999999999       
 0.7499999999999999       
 0.3749999999999997
```
"""
function laplacian_centrality(
    g::AbstractGraph; normalized::Bool=true, vs::Union{Nothing,Array{Int64,1}}=nothing
)
    nv(g) == 0 && error("null graph has no centrality defined")

    if vs !== nothing
        nodeset = Set([v for v in vs if has_vertex(g, v)])
        length(nodeset) !== length(vs) && error("vs has duplicate nodes or nodes not in g")
    else
        vs = vertices(g)
    end

    lap_matrix = laplacian_matrix(g; dir=is_directed(g) ? :out : :unspec)

    full_energy = energy(Matrix(lap_matrix))

    laplacian_centralities = zeros(length(vs))

    for (index, vertex_index) in enumerate(vs)
        laplacian_centralities[index] = _laplacian_vertex_centrality(
            lap_matrix, vertex_index, full_energy, normalized
        )
    end

    return laplacian_centralities
end
