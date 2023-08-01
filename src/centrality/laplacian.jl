using LinearAlgebra

"""
    energy(lap_matrix)

Computes the Laplacian energy of the graph which laplacian matrix is the given matrix.
"""
function energy(lap_matrix::Matrix)
    return sum(eigvals(lap_matrix) .^ 2)
end

function _laplacian_vertex_centrality(
    lap_matrix::SparseMatrixCSC,
    vertex_index::Integer,
    full_energy::Number,
    normalized::Bool,
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

julia> round.(laplacian_centrality(star_graph(3)), digits=1)
3-element Vector{Float64}:
 1.0
 0.6
 0.6

julia> round.(laplacian_centrality(path_graph(4)), digits=3)
4-element Vector{Float64}:
 0.375     
 0.75       
 0.75       
 0.375
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

    edge_weights = weights(g)

    adjacency_list = SimpleGraphs.adj(SimpleGraph(g))

    for (index, vertex_index) in enumerate(vs)
        # laplacian_centralities[index] = _laplacian_vertex_centrality(
        #     lap_matrix, vertex_index, full_energy, normalized
        # )
        laplacian_centralities[index] = _laplacian_vertex_centrality_v2(
            vertex_index, edge_weights, adjacency_list, full_energy, normalized
        )
    end

    return laplacian_centralities
end

"""
    energy_difference(g, weights, adjacency_list)

Computes the difference between Laplacian energy of a graph and that of its subgraph removed from a given vertex v.

The computation involves the count of 2-walks containing v:
Closed 2-walks nwc2, non-closed 2-walks containing vertex v as middle point nwm2 and as extreme point nwe2.

"""
function energy_difference(vertex_index::Integer, weights::AbstractMatrix, adjacency_list)
    v_neighbors = adjacency_list[vertex_index]

    num_neighbors = length(v_neighbors)

    nwc2 = sum([weights[vertex_index, i]^2 for i in v_neighbors])
    nwm2 = sum([
        weights[vertex_index, i] * weights[vertex_index, j] for i in
                                                                range(1, num_neighbors) for
        j in range(i + 1, num_neighbors)
    ])
    nwe2 = sum([
        sum([
            weights[vertex_index, y] * weights[y, z] for
            z in adjacency_list[y] if z !== vertex_index
        ]) for y in v_neighbors
    ])

    return 4 * nwc2 + 2 * nwe2 + 2 * nwm2
end

function _laplacian_vertex_centrality_v2(
    vertex_index::Integer,
    weights::AbstractMatrix,
    adjacency_list,
    full_energy::Number,
    normalized::Bool,
)
    lapl_cent = energy_difference(vertex_index, weights, adjacency_list)

    if normalized
        lapl_cent = lapl_cent / full_energy
    end

    return lapl_cent
end
