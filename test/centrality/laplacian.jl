@testset "Laplacian" begin
    g5 = path_graph(4)
    for g in test_generic_graphs(g5)
        z = @inferred(laplacian_centrality(g.g))
        @test round.(z, digits=2) == [0.38, 0.75, 0.75, 0.38]
        @test all(z .<= 1)
    end

    g5 = SimpleGraph(5)
    add_edge!(g5, 1, 2)
    for g in test_generic_graphs(g5)
        z = @inferred(laplacian_centrality(g.g))
        @test z[1] == z[2] == 1.0
        @test z[3] == z[4] == z[5] == 0.0
        @test all(z .<= 1)
    end

    g5 = SimpleGraph(5)
    add_edge!(g5, 1, 2)
    for g in test_generic_graphs(g5)
        z = @inferred(laplacian_centrality(g.g, normalized=false))
        @test z[1] == z[2] == 4.0
        @test z[3] == z[4] == z[5] == 0.0
    end

    g3 = SimpleGraph(3)
    for g in test_generic_graphs(g3)
        z = @inferred(laplacian_centrality(g.g))
        @test all(isnan.(z))
    end

    g0 = SimpleGraph(0)
    for g in test_generic_graphs(g0)
        @test_throws "null graph has no centrality defined" laplacian_centrality(g.g)
    end

    g_complete_directed = complete_digraph(3)
    for g in test_generic_graphs(g_complete_directed)
        z = @inferred(laplacian_centrality(g.g))
        @test round(z[1]; digits=2) ==
            round(z[2]; digits=2) ==
            round(z[3]; digits=2) ==
            0.78
    end

    g_complete = complete_graph(5)
    for g in test_generic_graphs(g_complete)
        z = @inferred(laplacian_centrality(g.g, vs=[1, 3, 5]))
        @test round(z[1]; digits=2) ==
            round(z[2]; digits=2) ==
            round(z[3]; digits=2) ==
            0.52
    end

    g_complete = complete_graph(5)
    for g in test_generic_graphs(g_complete)
        @test_throws "vs has duplicate nodes or nodes not in g" laplacian_centrality(
            g.g, vs=[1, 3, 3, 5]
        )
    end

    g_complete = complete_graph(5)
    for g in test_generic_graphs(g_complete)
        @test_throws "vs has duplicate nodes or nodes not in g" laplacian_centrality(
            g.g, vs=[1, 3, 5, 6]
        )
    end
end
