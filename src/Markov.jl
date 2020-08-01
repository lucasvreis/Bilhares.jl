include("Cylinders.jl")

module Markov

using LightGraphs, SimpleWeightedGraphs
using LinearAlgebra, NamedArrays, StaticArrays
using ColorSchemes, GraphRecipes, Plots
using .Cylinders

export Cylinder, intersect, univ, set, σ

estplot(g) = graphplot(
    g,
    arrow=true,
    marker=:circle,
    node_weights=prob,
    edgewidth= g./maximum(g),
    markercolors=get(colorschemes[:vik],prob./maximum(prob)),
)

function adjacencyMatrix(maps::AbstractMatrix) 
    P = copy(maps)
    for m in eachrow(maps)
        m ./= sum(m)
    end
    P
end

function adjacencyMatrix(maps::Vararg{Pair}) 
    l = length(maps)
    names = [m.first for m in maps]
    O = NamedArray(zeros(l,l), ( names, names ), ("I", "O"))
    for m in maps
        O["I"=> m.first] = m.second ./ sum(m.second)
    end
    O
end

randAdjacency(n) = adjacencyMatrix((i => rand(n) for i in 1:n)...)
randAdjacency(T,n) = adjacencyMatrix((i => rand(T,n) for i in 1:n)...)

p(g::SimpleWeightedDiGraph) = normalize(eigen(Array(adjacency_matrix(g))').vectors[:,end],1) .|> abs
p(g::AbstractMatrix) = normalize(eigen(g').vectors[:,end],1) .|> abs


