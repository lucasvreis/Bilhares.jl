module Bilhares

export main, markovmap_portion, adjacencyMatrix, estplot
export p

using DynamicalBilliards
using DynamicalBilliards: increment_counter
using StaticArrays, LinearAlgebra, Arpack, SparseArrays
using PyPlot, Plots, GraphRecipes, ColorSchemes
using ProgressMeter
const SV = SVector{2}


p(g::SparseMatrixCSC) = normalize(eigs(g';nev=1,which=:LR)[2][:,1],1) .|> abs

estplot(g) = (
  prob = p(g);
  graphplot(
    g,
    arrow=true,
    marker=:circle,
    node_weights=prob,
    edgewidth= g./maximum(g),
    markercolors=get(colorschemes[:vik],prob./maximum(prob)),
  )
)

function adjacencyMatrix(maps::AbstractMatrix) 
  P = copy(maps)
  for m in eachrow(P)
    if sum(m) != 0
      m ./= sum(m)
    else
      m .*= 0
    end
  end
  P
end

function markovmap_portion(bd::Billiard{T}, Nbounces::Int,
                           np::Int, Nξ, Nφ=Nξ;
                           intervals=arcintervals(bd)) where {T}
  δξ = totallength(bd) / Nξ # Tamanho dos compartimentos-arco
  δφ = 2 / Nφ # Tamanho dos compartimentos-ângulo
  # δφ = π / Nφ

  p = [randominside(bd) for _ in 1:np]
  # p = [Particle(from_bcoords(rξ,rφ,bd)[1:2]...)
  #      for (rξ,rφ) in 
  #           zip(totallength(bd)*rand(np),
  #            sin.(π*(rand(np) .- .5)))]

  pbounce!(x) = bounce!(x,bd)

  P = zeros(MMatrix{Nξ * Nφ,Nξ * Nφ,Int}) # Matriz estocástica
  # pP = [spzeros(Int, Nξ * Nφ, Nξ * Nφ) for i in 1:np]
  # d = Dict{SV{Int},Int}() 

  lastind = Vector{SV{Int}}() # Representa o último Índice

  vectoind(ind) = ind[2] * Nξ + ind[1] + 1
  # vectoind(ind) = ind[1] * Nφ + ind[2] + 1

  lastind = nothing
  # prog = Progress(Nbounces)
  @showprogress @distributed for N in 1:Nbounces
    i, = zip(pbounce!.(p)...)

    ξ, sφ = zip(to_bcoords.(p, [bd[k] for k in i])...)
    # φ = asin.(sφ)
    ξ = ξ .+ [intervals[k] for k in i]

    ind = SV{Int}.(floor.(Int, ξ ./ δξ), 
                   floor.(Int, (sφ .+ 1) ./ δφ))
                  #  floor.(Int, (φ .+ π/2) ./ δφ))
    
    if N > 1
      for k in 1:np
        P[vectoind(lastind[k]),vectoind(ind[k])] += 1
      end
    end

    lastind = ind
  end # collision number loop

  # total_boxes = Nξ * Nφ
  # ratio = length(keys(d)) / total_boxes

  return adjacencyMatrix(P./1)
end


function main()
  bd = billiard_bunimovich()

  pygui(true)


  n = 100 # how many particles to create
  t = 200 # how long to evolve each one


  bd = billiard_mushroom()

  markovmap_portion(bd, 2000, randominside(bd), 3, 3)

  # out = parallelize(markovmap_portion, bd, t, n)

  # colors = ["C$(rand(1:9))" for i in 1:n] # random colors

  # plot_boundarymap((bmap), arcs, color = colors) # (x->(x->[x[1],asin(x[2])]).(x)).

end


end
