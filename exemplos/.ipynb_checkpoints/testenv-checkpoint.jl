includet("../src/BilharesMarkov.jl")
using .BilharesMarkov
using DynamicalBilliards, PyCall
using StaticArrays, SparseArrays, Arpack
using LinearAlgebra, Statistics
using PyPlot, BenchmarkTools
const SV = SVector{2}
anim = pyimport("matplotlib.animation")

pygui(true)
bd = billiard_bunimovich()
# P = markovmap_portion(bd,1000,200,20,20)

function anim1()
    fig = figure("MyFigure",figsize=(5,5))
    Pt = copy(P)
    img = imshow(Pt)

    # Define the init function, which draws the first frame (empty, in this case)
    function init()
        img.set_data([[],[]])
        return (img,"")  # Union{} is the new word for None
    end

    # # Animate draws the i-th frame, where i starts at i=0 as in Python.
    function update(i)
        Pt *= P
        img.set_data(Pt ./ maximum(Pt))
        return (img,"")  # Union{} is the new word for None
    end

    ani = anim.FuncAnimation(fig, update, frames=60,
    init_func=init, blit=false)

    ani.save("circle.gif", writer="PillowWriter", fps=5)
    plt.show()
end


# Desenha a cadeia de Markov
# est_plot(g) = (
#   prob = p(g);
#   graphplot(
#     g,
#     arrow=true,
#     marker=:circle,
#     node_weights=prob,
#     edgewidth= g./maximum(g),
#     markercolors=get(colorschemes[:vik],prob./maximum(prob)),
#   )
# )