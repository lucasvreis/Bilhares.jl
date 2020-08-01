includet("Bilhares.jl")
using .Bilhares
using DynamicalBilliards, PyCall
using DynamicalBilliards: increment_counter
using StaticArrays, SparseArrays, Arpack
using LinearAlgebra, Statistics
using PyPlot, BenchmarkTools
const SV = SVector{2}
anim = pyimport("matplotlib.animation")

pygui(true)
bd = billiard_bunimovich()
# P = markovmap_portion(bd,1000,200,20,20)

function testplot()
    fig = figure("MyFigure",figsize=(20,20))
    ax = plt.axes(xlim=400,ylim=400)
    Pt = P
    img = ax.imshow(Pt)
    for i in 1:200
        img.set_data([[],[]])
        img = ax.imshow(Pt)
        Pt *= P
        # fig.colorbar(img)
        sleep(0.2)
    end
end


# # Define the init function, which draws the first frame (empty, in this case)
# function init()
#     global img
#     img.set_data([[],[]])
#     return (img,Union{})  # Union{} is the new word for None
# end

# # # Animate draws the i-th frame, where i starts at i=0 as in Python.
# function animate(i)
#     global img
#     img.set_data(P^i)
#     return (img,Union{})  # Union{} is the new word for None
# end

# bd = Billiard(Ellipse(SV{Float64}(0,0),3.,1.,false))