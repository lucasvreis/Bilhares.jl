

### Bilhar

function PyPlot.plot(bd::Billiard{T}, P::Type{<:AbstractParticle};
    ax = (PyPlot.figure(); PyPlot.gca())) where {T}
    PyPlot.sca(ax)
    for obst in bd; plot(obst, P); end
    # xmin, ymin, xmax, ymax = cellsize(bd)
    # dx = xmax - xmin; dy = ymax - ymin
    # ax.set_aspect("equal")
    # if !isinf(xmin) && !isinf(xmax)
    #     PyPlot.xlim(xmin - 0.1dx, xmax + 0.1dx)
    # end
    # if !isinf(ymin) && !isinf(ymax)
    #     PyPlot.ylim(ymin - 0.1dy, ymax + 0.1dy)
    # end
    return nothing
end


### Obstáculos

PyPlot.plot(o::Obstacle, ::Type{Particle}; kwargs...) = plot(o; kwargs...)

function _oc_plot(rtransf::Function, d::Circular; kwargs...)
    edgecolor = DynamicalBilliards.obcolor(d)
    facecolor = (edgecolor..., DynamicalBilliards.obalpha(d))
    circle1 = PyPlot.plt.Circle(d.c, d.r*rtransf(d.r^2);
        edgecolor = edgecolor, facecolor = facecolor,
        linestyle = DynamicalBilliards.obls(d), lw = 2.0, kwargs...)
    PyPlot.gca().add_artist(circle1)
end


### Partículas

function _p_plot(transf, ps, colrs, ax)
    x = []; y = []
    for p in ps
        x′, y′ = transf(p.pos...)
        push!(x,x′)
        push!(y,y′)
    end
    ax.scatter(x, y; color=colrs, s=20.0, zorder=99)
end

PyPlot.plot(ps::AbstractVector{<:AbstractParticle}, ::Type{Particle}, colrs; kwargs...) = 
    plot(ps, colrs; kwargs...)

PyPlot.plot(ps::AbstractVector{<:AbstractParticle}, colrs; ax=PyPlot.gca()) =
    _p_plot((x...) -> x, ps, colrs, ax)
