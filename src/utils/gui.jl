module BilliardGUI

using Interact, Blink

export billiard_gui

mush(s) = billiard_cmushroom(sw = .8, scale=s);
quad(s) = billiard_csq(s)
sinai(s) = billiard_csinai(s/4,s,s)
pol(s,n) = billiard_polygon(n,s)

function billiard_gui()
    Ns = spinbox(1:500,label="N")
    ts = spinbox(1:500,label="t")
    gcolrs = ["C$(rand(1:9))" for i in 1:500];

    curv = dropdown(Dict("Hiperbólico"=>HyperBKParticle, "Elíptico"=>EllipticParticle, "Plano"=>Particle))

    bils = dropdown(Dict(
        "cogumelo" => mush,
        "polígono" => pol,
        "sinai" => sinai
    ))

    bil_ui(bil,s) = (map(bil,s), nothing)
    function bil_ui(bil::typeof(pol),s) 
        n = spinbox(3:10,label="N lados")
        (map(bil,s,n),n)
    end
    b = button(label = "Outro!")

    pygui(true)
    gs = matplotlib.gridspec.GridSpec(1,2,width_ratios=[1,3])

    ph_plot =  map(curv) do Pt
        close("all")
        fig = plt.figure(figsize=(10,5))
        rcurv = togglebuttons(Dict('S'=>true, 'N'=>false), label="Mostra curvatura?")
        if Pt == HyperBKParticle
            scurv = spinbox(0.01:0.01:0.9, label = "Curvatura")
        elseif Pt == EllipticParticle
            scurv = spinbox(0.1:0.001:30., value=1., label = "Curvatura")
        else
            scurv = slider(0.01:0.1:30., label = "Curvatura")
        end
        map(rcurv, onchange(Ns), onchange(ts), bils, b) do r, N, t, bil, but
            Pp = r ? Pt : Particle
            colrs=gcolrs[1:N]
            bdo, bilwid = bil_ui(bil,onchange(scurv))
            map(bdo) do bd
                p = random_on_border(Pt,bd,N)
                bmap, arcs = parallelize(boundarymap, bd, t, p)
                fig.clear()
                ax0 = plt.subplot(get(gs,0))
                ax1 = plt.subplot(get(gs,1))
                plot_boundarymap((e1->(e2->convert.(Float64,e2)).(e1)).(bmap),convert.(Float64,arcs); color=colrs, ax=ax1)
                ax0.set_aspect("equal")
                plot(bd, Pp, ax=ax0)
                plot(p, Pp, colrs; ax=ax0)
                if r
                    xs, ys = timeseries!(p[1],bd,t)
                    ax0.plot(convert.(Float64,xs[2:end]),convert.(Float64,ys[2:end]))
                end
            end
            if Pt == Particle
                return bilwid
            else
                return vbox(rcurv, scurv,bilwid)
            end
        end
    end

    ui = vbox(hbox(Ns,ts),hbox(curv, bils, ph_plot),hline(),b)
    w = Window(Dict(:width=>600,:height=>200))
    body!(w, ui);
end

end