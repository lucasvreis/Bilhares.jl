function billiard_csq(s)
    o = typeof(s)(0)
    sp = [-s,s]; ep = [-s, -s]; n = [s,o]
    leftw = InfiniteWall(sp, ep, n, "Left wall")
    sp = [s,-s]; ep = [s, s]; n = [-s,o]
    rightw = InfiniteWall(sp, ep, n, "Right wall")
    sp = [s,s]; ep = [-s, s]; n = [o,-s]
    topw = InfiniteWall(sp, ep, n, "Top wall")
    sp = [-s,-s]; ep = [s, -s]; n = [o,s]
    botw = InfiniteWall(sp, ep, n, "Bottom wall")
    return Billiard(botw, rightw, topw, leftw)
end

function billiard_cmushroom(stem_length = 1.0, stem_width=0.2, cap_radius=1.0,
    stem_location = 0.0; sl = stem_length, sw = stem_width, cr = cap_radius,
    sloc = stem_location, door = false, scale=1.0)

    abs(sloc) + sw/2 > cr && error("Stem is outside the mushroom cap!")

    sl, sw, cr, sloc = scale .* promote(sl, sw, cr, sloc)
    T = eltype(sl)
    leftcorn = SV(-sw/2 + sloc, -sl)
    rightcorn = SV(sw/2 + sloc, -sl)
    upleftcorn = SV(-sw/2 + sloc, 0)
    uprightcorn = SV(sw/2 + sloc, 0)

    stembot = FiniteWall(leftcorn, rightcorn, SV(0, sw), door, "Stem bottom")
    stemleft = FiniteWall(upleftcorn, leftcorn, SV(sw, 0), false, "Stem left")
    stemright = FiniteWall(rightcorn, uprightcorn, SV(-sw, 0), false, "Stem right")

    farleft = SV(-cr, 0)
    farright = SV(cr, 0)

    capbotleft = FiniteWall(
    farleft, upleftcorn, SV(0, sw), false, "Cap bottom left")
    capbotright = FiniteWall(
    uprightcorn, farright, SV(0, sw), false, "Cap bottom right")

    cap = Semicircle([0, 0], cr, [0, -T(1.0)], "Mushroom cap")
    return Billiard(stembot, stemright, capbotright, cap, capbotleft, stemleft)
end

function billiard_crectangle(x′ = 1.0, y′ = 1.0;
    x = x′, y = y′, setting::String = "standard")

    x = convert(AbstractFloat, x)
    x, y = promote(x,y)
    ox, oy = -x/2, -y/2
    x, y = x/2, y/2
    if setting == "standard"
        sp = [ox,y]; ep = [ox, oy]; n = [1,0]
        leftw = InfiniteWall(sp, ep, n, "Left wall")
        sp = [x,oy]; ep = [x, y]; n = [-1,0]
        rightw = InfiniteWall(sp, ep, n, "Right wall")
        sp = [x,y]; ep = [ox, y]; n = [0,-1]
        topw = InfiniteWall(sp, ep, n, "Top wall")
        sp = [ox,oy]; ep = [x, oy]; n = [0,1]
        botw = InfiniteWall(sp, ep, n, "Bottom wall")
    elseif setting == "periodic"
        sp = [ox,y]; ep = [ox, oy]; n = [1,0]
        leftw = PeriodicWall(sp, ep, n, "Left periodic boundary")
        sp = [x,oy]; ep = [x, y]; n = [-1,0]
        rightw = PeriodicWall(sp, ep, n, "Right periodic boundary")
        sp = [x,y]; ep = [ox, y]; n = [0,-1]
        topw = PeriodicWall(sp, ep, n, "Top periodic boundary")
        sp = [ox,oy]; ep = [x, oy]; n = [0,1]
        botw = PeriodicWall(sp, ep, n, "Bottom periodic boundary")
    else
        throw(ArgumentError("The given setting=$setting is unknown."))
    end
    return Billiard(botw, rightw, topw, leftw)
end

"""
    billiard_sinai(r=0.25, x=1.0, y=1.0; setting = "standard")
Return a vector of obstacles that defines a Sinai billiard of size (`x`, `y`) with
a disk in its center, of radius `r`.

In the periodic case, the system is also known as "Lorentz Gas".

### Settings
* "standard" : Specular reflection occurs during collision.
* "periodic" : The walls are `PeriodicWall` type,
  enforcing periodicity at the boundaries
* "random" : The velocity is randomized upon collision.
* "ray-splitting" : All obstacles in the billiard allow for ray-splitting.
"""
function billiard_csinai(r′ = 0.25, x′ = 1.0, y′ = 1.0;
    r = r′, x = x′, y = y′, setting = "standard")

    if (setting == "periodic") && (r>=x/2 || r>=y/2)
        es = "Disk radius too big for a periodic Sinai billiard.\n"
        es*= "Obstacles must not overlap with `PeriodicWall`s."
        error(es)
    end
    r, x, y = promote(r,x,y)
    bdr = billiard_crectangle(x, y; setting = setting)

    centerdisk = Disk(typeof(x)[0.,0.], r, "Disk")

    return Billiard(centerdisk, bdr...)
end