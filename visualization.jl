function desaturate(px::RGB)
    g = gray(Gray(px))
    RGB(g, g, g)
end

function mark!(im::AbstractArray, center, radius::Integer=3, color=RGB(1., 0, 0))
    for dx in -radius:radius
        for dy in -radius:radius
            p = round.(Int, center .+ SVector(dx, dy))
            im[Tuple(p)...] = color
        end
    end
end

function imnormal(im)
    lb = minimum(x -> min(red(x), green(x), blue(x)), im)
    ub = maximum(x -> max(red(x), green(x), blue(x)), im)
    (im .- RGB(lb, lb, lb)) ./ (ub - lb)
end

function show_samples(c::CornerCamera)
    im = desaturate.(c.background)
    color = interpolate([RGB(1., 0, 0), RGB(0., 0, 1)], BSpline(Linear()), OnGrid())
    for (i, s) in enumerate(c.params.samples)
        pixel = c.homography(s)
        mark!(im, pixel, 1, color[1 + i / length(c.params.samples)])
    end
    im
end

