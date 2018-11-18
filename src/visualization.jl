function desaturate(px::RGB)
    g = gray(Gray(px))
    RGB(g, g, g)
end

function mark!(im::AbstractArray, center, radius::Integer=3, color=RGB(1., 0, 0))
    for dx in -radius:radius
        x = round(Int, center[2] + dx)
        if x ∈ axes(im, 2)
            for dy in -radius:radius
                y = round(Int, center[1] + dy)
                if y ∈ axes(im, 1)
                    im[y, x] = color
                end
            end
        end
    end
end

Base.show(io::IO, mime::MIME"image/png", stats::VideoStatistics) = show(io, mime, stats.mean_image)

function Base.show(io::IO, mime::MIME"image/png", source::StaticSource)
    im = copy(background(source))
    colors = [colorant"black", colorant"red", colorant"green", colorant"blue"]
    for (corner, color) in zip(DEFAULT_CORNERS, colors)
        mark!(im, source.homography(corner), 10, color)
    end
    show(io, mime, im)
end

function show_samples(c::EdgeCamera)
    im = desaturate.(background(c.source))
    cmap = x -> RGB(1 - x, 0.0, x)
    for (i, s) in enumerate(c.params.samples)
        pixel = c.source.homography(s)
        mark!(im, pixel, 1, cmap(i / length(c.params.samples)))
    end
    im
end

function Base.show(io::IO, mime::MIME"image/png", cam::EdgeCamera)
    show(io, mime, show_samples(cam))
end
