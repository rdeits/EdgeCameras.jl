polar_samples(θs, rs) = [CartesianFromPolar()(Polar(x[2], x[1])) for x in Iterators.product(θs, rs)]


function sample(cam::CornerCamera, im, blur)
    pixels = Array{eltype(im)}(size(cam.params.samples))
    sample!(pixels, cam, im, blur)
    pixels
end

function sample!(pixels::AbstractArray{<:Colorant},
                 cam::CornerCamera, im, blur)
    pixels .= sample_blurred.((im,), 
                              cam.source.homography.(cam.params.samples),
                              blur)
end

@generated function weights(blur::OffGridBlur{R, T}, y0) where {R, T}
    expr = Expr(:tuple, [:(exp(scale * ($dy + y_grid_offset)^2)) for dy in -R:R]...)
    quote
        scale = -1/(2 * blur.σ^2)
        ỹ = round(Int, y0)
        y_grid_offset = ỹ - y0
        $expr
    end
end

function sample_blurred(im::AbstractMatrix, 
                        pt::Union{Tuple, SVector}, 
                        blur::OffGridBlur{R, T}) where {R, T}
    y0, x0 = pt
    C = promote_type(eltype(im), RGB{T})
    weight_x::T = zero(T)
    total_x::C = zero(C)
    scale = -1/(2 * blur.σ^2)
    ỹ = round(Int, y0)
    x̃ = round(Int, x0)
    dx_range = max(-R, first(indices(im, 2)) - x̃):min(R, last(indices(im, 2)) - x̃)
    dy_range = max(-R, first(indices(im, 1)) - ỹ):min(R, last(indices(im, 1)) - ỹ)
    y_weights = weights(blur, y0)
    for dx in dx_range
        x = x̃ + dx
        weight_y::T = zero(T)
        total_y::C = zero(C)
        for dy in dy_range
            y = ỹ + dy
            w_y = y_weights[dy + R + 1]
            weight_y += w_y
            total_y += w_y * im[y, x]
        end
        w_x = exp(scale * (x - x0)^2)
        weight_x += w_x
        total_x += w_x * total_y / weight_y
    end
    total_x / weight_x
end

sample_blurred(im::AbstractMatrix, pt::CartesianIndex, args...) =
    sample_blurred(im, pt.I, args...)
