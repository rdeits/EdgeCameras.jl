polar_samples(θs, rs) = [CartesianFromPolar()(Polar(x[2], x[1])) for x in Iterators.product(θs, rs)]


function sample(cam::EdgeCamera, im, blur)
    pixels = Array{eltype(im)}(size(cam.params.samples))
    sample!(pixels, cam, im, blur)
    pixels
end

function sample!(pixels::AbstractArray{<:Colorant},
                 cam::EdgeCamera, im, blur)
    pixels .= sample_blurred.((im,), 
                              cam.source.homography.(cam.params.samples),
                              blur)
end

function sample_blurred(im::AbstractMatrix, 
                        pt::Union{Tuple, SVector}, 
                        blur::OffGridBlur{T}) where T
    y0, x0 = pt
    C = promote_type(eltype(im), RGB{T})
    ỹ = round(Int, y0)
    R = radius(blur)
    dy_range = max(-R, first(indices(im, 1)) - ỹ):min(R, last(indices(im, 1)) - ỹ)
    y_weights = weights(blur, ỹ - y0)
    weight_y::T = zero(T)
    for dy in dy_range
        weight_y += y_weights[dy + R + 1]
    end

    x̃ = round(Int, x0)
    dx_range = max(-R, first(indices(im, 2)) - x̃):min(R, last(indices(im, 2)) - x̃)
    x_weights = weights(blur, x̃ - x0)
    weight_x::T = zero(T)
    for dx in dx_range
        weight_x += x_weights[dx + R + 1]
    end

    total_x::C = zero(C)
    for dx in dx_range
        x = x̃ + dx
        total_y::C = zero(C)
        @inbounds for dy in dy_range
            y = ỹ + dy
            w_y = y_weights[dy + R + 1]
            total_y += w_y * im[y, x]
        end
        w_x = x_weights[dx + R + 1]
        total_x += (w_x / weight_y) * total_y
    end
    total_x / weight_x
end

sample_blurred(im::AbstractMatrix, pt::CartesianIndex, args...) =
    sample_blurred(im, pt.I, args...)
