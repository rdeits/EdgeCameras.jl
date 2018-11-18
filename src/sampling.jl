struct OffGridBlur{T, I <: AbstractInterpolation}
    σ::T
    radius::Int
    weights::I

    function OffGridBlur(σ::T, radius::Integer, samples=3) where {T}
        all_weights = [
            weights(σ, radius, convert(T, grid_offset)) for grid_offset in range(-0.5, stop=0.5, length=samples)
        ]
        itp = interpolate(all_weights, BSpline(Linear()))
        new{T, typeof(itp)}(σ, radius, itp)
    end
end

weights(blur::OffGridBlur, grid_offset) = blur.weights(2 * grid_offset + 2)
radius(blur::OffGridBlur) = blur.radius

function weights(σ, radius, grid_offset)
    scale = -1/(2 * σ^2)
    SVector{2 * radius + 1}([exp(scale * (dy + grid_offset)^2) for dy in -radius:radius])
end

"""
Sample from image `im` at non-integer point `pt`, using
a Gaussian blur described by `blur` to avoid aliasing.
"""
function sample_blurred(im::AbstractMatrix,
                        pt::Union{Tuple, SVector},
                        blur::OffGridBlur{T}) where T
    y0, x0 = pt
    C = promote_type(eltype(im), RGB{T})
    R = radius(blur)

    ỹ = round(Int, y0)
    dy_range = max(-R, first(axes(im, 1)) - ỹ):min(R, last(axes(im, 1)) - ỹ)
    y_weights = weights(blur, ỹ - y0)
    weight_y::T = zero(T)
    for dy in dy_range
        weight_y += y_weights[dy + R + 1]
    end

    x̃ = round(Int, x0)
    dx_range = max(-R, first(axes(im, 2)) - x̃):min(R, last(axes(im, 2)) - x̃)
    x_weights = weights(blur, x̃ - x0)
    weight_x::T = zero(T)
    for dx in dx_range
        weight_x += x_weights[dx + R + 1]
    end

    total_x::C = zero(C)
    for dx in dx_range
        total_y::C = zero(C)
        @inbounds for dy in dy_range
            w_y = y_weights[dy + R + 1]
            total_y += w_y * im[ỹ + dy, x̃ + dx]
        end
        w_x = x_weights[dx + R + 1]
        total_x += (w_x / weight_y) * total_y
    end
    total_x / weight_x
end

sample_blurred(im::AbstractMatrix, pt::CartesianIndex, args...) =
    sample_blurred(im, pt.I, args...)
