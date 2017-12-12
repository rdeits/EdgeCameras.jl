struct Homography2D{T, M <: AbstractMatrix{T}} <: Transformation
  H::M
end

Homography2D(H::M) where {T, M <: AbstractMatrix{T}} = Homography2D{T, M}(H)

function weights(σ, radius, grid_offset)
    scale = -1/(2 * σ^2)
    SVector{2 * radius + 1}([exp(scale * (dy + grid_offset)^2) for dy in -radius:radius]) 
end

struct OffGridBlur{T, I <: AbstractInterpolation}
    σ::T
    radius::Int
    weights::I 

    function OffGridBlur(σ::T, radius::Integer, samples=3) where {T}
        all_weights = [
            weights(σ, radius, convert(T, grid_offset)) for grid_offset in linspace(-0.5, 0.5, samples)
        ]
        itp = interpolate(all_weights, BSpline(Linear()), OnGrid())
        new{T, typeof(itp)}(σ, radius, itp)
    end
end

weights(blur::OffGridBlur, grid_offset) = blur.weights[2 * grid_offset + 2]
radius(blur::OffGridBlur) = blur.radius

@with_kw struct Params{Tθ, TS, B <: OffGridBlur}
    θs::Tθ = linspace(0, π/2, 200)
    samples::TS = polar_samples(linspace(0, π/2, 200), linspace(0.01, 1, 50))
    σ::Float64 = 0.085
    blur::B = OffGridBlur(3.0f0, 6)
end

struct StaticSource{R <: VideoReader, M <: AbstractMatrix{<:Colorant}, H <: Homography2D}
    video::R
    background::M
    homography::H
    λ::Float64
end

function StaticSource(video::R, 
                      corners_in_image::AbstractVector, 
                      im::M, 
                      λ=sqrt(2.7)) where {R <: VideoReader, M <: AbstractMatrix{<:Colorant}}
    H = inv(rectify(corners_in_image))
    StaticSource{R, M, typeof(H)}(video, im, H, λ)
end

struct EdgeCamera{S <: StaticSource, P <: Params}
    source::S
    params::P
    gain::Matrix{Float64}
end

function EdgeCamera(source::S, params::P) where {S <: StaticSource, P <: Params}
    A = visibility_gain(params.samples, params.θs)
    gain = edge_cam_gain(A, params.σ, source.λ)
    EdgeCamera{S, P}(source, params, gain)
end

