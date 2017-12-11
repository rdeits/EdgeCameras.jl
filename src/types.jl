struct Homography2D{T, M <: AbstractMatrix{T}} <: Transformation
  H::M
end

Homography2D(H::M) where {T, M <: AbstractMatrix{T}} = Homography2D{T, M}(H)

struct OffGridBlur{R, T}
    σ::T
    OffGridBlur{R}(σ::T) where {R, T} = new{R, T}(σ)
end

@with_kw struct Params{Tθ, TS, B <: OffGridBlur}
    θs::Tθ = linspace(0, π/2, 200)
    samples::TS = polar_samples(linspace(0, π/2, 200), linspace(0.01, 1, 50))
    σ::Float64 = 0.085
    blur::B = OffGridBlur{6}(3.0f0)
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

struct CornerCamera{S <: StaticSource, P <: Params}
    source::S
    params::P
    gain::Matrix{Float64}
end

function CornerCamera(source::S, params::P) where {S <: StaticSource, P <: Params}
    A = visibility_gain(params.samples, params.θs)
    gain = cornercam_gain(A, params.σ, source.λ)
    CornerCamera{S, P}(source, params, gain)
end

