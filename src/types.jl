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

struct CornerCamera{S <: VideoReader, H <: Homography2D, M <: AbstractArray, P <: Params}
    source::S
    homography::H
    background::M
    params::P
    gain::Matrix{Float64}
end

function CornerCamera(source::S, homography::H, background::M, params::P, λ) where {S, H, M, P}
    A = visibility_gain(params.samples, params.θs)
    gain = cornercam_gain(A, params.σ, λ)
    CornerCamera{S, H, M, P}(source, homography, background, params, gain)
end

