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
    σ::Float64 = 0.00033 # 0.085 in the original paper, scaled down by factor of 255
    blur::B = OffGridBlur(3.0f0, 6)
end

struct VideoStatistics{T, M <: AbstractArray{<:Colorant}}
    mean_image::M
    variance_image::M
    variance::T
end

function VideoStatistics(video::VideoReader, 
                         time_range::Tuple, 
                         target_rate::Quantity=framerate(video))
    mean_image = zeros(Float32, 3, video.height, video.width)
    variance_image = copy(mean_image) # initially accumulats sum of squared values
    times = eachframe(video, time_range, target_rate) do buffer
        mean_image .+= Float32.(channelview(buffer))
        variance_image .+= Float32.(channelview(buffer)) .^ 2
    end
    num_frames = length(times)
    mean_image ./= num_frames
    variance_image .= variance_image ./ num_frames .- mean_image .^ 2
    variance::Float32 = median(mean(variance_image, 1), (2, 3))[1]
    VideoStatistics(colorview(RGB, mean_image), 
                    colorview(RGB, variance_image), 
                    variance)
end

struct StaticSource{R <: VideoReader, H <: Homography2D, S <: VideoStatistics}
    video::R
    homography::H
    stats::S
end

function StaticSource(video::R, 
                      corners_in_image::AbstractVector, 
                      stats::S) where {R <: VideoReader, S <: VideoStatistics}
    H = inv(rectify(corners_in_image))
    StaticSource{R, typeof(H), S}(video, H, stats)
end

function StaticSource(video::VideoReader,
                      corners_in_image::AbstractVector,
                      time_range::Tuple,
                      target_rate::Quantity=framerate(video))
    stats = VideoStatistics(video, time_range, target_rate)
    StaticSource(video, corners_in_image, stats)
end

background(source::StaticSource) = source.stats.mean_image

struct EdgeCamera{S <: StaticSource, P <: Params}
    source::S
    params::P
    gain::Matrix{Float64}
end

function EdgeCamera(source::S, params::P=Params()) where {S <: StaticSource, P <: Params}
    A = visibility_gain(params.samples, params.θs)
    gain = edge_cam_gain(A, params.σ, sqrt(source.stats.variance))
    EdgeCamera{S, P}(source, params, gain)
end
