module Cornercam

export show_samples,
       imnormal,
       rectify

using ColorTypes
using FixedPointNumbers: N0f8
using VideoIO: VideoReader
using StaticArrays: SVector, SMatrix
using CoordinateTransformations
using Interpolations
using Unitful
using Parameters: @with_kw
using AxisArrays

include("types.jl")
include("visualization.jl")
include("homography.jl")
include("videos.jl")
include("sampling.jl")


function visibility_gain(samples, θs)
    N = length(θs)
    M = length(samples)
    A = zeros(N0f8, M, N)
    for (j, θ) in enumerate(θs)
        for (i, sample) in enumerate(samples)
            ρ = atan2(sample[2], sample[1])
            A[i, j] = ρ > θ
        end
    end
    A
end

function cornercam_gain(A::AbstractMatrix, σ, λ)
    Ã = hcat(ones(size(A, 1)), A)

    G = I - diagm(ones(size(Ã, 2) - 1), 1)
    G = G[1:end-1, :]
    G[1, :] = 0

    c = eye(size(G, 2))
    c[1, :] = 0

    R = Tridiagonal((G' * G + c)) / σ^2

    Σinv = Ã' * Ã / λ^2 + R
    gain = (Σinv \ (Ã' / λ^2));
end

function trace(cam::CornerCamera, time_range::Tuple, target_rate=framerate(cam.source))
    seek(cam.source, convert(Float64, time_range[1] / (1u"s")))

    background_samples = sample(cam, cam.background, cam.params.blur)
    pixels = copy(background_samples)
    frame = read(cam.source)

    frame_skip = max(1, round(Int, framerate(cam.source) / target_rate))
    closest_achievable_rate = framerate(cam.source) / frame_skip
    num_frames = round(Int, (time_range[2] - time_range[1]) * closest_achievable_rate)
    data = zeros(RGB{Float32}, length(cam.params.θs), num_frames)
    for i in 1:num_frames
        sample!(pixels, cam, frame, cam.params.blur)
        pixels .-= background_samples
        for k in 1:length(pixels)
            for j in 1:size(data, 1)
                data[j, i] += cam.gain[j + 1, k] * pixels[k]
            end
        end
        for j in 1:frame_skip
            read!(cam.source, frame)
        end
    end
    times = time_range[1]:(1/closest_achievable_rate):time_range[2]
    AxisArray(data, 
              Axis{:θ}(cam.params.θs), 
              Axis{:time}(times))
end


end
