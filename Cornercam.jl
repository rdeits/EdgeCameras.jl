module Cornercam

using Images
using VideoIO: openvideo, VideoReader
using StaticArrays: SVector, SMatrix
using CoordinateTransformations
using Unitful

function mean_frame(video::VideoReader, framerate, time_range)
    seek(video, convert(Float64, time_range[1] / (1u"s")))
    buf = read(video)
    total = RGB{Float32}.(copy(buf))
    num_frames = round(Int, (time_range[2] - time_range[1]) * framerate)
    for i in 2:num_frames
        read!(video, buf)
        total .+= buf
    end
    total ./ (num_frames)
end

struct Homography2D{T, M <: AbstractMatrix{T}} <: Transformation
  H::M
end

Homography2D(H::M) where {T, M <: AbstractMatrix{T}} = Homography2D{T, M}(H)

Base.inv(H::Homography2D) = Homography2D(inv(H.H))

@inline function (trans::Homography2D)(x)
    xprime = trans.H * SVector(x[1], x[2], 1)
    scale = 1 / xprime[end]
    SVector(xprime[1] * scale, xprime[2] * scale)
end

function rectify(original_corners, desired_corners)
    @assert length(original_corners) == length(desired_corners)
    N = length(original_corners)
    A = zeros(8, 8)
    for i in 1:N
        row = 2 * i - 1
        x1, y1 = original_corners[i]
        x2, y2 = desired_corners[i]
        A[row, :] = [x1, y1, 1, 0, 0, 0, -x1*x2, -y1*x2]
        A[row+1, :] = [0, 0, 0, x1, y1, 1, -x1*y2, -y1*y2]
    end
    b = vcat(desired_corners...)
    h = vcat(A \ b, 1)
    H = reshape(h, 3, 3)'
    Homography2D(SMatrix{3, 3}(H))
end

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

# function visibility_gain(frame_size, θs)
#     h, w = frame_size
#     M = h * w
#     N = length(θs)
#     A = zeros(N0f8, h, w, N)
#     for n in 1:N
#         θ = θs[n]
#         for j in 1:w
#             for i in 1:h
#                 ρ = atan2(j, i)
#                 A[i, j, n] = ρ >= θ
#             end
#         end
#     end
#     A
# end

# function estimate_scene(

end
