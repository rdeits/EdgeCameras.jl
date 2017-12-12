function visibility_gain(samples, θs)
    N = length(θs)
    M = length(samples)
    A = zeros(Float32, M, N)
    for (j, θ) in enumerate(θs)
        for (i, sample) in enumerate(samples)
            ρ = atan2(sample[2], sample[1])
            A[i, j] = ρ >= θ
        end
    end
    A
end

# Rather than computing G and then G' * G separately, we
# can just directly compute G'G as a Tridiagonal matrix. 
# See the tests for verification of this approach.
function G_times_G(m::Integer)
    off_diag = vcat(0.0, fill(-1.0, m - 2))
    Tridiagonal(off_diag,
                vcat(0.0, 1.0, fill(2.0, m - 3), 1.0),
                off_diag)
end

function regularizer(m::Integer, σ)
    G′G = G_times_G(m)
    c = Diagonal(vcat(0.0, fill(1.0, m - 1)))
    R = (G′G + c) / σ^2
end

function edge_cam_gain(A::AbstractMatrix, σ, λ)
    Ã = hcat(ones(size(A, 1)), A)
    R = regularizer(size(Ã, 2), σ)
    Σinv = Ã' * Ã / λ^2 + R
    gain = (Σinv \ (Ã' / λ^2));
end

function reconstruct(cam::EdgeCamera, time_range::Tuple, target_rate=framerate(cam.source))
    seek(cam.source.video, convert(Float64, time_range[1] / (1u"s")))

    background_samples = sample(cam, cam.source.background, cam.params.blur)
    pixels = copy(background_samples)
    frame = read(cam.source.video)

    frame_skip = max(1, round(Int, framerate(cam.source.video) / target_rate))
    closest_achievable_rate = framerate(cam.source.video) / frame_skip
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
            read!(cam.source.video, frame)
        end
    end
    times = time_range[1]:(1/closest_achievable_rate):time_range[2]
    AxisArray(data, 
              Axis{:θ}(cam.params.θs), 
              Axis{:time}(times))
end
