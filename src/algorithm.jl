struct Params{Tθ, TS, Tσ, B <: OffGridBlur}
    θs::Tθ
    samples::TS
    σ::Tσ
    blur::B
end

Params(;
       θs = range(0, stop=π/2, length=200),
       samples = polar_samples(range(0, stop=π/2, length=200), range(0.01, stop=1, length=50)),
       σ = 0.00033f0, # 0.085 in the original paper, scaled down by factor of 255,
       blur = OffGridBlur(3.0f0, 6)) =
    Params(θs, samples, σ, blur)

polar_samples(θs, rs) = [CartesianFromPolar()(Polar(x[2], x[1])) for x in Iterators.product(θs, rs)]

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

function sample(cam::EdgeCamera, im, blur)
    pixels = Array{eltype(im)}(undef, size(cam.params.samples))
    sample!(pixels, cam, im, blur)
    pixels
end

function sample!(pixels::AbstractArray{<:Colorant},
                 cam::EdgeCamera, im, blur)
    pixels .= sample_blurred.(Ref(im),
                              cam.source.homography.(cam.params.samples),
                              Ref(blur))
end

function visibility_gain(samples, θs)
    N = length(θs)
    M = length(samples)
    A = zeros(Float32, M, N)
    for (j, θ) in enumerate(θs)
        for (i, sample) in enumerate(samples)
            ρ = atan(sample[2], sample[1])
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
    background_samples = sample(cam, background(cam.source), cam.params.blur)
    pixels = copy(background_samples)
    times = frame_times(cam.source.video, time_range, target_rate)
    data = zeros(RGB{Float32}, length(cam.params.θs), length(times))
    enumerateframes(cam.source.video, time_range, target_rate) do i, buffer
        sample!(pixels, cam, buffer, cam.params.blur)
        pixels .-= background_samples
        for k in 1:length(pixels)
            for j in 1:size(data, 1)
                data[j, i] += cam.gain[j + 1, k] * pixels[k]
            end
        end
    end
    AxisArray(data,
              Axis{:θ}(cam.params.θs),
              Axis{:time}(times))
end
