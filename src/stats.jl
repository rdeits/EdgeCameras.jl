
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
    variance::Float32 = median(mean(variance_image, dims=1), dims=(2, 3))[1]
    VideoStatistics(colorview(RGB, mean_image),
                    colorview(RGB, variance_image),
                    variance)
end

