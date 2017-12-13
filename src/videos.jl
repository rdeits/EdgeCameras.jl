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


framerate(v::VideoReader) = v.stream_info.stream.avg_frame_rate.num / (v.stream_info.stream.avg_frame_rate.den * u"s")

const FRAME_TYPE = PermutedDimsArray{RGB{N0f8}, 2, (2, 1), (2, 1), Matrix{RGB{N0f8}}}

function enumerateframes(f::Function, video::VideoReader, time_range::Tuple, target_rate::Quantity)
    times = frame_times(video, time_range, target_rate)
    frames_per_sample = max(1, round(Int, framerate(video) / target_rate))
    buf::FRAME_TYPE = read(video)
    seek(video, convert(Float64, first(times) / u"s"))
    for (i, current_time) in enumerate(times)
        f(i, buf)
        skip_frames!(video, frames_per_sample - 1)
        read!(video, buf)
    end
    times
end

eachframe(f::Function, video::VideoReader, time_range::Tuple, target_rate::Quantity) = 
    enumerateframes((i, buf) -> f(buf), video, time_range, target_rate)

function frame_times(video::VideoReader, time_range::Tuple{Quantity, Quantity}, target_rate::Quantity=framerate(video))
    frames_per_sample = max(1, round(Int, framerate(video) / target_rate))
    closest_achievable_rate = framerate(video) / frames_per_sample
    time_range[1]:(1/closest_achievable_rate):time_range[2]
end

function skip_frame!(video::VideoReader)
    while !VideoIO.have_frame(video)
        idx = VideoIO.pump(video.avin)
        idx == video.stream_index0 && break
        idx == -1 && throw(EOFError())
    end
    VideoIO.reset_frame_flag!(video)
end

function skip_frames!(video::VideoReader, N::Integer)
    for i in 1:N
        skip_frame!(video)
    end
end
