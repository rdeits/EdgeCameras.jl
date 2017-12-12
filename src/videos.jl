framerate(v::VideoReader) = v.stream_info.stream.avg_frame_rate.num / (v.stream_info.stream.avg_frame_rate.den * u"s")

function frames_between(video::VideoReader, time_range::Tuple)
    rate = framerate(video)
    seek(video, convert(Float64, time_range[1] / (1u"s")))
    num_frames = round(Int, (time_range[2] - time_range[1]) * rate)
    return (read(video) for i in 1:num_frames)
end

const FRAME_TYPE = PermutedDimsArray{RGB{N0f8}, 2, (2, 1), (2, 1), Matrix{RGB{N0f8}}}

function enumerateframes(f::Function, video::VideoReader, times::AbstractVector{<:Quantity})
    buf::FRAME_TYPE = read(video)
    for (i, current_time) in enumerate(times)
        seek(video, convert(Float64, current_time / (1u"s")))
        read!(video, buf)
        f(i, buf)
    end
end

eachframe(f::Function, video::VideoReader, times::AbstractVector{<:Quantity}) = 
    enumerateframes((i, buf) -> f(buf), video, times)

frame_times(time_range::Tuple{Quantity, Quantity}, rate::Quantity) = 
    time_range[1]:(1/rate):time_range[2]
