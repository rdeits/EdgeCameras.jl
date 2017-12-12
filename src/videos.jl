framerate(v::VideoReader) = v.stream_info.stream.avg_frame_rate.num / (v.stream_info.stream.avg_frame_rate.den * u"s")

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

# struct FrameIterator{V <: VideoReader, Q1 <: Quantity, Q2 <: Quantity}
#     video::V
#     time_range::NTuple{2, Q1}
#     target_rate::Q2
# end

# function frames_between(video::VideoReader, time_range::Tuple, target_rate::Quantity=framerate(video))
#     FrameIterator(video, time_range, target_rate)
# end

# Base.start(it::FrameIterator) = (read(it.video)::FRAME_TYPE, it.time_range[1])
# function Base.done(it::FrameIterator, state)
#     buf, time = state
#     time > it.time_range[2]
# end
# function Base.next(it::FrameIterator, state)
#     buf, time = state
#     time += 1 / it.target_rate
#     seek(it.video, convert(Float64, time / (1u"s")))
#     read!(it.video, buf)
#     (buf, (buf, time))
# end
# Base.length(it::FrameIterator) = 1 + floor(Int, (it.time_range[2] - it.time_range[1]) * it.target_rate)
