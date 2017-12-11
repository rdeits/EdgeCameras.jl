framerate(v::VideoReader) = v.stream_info.stream.avg_frame_rate.num / (v.stream_info.stream.avg_frame_rate.den * u"s")

function frames_between(video::VideoReader, time_range::Tuple)
    rate = framerate(video)
    seek(video, convert(Float64, time_range[1] / (1u"s")))
    num_frames = round(Int, (time_range[2] - time_range[1]) * rate)
    return (read(video) for i in 1:num_frames)
end

background(video::VideoReader, time_range::Tuple) = 
    mean((RGB{Float32}.(im) for im in frames_between(video, time_range)))
