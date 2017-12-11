module CornerCameras

export show_samples,
       rectify,
       background,
       StaticSource,
       reconstruct

using Colors
import Images
using VideoIO: VideoReader
using StaticArrays: SVector, SMatrix
using CoordinateTransformations
using Unitful
using Parameters: @with_kw
using AxisArrays
using Interpolations: interpolate, BSpline, Linear, OnGrid, AbstractInterpolation

include("types.jl")
include("homography.jl")
include("videos.jl")
include("sampling.jl")
include("algorithm.jl")
include("visualization.jl")

end
