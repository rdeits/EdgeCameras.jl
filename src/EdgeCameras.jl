module EdgeCameras

export StaticSource,
       reconstruct,
       EdgeCamera

import Colors: N0f8, RGB, Colorant, @colorant_str, gray, Gray
import Images: colorview, channelview
import VideoIO
import VideoIO: VideoReader
import StaticArrays: SVector, SMatrix
import CoordinateTransformations: Transformation, CartesianFromPolar, Polar
import Unitful: Quantity, @u_str
import AxisArrays: AxisArray, Axis
import Interpolations: interpolate, BSpline, Linear, OnGrid, AbstractInterpolation

include("homography.jl")
include("stats.jl")
include("sampling.jl")
include("videos.jl")
include("algorithm.jl")
include("visualization.jl")

end
