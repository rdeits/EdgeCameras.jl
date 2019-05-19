__precompile__(false)  # https://github.com/kmsquire/VideoIO.jl/pull/114

module EdgeCameras

export StaticSource,
       reconstruct,
       EdgeCamera,
       Params,
       VideoStatistics

using Colors: N0f8, RGB, Colorant, @colorant_str, gray, Gray
using Images: colorview, channelview
using VideoIO: VideoIO, VideoReader
using StaticArrays: SVector, SMatrix
using CoordinateTransformations: Transformation, CartesianFromPolar, Polar
using Unitful: Quantity, @u_str
using AxisArrays: AxisArray, Axis
using Interpolations: interpolate, BSpline, Linear, OnGrid, AbstractInterpolation
using LinearAlgebra
using Statistics

include("homography.jl")
include("stats.jl")
include("sampling.jl")
include("videos.jl")
include("algorithm.jl")
include("visualization.jl")

end
