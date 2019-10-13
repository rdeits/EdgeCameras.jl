# EdgeCameras.jl: Edge Cameras in Julia

[![Build Status](https://travis-ci.org/rdeits/EdgeCameras.jl.svg?branch=master)](https://travis-ci.org/rdeits/EdgeCameras.jl)
[![codecov.io](https://codecov.io/github/rdeits/EdgeCameras.jl/coverage.svg?branch=master)](https://codecov.io/github/rdeits/EdgeCameras.jl?branch=master)

This project consists of an implementation of *edge cameras* based on the work of Bouman et al. [1]. An edge camera is formed when a sharp edge (such as the corner of a wall) creates a natural one-dimensonal pinhole camera, revealing the motions of objects which are completely obscured by the corner. More information from the original authors can be found at [people.csail.mit.edu](https://people.csail.mit.edu/klbouman/cornercameras.html).

This package consists of an entirely new implementation of the edge camera algorithm, based on the work presented in the paper, and done entirely in [Julia](https://julialang.org/) (except for the raw video I/O, which is ultimately handled by `ffmpeg`).

[1] Katherine L. Bouman, Vickie Ye, Gregory W. Wornell, Adam B. Yedidia, Antonio Torralba, William T. Freeman, and Frédo Durand. "Turning Corners into Cameras: Principles and Methods". ICCV 2017.

# Installation

You'll need to install Julia 1.0 or newer from <https://julialang.org/downloads/>

To install this package in Julia, you can do the following:

```julia
julia> using Pkg

julia> Pkg.add("EdgeCameras")
```

# Usage

* [Basic demo](notebooks/demo.ipynb)
* [Stereo reconstruction](notebooks/red_stereo.ipynb)
