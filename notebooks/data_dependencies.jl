using DataDeps

register(DataDep(
    "bouman_edge_cameras_data",
    "This dataset contains example videos provided with the original paper \"Turning Corners into Cameras\" by Bouman et al.",
    "https://people.csail.mit.edu/klbouman/pw/projects/cornercam/example_videos.zip",
    [],
    post_fetch_method=unpack
))