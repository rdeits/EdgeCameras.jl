using DataDeps

register(DataDep(
    "bouman_edge_cameras_data",
    "This dataset contains example videos provided with the original paper \"Turning Corners into Cameras\" by Bouman et al.",
    "https://people.csail.mit.edu/klbouman/pw/projects/cornercam/example_videos.zip",
    "ae7f1604dbf4ec3294a627bb2cd555dc5d6feef0762e37a0178329cac42e5258",
    post_fetch_method=unpack
))