using Base.Test
push!(LOAD_PATH, ".")

import Cornercam

@testset "homographies" begin
    original_corners = [[10, 22], 
        [12, 148], 
        [131, 155], 
        [123, 5]]

    desired_corners = [[0, 0], 
        [0, 1],
        [1, 1],
        [1, 0]]

    H1 = Cornercam.rectify(desired_corners, original_corners)
    @test all((H1.(desired_corners) .≈ original_corners))

    H2 = Cornercam.rectify(original_corners, desired_corners)
    @test all(norm.(H2.(original_corners) .- desired_corners) .< 1e-15)

    @test H1.H ≈ inv(H2).H ./ (inv(H2).H[3, 3])
end
