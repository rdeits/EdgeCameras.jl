using Base.Test
using CornerCameras

@testset "Corner Cameras" begin
    @testset "homographies" begin
        original_corners = [[10, 22], 
            [12, 148], 
            [131, 155], 
            [123, 5]]

        desired_corners = [[0, 0], 
            [0, 1],
            [1, 1],
            [1, 0]]

        H1 = rectify(desired_corners, original_corners)
        @test all((H1.(desired_corners) .≈ original_corners))

        H2 = rectify(original_corners, desired_corners)
        @test all(norm.(H2.(original_corners) .- desired_corners) .< 1e-15)

        @test H1.H ≈ inv(H2).H ./ (inv(H2).H[3, 3])
    end

    @testset "visibility gain" begin
        θs = linspace(0, π/2, 11)
        samples = CornerCameras.polar_samples(θs, [1.0])
        A = CornerCameras.visibility_gain(samples, θs)
        @test size(A, 1) == size(A, 2)
        @test A == LowerTriangular(ones(size(A, 1), size(A, 1)))
    end

    @testset "gradient matrix and regularizer" begin
        function G_reference(m)
            G = I - diagm(ones(m - 1), 1)
            G = G[1:end-1, :]
            G[1, :] = 0
            G
        end

        function R_reference(m, σ)
            G = G_reference(m)
            c = eye(m)
            c[1, :] = 0
            R = Tridiagonal((G' * G + c)) / σ^2
        end

        σ = 0.1
        for m in 5:10
            G′G = CornerCameras.G_times_G(m)
            G_expected = G_reference(m)
            @test G′G == G_expected' * G_expected

            R = CornerCameras.regularizer(m, σ)
            @test R == R_reference(m, σ)
        end
    end
end
