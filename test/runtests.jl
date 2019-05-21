using Test
using EdgeCameras
using Images
using LinearAlgebra
using Random

@testset "Edge Cameras" begin
    @testset "homographies" begin
        # Test homographies for image rectification.
        # Choose four corners in our warped image.
        original_corners = [[10, 22],
            [12, 148],
            [131, 155],
            [123, 5]]

        # For simplicity, we will transform the image
        # corners into a unit box
        desired_corners = [[0, 0],
            [0, 1],
            [1, 1],
            [1, 0]]

        # Compute a homography from desired to original:
        H1 = EdgeCameras.rectify(desired_corners, original_corners)
        # Verify that applying the homography gives the original coordinates
        @test all((H1.(desired_corners) .≈ original_corners))

        # Verify that the homography from original to desired
        # is just the inverse:
        H2 = EdgeCameras.rectify(original_corners, desired_corners)
        @test all(norm.(H2.(original_corners) .- desired_corners) .< 1e-15)
        @test H1.H ≈ inv(H2).H ./ (inv(H2).H[3, 3])
    end

    @testset "visibility gain" begin
        # Test the visibility gain matrix A for a toy example
        # (based on Fig. 2 in the paper).
        θs = range(0, stop=π/2, length=11)
        samples = EdgeCameras.polar_samples(θs, [1.0])
        A = EdgeCameras.visibility_gain(samples, θs)
        @test size(A, 1) == size(A, 2)
        @test A == LowerTriangular(ones(size(A, 1), size(A, 1)))
    end

    @testset "gradient matrix and regularizer" begin
        # Test that the efficient tri-diagonal computation of
        # the regularization term is correct.
        function G_reference(m)
            G = I - diagm(1 => ones(m - 1))
            G = G[1:end-1, :]
            G[1, :] .= 0
            G
        end

        function R_reference(m, σ)
            G = G_reference(m)
            c = Matrix(I, m, m)
            c[1, :] .= 0
            R = Tridiagonal((G' * G + c)) / σ^2
        end

        σ = 0.1
        for m in 5:10
            G′G = EdgeCameras.G_times_G(m)
            G_expected = G_reference(m)
            @test G′G == G_expected' * G_expected

            R = EdgeCameras.regularizer(m, σ)
            @test R == R_reference(m, σ)
        end
    end

    @testset "interpolated exponential" begin
        # Test that our linear interpolation for computing off-grid
        # gaussian weights is a close fit to the true values
        σ = 3.0
        radius = 6
        blur = EdgeCameras.OffGridBlur(σ, radius)
        for offset in range(-0.5, stop=0.5, length=11)
            @test isapprox(EdgeCameras.weights(blur, offset),
                           EdgeCameras.weights(σ, radius, offset),
                           rtol=1e-3,
                           atol=1e-2)
        end
    end

    @testset "blurred_sample" begin
        # Test that the lazy blurred samples produce the same
        # result as blurring the image and then sampling (at least
        # when on-grid).
        Random.seed!(1)

        # Generate a random image and blur it
        im = rand(RGB{Float32}, 100, 100)
        σ = 3.0
        kernel = Kernel.gaussian(σ)
        radius = (length(axes(kernel, 1)) - 1) ÷ 2
        filtered = imfilter(im, kernel)

        # Now generate what should be an identical image by
        # generating a lazy blurred sample at *every* coordinate:
        blur = EdgeCameras.OffGridBlur(σ, radius)
        sampled = EdgeCameras.sample_blurred.(
            (im,),
            collect(CartesianIndices((100, 100))),
            (blur,))

        # Verify that the images are the same (away from the boundary)
        @test maximum(abs.(sampled[10:90, 10:90] .- filtered[10:90, 10:90])) < 1e-6
    end
end

include("notebooks.jl")
