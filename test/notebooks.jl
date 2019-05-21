using Test
using NBInclude

notebookdir = joinpath(@__DIR__, "..", "notebooks")
for file in ["demo.ipynb", "red_stereo.ipynb"]
    @testset "Notebook: $file" begin
        @nbinclude(joinpath(notebookdir, file))
    end
end
