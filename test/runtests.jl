using GridapFSI
using Test

@testset "GridapFSI.jl" begin
    @testset "Analytical.jl" begin
      main(problemName="analytical")
    end
    @testset "ElasticFlag.jl" begin
      main()
    end
end
