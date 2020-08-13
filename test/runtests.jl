using GridapFSI
using Test

@testset "GridapFSI.jl" begin
  @testset "FSIDrivers.jl" begin include("FSIDriversTest.jl") end
  @testset "Analytical.jl" begin
    main(problemName="analytical",is_vtk=true)
  end
  @testset "ElasticFlag.jl" begin
    main()
  end
end
