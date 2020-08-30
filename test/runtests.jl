using GridapFSI
using Test

@testset "GridapFSI.jl" begin
  @testset "FSIDrivers.jl" begin include("FSIDriversTest.jl") end

  @testset "Analytical.jl" begin
    main(problemName="analytical",is_vtk=true)
    main(problemName="analytical",strategy="linearElasticity")
    #main(problemName="analytical",strategy="neoHookean")
    end
  @testset "ElasticFlag.jl" begin
    main()
    main(strategy="biharmonic")
    main(coupling="weak",strategy="biharmonic",is_vtk=true)
  end
end
