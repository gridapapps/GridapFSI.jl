module GridapFSI

using Gridap

export main

# Generic struct for multiple dispatch
struct Problem{Kind} end

# Include problem-specific implementations
include("ElasticFlag.jl")


function main(problemName::String="elasticFlag",modelName::String="elasticFlag_coarse")

    # Discrete model
    model = DiscreteModelFromFile("../models/"*modelName*".json")
    model_solid = DiscreteModel(model,"solid")
    model_fluid = DiscreteModel(model,"fluid")

    # Setup problem
    problem = Problem{Symbol(problemName)}()
    params = get_parameters(problem)
    show(params)

end

end # module
