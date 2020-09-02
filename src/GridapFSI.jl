module GridapFSI

using TimerOutputs

export main

include("WeakForms/WeakForms.jl")
include("FSIDrivers.jl")

using GridapFSI.FSIDrivers

function main(;problemName::String="elasticFlag",kwargs...)

    reset_timer!()

    # Setup problem
    problem = Problem{Symbol(problemName)}()

    # Execute driver
    outputs = execute(problem; kwargs...)

    print_timer()
    println()

    return outputs

end

end # module
