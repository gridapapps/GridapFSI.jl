module GridapFSI

using TimerOutputs

export main

include("WeakForms/WeakForms.jl")
include("Drivers/Drivers.jl")

using GridapFSI.Drivers

function main(;problemName::String="elasticFlag",kwargs...)

    reset_timer!()

    # Setup problem
    problem = get_problem(problemName,kwargs)

    # Execute driver
    outputs = execute(problem; kwargs...)

    print_timer()
    println()

    return outputs

end

end # module
