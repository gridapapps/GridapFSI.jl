module GridapFSI
using GridapPETSc
using MPI

using TimerOutputs

export main

include("WeakForms/WeakForms.jl")
include("FSIDrivers.jl")

using GridapFSI.FSIDrivers

function main(;problemName::String="elasticFlag",kwargs...)

    MPI.Init()
    GridapPETSc.Init()

    reset_timer!()

    # Setup problem
    problem = Problem{Symbol(problemName)}()

    # Execute driver
    outputs = execute(problem; kwargs...)

    print_timer()
    println()

    GridapPETSc.Finalize()
    MPI.Finalize()

    return outputs

end

end # module
