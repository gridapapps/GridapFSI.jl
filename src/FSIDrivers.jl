module FSIDrivers

using Gridap
using Gridap.Helpers
using Gridap.Geometry
using Gridap.Arrays
using Gridap.MultiField: ConsecutiveMultiFieldStyle
using GridapODEs.ODETools
using GridapODEs.TransientFETools
using GridapFSI.WeakForms
using TimerOutputs
using WriteVTK
using LineSearches: BackTracking, HagerZhang
using ForwardDiff
using Test

import GridapODEs.TransientFETools: ∂t

export Problem
export execute

struct Problem{Kind} end

include("ElasticFlag.jl")
include("Analytical.jl")
include("Oscillator.jl")

execute(problem::Problem; kwargs...) = @notimplemented("The driver for problem: $problem is not implemented")

# Output function
function writePVD(filePath::String, trian::Triangulation, sol; append=false)
    outfiles = paraview_collection(filePath, append=append) do pvd
        for (i, (xh, t)) in enumerate(sol)
            println("STEP: $i, TIME: $t")
            println("============================")
            uh = restrict(xh[1],trian)
            vh = restrict(xh[2],trian)
						ph = restrict(xh[3],trian)
            pvd[t] = createvtk(
                trian,
                filePath * "_$t.vtu",
                cellfields = ["uh" => uh, "vh" => vh, "ph" => ph]
            )
        end
    end
end

function _get_kwarg(kwarg,kwargs)
    try
        return kwargs[kwarg]
    catch
        s = "The key-word argument $(kwarg) is mandatory in the $problem driver"
        error(s)
    end
end

function _get_kwarg(kwarg,kwargs,value)
    try
        return kwargs[kwarg]
    catch
        return value
    end
end

end
