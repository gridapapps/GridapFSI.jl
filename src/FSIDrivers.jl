module FSIDrivers

using Gridap
using Gridap.Helpers
using Gridap.Geometry
using GridapODEs.ODETools
using GridapODEs.TransientFETools
using GridapFSI.WeakForms
using TimerOutputs
using WriteVTK
using LineSearches: BackTracking, HagerZhang

import GridapODEs.TransientFETools: ∂t

export Problem
export execute

struct Problem{Kind} end

### This is temporary. Need a better solution!!! ###
global u_in(x, t) = VectorValue(1.5 * 1.0 * x[2] * (0.41 - x[2]) / ((0.41 / 2)^2), 0.0)#VectorValue(0.0, 0.0)
global u_noSlip(x, t) = VectorValue(0.0, 0.0)
global ut_in(t::Real) = x -> u_in(x, t)
global ut_noSlip(t::Real) = x -> u_noSlip(x, t)
∂tu_in(t) = x -> VectorValue(0.0, 0.0)
∂tu_in(x, t) = ∂tu_in(t)(x)
∂t(::typeof(ut_in)) = ∂tu_in
∂t(::typeof(ut_noSlip)) = ∂tu_in
# ###############

include("ElasticFlag.jl")

execute(problem::Problem; kwargs...) = @notimplemented("The driver for problem: $problem is not implemented")

# Output function
function writePVD(filePath::String, trian::Triangulation, sol; append=false)
    outfiles = paraview_collection(filePath, append=append) do pvd
        for (i, (xh, t)) in enumerate(sol)
            uh = restrict(xh.blocks[1],trian)
            vh = restrict(xh.blocks[2],trian)
						ph = restrict(xh.blocks[3],trian)
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
