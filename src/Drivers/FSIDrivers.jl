module FSIDrivers

using Gridap
using Gridap.Helpers
using Gridap.Geometry
using Gridap.Arrays
using Gridap.MultiField: ConsecutiveMultiFieldStyle
using GridapODEs.ODETools
using GridapODEs.TransientFETools
using GridapFSI.WeakForms
using GridapEmbedded
using TimerOutputs
using WriteVTK
using LineSearches: BackTracking, HagerZhang
using ForwardDiff
using Test

import GridapODEs.TransientFETools: ∂t

export Problem
export execute

struct Problem{Kind} end

include("FSI_FESpaces.jl")
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

function get_FSI_triangulations(models,coupling)
  trian = Triangulation(models[:Ω])
  trian_s = Triangulation(models[:Ωs])
  trian_f = Triangulation(models[:Ωf])
  function Γi_triangulation(coupling)
    if typeof(coupling) == WeakForms.Coupling{:weak}
      InterfaceTriangulation(models[:Ωf],models[:Ωs])
    else
      BoundaryTriangulation(models[:Ωf],"interface")
    end
  end
  trian_Γi = Γi_triangulation(coupling)
  Dict(:Ω=>trian, :Ωs=>trian_s, :Ωf=>trian_f, :Γi=>trian_Γi)
end

function get_FSI_quadratures(triangulations,order)
  degree = 2*order
  bdegree = 2*order
  quad   = CellQuadrature(triangulations[:Ω],degree)
  quad_s = CellQuadrature(triangulations[:Ωs],degree)
  quad_f = CellQuadrature(triangulations[:Ωf],degree)
  quad_Γi = CellQuadrature(triangulations[:Γi],bdegree)
  Dict(:Ω=>quad, :Ωs=>quad_s, :Ωf=>quad_f, :Γi=>quad_Γi)
end

function get_Stokes_operator(FESpaces,strategy,trian,quad,μ,f)
  res_ST(x,y) = WeakForms.stokes_residual(strategy,x,y,μ,f)
  jac_ST(x,dx,y) = WeakForms.stokes_jacobian(strategy,dx,y,μ)
  t_ST_Ωf = FETerm(res_ST, jac_ST, trian, quad)
  X,Y = FESpaces
  op_ST = FEOperator(X,Y,t_ST_Ωf)
end

function get_FSI_operator(FESpaces,coupling,strategy,Tₕ,quads,params)
  X, Y = FESpaces
  m_params, f_params, s_params, Γi_params = params

  # Compute cell area (auxiliar quantity for mesh motion eq.)
  α_m = m_params[:α]
  if( m_params[:w_strategy] == "volume")
    volf = cell_measure(Tₕ[:Ωf],Tₕ[:Ω])
    vols = cell_measure(Tₕ[:Ωs],Tₕ[:Ω])
    α_Ωf = α_m * reindex(volf,Tₕ[:Ωf])
    α_Ωs = α_m * reindex(vols,Tₕ[:Ωs])
    if ( typeof(coupling) == WeakForms.Coupling{:weak} )
      α_Γi = α_m * reindex(volf,get_left_boundary(Tₕ[:Γi]))
    else
      α_Γi = α_m * reindex(volf,Tₕ[:Γi])
    end
  else
    α_Ωf = α_m; α_Ωs = α_m; α_Γi = α_m
  end

  # Compute interface element size (if weak coupling)
  if ( typeof(coupling) == WeakForms.Coupling{:weak} )
    trian_boundary_Γi = get_left_boundary(Tₕ[:Γi])
    hΓᵢ = reindex(cell_measure(trian_boundary_Γi,Tₕ[:Ω]),trian_boundary_Γi)
  else
    hΓᵢ = 0.0
  end

  # Interface normal vector
  n_Γi = get_normal_vector(Tₕ[:Γi])

  # Complete params
  push!(f_params, :α=>α_Ωf)
  push!(f_params, :E=>m_params[:E])
  push!(f_params, :ν=>m_params[:ν])
  push!(s_params, :α=>α_Ωs)
  push!(Γi_params, :α=>α_Γi)
  push!(Γi_params, :E=>m_params[:E])
  push!(Γi_params, :ν=>m_params[:ν])
  push!(Γi_params, :n=>n_Γi)
  push!(Γi_params, :h=>hΓᵢ)

  # Define operator
  res_FSI_Ωf(t,x,xt,y) = WeakForms.fluid_residual_Ω(strategy,coupling,t,x,xt,y,f_params)
  jac_FSI_Ωf(t,x,xt,dx,y) = WeakForms.fluid_jacobian_Ω(strategy,coupling,x,xt,dx,y,f_params)
  jac_t_FSI_Ωf(t,x,xt,dxt,y) = WeakForms.fluid_jacobian_t_Ω(strategy,coupling,x,xt,dxt,y,f_params)
  res_FSI_Ωs(t,x,xt,y) = WeakForms.solid_residual_Ω(strategy,coupling,t,x,xt,y,s_params)
  jac_FSI_Ωs(t,x,xt,dx,y) = WeakForms.solid_jacobian_Ω(strategy,coupling,x,xt,dx,y,s_params)
  jac_t_FSI_Ωs(t,x,xt,dxt,y) = WeakForms.solid_jacobian_t_Ω(strategy,coupling,x,xt,dxt,y,s_params)
  res_FSI_Γi(x,y) = WeakForms.fsi_residual_Γi(strategy,coupling,x,y,Γi_params)
  jac_FSI_Γi(x,dx,y) = WeakForms.fsi_jacobian_Γi(strategy,coupling,x,dx,y,Γi_params)
  t_FSI_Ωf = FETerm(res_FSI_Ωf, jac_FSI_Ωf, jac_t_FSI_Ωf, Tₕ[:Ωf], quads[:Ωf])
  t_FSI_Ωs = FETerm(res_FSI_Ωs, jac_FSI_Ωs, jac_t_FSI_Ωs, Tₕ[:Ωs], quads[:Ωs])
  t_FSI_Γi = FETerm(res_FSI_Γi,jac_FSI_Γi,Tₕ[:Γi],quads[:Γi])
  op_FSI = TransientFEOperator(X,Y,t_FSI_Ωf,t_FSI_Ωs,t_FSI_Γi)

end

end
